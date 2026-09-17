#!/usr/bin/env python
# run_pca_analysis.py

import argparse
import glob
import os
import sys
import cProfile
import pstats
import io
import json
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_MAX_THREADS", "128")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import Parallel, delayed

PTMPSI_PCA_DIR = Path("/scratch/path/to/scripts_perlmutter/PTMPSI/ptmpsi-PCA")
if str(PTMPSI_PCA_DIR) not in sys.path:
    sys.path.insert(0, str(PTMPSI_PCA_DIR))

import MDAnalysis as mda
from MDAnalysis import Universe, Writer
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.analysis import align

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from pca_config import time, STANDARD_NAMES
from pca_selection_features import (
    select_atoms_for_analysis,
    extract_features_mda,
    extract_and_save_features,
    write_pretty_csv,
    load_and_preprocess_subsystem_mda
)
from pca_clustering import (
    perform_pca,
    cluster_kmeans,
    build_transition_matrix,
    build_ktn,
    compute_stationary_distribution,
    compute_cluster_score,
    compute_mfpt,
    compute_residence_times,
    compute_local_score,
    overall_score
)

from pca_plotting import (
    plot_cluster_feature_ranges,
    plot_top_pcs,
    plot_clusters_pca,
    visualize_ktn,
    plot_cluster_frequency_3d,
    plot_ptm_importance,
    plot_top5_pcs_all_systems,
    plot_pairwise_top5_pcs,
    build_and_plot_residue_network,
)

from pca_config import time

import time
start_time = time.time()


def system_sort_key(name):
    if name == "REF":
        return (-1, -1, name)
    if name.isdigit():
        return (0, int(name), name)
    try:
        return (0, int(name.rsplit("_", 1)[-1]), name)
    except ValueError:
        return (1, 0, name)


def finite_mean(values):
    arr = np.asarray(values, dtype=float)
    finite = np.isfinite(arr)
    if not finite.any():
        return np.nan
    return float(arr[finite].mean())


def finite_fraction(values):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.nan
    return float(np.isfinite(arr).sum() / arr.size)


def write_feature_diagnostics(feat_names_dict, all_feat_names_raw, retained_names,
                              feature_variances, output_dir):
    rows = []
    retained = set(retained_names)
    for feature in all_feat_names_raw:
        rows.append({
            "feature": feature,
            "retained_for_pca": feature in retained,
            "variance": feature_variances.get(feature, np.nan),
            "present_in_systems": sum(feature in names for names in feat_names_dict.values()),
        })
    df_union = pd.DataFrame(rows)
    write_pretty_csv(df_union, os.path.join(output_dir, "feature_union_all.csv"))
    write_pretty_csv(
        df_union[df_union["retained_for_pca"]].reset_index(drop=True),
        os.path.join(output_dir, "feature_retained_for_pca.csv")
    )
    write_pretty_csv(
        df_union[~df_union["retained_for_pca"]].reset_index(drop=True),
        os.path.join(output_dir, "feature_dropped_zero_variance.csv")
    )

    count_rows = [
        {"system": name, "initial_feature_count": len(names)}
        for name, names in sorted(feat_names_dict.items(), key=lambda kv: system_sort_key(kv[0]))
    ]
    write_pretty_csv(
        pd.DataFrame(count_rows),
        os.path.join(output_dir, "per_system_feature_counts.csv")
    )


def n_components_for_variance(explained_variance_ratio, variance_threshold):
    ev = np.asarray(explained_variance_ratio, dtype=float)
    cum = np.cumsum(ev)
    n_comp = int(np.searchsorted(cum, variance_threshold)) + 1
    n_comp = min(n_comp, len(ev))
    return n_comp


def select_n_components_by_variance(explained_variance_ratio, variance_threshold):
    cum = np.cumsum(explained_variance_ratio)
    n_comp = n_components_for_variance(explained_variance_ratio, variance_threshold)
    print(
        f"Auto-selected n_comp = {n_comp} "
        f"({variance_threshold:.0%} cumulative variance threshold; "
        f"captured={cum[n_comp - 1]:.6f})"
    )
    return n_comp


def pca_distance_dataframes(sys_global_pca, explained_variance_abs,
                            explained_variance_ratio):
    ref_name = "REF"
    ref_mean = sys_global_pca[ref_name].mean(axis=0)
    pc_count = ref_mean.shape[0]
    var = np.asarray(explained_variance_abs, dtype=float)
    var = np.where(var > 0, var, np.nan)

    rows = []
    for name in sorted(sys_global_pca, key=system_sort_key):
        coords = sys_global_pca[name]
        mean_coords = coords.mean(axis=0)
        delta = mean_coords - ref_mean
        euclidean = float(np.linalg.norm(delta))
        variance_norm = float(np.sqrt(np.nansum((delta ** 2) / var)))
        row = {
            "system": name,
            "euclidean_distance_to_ref": euclidean,
            "variance_normalized_distance_to_ref": variance_norm,
            "n_pcs_used": pc_count,
            "cumulative_explained_variance_ratio": float(np.sum(explained_variance_ratio)),
        }
        for i, value in enumerate(mean_coords, start=1):
            row[f"mean_PC{i}"] = value
        for i, value in enumerate(delta, start=1):
            row[f"delta_PC{i}"] = value
        rows.append(row)

    df = pd.DataFrame(rows)
    ranked = df[df["system"] != ref_name].sort_values(
        "variance_normalized_distance_to_ref", ascending=False
    )
    return df, ranked.reset_index(drop=True)


def write_pca_distance_outputs(sys_global_pca, explained_variance_abs,
                               explained_variance_ratio, output_dir,
                               output_stem="pca_distance_to_ref"):
    df, ranked = pca_distance_dataframes(
        sys_global_pca,
        explained_variance_abs,
        explained_variance_ratio
    )
    write_pretty_csv(df, os.path.join(output_dir, f"{output_stem}.csv"))
    write_pretty_csv(
        ranked,
        os.path.join(output_dir, f"{output_stem}_ranked.csv")
    )
    return df, ranked


def write_full_scree_table(explained_variance_abs, explained_variance_ratio, output_dir):
    df_scree = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(explained_variance_ratio))],
        "explained_variance": explained_variance_abs,
        "explained_variance_ratio": explained_variance_ratio,
        "cumulative_explained_variance_ratio": np.cumsum(explained_variance_ratio),
    })
    fn = os.path.join(output_dir, "full_scree_variance.csv")
    write_pretty_csv(df_scree, fn)
    print(f"Full scree/variance table written to {fn}")
    return df_scree


def kneedle_pc(x, y):
    points = np.column_stack([x, y])
    start = points[0]
    end = points[-1]
    line = end - start
    denom = np.linalg.norm(line)
    if denom == 0:
        return int(x[0])
    distances = np.abs(
        line[0] * (start[1] - points[:, 1])
        - line[1] * (start[0] - points[:, 0])
    ) / denom
    return int(x[int(np.argmax(distances))])


def write_scree_plots(df_scree, output_dir, thresholds=None):
    if thresholds is None:
        thresholds = [0.80, 0.90, 0.95]

    pcs = np.arange(1, len(df_scree) + 1)
    evr = df_scree["explained_variance_ratio"].to_numpy(dtype=float)
    cum = df_scree["cumulative_explained_variance_ratio"].to_numpy(dtype=float)

    pc_all = kneedle_pc(pcs, evr)
    first_n = min(25, len(df_scree))
    pc_first = kneedle_pc(pcs[:first_n], evr[:first_n])
    pc_sub = None
    if len(df_scree) >= 25:
        pc_sub = kneedle_pc(pcs[2:25], evr[2:25])

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].plot(pcs, evr, marker="o", markersize=3, linewidth=1.2)
    axes[0].set_ylabel("Explained variance ratio")
    axes[0].set_title("PCA Scree Plot")
    axes[0].grid(alpha=0.25)

    axes[1].plot(pcs, cum, marker="o", markersize=3, linewidth=1.2, color="C1")
    axes[1].set_xlabel("Principal component")
    axes[1].set_ylabel("Cumulative explained variance")
    axes[1].set_ylim(0, 1.02)
    axes[1].grid(alpha=0.25)

    for threshold in thresholds:
        n_comp = n_components_for_variance(evr, threshold)
        for ax in axes:
            ax.axvline(n_comp, linestyle="--", linewidth=1, alpha=0.75)
        axes[1].axhline(threshold, linestyle=":", linewidth=1, alpha=0.75)
        axes[0].text(
            n_comp,
            axes[0].get_ylim()[1] * 0.92,
            f"{threshold:.0%}: PC{n_comp}",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
        )

    elbow_markers = [("all-curve", pc_all), ("first-25", pc_first)]
    if pc_sub is not None:
        elbow_markers.append(("PC3-PC25", pc_sub))
    seen = set()
    for label, pc in elbow_markers:
        if pc in seen:
            continue
        seen.add(pc)
        for ax in axes:
            ax.axvline(pc, color="black", linestyle="-.", linewidth=1, alpha=0.65)
        axes[0].text(
            pc,
            axes[0].get_ylim()[1] * 0.70,
            f"elbow {label}: PC{pc}",
            rotation=90,
            va="top",
            ha="left",
            fontsize=8,
            color="black",
        )

    fig.tight_layout()
    pdf = os.path.join(output_dir, "pca_scree_plot.pdf")
    png = os.path.join(output_dir, "pca_scree_plot.png")
    fig.savefig(pdf, dpi=300)
    fig.savefig(png, dpi=300)
    plt.close(fig)

    rows = [
        {"method": "kneedle_all_components", "elbow_pc": pc_all},
        {"method": "kneedle_first_25_components", "elbow_pc": pc_first},
    ]
    if pc_sub is not None:
        rows.append({"method": "kneedle_pc3_to_pc25", "elbow_pc": pc_sub})
    df_elbows = pd.DataFrame(rows)
    write_pretty_csv(df_elbows, os.path.join(output_dir, "pca_scree_elbow_candidates.csv"))
    print(f"Scree plots written to {pdf} and {png}")


def parse_sensitivity_thresholds(value):
    thresholds = []
    for raw in value.split(","):
        raw = raw.strip()
        if not raw:
            continue
        threshold = float(raw)
        if not 0 < threshold <= 1:
            raise ValueError(f"Sensitivity threshold must be in (0, 1]: {threshold}")
        thresholds.append(threshold)
    if not thresholds:
        raise ValueError("At least one PCA sensitivity threshold is required")
    return sorted(set(thresholds))


def threshold_label(threshold):
    pct = threshold * 100
    if abs(pct - round(pct)) < 1e-9:
        return f"{int(round(pct))}pct"
    return f"{pct:.1f}pct".replace(".", "p")


def write_pca_sensitivity_outputs(sys_global_pca_full, explained_variance_abs,
                                  explained_variance_ratio, thresholds,
                                  primary_threshold, output_dir):
    summary_rows = []
    rank_table = None
    cum = np.cumsum(explained_variance_ratio)

    for threshold in thresholds:
        n_comp = n_components_for_variance(explained_variance_ratio, threshold)
        label = threshold_label(threshold)
        sliced = {
            name: coords[:, :n_comp]
            for name, coords in sys_global_pca_full.items()
        }
        _, ranked = write_pca_distance_outputs(
            sliced,
            explained_variance_abs[:n_comp],
            explained_variance_ratio[:n_comp],
            output_dir,
            output_stem=f"pca_distance_to_ref_{label}"
        )
        ranked = ranked.copy()
        rank_col = f"rank_{label}"
        ranked[rank_col] = np.arange(1, len(ranked) + 1)
        ranks = ranked.set_index("system")[[rank_col]]
        rank_table = ranks if rank_table is None else rank_table.join(ranks, how="outer")
        summary_rows.append({
            "threshold": threshold,
            "label": label,
            "n_pcs_used": n_comp,
            "cumulative_explained_variance_ratio": float(cum[n_comp - 1]),
            "is_primary_threshold": abs(threshold - primary_threshold) < 1e-9,
            "ranking_file": f"pca_distance_to_ref_{label}_ranked.csv",
        })

    summary_df = pd.DataFrame(summary_rows)
    write_pretty_csv(summary_df, os.path.join(output_dir, "pca_sensitivity_summary.csv"))

    rank_table = rank_table.reset_index()
    write_pretty_csv(rank_table, os.path.join(output_dir, "pca_sensitivity_rank_matrix.csv"))

    stability_rows = []
    labels = [row["label"] for row in summary_rows]
    rank_indexed = rank_table.set_index("system")
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            ca = f"rank_{a}"
            cb = f"rank_{b}"
            pair = rank_indexed[[ca, cb]].dropna()
            spearman = pair[ca].corr(pair[cb], method="spearman")
            kendall = pair[ca].corr(pair[cb], method="kendall")
            top_n = min(25, len(pair))
            top100_n = min(100, len(pair))
            top_a = set(pair.nsmallest(top_n, ca).index) if top_n else set()
            top_b = set(pair.nsmallest(top_n, cb).index) if top_n else set()
            top100_a = set(pair.nsmallest(top100_n, ca).index) if top100_n else set()
            top100_b = set(pair.nsmallest(top100_n, cb).index) if top100_n else set()
            stability_rows.append({
                "threshold_a": a,
                "threshold_b": b,
                "n_systems_compared": len(pair),
                "spearman_rank_correlation": spearman,
                "kendall_rank_correlation": kendall,
                "top25_overlap_fraction": (
                    len(top_a & top_b) / top_n if top_n else np.nan
                ),
                "top100_overlap_fraction": (
                    len(top100_a & top100_b) / top100_n if top100_n else np.nan
                ),
            })

    stability_df = pd.DataFrame(stability_rows)
    write_pretty_csv(
        stability_df,
        os.path.join(output_dir, "pca_sensitivity_rank_stability.csv")
    )

    if stability_df.empty:
        min_spearman = np.nan
        min_top25 = np.nan
        strong_agreement = False
    else:
        min_spearman = stability_df["spearman_rank_correlation"].min()
        min_top25 = stability_df["top25_overlap_fraction"].min()
        strong_agreement = bool(min_spearman >= 0.90 and min_top25 >= 0.80)
    recommendation_fn = os.path.join(output_dir, "pca_sensitivity_recommendation.txt")
    with open(recommendation_fn, "w", encoding="utf-8") as handle:
        handle.write("PCA cutoff recommendation\n")
        handle.write("=========================\n\n")
        handle.write(
            f"Primary result: {primary_threshold:.0%} full PCA cumulative "
            "explained-variance cutoff.\n"
        )
        handle.write(
            "Sensitivity cutoffs tested: "
            + ", ".join(f"{t:.0%}" for t in thresholds)
            + ".\n"
        )
        handle.write(
            "Strong agreement heuristic: all pairwise Spearman rank "
            "correlations >= 0.90 and all top-25 overlap fractions >= 0.80.\n\n"
        )
        if strong_agreement:
            handle.write(
                "Recommendation: 80/90/95 rankings agree strongly by the "
                "script heuristic, so report the 90% full PCA result as the "
                "primary cutoff.\n"
            )
        else:
            handle.write(
                "Recommendation: 80/90/95 rankings diverge by the script "
                "heuristic. Inspect full_scree_variance.csv, PC_loadings.csv, "
                "and the sensitivity ranking files before choosing the final "
                "cutoff.\n"
            )
        handle.write(f"\nMinimum pairwise Spearman: {min_spearman:.6f}\n")
        handle.write(f"Minimum top-25 overlap: {min_top25:.6f}\n")
    print(f"PCA sensitivity recommendation written to {recommendation_fn}")


def frame_centroid_distances(pca_coords, labels, centroids):
    return np.array([
        np.linalg.norm(vec - centroids[label])
        for vec, label in zip(pca_coords, labels)
    ])


def append_kinetic_rows(system_name, labels, T, mfpt, taus, dt_ns,
                        system_rows, cluster_rows, transition_rows):
    labels = np.asarray(labels)
    finite_offdiag = mfpt[~np.eye(mfpt.shape[0], dtype=bool)]
    mfpt_frame = np.array([np.nanmean(mfpt[c, :]) for c in labels], dtype=float)
    tau_frame = np.array([taus[c] for c in labels], dtype=float)

    system_rows.append({
        "system": system_name,
        "n_frames": len(labels),
        "n_clusters": T.shape[0],
        "mean_mfpt_frames": finite_mean(mfpt_frame),
        "mean_mfpt_ns": finite_mean(mfpt_frame) * dt_ns,
        "mfpt_finite_fraction": finite_fraction(finite_offdiag),
        "mean_tau_frames": finite_mean(tau_frame),
        "mean_tau_ns": finite_mean(tau_frame) * dt_ns,
        "tau_finite_fraction": finite_fraction(taus),
    })

    for c in range(T.shape[0]):
        count = int(np.sum(labels == c))
        cluster_rows.append({
            "system": system_name,
            "cluster": c,
            "population_frames": count,
            "population_fraction": count / len(labels) if len(labels) else np.nan,
            "residence_time_frames": taus[c],
            "residence_time_ns": taus[c] * dt_ns if np.isfinite(taus[c]) else np.inf,
            "mean_outgoing_mfpt_frames": finite_mean(mfpt[c, :]),
            "mean_outgoing_mfpt_ns": finite_mean(mfpt[c, :]) * dt_ns,
        })

    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            transition_rows.append({
                "system": system_name,
                "from_cluster": i,
                "to_cluster": j,
                "transition_probability": T[i, j],
                "mfpt_frames": mfpt[i, j],
                "mfpt_ns": mfpt[i, j] * dt_ns if np.isfinite(mfpt[i, j]) else np.inf,
            })


def write_kinetic_outputs(system_rows, cluster_rows, transition_rows, output_dir):
    write_pretty_csv(
        pd.DataFrame(sorted(system_rows, key=lambda row: system_sort_key(row["system"]))),
        os.path.join(output_dir, "kinetic_observables_by_system.csv")
    )
    write_pretty_csv(
        pd.DataFrame(sorted(
            cluster_rows,
            key=lambda row: (system_sort_key(row["system"]), row["cluster"])
        )),
        os.path.join(output_dir, "cluster_observables.csv")
    )
    write_pretty_csv(
        pd.DataFrame(sorted(
            transition_rows,
            key=lambda row: (
                system_sort_key(row["system"]),
                row["from_cluster"],
                row["to_cluster"]
            )
        )),
        os.path.join(output_dir, "transition_observables.csv")
    )


def load_universe_for_output(info):
    u = mda.Universe(info["top"], *info["traj_list"])
    u.topology_filename = info["top"]
    return u


def process_system_outputs_worker(system_name, pca_coords, aligned_features,
                                  all_feat_names, local_feature_names, system_info,
                                  frame_indices, ptm_matrix, ptm_labels,
                                  explained_variance, n_clusters, dt_ns,
                                  fig_dir, pdb_dir, plot_score_hist,
                                  write_xtcs=False, xtc_dir=None):
    km = KMeans(n_clusters=n_clusters, n_init="auto", random_state=42)
    labels = km.fit_predict(pca_coords)
    T = build_transition_matrix(labels, n_clusters)
    G = build_ktn(T)
    mfpt = compute_mfpt(T)
    taus = compute_residence_times(T)
    scores = frame_centroid_distances(pca_coords, labels, km.cluster_centers_)

    local_system_rows = []
    local_cluster_rows = []
    local_transition_rows = []
    append_kinetic_rows(
        system_name, labels, T, mfpt, taus, dt_ns,
        local_system_rows, local_cluster_rows, local_transition_rows
    )

    plot_top_pcs(pca_coords, explained_variance, system_name, outdir=fig_dir)
    plot_clusters_pca(pca_coords, labels, km.cluster_centers_, system_name, outdir=fig_dir)
    visualize_ktn(G, T, mfpt, taus, system_name, dt_ns, labels, outdir=fig_dir)
    plot_cluster_frequency_3d(labels, scores, system_name, outdir=fig_dir)

    ranges = {}
    for c in range(n_clusters):
        idxs = np.where(labels == c)[0]
        if idxs.size:
            ranges[c] = {
                fname: (aligned_features[idxs, i].min(), aligned_features[idxs, i].max())
                for i, fname in enumerate(all_feat_names)
            }
    plot_cluster_feature_ranges(ranges, local_feature_names, system_name, outdir=fig_dir)

    if plot_score_hist:
        plt.figure()
        plt.hist(scores[np.isfinite(scores)], bins=30, alpha=0.7, color='C0')
        plt.xlabel("Score"); plt.ylabel("Frequency")
        plt.savefig(os.path.join(fig_dir, f"{system_name}_score_histogram.pdf"), dpi=300)
        plt.close('all')

    if ptm_labels:
        plot_ptm_importance(ptm_labels, ptm_matrix, system_name, outdir=fig_dir)

    u = load_universe_for_output(system_info)
    write_representative_structures(
        u, frame_indices, labels, scores, system_name, output_dir=pdb_dir
    )
    if write_xtcs:
        write_cluster_trajectories(
            u, frame_indices, labels, system_name, output_dir=xtc_dir
        )

    result_entry = {
        "labels": labels,
        "universe": None,
        "frames": frame_indices,
        "ptm_matrix": ptm_matrix,
        "ptm_labels": ptm_labels,
        "pc_coords": pca_coords,
    }
    return (
        system_name,
        result_entry,
        local_system_rows,
        local_cluster_rows,
        local_transition_rows,
    )


def write_and_plot_global_pca_ranking(sys_global_pca, outdir = "figures_pca"):
    """
    Rank each subsystem by its mean PC1 score.
    Writes:
      - PC_score_ranking.csv
      - PC_score_ranking.pdf
    """
    # Compute mean PC1 for each system (including REF)
    means = {name: np.mean(coords[:, 0])
             for name, coords in sys_global_pca.items()}
    ref_name = "REF"
    ref_mean = means[ref_name]

    # Sort only the subsystems by their numeric ID
    subs = sorted((s for s in means if s != ref_name),
                  key=lambda s: int(s.rsplit('_', 1)[-1]))
    ratios = [means[s] / ref_mean for s in subs]
    diffs  = [means[s] - ref_mean     for s in subs]

    csv_fn = os.path.join(outdir, "PC_score_ranking.csv")
    rows = [{"system": ref_name, "mean_PC1": ref_mean, "ratio_to_ref": 1.0, "diff_to_ref": 0.0}]
    rows += [
        {"system": s, "mean_PC1": means[s], "ratio_to_ref": ratios[i], "diff_to_ref": diffs[i]}
        for i, s in enumerate(subs)
    ]
    df_rank = pd.DataFrame(rows)
    write_pretty_csv(df_rank, csv_fn)

def write_representative_structures(u, frame_indices, labels, local_scores, system_name, output_dir):
    """
    Creates one PDB file containing the representative structure for each cluster.
    Before writing each representative frame, the structure is aligned to its respective
    reference topology (loaded from the original top file) using protein atoms.
    If writing a frame fails (e.g. due to out‐of‐range coordinates), the error and the full
    PDB dump for that frame are written to an error log file and also printed to the screen.
    
    Note: If alignment fails, a warning is printed and the unaligned structure is saved.
    :param u: MDAnalysis Universe.
    :param frame_indices: list of subsampled frame indices used for clustering.
    :param labels: cluster assignment array (same length as frame_indices).
    :param local_scores: per-frame local scores used to choose the representative.
    :param system_name: name prefix for the output file.
    :param output_dir: directory where the output file and error log will be written.
    """
    os.makedirs(output_dir, exist_ok=True)
    error_log_file = os.path.join(output_dir, "error_log.txt")
    
    # Determine the representative frame for each cluster.
    unique_clusters = np.unique(labels)
    reps = []
    for cl in unique_clusters:
        idxs = np.where(labels == cl)[0]
        best_idx = idxs[np.argmin(local_scores[idxs])]
        reps.append((cl, best_idx, local_scores[best_idx]))
    reps.sort(key=lambda x: x[2])
    
    out_filename = os.path.join(output_dir, f"{system_name}_representative.pdb")
    
    with Writer(out_filename, u.atoms.n_atoms) as W:
        for cl, best_idx, score in reps:
            try:
                # Select the representative frame for this cluster.
                u.trajectory[frame_indices[best_idx]]
                # Attempt to align the current frame to the reference topology.
                try:
                    # Load the reference structure using the original top file.
                    ref = Universe(u.topology_filename)
                    # Align the current frame's protein atoms to the reference protein atoms.
                    rmsd_val, _ = align.alignto(u.atoms, ref.atoms, select="name N CA C")
                    print(f"Frame {frame_indices[best_idx]} aligned to reference (RMSD = {rmsd_val:.2f}).")
                except Exception as align_err:
                    print(f"Warning: Alignment failed for frame {frame_indices[best_idx]}: {align_err}")
                # Write the (aligned or unaligned if alignment failed) representative structure.
                W.write(u.atoms)
            except Exception as e:
                err_msg = f"Error writing frame {frame_indices[best_idx]} for cluster {cl}: {e}"
                print(err_msg)
                pdb_buffer = io.StringIO()
                try:
                    with PDBWriter(pdb_buffer, u.atoms.n_atoms) as pdbW:
                        pdbW.write(u.atoms)
                    pdb_content = pdb_buffer.getvalue()
                except Exception as inner_e:
                    pdb_content = f"Could not dump PDB content: {inner_e}"
                print("Full PDB dump for problematic frame:")
                print(pdb_content)
                with open(error_log_file, "a") as f:
                    f.write("Full PDB dump for problematic frame:\n")
                    f.write(pdb_content + "\n")
    print(f"Saved representative structures to {out_filename}")

def write_cluster_trajectories(u, frame_indices, labels, system_name, output_dir="cluster_xtcs"):
    """
    Creates one .xtc file per cluster, containing all frames assigned to that cluster.

    :param u: MDAnalysis Universe.
    :param frame_indices: sub-sampled frame list used for clustering
    :param labels: cluster assignment array
    :param system_name: name prefix for output .xtc files.
    :param output_dir: directory for writing cluster_X.xtc files.
    """
    os.makedirs(output_dir, exist_ok=True)
    unique_clusters = np.unique(labels)
    
    for cl_id in unique_clusters:
        outname = f"{system_name}_cluster_{cl_id}.xtc"
        cl_outpath = os.path.join(output_dir, outname)
        
        idxs = np.where(labels == cl_id)[0]
        print(f"Cluster {cl_id}: {len(idxs)} frames -> {cl_outpath}")
        
        with mda.Writer(cl_outpath, n_atoms=u.atoms.n_atoms) as W:
            for idx in idxs:
                global_frame = frame_indices[idx]
                u.trajectory[global_frame]
                W.write(u.atoms)

def run_analysis_for_systems(
    systems_dict,
    output_dir,
    subsample_step,
    write_xtcs,
    n_comp=None,
    n_clusters=None,
    precomputed_feats=None,
    precomputed_names=None,
    plot_wheel=False,
    dt_ns=0.10,
    radius_ptm=3.0,
    radius_nad=3.0,
    radius_if=2.0,
    wheel_cutoff=3.0,
    dssp_window=25,
    contact_cutoff=3.0,
    variance_threshold=0.90,
    sensitivity_thresholds=None,
    per_system_jobs=32,
    plot_score_hist=False,
    plot_top5pcs=False,
    plot_pairwise=False,
    exclude_feature_prefixes=None
):
    """
    Runs the full pipeline on a dictionary of systems (including
    'REF' as reference) but *all* clustering/PCA is done
    in one GLOBAL embedding so that systems are directly comparable.

    Arguments:
      - systems_dict: mapping system_name -> info dict (must include 'no_PTMs')
      - output_dir: directory where all‐system CSVs and plots will be written
      - subsample_step: int, frames to skip when loading
      - n_comp: int, number of PCA components
      - n_clusters: int, KMeans clusters
      - write_xtcs: bool, whether to dump one XTC per cluster

    Assumes the helper functions from Cells 1–7 are already defined:
      load_and_preprocess_subsystem_mda,
      select_atoms_for_analysis,
      extract_features_mda,
      build_transition_matrix,
      build_ktn,
      compute_mfpt,
      compute_residence_times,
      plot_top_pcs,
      plot_clusters_pca,
      visualize_ktn,
      plot_cluster_frequency_3d,
      plot_cluster_feature_ranges,
      plot_ptm_importance,
      write_representative_structures,
      write_cluster_trajectories,
      build_and_plot_residue_network,
      plot_top5_pcs_all_systems,
      plot_pairwise_top5_pcs
    """
    if sensitivity_thresholds is None:
        sensitivity_thresholds = [0.80, 0.90, 0.95]
    exclude_feature_prefixes = tuple(exclude_feature_prefixes or [])

    # ─── Global dir for cross‐system outputs ───
    global_fig_dir = output_dir
    os.makedirs(global_fig_dir, exist_ok=True)
    fig_dir = global_fig_dir
    pdb_dir = os.path.join(output_dir, "representative_structures")
    os.makedirs(pdb_dir, exist_ok=True)

    # STEP 0: initialize bookkeeping dicts
    feats_dict       = {}  # raw feature arrays
    feat_names_dict  = {}  # column names for each system
    ptm_matrices     = {}  # per‐system PTM‐contact matrices
    ptm_labels_dict  = {}  # labels of PTM sites per system
    sel_info_dict    = {}  # atom‐selection dict (for PTM plots)
    frames_dict      = {}  # sub‐sampled frame lists

    # STEP 1: Extract (or load) features for every system into per‐system matrices
    feat_dir = os.path.join(output_dir, "feature_outputs")

    for sys_name, info in systems_dict.items():

        if precomputed_feats is not None and sys_name in precomputed_feats:
            # ─── Load precomputed features from disk ───
            feats      = precomputed_feats[sys_name]
            feat_names = precomputed_names[sys_name]
            ptm_mat    = np.load(os.path.join(feat_dir, f"{sys_name}_ptm_matrix.npy"))
            with open(os.path.join(feat_dir, f"{sys_name}_ptm_labels.json")) as f:
                ptm_labels = json.load(f)

            # Load Universe & subsampled frames
            u, frame_indices = load_and_preprocess_subsystem_mda(
                sys_name, info, subsample_step
            )
            # Build atom‐selection dict for later plotting
            sel_info = select_atoms_for_analysis(
                u,
                ptm_resnames      = info.get("ptm_resnames","AUTO"),
                radius_ptm        = radius_ptm,
                radius_nad        = radius_nad,
                radius_interfaces = radius_if
            )
        else:
            # ─── Original extraction path ───
            # Load Universe & subsampled frames
            u, frame_indices = load_and_preprocess_subsystem_mda(
                sys_name, info, subsample_step
            )
            # Build atom‐selection dict
            sel_info = select_atoms_for_analysis(
                u,
                ptm_resnames      = info.get("ptm_resnames","AUTO"),
                radius_ptm        = radius_ptm,
                radius_nad        = radius_nad,
                radius_interfaces = radius_if
            )
            # Extract features
            feats, feat_names, ptm_mat, ptm_labels = extract_features_mda(
                u,
                sel_info["final"],
                frame_indices,
                compute_pairwise_dist=False,
                native_contact_cutoff=contact_cutoff,
                dssp_window=dssp_window,
                ptm_resnames_for_contacts=info.get("ptm_resnames", [])
            )

        # PTM‐wheel network for this system (optional)
        if plot_wheel:
            build_and_plot_residue_network(
                u,
                sel_info,
                frame_indices,
                sys_name,
                group="ptm",
                cutoff=wheel_cutoff,
                outdir=fig_dir
            )

        # Common bookkeeping
        feats_dict[sys_name]      = feats
        feat_names_dict[sys_name] = feat_names
        ptm_matrices[sys_name]    = ptm_mat
        ptm_labels_dict[sys_name] = ptm_labels
        sel_info_dict[sys_name]   = sel_info
        frames_dict[sys_name]     = frame_indices

    # STEP 2: Build the full global feature-name union, align every system
    # to that union, then drop globally constant features once.
    all_feat_names_raw = sorted({
        name for names in feat_names_dict.values() for name in names
    })
    if exclude_feature_prefixes:
        included = [
            name for name in all_feat_names_raw
            if not name.startswith(exclude_feature_prefixes)
        ]
        excluded = [
            name for name in all_feat_names_raw
            if name.startswith(exclude_feature_prefixes)
        ]
        if not included:
            raise RuntimeError(
                "Feature-prefix exclusion removed all PCA features: "
                + ",".join(exclude_feature_prefixes)
            )
        print(
            "Excluding features from PCA by prefix "
            f"{exclude_feature_prefixes}: {len(excluded)} removed, {len(included)} remain"
        )
        write_pretty_csv(
            pd.DataFrame({"feature": excluded}),
            os.path.join(output_dir, "feature_excluded_by_prefix.csv")
        )
        all_feat_names_raw = included

    aligned_feats = {}
    for sys_name, feats in feats_dict.items():
        local_names = feat_names_dict[sys_name]
        n_frames    = feats.shape[0]
        A = np.zeros((n_frames, len(all_feat_names_raw)))
        local_index = {fname: idx for idx, fname in enumerate(local_names)}
        for j, fname in enumerate(all_feat_names_raw):
            if fname in local_index:
                idx = local_index[fname]
                A[:, j] = feats[:, idx]
        aligned_feats[sys_name] = A

    _all = np.vstack([aligned_feats[name] for name in systems_dict])
    _var = _all.var(axis=0)
    _mask = _var > 0
    feature_variances = {
        feature: variance
        for feature, variance in zip(all_feat_names_raw, _var)
    }
    all_feat_names = [f for f, keep in zip(all_feat_names_raw, _mask) if keep]
    if not all_feat_names:
        raise RuntimeError("No nonconstant features remain after global feature alignment.")
    for k in aligned_feats:
        aligned_feats[k] = aligned_feats[k][:, _mask]

    write_feature_diagnostics(
        feat_names_dict,
        all_feat_names_raw,
        all_feat_names,
        feature_variances,
        output_dir
    )

    # STEP 2.5: scale features, then auto-select n_comp if not provided.
    all_X = np.vstack([aligned_feats[name] for name in systems_dict])
    scaler = StandardScaler()
    scaler.fit(all_X)
    all_X_scaled = scaler.transform(all_X)

    full_pca_model = PCA()
    full_pca_model.fit(all_X_scaled)
    full_explained_variance = full_pca_model.explained_variance_ratio_
    full_explained_variance_abs = full_pca_model.explained_variance_
    write_full_scree_table(
        full_explained_variance_abs,
        full_explained_variance,
        output_dir
    )
    write_scree_plots(
        pd.read_csv(os.path.join(output_dir, "full_scree_variance.csv")),
        output_dir,
        sensitivity_thresholds
    )

    if n_comp is None:
        n_comp = select_n_components_by_variance(
            full_explained_variance,
            variance_threshold
        )
    else:
        n_comp = min(n_comp, len(all_feat_names))
        print(f"Using requested n_comp = {n_comp}")

    # auto-select number of clusters by max silhouette over k=2…10
    if n_clusters is None:
        proj = PCA(n_components=n_comp).fit_transform(all_X_scaled)
        def score_for_k(k):
            labs = KMeans(n_clusters=k, n_init="auto", random_state=42).fit_predict(proj)
            return silhouette_score(proj, labs)

        ks = list(range(2, min(11, len(proj))))
        sils = Parallel(n_jobs=-1)(
            delayed(score_for_k)(k) for k in ks
        )
        best_idx = int(np.argmax(sils))
        n_clusters = ks[best_idx]
        best_s = sils[best_idx]
        print(f"Auto-selected n_clusters = {n_clusters} (silhouette={best_s:.3f})")

    # STEP 3: Fit a GLOBAL PCA on the same scaled matrix used for PC selection.
    # ─── helper: only plot those global features that this system actually computed ───
    def local_features_for(sys_name):
        return feat_names_dict[sys_name]

    pca_model = full_pca_model

    # make pca_model our “global_pca_model” so the plotting calls work
    global_pca_model = pca_model
    explained_variance = full_explained_variance[:n_comp]

    # project each system into the full PCA space once, then slice for the
    # primary 90% result and sensitivity cutoffs.
    sys_global_pca_full = {}
    for sys_name, feats in aligned_feats.items():
        Xs = scaler.transform(feats)
        sys_global_pca_full[sys_name] = pca_model.transform(Xs)

    # ─── consistent PC1 sign (so Δ⟨PC1⟩ never flips between runs) ───
    ref_name = "REF"
    if sys_global_pca_full[ref_name][:, 0].mean() < 0:
        for name in sys_global_pca_full:
            sys_global_pca_full[name][:, 0] *= -1
        global_pca_model.components_[0] *= -1

    # save loadings (tab‑sep with proper index name)
    df_loadings = pd.DataFrame(
        pca_model.components_[:n_comp].T,
        index=all_feat_names,
        columns=[f"PC{i+1}" for i in range(n_comp)]
    )
    df_loadings.index.name = "Feature" 
    loadings_fn = os.path.join(output_dir, "PC_loadings.csv")
    write_pretty_csv(df_loadings.reset_index().rename(columns={"index":"Feature"}), loadings_fn)
    print(f"Global PCA loadings written to {loadings_fn}")

    explained_var_fn = os.path.join(output_dir, "explained_variance.csv")
    df_exp_var = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(explained_variance))],
        "explained_variance": full_explained_variance_abs[:n_comp],
        "explained_variance_ratio": explained_variance,
        "cumulative_explained_variance_ratio": np.cumsum(explained_variance),
    })
    write_pretty_csv(df_exp_var, explained_var_fn)
    print(f"Explained variance ratios written to {explained_var_fn}")

    # ——— Compute & write top5 PC1 features per system ———
    df_mean_feats = pd.DataFrame(
        { sys: feats.mean(axis=0) for sys, feats in aligned_feats.items() },
        index=all_feat_names
    ).T
    df_mean_feats.index.name  = "system"
    df_mean_feats.columns.name = "feature"

    contrib_rows = []
    for sys in df_mean_feats.index:
        for feat in df_mean_feats.columns:
            mean_val = df_mean_feats.at[sys, feat]
            loading  = df_loadings.at[feat, "PC1"]
            contrib_rows.append((sys, feat, mean_val * loading))
    df_contrib = pd.DataFrame(contrib_rows, columns=["system","feature","contribution"])
    df_contrib["abs_c"] = df_contrib["contribution"].abs()

    top5 = (
        df_contrib
        .sort_values(["system","abs_c"], ascending=[True, False])
        .groupby("system", group_keys=False)
        .head(5)
    )

    rows = []
    for sys, grp in top5.groupby("system"):
        grp = grp.sort_values("abs_c", ascending=False)
        feats    = grp["feature"].tolist()
        contribs = grp["contribution"].tolist()
        while len(feats) < 5:
            feats.append("")
            contribs.append(0.0)
        cells = [f"{feats[i]}\n{contribs[i]:.6f}" for i in range(5)]
        rows.append({
            "system":    sys,
            "PC":        "PC1",
            "feature_1": cells[0],
            "feature_2": cells[1],
            "feature_3": cells[2],
            "feature_4": cells[3],
            "feature_5": cells[4],
        })

    df_top5_pc1 = pd.DataFrame(rows, columns=[
        "system","PC","feature_1","feature_2","feature_3","feature_4","feature_5"
    ])
    top5_fn = os.path.join(output_dir, "system_top5-PC1_contributions.csv")
    write_pretty_csv(df_top5_pc1, top5_fn)
    print(f"Top 5 PC1 contributions CSV written to {top5_fn}")

    sys_global_pca = {
        name: coords[:, :n_comp]
        for name, coords in sys_global_pca_full.items()
    }

    write_pca_distance_outputs(
        sys_global_pca,
        full_explained_variance_abs[:n_comp],
        explained_variance,
        output_dir
    )

    write_pca_sensitivity_outputs(
        sys_global_pca_full,
        full_explained_variance_abs,
        full_explained_variance,
        sensitivity_thresholds,
        variance_threshold,
        output_dir
    )

    # Preserve the legacy PC1-only ranking as a diagnostic, not the primary ranking.
    write_and_plot_global_pca_ranking(sys_global_pca, outdir=output_dir)

    # STEP 4: For each system (ref first), cluster and compute KTN observables.
    ref_name = "REF"
    system_kinetic_rows = []
    cluster_kinetic_rows = []
    transition_kinetic_rows = []
    results_dict    = {}
    #dt_ns = 0.10

    def _process_system(system_name, pca_coords):
        """Cluster, build KTN, compute MFPT/tau, and centroid-distance frame scores."""
        km    = KMeans(n_clusters=n_clusters, n_init="auto", random_state=42)
        labels= km.fit_predict(pca_coords)
        T     = build_transition_matrix(labels, n_clusters)
        G     = build_ktn(T)
        mfpt  = compute_mfpt(T)
        taus  = compute_residence_times(T)
        frame_scores = frame_centroid_distances(pca_coords, labels, km.cluster_centers_)
        local_system_rows = []
        local_cluster_rows = []
        local_transition_rows = []
        append_kinetic_rows(
            system_name, labels, T, mfpt, taus, dt_ns,
            local_system_rows, local_cluster_rows, local_transition_rows
        )
        return (
            labels, km, G, T, mfpt, taus, frame_scores,
            local_system_rows, local_cluster_rows, local_transition_rows
        )

    # ---- Reference system ----
    ref_coords = sys_global_pca[ref_name]
    (ref_labels, ref_km, ref_G, ref_T,
     ref_mfpt, ref_taus, ref_scores,
     ref_system_rows, ref_cluster_rows, ref_transition_rows) = _process_system(ref_name, ref_coords)
    system_kinetic_rows.extend(ref_system_rows)
    cluster_kinetic_rows.extend(ref_cluster_rows)
    transition_kinetic_rows.extend(ref_transition_rows)

    results_dict[ref_name] = dict(
        labels=ref_labels,
        # reload Universe only when needed
        universe=load_universe_for_output(systems_dict[ref_name]),
        frames=frames_dict[ref_name],
        ptm_matrix=ptm_matrices[ref_name],
        ptm_labels=ptm_labels_dict[ref_name],
        pc_coords=ref_coords
    )

    # Reference plots
    plot_top_pcs(ref_coords, explained_variance, ref_name, outdir=fig_dir)
    plot_clusters_pca(ref_coords, ref_labels, ref_km.cluster_centers_, ref_name, outdir=fig_dir)
    visualize_ktn(ref_G, ref_T, ref_mfpt, ref_taus, ref_name, dt_ns, ref_labels, outdir=fig_dir)
    plot_cluster_frequency_3d(ref_labels, ref_scores, ref_name, outdir=fig_dir)
    # ─── compute feature-ranges per cluster for the reference ───
    ranges_ref = {}
    feats_ref  = aligned_feats[ref_name]
    for c in range(n_clusters):
        idxs = np.where(ref_labels == c)[0]
        if idxs.size:
            ranges_ref[c] = {
                fname: (feats_ref[idxs, i].min(), feats_ref[idxs, i].max())
                for i, fname in enumerate(all_feat_names)
            }

    plot_cluster_feature_ranges(ranges_ref, local_features_for(ref_name), ref_name, outdir=fig_dir)

    # histogram
    if plot_score_hist:
        plt.figure()
        plt.hist(ref_scores[np.isfinite(ref_scores)], bins=30, alpha=0.7, color='C0')
        plt.xlabel("Score"); plt.ylabel("Frequency")
        plt.savefig(os.path.join(fig_dir, f"{ref_name}_score_histogram.pdf"), dpi=300)
        plt.close('all')
    # PTM importance & representatives
    if ptm_labels_dict[ref_name]:
        plot_ptm_importance(ptm_labels_dict[ref_name], ptm_matrices[ref_name], ref_name, outdir=fig_dir)
    write_representative_structures(
        results_dict[ref_name]['universe'], frames_dict[ref_name],
        ref_labels, ref_scores, ref_name,
        output_dir=pdb_dir
    )
    if write_xtcs:
        write_cluster_trajectories(
            results_dict[ref_name]["universe"], frames_dict[ref_name],
            ref_labels, ref_name,
            output_dir=xtc_dir
        )

    # ---- Other systems ----
    nonref_systems = [
        name for name in sorted(systems_dict, key=system_sort_key)
        if name != ref_name
    ]
    if per_system_jobs and per_system_jobs > 1:
        worker_count = min(per_system_jobs, len(nonref_systems))
        print(f"Processing {len(nonref_systems)} non-reference systems with {worker_count} parallel workers")
        nonref_results = Parallel(n_jobs=worker_count)(
            delayed(process_system_outputs_worker)(
                sys_name,
                sys_global_pca[sys_name],
                aligned_feats[sys_name],
                all_feat_names,
                local_features_for(sys_name),
                systems_dict[sys_name],
                frames_dict[sys_name],
                ptm_matrices[sys_name],
                ptm_labels_dict[sys_name],
                explained_variance,
                n_clusters,
                dt_ns,
                fig_dir,
                pdb_dir,
                plot_score_hist,
                write_xtcs,
                xtc_dir if write_xtcs else None,
            )
            for sys_name in nonref_systems
        )
    else:
        print(f"Processing {len(nonref_systems)} non-reference systems serially")
        nonref_results = [
            process_system_outputs_worker(
                sys_name,
                sys_global_pca[sys_name],
                aligned_feats[sys_name],
                all_feat_names,
                local_features_for(sys_name),
                systems_dict[sys_name],
                frames_dict[sys_name],
                ptm_matrices[sys_name],
                ptm_labels_dict[sys_name],
                explained_variance,
                n_clusters,
                dt_ns,
                fig_dir,
                pdb_dir,
                plot_score_hist,
                write_xtcs,
                xtc_dir if write_xtcs else None,
            )
            for sys_name in nonref_systems
        ]

    for (sys_name, result_entry,
         local_system_rows, local_cluster_rows, local_transition_rows) in nonref_results:
        results_dict[sys_name] = result_entry
        system_kinetic_rows.extend(local_system_rows)
        cluster_kinetic_rows.extend(local_cluster_rows)
        transition_kinetic_rows.extend(local_transition_rows)

    # STEP 5: Write separate kinetic/KTN observables. Combined score output is
    # intentionally not produced by this fixed script.
    write_kinetic_outputs(
        system_kinetic_rows,
        cluster_kinetic_rows,
        transition_kinetic_rows,
        output_dir
    )

    # STEP 6: Global summary plots across all systems
    # prepare input dicts
    pca_data_dict = {
        name: (sys_global_pca[name], explained_variance)
        for name in systems_dict
    }
    cluster_labels_dict = {
        name: results_dict[name]["labels"]
        for name in systems_dict
    }

    if plot_top5pcs:
        plot_top5_pcs_all_systems(
            pca_data_dict,
            outdir=output_dir,
            filename="all_systems_top5_PCs.pdf"
        )
    if plot_pairwise:
        plot_pairwise_top5_pcs(
            pca_data_dict,
            cluster_labels_dict,
            outdir=output_dir,
            filename="pairwise_top5_PCs.pdf"
        )

    return sys_global_pca, results_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PTM pipeline over a set of subsystem folders")
    parser.add_argument("--base_dir", required=True,
                        help="Base directory containing all subsystem folders and the REF folder")
    parser.add_argument("--ref_dir", required=True,
                        help="Directory of your reference subsystem (e.g. POST_BASE/REF)")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write all‐system CSVs and global plots")
    parser.add_argument("--topology", required=True,
                        help="Topology filename in each subsystem folder (e.g. solute_fixed.pdb)")
    parser.add_argument("--trajectories", required=True, nargs="+",
                        help="Trajectory file patterns in each subsystem folder (e.g. solute_fit.xtc)")
    parser.add_argument("--step",        type=int,   default=10,
                        help="Frame subsampling step (default: %(default)s)")
    parser.add_argument("--n_comp",     type=int,   default=None,
                        help="Number of PCA components (omit to auto-select)")
    parser.add_argument("--pca_variance_threshold", type=float, default=0.90,
                        help="Cumulative variance threshold for auto-selected PCs")
    parser.add_argument("--pca_sensitivity_thresholds", default="0.80,0.90,0.95",
                        help="Comma-separated cumulative variance cutoffs for sensitivity rankings")
    parser.add_argument("--exclude_feature_prefixes", default="",
                        help="Comma-separated feature prefixes to exclude from PCA after feature loading/alignment")
    parser.add_argument("--n_clusters", type=int,   default=None,
                        help="Number of KMeans clusters (omit to auto-select)")
    parser.add_argument("--per_system_jobs", type=int, default=32,
                        help="Parallel workers for independent per-system post-PCA analysis")
    parser.add_argument("--write-xtcs",  action="store_true",
                        help="Write one .xtc per cluster")
    # ─── flags for all of your magic numbers ───
    parser.add_argument("--dt_ns",           type=float, default=0.10,
                        help="Time between frames in ns")
    parser.add_argument("--radius_ptm",      type=float, default=3.0,
                        help="Radius around PTM (Å)")
    parser.add_argument("--radius_nad",      type=float, default=3.0,
                        help="Radius around NAD (Å)")
    parser.add_argument("--radius_if",       type=float, default=2.0,
                        help="Interface‐atom radius (Å)")
    parser.add_argument("--dssp_window",     type=int,   default=25,
                        help="Half‐window for sliding B-factor")
    parser.add_argument("--contact_cutoff",  type=float, default=3.0,
                        help="Native contact cutoff (Å)")
    parser.add_argument("--wheel_cutoff",    type=float, default=3.0,
                        help="PTM‐wheel network cutoff (Å)")
    parser.add_argument("--plot_score_hist", action="store_true",
                        help="Enable per-system score‐histogram plots")
    parser.add_argument("--plot_top5pcs",    action="store_true",
                        help="Enable combined Top-5-PCs plot")
    parser.add_argument("--plot_pairwise",   action="store_true",
                        help="Enable pairwise Top-5 PCs scatter plots")
    parser.add_argument("--plot-wheel",      action="store_true",
                        help="Enable PTM-wheel network plots")
    parser.add_argument("--profile", action="store_true",
                        help="Profile the main analysis and print top-20 slowest calls")
    args = parser.parse_args()
    sensitivity_thresholds = parse_sensitivity_thresholds(args.pca_sensitivity_thresholds)
    exclude_feature_prefixes = [
        x.strip() for x in args.exclude_feature_prefixes.split(",")
        if x.strip()
    ]

    # ── Build systems dict from CLI args ──
    systems = {}
    # reference system, name it after the basename of --ref_dir
    ref_name = os.path.basename(os.path.normpath(args.ref_dir))
    ref_trajs = []
    for pat in args.trajectories:
        ref_trajs.extend(glob.glob(os.path.join(args.ref_dir, pat)))
    systems[ref_name] = {
        "category":     ref_name,
        "ptm_resnames": [],
        "ptm_type":     ref_name,
        "top":          os.path.join(args.ref_dir, args.topology),
        "traj_list":    ref_trajs
    }

    for sys_dir in sorted(glob.glob(os.path.join(args.base_dir, "*"))):
        # skip whichever folder was passed as --ref_dir
        if os.path.abspath(sys_dir) == os.path.abspath(args.ref_dir):
            continue
        name    = os.path.basename(sys_dir)
        topo_path = os.path.join(sys_dir, args.topology)
        trajs = []
        for pat in args.trajectories:
            trajs.extend(glob.glob(os.path.join(sys_dir, pat)))
        if os.path.isfile(topo_path) and trajs:
            systems[name] = {
                "category":     name,
                "ptm_resnames": "AUTO",
                "ptm_type":     name,
                "top":          topo_path,
                "traj_list":    trajs
            }

    # directory to dump features
    feat_dir = os.path.join(args.output_dir, "feature_outputs")
    os.makedirs(feat_dir, exist_ok=True)

    # parallel feature dumps
    def safe_extract(name, info,
                     step, feat_dir,
                     radius_ptm, radius_nad, radius_if,
                     dssp_window, contact_cutoff):
        out_np = os.path.join(feat_dir, f"{name}_feats.npy")
        if os.path.exists(out_np):
            return name
        return extract_and_save_features(
            name, info, step, feat_dir,
            radius_ptm, radius_nad,
            radius_if, dssp_window,
            contact_cutoff
        )

    feat_names = Parallel(n_jobs=-1)(
        delayed(safe_extract)(
            name, info,
            args.step,
            feat_dir,
            args.radius_ptm,
            args.radius_nad,
            args.radius_if,
            args.dssp_window,
            args.contact_cutoff
        )
        for name, info in systems.items()
    )

    # load back feature arrays & names
    precomp_feats = {}
    precomp_names = {}
    for name in feat_names:
        precomp_feats[name]  = np.load(os.path.join(feat_dir, f"{name}_feats.npy"))
        precomp_names[name]  = json.load(open(os.path.join(feat_dir, f"{name}_featnames.json")))

    # run main analysis (optionally under cProfile)
    if args.profile:
        profiler = cProfile.Profile()
        profiler.enable()

        # your original call, unmodified:
        run_analysis_for_systems(
            systems_dict      = systems,
            output_dir        = args.output_dir,
            subsample_step    = args.step,
            n_comp            = args.n_comp,
            n_clusters        = args.n_clusters,
            write_xtcs        = args.write_xtcs,
            precomputed_feats = precomp_feats,
            precomputed_names = precomp_names,
            plot_wheel        = args.plot_wheel,
            dt_ns             = args.dt_ns,
            radius_ptm        = args.radius_ptm,
            radius_nad        = args.radius_nad,
            radius_if         = args.radius_if,
            wheel_cutoff      = args.wheel_cutoff,
            dssp_window       = args.dssp_window,
            contact_cutoff    = args.contact_cutoff,
            variance_threshold= args.pca_variance_threshold,
            sensitivity_thresholds=sensitivity_thresholds,
            exclude_feature_prefixes=exclude_feature_prefixes,
            per_system_jobs   = args.per_system_jobs,
            plot_score_hist   = args.plot_score_hist,
            plot_top5pcs      = args.plot_top5pcs,
            plot_pairwise     = args.plot_pairwise
        )

        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats("cumtime")
        stats.print_stats(20)   # top 20 hotspots
    else:
        # same call when not profiling
        run_analysis_for_systems(
            systems_dict      = systems,
            output_dir        = args.output_dir,
            subsample_step    = args.step,
            n_comp            = args.n_comp,
            n_clusters        = args.n_clusters,
            write_xtcs        = args.write_xtcs,
            precomputed_feats = precomp_feats,
            precomputed_names = precomp_names,
            plot_wheel        = args.plot_wheel,
            dt_ns             = args.dt_ns,
            radius_ptm        = args.radius_ptm,
            radius_nad        = args.radius_nad,
            radius_if         = args.radius_if,
            wheel_cutoff      = args.wheel_cutoff,
            dssp_window       = args.dssp_window,
            contact_cutoff    = args.contact_cutoff,
            variance_threshold= args.pca_variance_threshold,
            sensitivity_thresholds=sensitivity_thresholds,
            exclude_feature_prefixes=exclude_feature_prefixes,
            per_system_jobs   = args.per_system_jobs,
            plot_score_hist   = args.plot_score_hist,
            plot_top5pcs      = args.plot_top5pcs,
            plot_pairwise     = args.plot_pairwise
        )

    # total runtime
    end_time = time.time()
    print(f"Total script runtime: {end_time - start_time:.2f} seconds.")
