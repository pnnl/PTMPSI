#!/usr/bin/env python
# run_pca_analysis.py

import argparse
import glob
import os
import cProfile
import pstats
import io
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import Parallel, delayed

import MDAnalysis as mda
from MDAnalysis import Universe, Writer
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.analysis import align

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, IncrementalPCA
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
    overall_score,
    compute_combined_score
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
    plot_score_hist=False,
    plot_top5pcs=False,
    plot_pairwise=False
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
      compute_combined_score,
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

    # ─── Global dir for cross‐system outputs ───
    global_fig_dir = output_dir
    os.makedirs(global_fig_dir, exist_ok=True)
    fig_dir = global_fig_dir
    pdb_dir = os.path.join(os.path.dirname(output_dir), "representative_structures")
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

    # STEP 2: Build the global feature‐name list & align each system’s features
    all_feat_names = sorted({
        name for names in feat_names_dict.values() for name in names
    })

    aligned_feats = {}
    for sys_name, feats in feats_dict.items():
        local_names = feat_names_dict[sys_name]
        n_frames    = feats.shape[0]
        A = np.zeros((n_frames, len(all_feat_names)))
        for j, fname in enumerate(all_feat_names):
            if fname in local_names:
                idx = local_names.index(fname)
                A[:, j] = feats[:, idx]
        aligned_feats[sys_name] = A
        # ── Drop any constant (zero-variance) features ──
        # (stack them to compute variance, then mask)
        _all = np.vstack(list(aligned_feats.values()))
        _var = _all.var(axis=0)
        _mask = _var > 0
        # filter feature names & data
        all_feat_names = [f for f, keep in zip(all_feat_names, _mask) if keep]
        for k in aligned_feats:
            aligned_feats[k] = aligned_feats[k][:, _mask] 

    # STEP 2.5: auto-select n_comp and n_clusters if not provided
    # stack all aligned features into one matrix
    all_X = np.vstack([aligned_feats[name] for name in systems_dict])

    # auto-select number of PCA components to reach 90% explained variance
    if n_comp is None:
        tmp = PCA().fit(all_X)
        cum = np.cumsum(tmp.explained_variance_ratio_)
        n_comp = int(np.searchsorted(cum, 0.90)) + 1
        print(f"Auto-selected n_comp = {n_comp} (90% variance)")

    # auto-select number of clusters by max silhouette over k=2…10
    if n_clusters is None:
        proj = PCA(n_components=n_comp).fit_transform(all_X)
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

    # STEP 3: Fit a GLOBAL IncrementalPCA on the aligned features
    scaler = StandardScaler()

    # ─── helper: only plot those global features that this system actually computed ───
    def local_features_for(sys_name):
        return feat_names_dict[sys_name]
        
    # accumulate mean/var
    for feats in aligned_feats.values():
        scaler.partial_fit(feats)

    # build components incrementally
    ipca = IncrementalPCA(n_components=n_comp)
    for feats in aligned_feats.values():
        Xs = scaler.transform(feats)
        ipca.partial_fit(Xs)

    # make ipca our “global_pca_model” so the plotting calls work
    global_pca_model = ipca
    explained_variance = ipca.explained_variance_ratio_

    # save loadings (tab‑sep with proper index name)
    df_loadings = pd.DataFrame(
        ipca.components_.T,
        index=all_feat_names,
        columns=[f"PC{i+1}" for i in range(n_comp)]
    )
    df_loadings.index.name = "Feature" 
    loadings_fn = os.path.join(output_dir, "PC_loadings.csv")
    write_pretty_csv(df_loadings.reset_index().rename(columns={"index":"Feature"}), loadings_fn)
    print(f"Global PCA loadings written to {loadings_fn}")

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


    # project each system
    sys_global_pca = {}
    for sys_name, feats in aligned_feats.items():
        Xs = scaler.transform(feats)
        sys_global_pca[sys_name] = ipca.transform(Xs)
    # ─── consistent PC1 sign (so Δ⟨PC1⟩ never flips between runs) ───
    ref_name = "REF"
    if sys_global_pca[ref_name][:, 0].mean() < 0:
        for name in sys_global_pca:
            sys_global_pca[name][:, 0] *= -1
        global_pca_model.components_[0] *= -1

    # ─── Build global z‐scaler for PC1/MFPT/τ ───
    all_blocks = []
    for name, coords in sys_global_pca.items():
        # just cluster → build T → get MFPT & τ (no scoring)
        km      = KMeans(n_clusters=n_clusters, n_init="auto", random_state=42)
        labels  = km.fit_predict(coords)
        T       = build_transition_matrix(labels, n_clusters)
        mfpt    = compute_mfpt(T)
        taus    = compute_residence_times(T)
        mfpt_fr = np.array([np.nanmean(mfpt[c, :]) for c in labels])
        tau_fr  = np.array([taus[c]                for c in labels])
        block   = np.vstack((coords[:, 0], mfpt_fr, tau_fr)).T
        all_blocks.append(block)
    global_score_X = np.vstack(all_blocks)

    # ——— remove infinities/NaNs before scaling ———
    finite = np.isfinite(global_score_X)
    if not finite.all():
        max_f = np.max(global_score_X[finite])
        min_f = np.min(global_score_X[finite])
        global_score_X = np.nan_to_num(
            global_score_X,
            nan=0.0,
            posinf=max_f,
            neginf=min_f
        )

    global_scaler = StandardScaler().fit(global_score_X)

    # write & plot global-PCA ranking
    write_and_plot_global_pca_ranking(sys_global_pca, outdir=output_dir)

    # STEP 4: For each system (ref first), cluster & score
    ref_name = "REF"
    all_scores      = {}
    relative_scores = {}
    results_dict    = {}
    #dt_ns = 0.10

    def _process_system(pca_coords, global_scaler=None):
        """Cluster, build KTN, compute MFPT/tau, compute combined scores."""
        km    = KMeans(n_clusters=n_clusters, n_init="auto", random_state=42)
        labels= km.fit_predict(pca_coords)
        T     = build_transition_matrix(labels, n_clusters)
        G     = build_ktn(T)
        mfpt  = compute_mfpt(T)
        taus  = compute_residence_times(T)
        scores = compute_combined_score(
            pc1_values     = pca_coords[:,0],
            cluster_labels = labels,
            mfpt           = mfpt,
            taus           = taus,
            global_scaler  = global_scaler,
            alpha=1, beta=1, gamma=1
        )
        mean_score = np.nanmean(scores)
        return labels, km, G, T, mfpt, taus, scores, mean_score

    # ---- Reference system ----
    ref_coords = sys_global_pca[ref_name]
    (ref_labels, ref_km, ref_G, ref_T,
     ref_mfpt, ref_taus, ref_scores,
     ref_mean) = _process_system(ref_coords, global_scaler)

    all_scores[ref_name] = ref_mean
    results_dict[ref_name] = dict(
        labels=ref_labels,
        # reload Universe only when needed
        universe=mda.Universe(
            systems_dict[ref_name]["top"],
            *systems_dict[ref_name]["traj_list"]
        ),
        frames=frames_dict[ref_name],
        ptm_matrix=ptm_matrices[ref_name],
        ptm_labels=ptm_labels_dict[ref_name],
        pc_coords=ref_coords
    )

    # Reference plots
    plot_top_pcs(ref_coords, explained_variance, ref_name)
    plot_clusters_pca(ref_coords, ref_labels, ref_km.cluster_centers_, ref_name)
    visualize_ktn(ref_G, ref_T, ref_mfpt, ref_taus, ref_name, dt_ns, ref_labels)
    plot_cluster_frequency_3d(ref_labels, ref_scores, ref_name)
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

    plot_cluster_feature_ranges(ranges_ref, local_features_for(ref_name), ref_name)

    # histogram
    if plot_score_hist:
        plt.figure()
        plt.hist(ref_scores[np.isfinite(ref_scores)], bins=30, alpha=0.7, color='C0')
        plt.xlabel("Score"); plt.ylabel("Frequency")
        plt.savefig(os.path.join(fig_dir, f"{ref_name}_score_histogram.pdf"), dpi=300)
        plt.close('all')
    # PTM importance & representatives
    if ptm_labels_dict[ref_name]:
        plot_ptm_importance(ptm_labels_dict[ref_name], ptm_matrices[ref_name], ref_name)
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
    for sys_name in sorted(systems_dict):
        if sys_name == ref_name:
            continue
        coords = sys_global_pca[sys_name]
        (labels, km, G, T, mfpt, taus, scores, mean_score
        ) = _process_system(coords, global_scaler)

        all_scores[sys_name]      = mean_score
        relative_scores[sys_name] = {
            "ratio": mean_score / ref_mean,
            "diff":  mean_score - ref_mean
        }
        results_dict[sys_name] = dict(
            labels=labels,
            universe=mda.Universe(
                systems_dict[sys_name]["top"],
                *systems_dict[sys_name]["traj_list"]
            ),
            frames=frames_dict[sys_name],
            ptm_matrix=ptm_matrices[sys_name],
            ptm_labels=ptm_labels_dict[sys_name],
            pc_coords=coords
        )

        # per‐system plots
        plot_top_pcs(coords, explained_variance, sys_name)
        plot_clusters_pca(coords, labels, km.cluster_centers_, sys_name)
        visualize_ktn(G, T, mfpt, taus, sys_name, dt_ns, labels)
        plot_cluster_frequency_3d(labels, scores, sys_name)
        # feature‐range
        feats_sys = aligned_feats[sys_name]
        ranges_sys = {}
        for c in range(n_clusters):
            idxs = np.where(labels==c)[0]
            if idxs.size:
                ranges_sys[c] = {
                    fname: (feats_sys[idxs,i].min(), feats_sys[idxs,i].max())
                    for i,fname in enumerate(all_feat_names)
                }
        plot_cluster_feature_ranges(ranges_sys, local_features_for(sys_name), sys_name)

        # histogram
        plt.figure()
        plt.hist(scores[np.isfinite(scores)], bins=30, alpha=0.7, color='C0')
        plt.xlabel("Score"); plt.ylabel("Frequency")
        plt.savefig(os.path.join(fig_dir, f"{sys_name}_score_histogram.pdf"), dpi=300)
        plt.close('all')
        # PTM importance & representatives
        if ptm_labels_dict[sys_name]:
            plot_ptm_importance(ptm_labels_dict[sys_name], ptm_matrices[sys_name], sys_name)
        write_representative_structures(
            results_dict[sys_name]['universe'], frames_dict[sys_name],
            labels, scores, sys_name,
            output_dir=pdb_dir
        )
        if write_xtcs:         
            write_cluster_trajectories(
                results_dict[sys_name]["universe"], frames_dict[sys_name],
                labels, sys_name,
                output_dir=xtc_dir
            )    

    # STEP 5: Write final CSV summary and stacked‐bar of Δ(Combined Score)
    # Compute the reference‐system mean contributions (PC1, MFPT, τ)
    ref_pc1_mean   = np.mean(sys_global_pca[ref_name][:, 0])
    ref_mfpt_means = np.mean([np.nanmean(ref_mfpt[c, :]) for c in ref_labels])
    ref_tau_means  = np.mean([ref_taus[c] for c in ref_labels])

    # Define the list of subsystems (excluding the reference)
    subs = sorted(
        (s for s in all_scores if s != ref_name),
        key=lambda s: int(s.rsplit('_', 1)[-1])
    )

    # Gather per‐system mean contributions
    pc1_means  = []
    mfpt_means = []
    tau_means  = []
    for s in subs:
        # PC1 mean
        coords   = sys_global_pca[s]
        labels_s = results_dict[s]['labels']

        # build a fresh transition matrix for this system
        T_s      = build_transition_matrix(labels_s, n_clusters)
        mfpt_s   = compute_mfpt(T_s)
        taus_s   = compute_residence_times(T_s)

        pc1_means.append(np.mean(coords[:, 0]))
        mfpt_by_frame = [np.nanmean(mfpt_s[c, :]) for c in labels_s]
        tau_by_frame  = [taus_s[c] for c in labels_s]

        mfpt_means.append(np.mean(mfpt_by_frame))
        tau_means.append(np.mean(tau_by_frame))

    pc1_means  = np.array(pc1_means)
    mfpt_means = np.array(mfpt_means)
    tau_means  = np.array(tau_means)

    # Compute deltas relative to reference, weighted by α, β, γ = 1,1,1
    del_pc1  = (pc1_means  - ref_pc1_mean)  * 1.0
    del_mfpt = (mfpt_means - ref_mfpt_means) * 1.0
    del_tau  = (tau_means  - ref_tau_means)  * 1.0

    # Write CSV with Δ-columns
    del_pc1_vals  = dict(zip(subs, del_pc1))
    del_mfpt_vals = dict(zip(subs, del_mfpt))
    del_tau_vals  = dict(zip(subs, del_tau))

    csv_fn = os.path.join(output_dir, "combined_score_ranking.csv")
    rows = [{
        "system":               ref_name,
        "mean_combined_score":  all_scores[ref_name],
        "ratio_to_ref":         1.0,
        "diff_to_ref":          0.0,
        "delta_PC1":            0.0,
        "delta_MFPT":           0.0,
        "delta_tau":            0.0
    }] + [
        {
            "system":               s,
            "mean_combined_score":  all_scores[s],
            "ratio_to_ref":         relative_scores[s]["ratio"],
            "diff_to_ref":          relative_scores[s]["diff"],
            "delta_PC1":            del_pc1_vals[s],
            "delta_MFPT":           del_mfpt_vals[s],
            "delta_tau":            del_tau_vals[s]
        }
        for s in subs
    ]
    df_comb = pd.DataFrame(rows)
    write_pretty_csv(df_comb, csv_fn)
    print(f"Detailed Combined‑score ranking CSV written to {csv_fn}")

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

    return all_scores, relative_scores, results_dict


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
    parser.add_argument("--n_clusters", type=int,   default=None,
                        help="Number of KMeans clusters (omit to auto-select)")
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
            plot_score_hist   = args.plot_score_hist,
            plot_top5pcs      = args.plot_top5pcs,
            plot_pairwise     = args.plot_pairwise
        )

    # total runtime
    end_time = time.time()
    print(f"Total script runtime: {end_time - start_time:.2f} seconds.")
