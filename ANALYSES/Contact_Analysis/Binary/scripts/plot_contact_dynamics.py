#!/usr/bin/env python3
"""
Plots for binary contact dynamics.

Reads analyze_contact_dynamics.py outputs in `contact_dynamics/`:
  per_pair_stats.csv.gz
  per_pair_hmm_aggregated.csv
  per_class_summary.csv
  condition_summary.csv

Produces:
  contact_dynamics_plots/
    fig_persistence_heatmap_<partner>.{pdf,png}      mean dwell time per HMM col
    fig_kinetic_heatmap_<partner>.{pdf,png}          k_off vs k_on per HMM col
    fig_class_violins.{pdf,png}                       dwell-time distribution by class
    fig_class_kinetics_strip.{pdf,png}                turnover & k_off per class
    fig_host_redox_swarm_<partner>.{pdf,png}          per-HMM dwell points across 4 conditions
    fig_survival_curves_<partner>.{pdf,png}           median S(t) by condition
    contact_dynamics_text_summary.md                  numbers for the manuscript
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


CONDITION_ORDER = ["med4_ox", "med4_red", "hm2_ox", "hm2_red"]
CONDITION_COLORS = {
    "med4_ox":  "#1f77b4",
    "med4_red": "#aec7e8",
    "hm2_ox":   "#d62728",
    "hm2_red":  "#ff9896",
}
CLASS_ORDER = ["hydrophobic", "polar_hbond", "salt_bridge", "cation_pi", "mixed"]
CLASS_COLORS = {
    "hydrophobic": "#666666",
    "polar_hbond": "#2ca02c",
    "salt_bridge": "#9467bd",
    "cation_pi":   "#ff7f0e",
    "mixed":       "#bbbbbb",
}


def savefig(fig: plt.Figure, base: Path, dpi: int = 300) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=dpi)
    plt.close(fig)


def make_persistence_heatmap(per_hmm: pd.DataFrame, partner: str, outdir: Path) -> None:
    df = per_hmm[per_hmm["partner"] == partner].copy()
    if df.empty:
        return
    # mean dwell per (condition, hmm_col) — averaged across contact classes within position
    pivot_mean = df.groupby(["condition", "hmm_col"], as_index=False).agg(
        mean_on_ps=("mean_on_ps_mean", "mean")
    )
    mat = pivot_mean.pivot(index="condition", columns="hmm_col", values="mean_on_ps")
    mat = mat.reindex(CONDITION_ORDER)
    all_cols = list(range(1, int(per_hmm["hmm_col"].max()) + 1))
    mat = mat.reindex(columns=all_cols)

    fig, ax = plt.subplots(figsize=(14, 2.8))
    # log-scale color since dwell distribution spans 2-3 decades
    arr = mat.to_numpy()
    arr_log = np.log10(np.clip(arr, 1, None))  # mask <1ps to 0
    im = ax.imshow(arr_log, aspect="auto", cmap="viridis",
                   interpolation="none", origin="upper",
                   extent=(0.5, len(all_cols) + 0.5, 4 - 0.5, -0.5))
    ax.set_yticks(range(len(CONDITION_ORDER)))
    ax.set_yticklabels(CONDITION_ORDER)
    ax.set_xticks(list(range(5, len(all_cols) + 1, 5)))
    ax.set_xlabel("CP12 HMM alignment position")
    ax.set_title(f"Mean contact dwell time per HMM column — CP12 vs {partner.upper()} (log10 ps)")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("log10 mean τ_on (ps)")
    fig.tight_layout()
    savefig(fig, outdir / f"fig_persistence_heatmap_{partner}")


def make_kinetic_scatter(per_hmm: pd.DataFrame, partner: str, outdir: Path) -> None:
    """Scatter of k_off vs k_on per HMM column, colored by condition."""
    df = per_hmm[per_hmm["partner"] == partner].copy()
    if df.empty:
        return
    # collapse multi-class rows to per (condition, hmm_col)
    df = df.groupby(["condition", "hmm_col"], as_index=False).agg(
        k_off=("k_off_per_ns_mean", "mean"),
        k_on=("k_on_per_ns_mean", "mean"),
        occ=("occupancy_mean", "mean"),
    )
    fig, ax = plt.subplots(figsize=(6, 5.5))
    for cond in CONDITION_ORDER:
        sub = df[df["condition"] == cond]
        ax.scatter(sub["k_on"], sub["k_off"], s=24, alpha=0.7,
                   color=CONDITION_COLORS[cond], edgecolor="white", linewidth=0.4,
                   label=cond)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k_on (per ns)")
    ax.set_ylabel("k_off (per ns) — lower = stronger contact")
    ax.set_title(f"Two-state kinetics per HMM column — CP12 vs {partner.upper()}")
    ax.legend(frameon=False, fontsize=8, loc="best")
    # reference: occupancy = k_on / (k_on + k_off) -> diagonal of occupancy=0.5
    xs = np.array(ax.get_xlim())
    ax.plot(xs, xs, "k--", lw=0.5, alpha=0.4, label="occ=0.5")
    fig.tight_layout()
    savefig(fig, outdir / f"fig_kinetic_scatter_{partner}")


def make_class_violins(per_pair: pd.DataFrame, outdir: Path) -> None:
    """Dwell-time distribution per chemistry class per condition."""
    df = per_pair.copy()
    df = df.dropna(subset=["mean_on_ps", "contact_class"])
    df = df[df["mean_on_ps"] > 0]
    df["log_dwell"] = np.log10(df["mean_on_ps"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    for ax, partner in zip(axes, ["prk", "gap2"]):
        sub = df[df["partner"] == partner]
        if sub.empty:
            ax.set_visible(False)
            continue
        pos = 0
        xticks, xlabels = [], []
        for cls in CLASS_ORDER:
            for cond in CONDITION_ORDER:
                vals = sub.loc[(sub["contact_class"] == cls) & (sub["condition"] == cond),
                               "log_dwell"]
                if len(vals) < 3:
                    pos += 1
                    continue
                parts = ax.violinplot(vals, positions=[pos], widths=0.7, showmeans=False,
                                      showmedians=True, showextrema=False)
                for pc in parts["bodies"]:
                    pc.set_facecolor(CLASS_COLORS[cls])
                    pc.set_alpha(0.45)
                    pc.set_edgecolor(CONDITION_COLORS[cond])
                    pc.set_linewidth(0.8)
                parts["cmedians"].set_color(CONDITION_COLORS[cond])
                pos += 1
            xticks.append(pos - 2)  # roughly center of 4-cond block
            xlabels.append(cls.replace("_", "\n"))
            pos += 1  # gap between classes
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, fontsize=9)
        ax.set_ylabel("log10 mean dwell time (ps)")
        ax.set_title(f"CP12 vs {partner.upper()}")
        ax.grid(axis="y", lw=0.3, alpha=0.4)
    legend = [Patch(facecolor=CONDITION_COLORS[c], edgecolor="k", label=c) for c in CONDITION_ORDER]
    legend.append(Patch(facecolor="white", edgecolor="none", label=""))
    legend += [Patch(facecolor=CLASS_COLORS[cl], edgecolor="k", label=cl) for cl in CLASS_ORDER]
    fig.legend(handles=legend, loc="lower center", ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    savefig(fig, outdir / "fig_class_violins")


def make_class_kinetics_strip(per_pair: pd.DataFrame, outdir: Path) -> None:
    """Median k_off + turnover stripped per class & condition (with bootstrap CI)."""
    rng = np.random.default_rng(0)
    df = per_pair.dropna(subset=["k_off_per_ns", "turnover_per_ns", "contact_class"]).copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharex=True)
    for ax, metric, title in zip(axes,
                                  ["k_off_per_ns", "turnover_per_ns"],
                                  ["k_off (per ns) — lower = stronger",
                                   "turnover rate (per ns) — higher = more sampling"]):
        x_positions = []
        x_labels = []
        for partner_idx, partner in enumerate(["prk", "gap2"]):
            sub = df[df["partner"] == partner]
            for cls_idx, cls in enumerate(CLASS_ORDER):
                for cond_idx, cond in enumerate(CONDITION_ORDER):
                    vals = sub.loc[(sub["contact_class"] == cls) & (sub["condition"] == cond),
                                   metric].to_numpy()
                    vals = vals[np.isfinite(vals)]
                    if vals.size < 3:
                        continue
                    pos = partner_idx * (len(CLASS_ORDER) * 5) + cls_idx * 5 + cond_idx
                    med = np.median(vals)
                    boot_meds = np.array([np.median(rng.choice(vals, size=vals.size))
                                          for _ in range(500)])
                    lo, hi = np.percentile(boot_meds, [2.5, 97.5])
                    ax.errorbar(pos, med, yerr=[[med - lo], [hi - med]],
                                fmt="o", color=CONDITION_COLORS[cond],
                                ecolor=CONDITION_COLORS[cond], capsize=2, ms=4, alpha=0.85)
                    x_positions.append(pos)
                    x_labels.append(f"{partner}/{cls}/{cond}")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(axis="y", lw=0.3, alpha=0.4)
    axes[0].set_ylabel("k_off (1/ns)")
    axes[1].set_ylabel("turnover (per ns)")
    # x-axis: label by class blocks
    block_centers = []
    for partner_idx, partner in enumerate(["prk", "gap2"]):
        for cls_idx, cls in enumerate(CLASS_ORDER):
            block_centers.append((partner_idx * (len(CLASS_ORDER) * 5) + cls_idx * 5 + 1.5,
                                  f"{partner}\n{cls.replace('_','-')}"))
    for ax in axes:
        ax.set_xticks([p for p, _ in block_centers])
        ax.set_xticklabels([t for _, t in block_centers], fontsize=7, rotation=45,
                           ha="right")
    legend = [Patch(facecolor=CONDITION_COLORS[c], edgecolor="k", label=c) for c in CONDITION_ORDER]
    fig.legend(handles=legend, loc="upper center", ncol=4, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    savefig(fig, outdir / "fig_class_kinetics_strip")


def make_host_redox_dwell_swarm(per_hmm: pd.DataFrame, partner: str, outdir: Path) -> None:
    """For each HMM column, scatter mean dwell across 4 conditions."""
    df = per_hmm[per_hmm["partner"] == partner].copy()
    if df.empty:
        return
    df = df.groupby(["condition", "hmm_col"], as_index=False).agg(
        mean_on_ps=("mean_on_ps_mean", "mean"),
        occ=("occupancy_mean", "mean"),
    )

    hmm_cols = sorted(df["hmm_col"].unique())
    fig, ax = plt.subplots(figsize=(15, 4.0))
    for cond in CONDITION_ORDER:
        sub = df[df["condition"] == cond].copy()
        ax.scatter(sub["hmm_col"], sub["mean_on_ps"], s=18,
                   color=CONDITION_COLORS[cond], alpha=0.85,
                   edgecolor="white", linewidth=0.4, label=cond)
    ax.set_yscale("log")
    ax.set_xlabel("CP12 HMM alignment position")
    ax.set_ylabel("Mean dwell time τ_on (ps)")
    ax.set_title(f"Per-position contact persistence — CP12 vs {partner.upper()}")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="y", lw=0.3, alpha=0.4)
    fig.tight_layout()
    savefig(fig, outdir / f"fig_host_redox_dwell_{partner}")


def make_survival_overview(per_pair: pd.DataFrame, outdir: Path) -> None:
    """Median survival fractions per condition for each partner."""
    df = per_pair.copy()
    rows = []
    for partner in ["prk", "gap2"]:
        for cond in CONDITION_ORDER:
            sub = df[(df["partner"] == partner) & (df["condition"] == cond)]
            if sub.empty:
                continue
            rows.append({
                "partner": partner, "condition": cond,
                "S1ns":  float(sub["survival_1ns"].median()),
                "S5ns":  float(sub["survival_5ns"].median()),
                "S20ns": float(sub["survival_20ns"].median()),
                "n_pairs": len(sub),
            })
    if not rows:
        return
    surv = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, partner in zip(axes, ["prk", "gap2"]):
        sub = surv[surv["partner"] == partner]
        x = [0, 1, 2]
        for _, row in sub.iterrows():
            ax.plot(x, [row["S1ns"], row["S5ns"], row["S20ns"]],
                    "-o", color=CONDITION_COLORS[row["condition"]], lw=2, ms=6,
                    label=row["condition"])
        ax.set_xticks(x)
        ax.set_xticklabels(["S(1 ns)", "S(5 ns)", "S(20 ns)"])
        ax.set_ylim(0, 1)
        ax.set_title(f"CP12 vs {partner.upper()}")
        ax.grid(lw=0.3, alpha=0.4)
    axes[0].set_ylabel("Median survival (across pairs)")
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    savefig(fig, outdir / "fig_survival_overview")


def write_text_summary(per_pair: pd.DataFrame, per_class: pd.DataFrame,
                       cond_summary: pd.DataFrame, outdir: Path) -> None:
    """Markdown summary with the numbers the manuscript will quote."""
    lines = ["# Contact dynamics summary (binary 5.0 Å heavy-atom + 3-frame min-run filter)",
             "", "## Per-condition overview", ""]
    if not cond_summary.empty:
        cols = ["condition", "partner", "n_pairs",
                "occupancy_mean", "mean_on_ps_mean", "p90_on_ps_mean",
                "turnover_per_ns_mean", "k_off_per_ns_mean",
                "survival_5ns_mean", "survival_20ns_mean"]
        keep = [c for c in cols if c in cond_summary.columns]
        # round numerics
        cs = cond_summary[keep].copy()
        for c in keep:
            if cs[c].dtype.kind == "f":
                cs[c] = cs[c].round(3)
        lines.append(cs.to_markdown(index=False))
    else:
        lines.append("(no data)")

    lines += ["", "## Class-stratified summary (median across pairs)", ""]
    if not per_class.empty:
        keep = ["condition", "partner", "contact_class", "n_pairs",
                "occupancy_median", "mean_on_ps_median", "turnover_per_ns_median",
                "k_off_per_ns_median", "survival_5ns_median"]
        keep = [c for c in keep if c in per_class.columns]
        pc = per_class[keep].copy()
        for c in keep:
            if pc[c].dtype.kind == "f":
                pc[c] = pc[c].round(3)
        lines.append(pc.to_markdown(index=False))
    else:
        lines.append("(no data)")

    # Host-by-redox interaction: dwell distribution test per partner
    lines += ["", "## Host & redox effect on mean dwell (median bootstrap CI)", ""]
    rng = np.random.default_rng(0)
    for partner in ["prk", "gap2"]:
        sub = per_pair[per_pair["partner"] == partner].dropna(subset=["mean_on_ps"])
        if sub.empty:
            continue
        lines.append(f"### {partner.upper()}")
        lines.append("| condition | n pairs | median τ_on (ps) | [95% CI] |")
        lines.append("|---|---|---|---|")
        for cond in CONDITION_ORDER:
            vals = sub.loc[sub["condition"] == cond, "mean_on_ps"].to_numpy()
            vals = vals[np.isfinite(vals)]
            if vals.size < 5:
                continue
            med = np.median(vals)
            boots = np.array([np.median(rng.choice(vals, size=vals.size)) for _ in range(2000)])
            lo, hi = np.percentile(boots, [2.5, 97.5])
            lines.append(f"| {cond} | {vals.size} | {med:.0f} | [{lo:.0f}, {hi:.0f}] |")
        lines.append("")

    (outdir / "contact_dynamics_text_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="contact_dynamics")
    ap.add_argument("--outdir", default="contact_dynamics_plots")
    args = ap.parse_args()

    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    per_pair = pd.read_csv(indir / "per_pair_stats.csv.gz")
    per_hmm = pd.read_csv(indir / "per_pair_hmm_aggregated.csv")
    per_class = pd.read_csv(indir / "per_class_summary.csv")
    cond_summary = pd.read_csv(indir / "condition_summary.csv")

    print(f"loaded {len(per_pair):,} per-pair rows, "
          f"{len(per_hmm):,} per-(hmm,class) rows", flush=True)

    for partner in ["prk", "gap2"]:
        make_persistence_heatmap(per_hmm, partner, outdir)
        make_kinetic_scatter(per_hmm, partner, outdir)
        make_host_redox_dwell_swarm(per_hmm, partner, outdir)
    make_class_violins(per_pair, outdir)
    make_class_kinetics_strip(per_pair, outdir)
    make_survival_overview(per_pair, outdir)
    write_text_summary(per_pair, per_class, cond_summary, outdir)

    print(f"wrote plots and summary to {outdir}")


if __name__ == "__main__":
    main()
