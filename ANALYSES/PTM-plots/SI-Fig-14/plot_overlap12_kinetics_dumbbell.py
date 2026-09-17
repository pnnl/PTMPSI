#!/usr/bin/env python3
"""Overlap-12 kinetics dumbbell: MED4 vs P-HM2 matched pairs.

The 12 systems in both ensembles' top-5% carry identical GAP2+PRK PTM sets in
both CP12 backgrounds, so they are matched pairs differing only in CP12 origin.
This figure shows, per system, log2(metric/REF) for MED4 (host) and P-HM2
(viral) as paired dots joined by a connector, in two panels:

  left  — mean first-passage time (MFPT)
  right — intra-cluster residence time (tau)

The host->viral shift is a left-pointing connector in the tau panel (P-HM2
accelerates) but a short connector in the MFPT panel (little change) — the
tau-specific effect, read per system.

Data: CANONICAL kinetic_summary_structural_only_90pct. REF-relative log2 ratios precomputed in
combined_kinetics_ref_relative.csv.


Outputs (kinetics_vs_pca/):
  fig_overlap12_kinetics_dumbbell.{pdf,png,tiff}
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MM = 1 / 25.4  # mm -> inch

SCALE = 180 / 150
FIG_W_MM, FIG_H_MM = 150 * SCALE, 92 * SCALE  # = 180.0 x 110.4 mm
MAX_H_MM = 210  # 2-column
assert FIG_H_MM <= MAX_H_MM

mpl.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6 * SCALE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

SCRATCH = Path("/scratch/path/to")
KIN = (SCRATCH / "/kinetic_summary_structural_only_90pct"
       / "combined_kinetics_ref_relative.csv")
OUT = SCRATCH / "/kinetics_vs_pca"
OUT.mkdir(parents=True, exist_ok=True)

OVERLAP = ["0273", "0332", "0334", "0620", "0925", "0928",
           "0929", "0935", "0956", "0997", "1005", "1006"]
ENS_COLORS = {"MED4": "#428BBE", "HM2": "#FF9333"}
ENS_DISPLAY = {"MED4": "MED4", "HM2": "P-HM2"}
LM = "log2_mean_mfpt_ns_ratio_vs_ref"
LT = "log2_mean_tau_ns_ratio_vs_ref"
PANELS = [(LM, "log$_2$(MFPT / REF)", "MFPT"),
          (LT, "log$_2$(τ / REF)", "τ")]


def save_all(fig, base):
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=200)
    fig.savefig(base.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def load():
    df = pd.read_csv(KIN, dtype={"system": str})
    out = {}
    for ens in ("MED4", "HM2"):
        d = df[df.ensemble == ens].set_index("system")
        out[ens] = d
    return out


def main():
    data = load()
    # order systems top-to-bottom by P-HM2 tau shift (most accelerated at top)
    hm2_tau = data["HM2"].reindex(OVERLAP)[LT]
    order = hm2_tau.sort_values().index.tolist()  # most negative first
    y = {s: i for i, s in enumerate(order)}

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W_MM * MM, FIG_H_MM * MM), sharey=True)
    for ax, (col, xlabel, title) in zip(axes, PANELS):
        for s in order:
            ym = data["MED4"].at[s, col] if s in data["MED4"].index else np.nan
            yh = data["HM2"].at[s, col] if s in data["HM2"].index else np.nan
            yi = y[s]
            if not (np.isnan(ym) or np.isnan(yh)):
                ax.plot([ym, yh], [yi, yi], color="#9a9a9a", lw=0.9 * SCALE, zorder=1)
            ax.plot(ym, yi, "o", color=ENS_COLORS["MED4"], ms=5 * SCALE,
                    markeredgecolor="black", markeredgewidth=0.4 * SCALE, zorder=3)
            ax.plot(yh, yi, "o", color=ENS_COLORS["HM2"], ms=5 * SCALE,
                    markeredgecolor="black", markeredgewidth=0.4 * SCALE, zorder=3)
        ax.axvline(0, color="black", lw=0.6 * SCALE, ls="--", zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
    axes[0].set_yticks(list(y.values()))
    axes[0].set_yticklabels(order, fontsize=6)
    axes[0].set_ylabel("System")
    axes[0].set_ylim(-0.6, len(order) - 0.4)

    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=ENS_COLORS[e], markeredgecolor="black",
                      markeredgewidth=0.4 * SCALE, markersize=5 * SCALE, label=ENS_DISPLAY[e])
               for e in ("MED4", "HM2")]
    fig.legend(handles=handles, frameon=False, loc="lower center",
               ncol=2, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    save_all(fig, OUT / "fig_overlap12_kinetics_dumbbell")
    print(f"Saved to {OUT}")


if __name__ == "__main__":
    main()
