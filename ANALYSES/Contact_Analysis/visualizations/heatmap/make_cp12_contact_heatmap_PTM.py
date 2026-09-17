#!/usr/bin/env python3
"""CP12 contact-frequency heatmap — full PTM condition range.

Layout:
  Grid of N_COLS super-panels per row (default 4), one per numbered PTM condition.
  Each super-panel has 6 narrow columns:
    1. MED4 contact-frequency heatmap (receptor = PRK+GAP2)
    2. HM2  contact-frequency heatmap (receptor)
    3. MED4 partner-classification annotation strip
    4. HM2  partner-classification annotation strip
    5. MED4 # contacts bar (PRK in blue, GAP2 in green, overlaid)
    6. HM2  # contacts bar

y-axis (shared across the leftmost panel of each row): HMM cols 1-70.

Inputs:
  FINAL_FIGURES/cp12_contact_freq_heatmap_data_master_PTM.csv
"""
from __future__ import annotations
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

BASE       = Path(".")
MASTER_CSV = BASE / "FINAL_FIGURES/cp12_contact_freq_heatmap_data_master_PTM.csv"
OUT        = BASE / "FINAL_FIGURES"

# Heatmap and annotation palettes
FREQ_CMAP = LinearSegmentedColormap.from_list(
    "freq", ["#ffffff", "#fff5b1", "#fdae61", "#d7301f", "#7f0000"]
)
ANNOT_COLORS = {
    "PRK":     "violet",
    "GAP2":    "lightsteelblue",
    "Both":    "#ff7f0e",
    "Neither": "#e5e5e5",
    "Gap":     "#7f7f7f",
}

PALETTES = {
    "navy_forest":   {"prk_bar": "#08306b", "gap2_bar": "#00441b"},
    "purple_orange": {"prk_bar": "#762a83", "gap2_bar": "#e08214"},
    "teal_magenta":  {"prk_bar": "#01665e", "gap2_bar": "#c51b7d"},
    "charcoal_gold": {"prk_bar": "#525252", "gap2_bar": "#bf812d"},
    "blue_green":    {"prk_bar": "#1f77b4", "gap2_bar": "#2ca02c"},
    "lavender_pink": {
        "prk_interface":  "#E691E8",
        "gap2_interface": "#D3D3EC",
        "prk_bar":        "#8E1F90",
        "gap2_bar":       "lightsteelblue",
    },
}
DEFAULT_PALETTE = "lavender_pink"

# Layout
N_COLS    = 4         # super-panels per row; rows are computed automatically
N_HMM     = 70
PAIR_THRESH      = 0.10
INTERFACE_THRESH = 0.50

mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype":     42,
    "axes.linewidth":   0.6,
})


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def load_master() -> pd.DataFrame:
    """Load the PTM master CSV."""
    df = pd.read_csv(MASTER_CSV)
    df["aln_col_1based"] = df["aln_col_1based"].astype(int)
    df["condition"]      = df["condition"].astype(int)
    return df


def column_arrays(df: pd.DataFrame, host: str, condition: int):
    """Pull the six arrays needed for one super-panel column from the master CSV."""
    sub = (df[(df.host == host) & (df.condition == condition)]
           .set_index("aln_col_1based")
           .sort_index())
    cols = range(1, N_HMM + 1)
    rec_freq  = np.array([sub.loc[c, "receptor_freq_chain_rep_avg"] if c in sub.index else np.nan
                          for c in cols], dtype=float)
    prk_freq  = np.array([sub.loc[c, "prk_freq_chain_rep_avg"]      if c in sub.index else np.nan
                          for c in cols], dtype=float)
    gap_freq  = np.array([sub.loc[c, "gap2_freq_chain_rep_avg"]     if c in sub.index else np.nan
                          for c in cols], dtype=float)
    iface_cls = np.array([sub.loc[c, "interface_class"]             if c in sub.index else "Gap"
                          for c in cols], dtype=object)
    prk_cnt   = np.array([int(sub.loc[c, "prk_partner_count"])      if c in sub.index else 0
                          for c in cols], dtype=int)
    gap_cnt   = np.array([int(sub.loc[c, "gap2_partner_count"])     if c in sub.index else 0
                          for c in cols], dtype=int)
    rec_freq[iface_cls == "Gap"] = np.nan
    return rec_freq, prk_freq, gap_freq, iface_cls, prk_cnt, gap_cnt


def yaxis_labels_from_master(df: pd.DataFrame) -> list[str]:
    """y-axis labels: 'HMMcol  MED4aa/HM2aa'  (uses the first available condition)."""
    first_cond = sorted(df["condition"].unique())[0]
    sub = (df[(df.host == "med4") & (df.condition == first_cond)]
           .set_index("aln_col_1based")
           .sort_index())
    return [f"{c:>2}  {sub.loc[c, 'med4_aa']}/{sub.loc[c, 'hm2_aa']}"
            for c in range(1, N_HMM + 1)]


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def draw_freq_col(ax, arr, title: str):
    ax.set_facecolor(ANNOT_COLORS["Gap"])
    im = ax.imshow(arr.reshape(-1, 1), aspect="auto", cmap=FREQ_CMAP, vmin=0.0, vmax=1.0,
                   interpolation="nearest", extent=[0, 1, N_HMM + 0.5, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=7, fontweight="bold", pad=2)
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.7); s.set_color("black")
    return im


def draw_annot_col(ax, labels, title: str):
    code = {"PRK": 0, "GAP2": 1, "Both": 2, "Neither": 3, "Gap": 4}
    cmap = ListedColormap([ANNOT_COLORS[k] for k in ["PRK", "GAP2", "Both", "Neither", "Gap"]])
    img  = np.array([code[x] for x in labels]).reshape(-1, 1)
    ax.imshow(img, aspect="auto", cmap=cmap, vmin=0, vmax=4,
              interpolation="nearest", extent=[0, 1, N_HMM + 0.5, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=7, fontweight="bold", pad=2)
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.7); s.set_color("black")


def draw_bar_col(ax, prk_counts, gap_counts, title: str, xmax: int, side: str = "left"):
    """Mirrored horizontal bars: PRK right, GAP2 left."""
    y = np.arange(1, N_HMM + 1)
    ax.barh(y,  prk_counts, color=PRK_BAR, edgecolor="none", height=0.85, alpha=0.85)
    ax.barh(y, -gap_counts, color=GAP_BAR, edgecolor="none", height=0.85, alpha=0.85)
    ax.axvline(0, color="black", lw=0.4)
    ax.set_yticks([])
    ax.set_ylim(N_HMM + 0.5, 0.5)
    ax.set_xlim(-xmax, xmax)
    if side == "left":
        ax.set_xticks([-xmax, 0])
        ax.set_xticklabels([str(xmax), "0"], fontsize=5.5)
    else:
        ax.set_xticks([0, xmax])
        ax.set_xticklabels(["0", str(xmax)], fontsize=5.5)
    ax.tick_params(axis="x", length=2, pad=1)
    ax.set_title(title, fontsize=7, fontweight="bold", pad=2)
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.7); s.set_color("black")


def super_panel(fig, gs, df: pd.DataFrame, condition: int,
                show_yticks: bool = False, xmax_bars: int = 30):
    """Build one super-panel for a single PTM condition. Returns the freq imshow handle."""
    med4_rec, med4_prk, med4_gap, med4_class, med4_prk_cnt, med4_gap_cnt = \
        column_arrays(df, "med4", condition)
    hm2_rec,  hm2_prk,  hm2_gap,  hm2_class,  hm2_prk_cnt,  hm2_gap_cnt  = \
        column_arrays(df, "hm2",  condition)

    inner = gs.subgridspec(1, 6, width_ratios=[1.1, 1.1, 0.5, 0.5, 1.6, 1.6], wspace=0.05)
    ax_m_freq = fig.add_subplot(inner[0, 0])
    ax_h_freq = fig.add_subplot(inner[0, 1])
    ax_m_anno = fig.add_subplot(inner[0, 2])
    ax_h_anno = fig.add_subplot(inner[0, 3])
    ax_m_bar  = fig.add_subplot(inner[0, 4])
    ax_h_bar  = fig.add_subplot(inner[0, 5])

    im = draw_freq_col(ax_m_freq, med4_rec, "MED4")
    draw_freq_col(ax_h_freq, hm2_rec,  "HM2")
    draw_annot_col(ax_m_anno, med4_class, "M")
    draw_annot_col(ax_h_anno, hm2_class,  "H")
    draw_bar_col(ax_m_bar, med4_prk_cnt, med4_gap_cnt, "MED4 #", xmax_bars, side="left")
    draw_bar_col(ax_h_bar,  hm2_prk_cnt,  hm2_gap_cnt, "HM2 #",  xmax_bars, side="right")

    if show_yticks:
        labels = yaxis_labels_from_master(df)
        ax_m_freq.set_yticks(np.arange(1, N_HMM + 1))
        ax_m_freq.set_yticklabels(labels, fontsize=7, family="monospace")
        ax_m_freq.tick_params(axis="y", length=0, pad=2)

    # Cysteine markers
    for ax in (ax_m_freq, ax_h_freq, ax_m_anno, ax_h_anno):
        for c in (57, 66):
            ax.axhline(c, color="black", lw=0.4, ls=":", alpha=0.7)

    return im


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    df         = load_master()
    conditions = sorted(df["condition"].unique().tolist())
    n_conds    = len(conditions)
    n_cols     = min(N_COLS, n_conds)
    n_rows     = math.ceil(n_conds / n_cols)

    # Figure size: each super-panel is roughly 3 inches wide × 14.5 inches tall
    fig_w = 0.95 + n_cols * 3.2   # left margin + cols
    fig_h = n_rows * 14.5 / 2     # scale height with rows

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Outer grid: n_rows × n_cols, each cell holds one super-panel
    outer = fig.add_gridspec(
        n_rows, n_cols,
        hspace=0.18, wspace=0.14,
        left=0.07, right=0.93,
        top=0.93, bottom=0.06,
    )

    first_im = None
    for idx, cond in enumerate(conditions):
        row = idx // n_cols
        col = idx % n_cols
        gs  = outer[row, col]

        # Show y-axis tick labels only on the leftmost column of each row
        show_y = (col == 0)

        im = super_panel(fig, gs, df, condition=cond,
                         show_yticks=show_y, xmax_bars=30)
        if first_im is None:
            first_im = im

        # Condition label above each super-panel
        # Compute figure-space centre of this grid cell
        ss   = gs.get_position(fig)  # SubplotSpec → bbox in figure coords
        cx   = (ss.x0 + ss.x1) / 2
        top  = ss.y1 + 0.005
        fig.text(cx, top, f"Condition {cond}", ha="center", va="bottom",
                 fontsize=10, fontweight="bold")

    # Shared colorbar on the far right
    cax = fig.add_axes([0.945, 0.30, 0.010, 0.40])
    cb  = fig.colorbar(first_im, cax=cax)
    cb.set_label("Avg. contact frequency", fontsize=11)
    cb.ax.tick_params(labelsize=7)

    # y-axis label
    fig.text(0.012, 0.50, "CP12 HMM alignment position",
             rotation=90, va="center", ha="center", fontsize=13, fontweight="bold")

    # Figure title
    fig.suptitle("CP12 contact frequency — PTM conditions", fontsize=14,
                 fontweight="bold", y=0.975)

    # Legend
    handles = [
        Patch(facecolor=ANNOT_COLORS["PRK"],     edgecolor="none",  label="PRK interface"),
        Patch(facecolor=ANNOT_COLORS["GAP2"],    edgecolor="none",  label="GAP2 interface"),
        Patch(facecolor=ANNOT_COLORS["Both"],    edgecolor="none",  label="PRK & GAP2 interface"),
        Patch(facecolor=ANNOT_COLORS["Neither"], edgecolor="black", linewidth=0.5,
              label="Non-interface"),
        Patch(facecolor=ANNOT_COLORS["Gap"],     edgecolor="none",  label="Alignment gap"),
        Patch(facecolor=PRK_BAR, edgecolor="none", label="PRK contact count"),
        Patch(facecolor=GAP_BAR, edgecolor="none", label="GAP2 contact count"),
        Line2D([0], [0], color="black", lw=1.0, ls=":", label="Regulatory cysteine"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.002),
               ncol=4, fontsize=9, frameon=True, columnspacing=1.2, handlelength=1.6,
               handleheight=1.1, borderpad=0.5)

    OUT.mkdir(exist_ok=True)
    base = f"cp12_contact_freq_heatmap_PTM_{BAR_COLOR_TAG}"
    for ext in ("pdf", "png", "tiff"):
        path = OUT / f"{base}.{ext}"
        fig.savefig(path,
                    dpi=600 if ext == "tiff" else 300,
                    format=ext if ext == "tiff" else None,
                    bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="CP12 contact-freq heatmap for all numbered PTM conditions."
    )
    ap.add_argument("--palette", choices=list(PALETTES.keys()) + ["all"],
                    default=DEFAULT_PALETTE,
                    help="Bar-color palette. Use 'all' to render every variant.")
    args = ap.parse_args()

    palettes_to_run = list(PALETTES.keys()) if args.palette == "all" else [args.palette]
    _DEFAULT_PRK_IFACE  = ANNOT_COLORS["PRK"]
    _DEFAULT_GAP2_IFACE = ANNOT_COLORS["GAP2"]

    for tag in palettes_to_run:
        spec = PALETTES[tag]
        PRK_BAR       = spec["prk_bar"]
        GAP_BAR       = spec["gap2_bar"]
        BAR_COLOR_TAG = tag
        ANNOT_COLORS["PRK"]  = spec.get("prk_interface",  _DEFAULT_PRK_IFACE)
        ANNOT_COLORS["GAP2"] = spec.get("gap2_interface", _DEFAULT_GAP2_IFACE)
        print(f"\n>>> Rendering palette '{tag}'  "
              f"(bar PRK={PRK_BAR}, bar GAP2={GAP_BAR}, "
              f"iface PRK={ANNOT_COLORS['PRK']}, iface GAP2={ANNOT_COLORS['GAP2']})")
        main()
