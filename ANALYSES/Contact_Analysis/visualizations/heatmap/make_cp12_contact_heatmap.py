#!/usr/bin/env python3
"""vertical CP12 contact-frequency heatmap.

Layout:
  Two side-by-side super-panels (left = reduced/system 0000, right = oxidized/system 0003).
  Each super-panel has 6 narrow columns:
    1. MED4 contact-frequency heatmap (receptor = PRK+GAP2)
    2. HM2  contact-frequency heatmap (receptor)
    3. MED4 partner-classification annotation strip
    4. HM2  partner-classification annotation strip
    5. MED4 # contacts bar (PRK in blue, GAP2 in green, overlaid)
    6. HM2  # contacts bar (PRK in blue, GAP2 in green, overlaid)

y-axis (shared across both super-panels): HMM cols 1-70. Residue labels show paired
MED4-residue and HM2-residue numbers, with '-' for gaps.

Inputs:
  interface_results/redox_cutoff_5p0A/cp12_contact_freq_chain_and_rep_averaged.csv
  cp12_contact_pairs_5p0A/differential_pairs_{prk,gap2}.csv
  CP12_MED4_vs_HM2_alignment_annotation.csv
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch

BASE = Path(".")
# Single master CSV with every value needed for this figure (built from the underlying
# per-residue contact-freq + pairwise-contact CSVs; see FINAL_FIGURES/README).
MASTER_CSV = BASE / "FINAL_FIGURES/cp12_contact_freq_heatmap_data_master.csv"
OUT        = BASE / "FINAL_FIGURES"

# Heatmap and annotation palettes
FREQ_CMAP = LinearSegmentedColormap.from_list(
    "freq", ["#ffffff", "#fff5b1", "#fdae61", "#d7301f", "#7f0000"]
)
ANNOT_COLORS = {
    "PRK":     "violet",   # blue
    "GAP2":    "lightsteelblue",   # green
    "Both":    "#ff7f0e",   # orange
    "Neither": "#e5e5e5",   # light gray
    "Gap":     "#7f7f7f",   # mid gray
}
# Bar palettes deliberately distinct from the interface annotation blue/green.
# Pick one via CLI --palette ; the chosen tag is encoded into the output filename so
# multiple palette renders don't overwrite each other.
#
# Each palette dict has two required keys (PRK / GAP2 bar colors) and two optional
# keys (PRK / GAP2 interface annotation colors). If the interface colors are omitted,
# the figure uses the canonical blue/green annotation strip; if supplied, they override
# the global ANNOT_COLORS — useful for palettes where the interface strip and the bars
# are meant to share a hue family (e.g., "lavender_pink": light shade on the strip,
# dark shade on the bars).
PALETTES = {
    "navy_forest":    {"prk_bar": "#08306b", "gap2_bar": "#00441b"},
    "purple_orange":  {"prk_bar": "#762a83", "gap2_bar": "#e08214"},
    "teal_magenta":   {"prk_bar": "#01665e", "gap2_bar": "#c51b7d"},
    "charcoal_gold":  {"prk_bar": "#525252", "gap2_bar": "#bf812d"},
    "blue_green":     {"prk_bar": "#1f77b4", "gap2_bar": "#2ca02c"},
    "lavender_pink": {
        "prk_interface":  "#E691E8",  # light pink — PRK interface annotation
        "gap2_interface": "#D3D3EC",  # light lavender — GAP2 interface annotation
        "prk_bar":        "#8E1F90",  # deep magenta — PRK # contacts
        "gap2_bar":       "lightsteelblue",  # deep indigo — GAP2 # contacts
    },
}
DEFAULT_PALETTE = "lavender_pink"

N_HMM = 70
PAIR_THRESH = 0.10       # for partner count
INTERFACE_THRESH = 0.50  # for "PRK / GAP2 / Both" classification

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})


def load_master():
    """Load the single master CSV and pivot per (host, redox) into N_HMM-length arrays."""
    df = pd.read_csv(MASTER_CSV)
    df.aln_col_1based = df.aln_col_1based.astype(int)
    return df


def column_arrays(df, host, redox):
    """Pull the four arrays needed for one super-panel column from the master CSV."""
    sub = df[(df.host == host) & (df.redox == redox)].set_index("aln_col_1based").sort_index()
    rec_freq   = np.array([sub.loc[c, "receptor_freq_chain_rep_avg"] if c in sub.index else np.nan
                           for c in range(1, N_HMM + 1)], dtype=float)
    prk_freq   = np.array([sub.loc[c, "prk_freq_chain_rep_avg"] if c in sub.index else np.nan
                           for c in range(1, N_HMM + 1)], dtype=float)
    gap_freq   = np.array([sub.loc[c, "gap2_freq_chain_rep_avg"] if c in sub.index else np.nan
                           for c in range(1, N_HMM + 1)], dtype=float)
    iface_cls  = np.array([sub.loc[c, "interface_class"] if c in sub.index else "Gap"
                           for c in range(1, N_HMM + 1)], dtype=object)
    prk_cnt    = np.array([int(sub.loc[c, "prk_partner_count"]) if c in sub.index else 0
                           for c in range(1, N_HMM + 1)], dtype=int)
    gap_cnt    = np.array([int(sub.loc[c, "gap2_partner_count"]) if c in sub.index else 0
                           for c in range(1, N_HMM + 1)], dtype=int)
    # If a residue is absent, the master CSV stores "" → pandas reads NaN
    # but the receptor_freq_chain_rep_avg also reads NaN for those rows. NaN cells will
    # display as gap-gray in draw_freq_col.
    rec_freq[iface_cls == "Gap"] = np.nan
    return rec_freq, prk_freq, gap_freq, iface_cls, prk_cnt, gap_cnt


def yaxis_labels_from_master(df):
    """y-axis labels = 'HMMcol  MED4aa/HM2aa', e.g., '57 C/C'. Reads from master CSV."""
    # The aa info is in the master too; use any single host's rows to read it
    sub = df[(df.host == "med4") & (df.redox == "ox")].set_index("aln_col_1based").sort_index()
    return [f"{c:>2}  {sub.loc[c, 'med4_aa']}/{sub.loc[c, 'hm2_aa']}"
            for c in range(1, N_HMM + 1)]


def yaxis_labels(med4_at, hm2_at):
    """y-axis labels = 'HMMcol  MED4aa / HM2aa', e.g., ' 57  C / C', ' 19  S / -'.
    Uses HMM column 1-70 numbering with the single-letter MED4/HM2 amino acids appended."""
    labels = []
    for c in range(1, N_HMM + 1):
        m_aa = med4_at.get(c, (None, None))[1] or "-"
        h_aa = hm2_at.get(c, (None, None))[1]  or "-"
        labels.append(f"{c:>2}  {m_aa}/{h_aa}")
    return labels


def draw_freq_col(ax, arr, title):
    """Draw a single-column heatmap (height = N_HMM)."""
    ax.set_facecolor(ANNOT_COLORS["Gap"])  # NaN cells show through as gap-gray
    im = ax.imshow(arr.reshape(-1, 1), aspect="auto", cmap=FREQ_CMAP, vmin=0.0, vmax=1.0,
                   interpolation="nearest", extent=[0, 1, N_HMM + 0.5, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=3)
    # Full border on all 4 sides
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.7)
        s.set_color("black")
    return im


def draw_annot_col(ax, labels, title):
    """Annotation strip with per-row partner classification color."""
    code = {"PRK": 0, "GAP2": 1, "Both": 2, "Neither": 3, "Gap": 4}
    cmap = ListedColormap([ANNOT_COLORS[k] for k in ["PRK", "GAP2", "Both", "Neither", "Gap"]])
    img = np.array([code[x] for x in labels]).reshape(-1, 1)
    ax.imshow(img, aspect="auto", cmap=cmap, vmin=0, vmax=4, interpolation="nearest",
              extent=[0, 1, N_HMM + 0.5, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=3)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.7)
        s.set_color("black")


def draw_bar_col(ax, prk_counts, gap_counts, title, xmax, side="left"):
    """Horizontal mirrored bars: PRK (blue, right), GAP2 (green, left) at each row.

    `side` ∈ {"left", "right"} controls which inner tick label is dropped so that, when two
    bar columns are placed adjacent to each other, only ONE "<xmax>" label sits at the
    inner boundary instead of two overlapping ones.
        side="left"  -> bar is on the LEFT of the pair -> drop the right-edge label
        side="right" -> bar is on the RIGHT of the pair -> drop the left-edge label
    """
    y = np.arange(1, N_HMM + 1)
    ax.barh(y, prk_counts, color=PRK_BAR, edgecolor="none", height=0.85, alpha=0.85, label="PRK")
    ax.barh(y, -gap_counts, color=GAP_BAR, edgecolor="none", height=0.85, alpha=0.85, label="GAP2")
    ax.axvline(0, color="black", lw=0.4)
    ax.set_yticks([])
    ax.set_ylim(N_HMM + 0.5, 0.5)
    ax.set_xlim(-xmax, xmax)
    if side == "left":
        ax.set_xticks([-xmax, 0])
        ax.set_xticklabels([str(xmax), "0"], fontsize=6.5)
    else:  # right
        ax.set_xticks([0, xmax])
        ax.set_xticklabels(["0", str(xmax)], fontsize=6.5)
    ax.tick_params(axis="x", length=2, pad=1)
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=3)
    # Full border on all 4 sides
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.7)
        s.set_color("black")


def super_panel(fig, gs, df, redox_label,
                show_yticks=False, xmax_bars=30):
    """Build one super-panel (one redox state). Returns the freq imshow handle."""
    # MED4 and HM2 arrays for this redox state, sourced from the single master CSV
    med4_rec, med4_prk, med4_gap, med4_class, med4_prk_cnt, med4_gap_cnt = \
        column_arrays(df, "med4", redox_label)
    hm2_rec,  hm2_prk,  hm2_gap,  hm2_class,  hm2_prk_cnt,  hm2_gap_cnt  = \
        column_arrays(df, "hm2",  redox_label)

    # Six columns
    inner = gs.subgridspec(1, 6, width_ratios=[1.1, 1.1, 0.5, 0.5, 1.6, 1.6], wspace=0.05)
    ax_m_freq  = fig.add_subplot(inner[0, 0])
    ax_h_freq  = fig.add_subplot(inner[0, 1])
    ax_m_anno  = fig.add_subplot(inner[0, 2])
    ax_h_anno  = fig.add_subplot(inner[0, 3])
    ax_m_bar   = fig.add_subplot(inner[0, 4])
    ax_h_bar   = fig.add_subplot(inner[0, 5])

    im = draw_freq_col(ax_m_freq, med4_rec, "MED4")
    draw_freq_col(ax_h_freq, hm2_rec, "PHM2")
    draw_annot_col(ax_m_anno, med4_class, "MED4")
    draw_annot_col(ax_h_anno, hm2_class,  "PHM2")
    draw_bar_col(ax_m_bar, med4_prk_cnt, med4_gap_cnt, "MED4 # contacts", xmax_bars, side="left")
    draw_bar_col(ax_h_bar, hm2_prk_cnt,  hm2_gap_cnt,  "HM2 # contacts",  xmax_bars, side="right")

    # y-tick labels go on the leftmost axis of the leftmost super-panel
    if show_yticks:
        labels = yaxis_labels_from_master(df)
        ax_m_freq.set_yticks(np.arange(1, N_HMM + 1))
        ax_m_freq.set_yticklabels(labels, fontsize=8, family="monospace")
        ax_m_freq.tick_params(axis="y", length=0, pad=2)

    # Cysteine markers (HMM cols 57 and 66) — dotted line across each freq/annot column
    for ax in (ax_m_freq, ax_h_freq, ax_m_anno, ax_h_anno):
        for c in (57, 66):
            ax.axhline(c, color="black", lw=0.4, ls=":", alpha=0.7)

    # Super-panel title (overall)
    # We'll add via fig.text after layout

    return im


def main():
    df = load_master()

    fig = plt.figure(figsize=(13.5, 14.5))
    outer = fig.add_gridspec(1, 2, wspace=0.10, left=0.07, right=0.92, top=0.91, bottom=0.07)
    # Oxidized (system 0003) on LEFT, Reduced (system 0000) on RIGHT
    im_l = super_panel(fig, outer[0, 0], df, redox_label="ox",
                       show_yticks=True,  xmax_bars=30)
    im_r = super_panel(fig, outer[0, 1], df, redox_label="red",
                       show_yticks=False, xmax_bars=30)

    # Super-panel banner titles — sit above the per-column "MED4 / HM2 / M / H / ..." titles
    fig.text(0.245, 0.955, "Oxidized", fontsize=16, fontweight="bold",
             ha="center", va="bottom")
    fig.text(0.745, 0.955, "Reduced", fontsize=16, fontweight="bold",
             ha="center", va="bottom")

    # Shared colorbar (frequency) — title cleaned up
    cax = fig.add_axes([0.94, 0.30, 0.011, 0.40])
    cb = fig.colorbar(im_l, cax=cax)
    cb.set_label("Avg. contact frequency", fontsize=14)
    cb.ax.tick_params(labelsize=7)

    # y-axis label
    fig.text(0.012, 0.49, "CP12 HMM alignment position",
             rotation=90, va="center", ha="center", fontsize=15, fontweight="bold")

    # Legend (annotation strip + bar colors + cysteine dashed-line key)
    # No "Both" since the data never reaches both partners simultaneously at the threshold;
    # "(annot)" suffix dropped.
    from matplotlib.lines import Line2D
    handles = [
        Patch(facecolor=ANNOT_COLORS["PRK"],     edgecolor="none",  label="PRK interface"),
        Patch(facecolor=ANNOT_COLORS["GAP2"],    edgecolor="none",  label="GAP2 interface"),
        Patch(facecolor=ANNOT_COLORS["Neither"], edgecolor="black", linewidth=0.5, label="Non-interface"),
        Patch(facecolor=ANNOT_COLORS["Gap"],     edgecolor="none",  label="Alignment gap"),
        Patch(facecolor=PRK_BAR, edgecolor="none", label="PRK contact count"),
        Patch(facecolor=GAP_BAR, edgecolor="none", label="GAP2 contact count"),
        Line2D([0], [0], color="black", lw=1.0, ls=":",
               label="Regulatory cysteine"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=4, fontsize=10, frameon=True, columnspacing=1.5, handlelength=1.8,
               handleheight=1.2, borderpad=0.6)

    OUT.mkdir(exist_ok=True)
    base = f"cp12_contact_freq_heatmap_{BAR_COLOR_TAG}"
    for ext in ("pdf", "png", "tiff"):
        path = OUT / f"{base}.{ext}"
        fig.savefig(path, dpi=600 if ext == "tiff" else 300,
                    format=ext if ext == "tiff" else None,
                    bbox_inches="tight")
        print(f"Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--palette", choices=list(PALETTES.keys()) + ["all"],
                    default=DEFAULT_PALETTE,
                    help="Bar-color palette. Use 'all' to render every variant in one go.")
    args = ap.parse_args()
    palettes_to_run = list(PALETTES.keys()) if args.palette == "all" else [args.palette]
    # Snapshot the canonical interface colors so we can restore between palette runs
    _DEFAULT_PRK_IFACE  = ANNOT_COLORS["PRK"]
    _DEFAULT_GAP2_IFACE = ANNOT_COLORS["GAP2"]
    for tag in palettes_to_run:
        spec = PALETTES[tag]
        PRK_BAR = spec["prk_bar"]
        GAP_BAR = spec["gap2_bar"]
        BAR_COLOR_TAG = tag
        # Apply / restore interface-color overrides if the palette specifies them
        ANNOT_COLORS["PRK"]  = spec.get("prk_interface",  _DEFAULT_PRK_IFACE)
        ANNOT_COLORS["GAP2"] = spec.get("gap2_interface", _DEFAULT_GAP2_IFACE)
        print(f"\n>>> Rendering palette '{tag}'  (bar PRK={PRK_BAR}, bar GAP2={GAP_BAR}, "
              f"iface PRK={ANNOT_COLORS['PRK']}, iface GAP2={ANNOT_COLORS['GAP2']})")
        main()
