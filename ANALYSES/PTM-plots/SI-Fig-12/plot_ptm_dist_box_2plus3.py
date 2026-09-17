#!/usr/bin/env python3
"""Combined figure: total-PTM-count box (a, top, full width) over the
per-subunit 3-panel box (b, bottom, CP12 | PRK | GAP2). Reuses load +
drawing helpers from plot_ptm_chaincount_vs_distance.py so the panels match
the standalone figures fig_ptm_chaincount_vs_distance_box_combined and
..._box_combined_3panel.

"""
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
sys.path.insert(0, "/scratch/path/to/PRK_SSG_RERUN/ensemble")
import plot_ptm_chaincount_vs_distance as P  # noqa: E402

MM = 1 / 25.4
WIDTH_MM = 180.0   # 2-column fixed width
MAX_H_MM = 185.0   # 2-column cap for a <300-word caption (176 words)
TOP_H_MM = 62.0    # panel (a): total-PTM-count box, full width
BOT_H_MM = 88.0    # panel (b): 3-panel subunit box (CP12 | PRK | GAP2)
FIG_H_MM = TOP_H_MM + BOT_H_MM + 6.0  # + inter-row gap
assert FIG_H_MM <= MAX_H_MM, f"{FIG_H_MM} exceeds {MAX_H_MM} mm cap"

# P's own import sets 11/12pt fonts and savefig.bbox='tight'; override for
# max 7pt sans-serif; exact canvas size, no auto-crop.
mpl.rcParams.update({
    "font.size": 6.5, "axes.labelsize": 7, "axes.titlesize": 7,
    "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "savefig.bbox": "standard",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    # _draw_combined_box doesn't pass box/whisker/cap linewidths explicitly
    # (only medianprops), so it inherits these rcParams -- thin them down
    # from the default ~1.0pt, which reads as too heavy at this panel size.
    "boxplot.boxprops.linewidth": 0.6,
    "boxplot.whiskerprops.linewidth": 0.6,
    "boxplot.capprops.linewidth": 0.6,
})


def _add_corr_bottom_right(ax, md, hd, col):
    """Place compact per-ensemble Spearman text in the lower-right corner."""
    data = {"MED4": md, "HM2": hd}
    for k, ens in enumerate(("MED4", "HM2")):
        sub = data[ens].dropna(subset=[col, P.DIST_COL])
        rho, p = P.stats.spearmanr(sub[col], sub[P.DIST_COL])
        ax.text(0.985, 0.012 + (1 - k) * 0.038,
                f"{P.ENS_DISPLAY[ens]}: Spearman ρ = {rho:+.3f}, p = {p:.2g}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8,
                          pad=0.6))


def _fix_box_alpha(ax):
    """_draw_combined_box sets a single patch.set_alpha(0.7), which also
    makes the (already too-thick) box edge translucent -- overlapping
    translucent edges read as fuzzy/blurry. Bake the alpha into the face
    color only and keep the edge solid opaque black."""
    for patch in ax.patches:
        if patch.get_alpha() is not None:
            fc = patch.get_facecolor()
            patch.set_facecolor((fc[0], fc[1], fc[2], 0.7))
            patch.set_edgecolor((0, 0, 0, 1))
            patch.set_alpha(None)

    for line in ax.lines:
        x = line.get_xdata()
        y = line.get_ydata()
        if len(x) == 2 and len(y) == 2 and y[0] == y[1]:
            line.set_linewidth(0.65)


def _shrink_texts(ax, size=6):
    """P's helper functions hardcode fontsizes; clamp
    after the fact rather than editing the shared module."""
    for t in ax.texts:
        t.set_fontsize(size)
    plt.setp(ax.get_xticklabels(), fontsize=6)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    leg = ax.get_legend()
    if leg is not None:
        for t in leg.get_texts():
            t.set_fontsize(6)


def _set_zero_baseline(*axes, top_pad=0.0):
    """Use a true zero baseline so 0-1 has the same scale as other y ticks."""
    ymax = max(ax.get_ylim()[1] for ax in axes) + top_pad
    yticks = range(0, int(ymax) + 1)
    for ax in axes:
        ax.set_ylim(0, ymax)
        ax.set_yticks(yticks)


def main():
    md, hd = P.load("MED4"), P.load("HM2")
    handles = [Patch(facecolor=P.ENS_COLORS[e], alpha=0.7, label=P.ENS_DISPLAY[e])
               for e in ("MED4", "HM2")]

    fig = plt.figure(figsize=(WIDTH_MM * MM, FIG_H_MM * MM))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[TOP_H_MM, BOT_H_MM],
                 hspace=0.34, wspace=0.10, top=0.95, bottom=0.09,
                 left=0.09, right=0.98)

    # --- (a) top: total PTM count, full width ---
    axT = fig.add_subplot(gs[0, :])
    vals = P.count_vals(md, hd, "count_total")
    P._draw_combined_box(axT, md, hd, "count_total", vals)
    _fix_box_alpha(axT)
    axT.set_xlabel("Total PTMs"); axT.set_ylabel("PCA distance")
    axT.legend(handles=handles, frameon=False, loc="upper left")
    P._add_corr_text(axT, md, hd, "count_total")
    _set_zero_baseline(axT)
    _shrink_texts(axT)

    # --- (b) bottom: per-subunit 3-panel (CP12 | PRK | GAP2), shared y ---
    axes = []
    for k, (col, xlabel, title) in enumerate(P.PANELS):
        ax = fig.add_subplot(gs[1, k], sharey=axes[0] if axes else None)
        P._draw_combined_box(ax, md, hd, col, P.count_vals(md, hd, col))
        _fix_box_alpha(ax)
        ax.set_xlabel(xlabel); ax.set_title(title)
        if k == 0:
            # color legend already shown once in panel (a); avoid overlap here
            ax.set_ylabel("PCA distance")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        axes.append(ax)

    for ax, (col, xlabel, title) in zip(axes, P.PANELS):
        _add_corr_bottom_right(ax, md, hd, col)
        _shrink_texts(ax)
    _set_zero_baseline(*axes, top_pad=0.8)

    # color legend for panel (b): keep it inside CP12, using the shared
    # y-axis headroom above the tallest CP12 whisker as a buffer.
    axes[0].legend(handles=handles, frameon=False, loc="upper left")
    _shrink_texts(axes[0])

    fig.text(0.01, 0.965, "a", fontsize=8, fontweight="bold")
    fig.text(0.01, 0.585, "b", fontsize=8, fontweight="bold")

    base = P.OUT / "fig_ptm_chaincount_vs_distance_box_combined_2plus3"
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=450)
    fig.savefig(base.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"done: {base.name} ({WIDTH_MM:.1f} x {FIG_H_MM:.1f} mm)")


if __name__ == "__main__":
    main()
