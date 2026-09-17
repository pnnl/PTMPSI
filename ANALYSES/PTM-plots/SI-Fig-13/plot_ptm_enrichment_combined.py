#!/usr/bin/env python3
"""Combined PTM enrichment figure: fold enrichment by type (top) over the
per-site enrichment dotplot (bottom).

(a) fig_ptmtype_fold_enrichment content: log2(mean count top-5% / mean count
    background) per PTM type (SSG/SOH/SNO), faceted MED4 | P-HM2.
(b) fig_ptm_enrichment_dotplot content: per-site PTM ratio (x) vs feature (y),
    dot size = count in top-5%, dot color = -log10(BH q), faceted MED4 | P-HM2.


Data: {enrichment_by_type,enrichment_per_feature}.csv (same source as the two standalone figures).

Output: fig_ptm_enrichment_combined.{pdf,png,tiff}
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

MM = 1 / 25.4  # mm -> inch
WIDTH_MM = 180.0          # journal 2-column fixed width
MAX_H_MM = 185.0          # journal 2-column cap for a <300-word caption
TOP_H_MM = 55.0           # panel (a): fold enrichment
BOT_H_MM = 118.0          # panel (b): dotplot (18 feature rows)
FIG_H_MM = TOP_H_MM + BOT_H_MM + 4.0  # + inter-row gap
assert FIG_H_MM <= MAX_H_MM, f"{FIG_H_MM} exceeds {MAX_H_MM} mm cap"

mpl.rcParams.update({
    "font.size": 6.5, "axes.labelsize": 7, "axes.titlesize": 7,
    "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

SCRATCH = Path("/scratch/path/to")
OUT = SCRATCH / "ptm_enrichment"
OUT.mkdir(parents=True, exist_ok=True)

N_TOP = 51
ENS_DISPLAY = {"MED4": "MED4", "HM2": "P-HM2"}
ENS_COLORS = {"MED4": "#428BBE", "HM2": "#FF9333"}
PTM_DISPLAY = {"CGL": "SSG", "CSO": "SOH", "SNC": "SNO"}
PTM_ORDER = ["CGL", "CSO", "SNC"]
CMAP = "viridis_r"


def sig_stars(q):
    if q < 1e-4: return "****"
    if q < 1e-3: return "***"
    if q < 1e-2: return "**"
    if q < 5e-2: return "*"
    return ""


def save_all(fig, base):
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    fig.savefig(base.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


# ---------- panel (a): fold enrichment by PTM type ----------
def load_occ(ens):
    occ = pd.read_csv(
        SCRATCH / f"{ens}_NAD/postprocess_final/ptm_enrichment/ptm_occupancy_binary.csv",
        dtype={"system": str})
    occ = occ[occ["in_top5pct"].isin([True, False, "True", "False"])].copy()
    occ["top"] = occ["in_top5pct"].astype(str) == "True"
    for t in PTM_ORDER:
        cols = [c for c in occ.columns if c.startswith(f"PTM_{t}-")]
        occ[t] = occ[cols].sum(axis=1)
    return occ


def load_bhq_type(ens):
    bt = pd.read_csv(SCRATCH / f"{ens}_NAD/postprocess_final/ptm_enrichment/enrichment_by_type.csv")
    return dict(zip(bt["group"], bt["mw_count_q_bh"]))


def fold_enrichment_table():
    rows = []
    for ens in ("MED4", "HM2"):
        occ = load_occ(ens)
        bhq = load_bhq_type(ens)
        for t in PTM_ORDER:
            top_mean = occ.loc[occ["top"], t].mean()
            bg_mean = occ.loc[~occ["top"], t].mean()
            log2fe = np.log2(top_mean / bg_mean) if bg_mean > 0 else np.nan
            rows.append({"ensemble": ens, "ptm": PTM_DISPLAY[t],
                         "log2fe": log2fe, "q": bhq[t]})
    return pd.DataFrame(rows)


def draw_fold_enrichment(fig, gs_row):
    fe = fold_enrichment_table()
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_row, wspace=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
    for ax, ens in zip(axes, ("MED4", "HM2")):
        sub = fe[fe["ensemble"] == ens].set_index("ptm").loc[[PTM_DISPLAY[t] for t in PTM_ORDER]]
        x = np.arange(len(sub))
        ax.bar(x, sub["log2fe"], color=ENS_COLORS[ens], edgecolor="black",
               linewidth=0.5, width=0.6)
        for xi, (lf, q) in enumerate(zip(sub["log2fe"], sub["q"])):
            # avoid a misleading "-0.00" when the rounded value is exactly zero
            label = "0.00" if round(lf, 2) == 0 else f"{lf:+.2f}"
            ax.text(xi, lf + 0.03 if lf >= 0 else lf - 0.03, label,
                    ha="center", va="bottom" if lf >= 0 else "top", fontsize=6)
            stars = sig_stars(q)
            if stars:
                ax.text(xi, (lf if lf > 0 else 0) + 0.14, stars,
                        ha="center", va="bottom", fontsize=7, fontweight="bold")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(x); ax.set_xticklabels(sub.index)
        ax.set_title(ENS_DISPLAY[ens])
        ax.set_xlabel("PTM type")
        ax.set_ylim(-0.45, 1.75)
        if ens == "HM2":
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel(r"log$_2$ fold enrichment" "\n(top 5% / background)")
    return axes


# ---------- panel (b): per-site dotplot ----------
def load_feat(ens):
    df = pd.read_csv(SCRATCH / f"{ens}_NAD/postprocess_final/ptm_enrichment/enrichment_per_feature.csv")
    df["ptm_ratio"] = df["n_present_top"] / N_TOP
    df["count"] = df["n_present_top"]
    df["neglog10_padj"] = -np.log10(df["fisher_q_bh"].clip(lower=1e-10))
    df["label"] = df.apply(
        lambda r: f"{PTM_DISPLAY.get(r['ptm_type'], r['ptm_type'])}-"
                  f"{str(r['subunit']).replace('GAPDH','GAP2')} Cys{int(r['residue'])}",
        axis=1)
    df["ensemble"] = ens
    return df


def draw_dotplot(fig, gs_row):
    combined = pd.concat([load_feat(e) for e in ("MED4", "HM2")], ignore_index=True)
    order = combined.groupby("label")["ptm_ratio"].max().sort_values().index.tolist()
    y_for = {lab: i for i, lab in enumerate(order)}
    qmin, qmax = combined["neglog10_padj"].min(), combined["neglog10_padj"].max()
    norm = Normalize(vmin=qmin, vmax=qmax)
    cmap = mpl.colormaps[CMAP]

    # 3 columns: MED4 panel | P-HM2 panel | colorbar
    gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_row,
                                          width_ratios=[1, 1, 0.05], wspace=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
    for ax, ens in zip(axes, ("MED4", "HM2")):
        d = combined[combined["ensemble"] == ens]
        sizes = d["count"] * 4.5
        ax.scatter(d["ptm_ratio"], [y_for[l] for l in d["label"]],
                   s=sizes, c=d["neglog10_padj"], cmap=cmap, norm=norm,
                   edgecolor="none", alpha=0.95)
        ax.set_title(ENS_DISPLAY[ens])
        ax.set_xlabel("PTM ratio")
        ax.grid(axis="both", color="#e8e8e8", lw=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.02, combined["ptm_ratio"].max() * 1.12)
        if ens == "HM2":
            ax.tick_params(labelleft=False)
    axes[0].set_yticks(list(y_for.values()))
    axes[0].set_yticklabels(order, fontsize=6)

    # reserve headroom above the colorbar for the size legend (top-right of panel b)
    cb_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 2],
                                             height_ratios=[0.30, 0.70], hspace=0.0)
    cax = fig.add_subplot(cb_gs[1, 0])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cbar.set_label(r"$-\log_{10}$ BH q", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.5)

    size_vals = [5, 10, 15, 19]
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="grey", markeredgecolor="none",
                      markersize=np.sqrt(v * 4.5) / 1.1, label=str(v))
               for v in size_vals]
    legend_ax = fig.add_subplot(cb_gs[0, 0]); legend_ax.axis("off")
    legend_ax.legend(handles=handles, title="Count", loc="center",
                     frameon=False, labelspacing=1.0, fontsize=6, title_fontsize=6.5)
    return axes


def main():
    fig = plt.figure(figsize=(WIDTH_MM * MM, FIG_H_MM * MM))
    outer = gridspec.GridSpec(2, 1, height_ratios=[TOP_H_MM, BOT_H_MM],
                              hspace=0.32, top=0.94, bottom=0.09,
                              left=0.16, right=0.93)
    draw_fold_enrichment(fig, outer[0])
    draw_dotplot(fig, outer[1])

    # panel letters
    fig.text(0.01, 0.965, "a", fontsize=8, fontweight="bold")
    fig.text(0.01, 0.66, "b", fontsize=8, fontweight="bold")

    base = OUT / "fig_ptm_enrichment_combined"
    save_all(fig, base)
    print(f"Saved to {base} ({WIDTH_MM:.1f} x {FIG_H_MM:.1f} mm)")


if __name__ == "__main__":
    main()
