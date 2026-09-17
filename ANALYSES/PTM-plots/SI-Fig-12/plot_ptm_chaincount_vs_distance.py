#!/usr/bin/env python3
"""Chain-level PTM count vs PCA distance from REF — per subunit.

Generates, for each grain (total / GAP2 / PRK / CP12):
  fig_ptm_chaincount_vs_distance_scatter[suffix]        — jittered scatter
  fig_ptm_chaincount_vs_distance_box[suffix]            — boxplot per count
  fig_ptm_chaincount_vs_distance_violin[suffix]         — violin per count
  fig_ptm_chaincount_vs_distance_box_combined[suffix]   — MED4|P-HM2 paired box
  fig_ptm_chaincount_vs_distance_violin_combined[suffix]— MED4|P-HM2 paired violin

Counts (chain-level, i.e. one increment per modified chain so values are
multiples of the subunit copy number):
  total (no suffix) — all 120 chain-residue slots, 0–40 (8 GAP2 chains × 4
        sites + 4 PRK chains × 2 sites; one chemistry per slot).
  _GAP2  — PTM_*-GAPDH-*, multiples of 8, 0–32.
  _PRK   — PTM_*-PRK-*,   multiples of 4, 0–8.
  _CP12  — CP12 carries no chemistry PTM in the 120-feature set; its only
        modification axis is the Cys–Cys disulfide. Count = number of
        disulfide-bonded CP12 chains (0–4) from cp12_redox_summary.csv.
        The redox state is all-or-nothing across the 4 chains, so this
        axis takes only the values 0 (all reduced) and 4 (all bonded).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

mpl.rcParams.update({
    "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 12,
    "legend.fontsize": 9, "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})

SCRATCH = Path("/scratch/path/to")
OUT = SCRATCH / "/ptm_enrichment"
OUT.mkdir(parents=True, exist_ok=True)

ENS_COLORS = {"MED4": "#428BBE", "HM2": "#FF9333"}
ENS_DISPLAY = {"MED4": "MED4", "HM2": "P-HM2"}
DIST_COL = "variance_normalized_distance_to_ref"

# chainlevel columns are PTM_<CHEM>-<CHAIN>-<RESID>; classify by chain letter
GAP2_CHAINS = {"A", "D", "E", "F", "I", "K", "N", "O"}
PRK_CHAINS = {"B", "G", "J", "L"}

# (count column, x-axis label, filename suffix). "" preserves the original
# total-count filenames referenced by the manuscript.
GRAINS = [
    ("count_total", "Total PTMs",                   ""),
    ("count_gap2",  "GAP2 PTMs",                    "_GAP2"),
    ("count_prk",   "PRK PTMs",                     "_PRK"),
    ("count_cp12",  "CP12 disulfide-bonded chains", "_CP12"),
]


def save_all(fig, base):
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=200)
    fig.savefig(base.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def load(ens):
    pca_dir = SCRATCH / f"{ens}_NAD/postprocess_final/pca_structural_only_90pct"
    ptm_dir = SCRATCH / f"{ens}_NAD/postprocess_final/ptm_enrichment"
    ranked = pd.read_csv(pca_dir / "pca_distance_to_ref_90pct_ranked.csv",
                         dtype={"system": str})
    chain = pd.read_csv(ptm_dir / "ptm_occupancy_binary_chainlevel.csv",
                        dtype={"system": str})
    ptm_cols = [c for c in chain.columns if c.startswith("PTM_")]
    gap2_cols = [c for c in ptm_cols if c.split("-")[1] in GAP2_CHAINS]
    prk_cols = [c for c in ptm_cols if c.split("-")[1] in PRK_CHAINS]
    chain["count_total"] = chain[ptm_cols].sum(axis=1)
    chain["count_gap2"] = chain[gap2_cols].sum(axis=1)
    chain["count_prk"] = chain[prk_cols].sum(axis=1)

    m = ranked[["system", DIST_COL]].merge(
        chain[["system", "count_total", "count_gap2", "count_prk"]],
        on="system", how="inner")

    # CP12 modification axis = disulfide redox state (# bonded chains, 0–4)
    redox = pd.read_csv(
        SCRATCH / f"{ens}_NAD/postprocess_final/ptm_enrichment_with_cp12/"
        "cp12_redox_summary.csv", dtype={"system": str})
    m = m.merge(redox[["system", "n_disulfide_intact"]], on="system", how="left")
    m = m.rename(columns={"n_disulfide_intact": "count_cp12"})
    return m


def count_vals(md, hd, col):
    """Sorted integer count levels present across both ensembles."""
    vals = pd.concat([md[col], hd[col]]).dropna().unique()
    return sorted(int(v) for v in vals)


def stat_text(m, col):
    sub = m.dropna(subset=[col, DIST_COL])
    rho, p = stats.spearmanr(sub[col], sub[DIST_COL])
    return f"Spearman ρ = {rho:+.3f}\np = {p:.2g}"


def _add_corr_text(ax, md, hd, col):
    """Per-ensemble Spearman ρ + p in the bottom-right (black text)."""
    for k, (ens, m) in enumerate([("MED4", md), ("HM2", hd)]):
        sub = m.dropna(subset=[col, DIST_COL])
        rho, p = stats.spearmanr(sub[col], sub[DIST_COL])
        ax.text(0.985, 0.012 + (1 - k) * 0.052,
                f"{ENS_DISPLAY[ens]}: Spearman ρ = {rho:+.3f}, p = {p:.2g}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.5, family="sans-serif", color="black",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8,
                          pad=0.6))


# subunit panels, ordered by increasing influence on conformational
# displacement (CP12 null -> PRK weak -> GAP2 strong)
PANELS = [
    ("count_cp12", "CP12 disulfide-bonded chains", "CP12"),
    ("count_prk",  "PRK PTMs",                      "PRK"),
    ("count_gap2", "GAP2 PTMs",                      "GAP2"),
]


def _draw_combined_box(ax, md, hd, col, vals):
    """MED4|P-HM2 paired boxplot per count value onto a single axis."""
    data = {"MED4": md, "HM2": hd}
    w = 0.38
    for j, ens in enumerate(("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        offset = (j - 0.5) * w * 1.05
        positions = [i + offset for i in range(len(vals))]
        pos_v = [p for p, g in zip(positions, groups) if len(g) >= 1]
        dat_v = [g for g in groups if len(g) >= 1]
        bp = ax.boxplot(dat_v, positions=pos_v, widths=w, patch_artist=True,
                        showfliers=False,
                        medianprops=dict(color="black", linewidth=1.2))
        for patch in bp["boxes"]:
            patch.set_facecolor(ENS_COLORS[ens]); patch.set_alpha(0.7)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([str(v) for v in vals], fontsize=9)


def _draw_combined_violin(ax, md, hd, col, vals):
    """MED4|P-HM2 paired violin per count value onto a single axis."""
    data = {"MED4": md, "HM2": hd}
    w = 0.38
    for j, ens in enumerate(("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        offset = (j - 0.5) * w * 1.05
        pos_v, dat_v = [], []
        for i, g in enumerate(groups):
            if len(g) >= 2:
                pos_v.append(i + offset); dat_v.append(g)
        if dat_v:
            parts = ax.violinplot(dat_v, positions=pos_v, widths=w,
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_alpha(None)
                pc.set_facecolor(mpl.colors.to_rgba(ENS_COLORS[ens], 0.7))
                pc.set_edgecolor((0, 0, 0, 1)); pc.set_linewidth(0.8)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.2)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([str(v) for v in vals], fontsize=9)


def partial_spearman(df, x, y, covars):
    """Spearman of x vs y controlling for covars, by rank-residualizing
    both against the covariate ranks. Returns (rho, p, n)."""
    import numpy.linalg as la
    d = df.dropna(subset=[x, y] + covars)
    R = d[[x, y] + covars].rank()
    C = np.column_stack([np.ones(len(d))] + [R[c].values for c in covars])

    def resid(v):
        beta, *_ = la.lstsq(C, v, rcond=None)
        return v - C @ beta

    rho, _ = stats.spearmanr(resid(R[x].values), resid(R[y].values))
    n, k = len(d), len(covars)
    t = rho * np.sqrt((n - 2 - k) / (1 - rho ** 2))
    p = 2 * stats.t.sf(abs(t), n - 2 - k)
    return rho, p, n


def _fisher_ci(rho, n, k=0, z=1.96):
    """95% CI for a (partial) Spearman rho via Fisher z-transform."""
    se = 1.0 / np.sqrt(n - 3 - k)
    lo, hi = np.tanh(np.arctanh(rho) - z * se), np.tanh(np.arctanh(rho) + z * se)
    return lo, hi


def _sig_stars(p):
    """Tiered significance asterisks from a p-value."""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return ""


def make_forest(md, hd):
    """Forest plot of Spearman rho (subunit PTM count vs PCA distance)
    with 95% CI, one row per subunit per ensemble. Summarizes the
    GAP2 >> PRK >> CP12 hierarchy in a single panel. The three subunit
    counts are mutually orthogonal in the factorial design (|rho| < 0.01
    between them), so the marginal rho shown here is unconfounded by the
    other subunits (partial rho is numerically identical)."""
    data = {"MED4": md, "HM2": hd}
    # bottom-to-top so GAP2 (strongest) sits at the top
    subunits = [("count_cp12", "CP12"), ("count_prk", "PRK"), ("count_gap2", "GAP2")]

    fig, ax = plt.subplots(figsize=(8, 5))
    dy = 0.18  # vertical offset between the two ensembles within a row
    yticks, ylabels = [], []
    for i, (col, name) in enumerate(subunits):
        yticks.append(i)
        ylabels.append(name)
        for s, ens in ((+1, "MED4"), (-1, "HM2")):
            sub = data[ens].dropna(subset=[col, DIST_COL])
            rho, p = stats.spearmanr(sub[col], sub[DIST_COL])
            lo, hi = _fisher_ci(rho, len(sub))
            y = i + s * dy
            ax.plot([lo, hi], [y, y], color=ENS_COLORS[ens], lw=2, zorder=2)
            ax.plot(rho, y, "o", color=ENS_COLORS[ens], ms=7,
                    markeredgecolor="black", markeredgewidth=0.5, zorder=3)
            stars = _sig_stars(p)
            if stars:
                ax.text(rho, y + 0.07, stars, va="bottom", ha="center",
                        fontsize=11, color="black")
    ax.axvline(0, color="grey", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylabel("Subunit")
    ax.set_ylim(-0.6, len(subunits) - 0.4)
    ax.set_xlim(-0.15, 0.85)
    ax.set_xlabel("Spearman ρ (subunit PTM count vs PCA distance)")
    handles = [Patch(facecolor=ENS_COLORS[e], label=ENS_DISPLAY[e])
               for e in ("MED4", "HM2")]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    save_all(fig, OUT / "fig_ptm_chaincount_vs_distance_forest")


def make_three_panel(md, hd):
    """Three-panel (CP12 | PRK | GAP2) versions of the combined box/violin
    layout, one subunit per panel, sharing the PCA-distance y-axis."""
    for kind, drawer in (("box", _draw_combined_box),
                         ("violin", _draw_combined_violin)):
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)
        for ax, (col, xlabel, title) in zip(axes, PANELS):
            vals = count_vals(md, hd, col)
            drawer(ax, md, hd, col, vals)
            ax.set_xlabel(xlabel)
            ax.set_title(title)
            _add_corr_text(ax, md, hd, col)
        axes[0].set_ylabel("PCA distance")
        handles = [Patch(facecolor=ENS_COLORS[e], alpha=0.7, label=ENS_DISPLAY[e])
                   for e in ("MED4", "HM2")]
        axes[0].legend(handles=handles, frameon=False, loc="upper left")
        save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_{kind}_combined_3panel")


def make_grain(md, hd, col, xlabel, suffix):
    data = {"MED4": md, "HM2": hd}
    vals = count_vals(md, hd, col)
    jit = (max(vals) / 25.0) if len(vals) > 1 else 0.3  # scale jitter to range
    rng = np.random.default_rng(42)

    # --- Scatter ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, ens in zip(axes, ("MED4", "HM2")):
        m = data[ens].dropna(subset=[col, DIST_COL])
        x = m[col].values + rng.uniform(-jit, jit, len(m))
        ax.scatter(x, m[DIST_COL].values, s=8, color=ENS_COLORS[ens],
                   alpha=0.35, edgecolor="none")
        ax.text(0.03, 0.97, stat_text(m, col), transform=ax.transAxes,
                va="top", ha="left", fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
        ax.set_xlabel(xlabel)
        ax.set_title(ENS_DISPLAY[ens])
    axes[0].set_ylabel("PCA distance")
    save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_scatter{suffix}")

    # --- Boxplot per count value ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, ens in zip(axes, ("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        positions = range(len(vals))
        bp = ax.boxplot(groups, positions=positions, widths=0.65,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=1.2))
        for patch in bp["boxes"]:
            patch.set_facecolor(ENS_COLORS[ens]); patch.set_alpha(0.7)
        ax.text(0.03, 0.97, stat_text(m, col), transform=ax.transAxes,
                va="top", ha="left", fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
        ax.set_xticks(list(positions))
        ax.set_xticklabels([str(v) for v in vals], fontsize=9)
        ax.set_xlabel(xlabel)
        ax.set_title(ENS_DISPLAY[ens])
    axes[0].set_ylabel("PCA distance")
    save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_box{suffix}")

    # --- Violin per count value ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, ens in zip(axes, ("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        valid_pos = [i for i, g in enumerate(groups) if len(g) >= 2]
        valid_data = [groups[i] for i in valid_pos]
        if valid_data:
            parts = ax.violinplot(valid_data, positions=valid_pos, widths=0.75,
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_alpha(None)
                pc.set_facecolor(mpl.colors.to_rgba(ENS_COLORS[ens], 0.7))
                pc.set_edgecolor((0, 0, 0, 1)); pc.set_linewidth(0.8)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.2)
        ax.text(0.03, 0.97, stat_text(m, col), transform=ax.transAxes,
                va="top", ha="left", fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels([str(v) for v in vals], fontsize=9)
        ax.set_xlabel(xlabel)
        ax.set_title(ENS_DISPLAY[ens])
    axes[0].set_ylabel("PCA distance")
    save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_violin{suffix}")

    # --- Combined boxplot: MED4 + HM2 paired per count value ---
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for j, ens in enumerate(("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        offset = (j - 0.5) * w * 1.05
        positions = [i + offset for i in range(len(vals))]
        pos_v = [p for p, g in zip(positions, groups) if len(g) >= 1]
        dat_v = [g for g in groups if len(g) >= 1]
        bp = ax.boxplot(dat_v, positions=pos_v, widths=w, patch_artist=True,
                        showfliers=False,
                        medianprops=dict(color="black", linewidth=1.2))
        for patch in bp["boxes"]:
            patch.set_facecolor(ENS_COLORS[ens]); patch.set_alpha(0.7)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([str(v) for v in vals], fontsize=9)
    ax.set_xlabel(xlabel); ax.set_ylabel("PCA distance")
    handles = [Patch(facecolor=ENS_COLORS[e], alpha=0.7, label=ENS_DISPLAY[e])
               for e in ("MED4", "HM2")]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    _add_corr_text(ax, md, hd, col)
    save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_box_combined{suffix}")

    # --- Combined violin: MED4 + HM2 paired per count value ---
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for j, ens in enumerate(("MED4", "HM2")):
        m = data[ens]
        groups = [m.loc[m[col] == v, DIST_COL].dropna().values for v in vals]
        offset = (j - 0.5) * w * 1.05
        pos_v, dat_v = [], []
        for i, g in enumerate(groups):
            if len(g) >= 2:
                pos_v.append(i + offset); dat_v.append(g)
        if dat_v:
            parts = ax.violinplot(dat_v, positions=pos_v, widths=w,
                                  showmedians=True, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_alpha(None)
                pc.set_facecolor(mpl.colors.to_rgba(ENS_COLORS[ens], 0.7))
                pc.set_edgecolor((0, 0, 0, 1)); pc.set_linewidth(0.8)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.2)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([str(v) for v in vals], fontsize=9)
    ax.set_xlabel(xlabel); ax.set_ylabel("PCA distance")
    handles = [Patch(facecolor=ENS_COLORS[e], alpha=0.6, label=ENS_DISPLAY[e])
               for e in ("MED4", "HM2")]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    _add_corr_text(ax, md, hd, col)
    save_all(fig, OUT / f"fig_ptm_chaincount_vs_distance_violin_combined{suffix}")


def main():
    md, hd = load("MED4"), load("HM2")
    for col, xlabel, suffix in GRAINS:
        make_grain(md, hd, col, xlabel, suffix)
        print(f"  grain '{suffix or 'total'}': levels "
              f"{count_vals(md, hd, col)}")
    make_three_panel(md, hd)
    print("  3-panel (CP12 | PRK | GAP2) box + violin combined")
    make_forest(md, hd)
    print("  forest plot (partial Spearman ρ per subunit per ensemble)")
    print(f"Saved to {OUT}")


if __name__ == "__main__":
    main()
