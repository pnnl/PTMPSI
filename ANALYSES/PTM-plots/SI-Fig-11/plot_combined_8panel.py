#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_combined_8panel.py — publication-quality 8-panel combined figure (MED4 + HM2).

Layout (2 cols × 4 rows):
         MED4            HM2
  (A) RMSD         (B) RMSD
  (C) RMSF         (D) RMSF
  (E) Rg           (F) Rg
  (G) SecStruct    (H) SecStruct

Reads pre-computed cache .npz files only (no trajectory access required).

Usage:
  python combined_8panel/plot_combined_8panel.py \
    --med4_cache MED4_NAD/combined_analysis/all_systems/cache_matrices \
    --hm2_cache  HM2_NAD/combined_analysis/all_systems/cache_matrices
"""
import os
import glob
import argparse
from collections import defaultdict

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.patches import Patch

# ── Constants ────────────────────────────────────────────────────────────────
STRAIN_COLORS = {"MED4": "#1f4e79", "HM2": "#c55a11"}
ALL_CODES = ["H", "G", "I", "E", "B", "T", "S", "C"]
SUBUNIT_BY_SEGID = {
    "A": "GAP2", "D": "GAP2", "E": "GAP2", "F": "GAP2",
    "I": "GAP2", "K": "GAP2", "N": "GAP2", "O": "GAP2",
    "B": "PRK",  "G": "PRK",  "J": "PRK",  "L": "PRK",
    "C": "CP12", "H": "CP12", "M": "CP12", "P": "CP12",
    "GAP": "GAP2", "GAPDH": "GAP2", "GA": "GAP2", "GB": "GAP2",
    "GAP_A": "GAP2", "GAP_B": "GAP2",
    "PRK": "PRK", "PR": "PRK", "PK": "PRK", "PRK_A": "PRK", "PRK_B": "PRK",
    "CP12": "CP12", "CP": "CP12", "CP_A": "CP12", "CP_B": "CP12",
}
PANEL_LETTERS = list("abcd")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "legend.fontsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.labelpad": 3,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
})

# ── Cache loading ─────────────────────────────────────────────────────────────

def _load_npz(cache_dir, metric):
    matches = sorted(glob.glob(os.path.join(cache_dir, f"{metric}_*.npz")))
    if not matches:
        raise FileNotFoundError(f"No {metric}_*.npz in {cache_dir}")
    print(f"[CACHE] {metric}: {os.path.basename(matches[0])}")
    return np.load(matches[0], allow_pickle=True)


def load_line(cache_dir, metric):
    C = _load_npz(cache_dir, metric)
    systems = C["systems"].tolist()
    return {s: (np.asarray(x, float), np.asarray(y, float))
            for s, x, y in zip(systems, C["x_list"].tolist(), C["y_list"].tolist())}


def load_rmsf(cache_dir):
    C = _load_npz(cache_dir, "rmsf")
    return C["arr"].copy(), C["labs0"].tolist()


def load_secstruct(cache_dir):
    C = _load_npz(cache_dir, "secstruct")
    systems = C["valid_systems"].tolist()
    frac_arr = C["frac_arr"]
    fracs = {s: {code: float(frac_arr[i, j]) for j, code in enumerate(ALL_CODES)}
             for i, s in enumerate(systems)}
    return systems, fracs

# ── Draw: RMSD / Rg line panels ───────────────────────────────────────────────

def draw_line(ax, datasets, strain_color, ylim=None,
              nonref_alpha=0.25, nonref_lw=0.8, ref_lw=1.8,
              ylabel=None, xlabel=None, yticks=True, title=None):
    cmap = plt.get_cmap("tab20")
    for sid, (x, y) in datasets.items():
        if sid == "REF":
            continue
        try:
            c = cmap(int(sid) % 20)
        except (ValueError, TypeError):
            c = cmap(abs(hash(sid)) % 20)
        ax.plot(x, y, color=c, alpha=nonref_alpha, lw=nonref_lw, zorder=1)
    if "REF" in datasets:
        xr, yr = datasets["REF"]
        ax.plot(xr, yr, color="black", lw=ref_lw, zorder=10)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if not yticks:
        ax.tick_params(axis="y", labelleft=False)
    ax.margins(x=0)
    ax.tick_params(axis="both", labelsize=6)
    if title:
        ax.set_title(title, fontsize=7, fontweight="bold", color="black", pad=6)
    return ax

# ── Draw: RMSF heatmap ────────────────────────────────────────────────────────

def _col_order(labs0):
    chain_ids = [lbl.split("-")[-1] for lbl in labs0]
    chain_to_idx = defaultdict(list)
    for i, ch in enumerate(chain_ids):
        chain_to_idx[ch].append(i)
    def rank(ch):
        sub = SUBUNIT_BY_SEGID.get(ch, ch)
        return ({"GAP2": 0, "PRK": 1, "CP12": 2}.get(sub, 3), ch)
    order = sorted(chain_to_idx.keys(), key=rank)
    return np.array([i for ch in order for i in chain_to_idx[ch]], dtype=int)


def draw_rmsf(fig, ss, arr, labs0, vmin, vmax,
              show_ylabel=True, show_xlabel=True, show_colorbar=True, title=None):
    """Draw RMSF heatmap into SubplotSpec ss. Returns primary ax for panel letter."""
    col_ord = _col_order(labs0)
    if col_ord.size == arr.shape[1]:
        arr = arr[:, col_ord]
        labs0 = [labs0[i] for i in col_ord]

    n_nonref = max(1, arr.shape[0] - 1)
    ref_ratio = max(18, int(0.04 * n_nonref))

    Z = None
    row_order = np.arange(n_nonref)
    if n_nonref > 2:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            Z = linkage(arr[1:], method="average", metric="euclidean")
            row_order = leaves_list(Z)
        except Exception:
            Z = None
    arr_plot = np.vstack([arr[0:1], arr[1:][row_order]])

    if show_colorbar:
        # 4 cols: dendrogram | heatmap | gap | colorbar
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 4, subplot_spec=ss,
            width_ratios=[1.1, 20.6, 0.5, 0.7],
            height_ratios=[ref_ratio, n_nonref],
            hspace=0.05, wspace=0.02,
        )
        cax = fig.add_subplot(inner[:, 3])
    else:
        # 2 cols: dendrogram | heatmap (no colorbar)
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=ss,
            width_ratios=[1.1, 20.6],
            height_ratios=[ref_ratio, n_nonref],
            hspace=0.05, wspace=0.02,
        )
        cax = None

    axd0 = fig.add_subplot(inner[0, 0])
    ax0  = fig.add_subplot(inner[0, 1])
    axd1 = fig.add_subplot(inner[1, 0])
    ax1  = fig.add_subplot(inner[1, 1], sharex=ax0)

    ax0.imshow(arr_plot[0:1], aspect="auto", origin="upper",
               cmap="YlGn", vmin=vmin, vmax=vmax, interpolation="nearest")
    im = ax1.imshow(arr_plot[1:], aspect="auto", origin="upper",
                    cmap="YlGn", vmin=vmin, vmax=vmax, interpolation="nearest")

    axd0.axis("off")
    if Z is not None and n_nonref > 2:
        try:
            from scipy.cluster.hierarchy import dendrogram
            dendrogram(Z, orientation="left", no_labels=True,
                       truncate_mode="lastp", p=min(50, n_nonref),
                       color_threshold=0, above_threshold_color="k",
                       link_color_func=lambda _: "k", ax=axd1)
            for coll in axd1.collections:
                coll.set_linewidth(0.35)
            axd1.set_xticks([]); axd1.set_yticks([])
            for side in ("top", "right", "bottom", "left"):
                axd1.spines[side].set_visible(False)
            axd1.margins(x=0, y=0)
        except Exception:
            axd1.axis("off")
    else:
        axd1.axis("off")

    chain_ids = [lbl.split("-")[-1] for lbl in labs0]
    bounds = [i for i in range(1, len(chain_ids)) if chain_ids[i] != chain_ids[i - 1]]
    for b in bounds:
        ax0.axvline(b - 0.5, ls="-", c="k", lw=0.35, alpha=0.9)
        ax1.axvline(b - 0.5, ls="-", c="k", lw=0.35, alpha=0.9)

    ax0.set_yticks([0]); ax0.set_yticklabels(["REF"], fontsize=6)
    ax1.set_yticks([]); ax1.set_yticklabels([])
    ax0.xaxis.set_tick_params(bottom=False, labelbottom=False)
    ax1.set_xticks([]); ax1.set_xticklabels([])
    ax1.set_xlim(0, len(labs0) - 1)

    if show_xlabel:
        starts = [0] + bounds
        ends   = bounds + [len(chain_ids)]
        grouped = []
        for s, e in zip(starts, ends):
            sub = SUBUNIT_BY_SEGID.get(chain_ids[s], chain_ids[s])
            if grouped and grouped[-1][2] == sub:
                grouped[-1] = (grouped[-1][0], e, sub)
            else:
                grouped.append((s, e, sub))
        trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
        for s, e, sub in grouped:
            x0, x1 = s - 0.5, e - 0.5
            xm = (x0 + x1) / 2.0
            ax1.plot([x0, x1], [-0.07, -0.07], transform=trans, c="k", lw=0.5, clip_on=False)
            ax1.plot([x0, x0], [-0.07, -0.045], transform=trans, c="k", lw=0.5, clip_on=False)
            ax1.plot([x1, x1], [-0.07, -0.045], transform=trans, c="k", lw=0.5, clip_on=False)
            ax1.text(xm, -0.10, sub, transform=trans, ha="center", va="top",
                     fontsize=7, clip_on=False)

    if cax is not None:
        cbar = fig.colorbar(im, cax=cax)
        if show_ylabel:
            cbar.set_label(r"RMSF ($\mathrm{\AA}$)", fontsize=7, labelpad=2)
        cbar.ax.tick_params(labelsize=6)

    if title:
        ax0.set_title(title, fontsize=7, fontweight="bold", color="black", pad=6)

    return ax0  # primary ax for panel letter

# ── Draw: Secondary structure ─────────────────────────────────────────────────

def draw_secstruct(fig, ss, valid_systems, fracs, show_ylabel=True, show_legend=True,
                   title=None):
    """Draw stacked bar secondary structure into SubplotSpec ss. Returns primary ax."""
    nonref = [s for s in valid_systems if s != "REF"]
    has_ref = "REF" in valid_systems

    Z = None
    if len(nonref) > 2:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            M = np.array([[fracs[s][c] for c in ALL_CODES] for s in nonref], float)
            Z = linkage(M, method="average", metric="euclidean")
            nonref = [nonref[i] for i in leaves_list(Z)]
        except Exception:
            Z = None

    ordered = (["REF"] + nonref) if has_ref else nonref
    do_dendro = Z is not None and len(nonref) > 2
    cmap = plt.get_cmap("tab10")
    # Only show legend entries for codes that actually occur (drop all-zero codes).
    present_idx = [i for i, c in enumerate(ALL_CODES)
                   if max((fracs[s][c] for s in valid_systems), default=0.0) > 1e-9]
    legend_handles = [Patch(facecolor=cmap(i), edgecolor="none", label=ALL_CODES[i])
                      for i in present_idx]

    if has_ref:
        if do_dendro:
            inner = gridspec.GridSpecFromSubplotSpec(
                2, 2, subplot_spec=ss,
                width_ratios=[1.6, 8.4], height_ratios=[0.45, 5.55],
                wspace=0.08, hspace=0.00)
            ax_blank = fig.add_subplot(inner[0, 0]); ax_blank.axis("off")
            ax_den = fig.add_subplot(inner[0, 1])
            ax_ref = fig.add_subplot(inner[1, 0])
            ax     = fig.add_subplot(inner[1, 1], sharey=ax_ref)
        else:
            inner = gridspec.GridSpecFromSubplotSpec(
                1, 2, subplot_spec=ss, width_ratios=[1.6, 8.4], wspace=0.08)
            ax_ref = fig.add_subplot(inner[0, 0])
            ax_den = None
            ax     = fig.add_subplot(inner[0, 1], sharey=ax_ref)
    else:
        inner  = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=ss)
        ax     = fig.add_subplot(inner[0, 0])
        ax_ref = None
        ax_den = None

    x_ref    = np.array([0.0]) if has_ref else np.array([])
    x_nonref = ((5.0 + 10.0 * np.arange(len(nonref), dtype=float))
                if do_dendro else np.arange(len(nonref), dtype=float))
    w_ref    = 0.475
    w_nonref = 7.2 if do_dendro else 0.78
    bot_ref    = np.zeros(len(x_ref))
    bot_nonref = np.zeros(len(x_nonref))

    for i, code in enumerate(ALL_CODES):
        vals = np.array([fracs[s][code] for s in ordered], float)
        if has_ref and ax_ref is not None:
            ax_ref.bar(x_ref, vals[:1], bottom=bot_ref, width=w_ref,
                       color=cmap(i), align="center", edgecolor="black", linewidth=0.9)
            bot_ref += vals[:1]
            vals_nr = vals[1:]
        else:
            vals_nr = vals
        if len(vals_nr) and len(x_nonref):
            ax.bar(x_nonref, vals_nr, bottom=bot_nonref, width=w_nonref,
                   color=cmap(i), align="center")
            bot_nonref += vals_nr

    if has_ref and ax_ref is not None:
        ax_ref.set_xticks([0.0]); ax_ref.set_xticklabels(["REF"], fontsize=6)
        ax_ref.set_ylim(0, 1); ax_ref.set_xlim(-0.85, 0.85)
        ax_ref.tick_params(axis="y", labelsize=6)
        if show_ylabel:
            ax_ref.set_ylabel("Fraction", fontsize=7, labelpad=2)
        ax.set_xticks([]); ax.set_xticklabels([])
        ax.set_ylim(0, 1)
        ax.tick_params(axis="y", left=False, labelleft=False)
        if do_dendro and ax_den is not None:
            try:
                from scipy.cluster.hierarchy import dendrogram
                dendrogram(Z, orientation="top", no_labels=True,
                           truncate_mode="lastp", p=min(80, len(nonref)),
                           color_threshold=0, above_threshold_color="k",
                           link_color_func=lambda _: "k", ax=ax_den)
                for coll in ax_den.collections:
                    coll.set_linewidth(0.35)
                ax_den.set_xticks([]); ax_den.set_yticks([])
                for side in ("top", "right", "bottom", "left"):
                    ax_den.spines[side].set_visible(False)
                ax_den.tick_params(axis="both", which="both", length=0)
                ax.set_xlim(0, 10.0 * len(nonref))
            except Exception:
                pass
        else:
            ax.set_xlim(-0.6, max(0.6, len(x_nonref) - 0.4))
        if show_legend:
            ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 0.5),
                      loc="center left", ncol=1, fontsize=7, frameon=False)
        primary = ax_ref
    else:
        ax.set_ylim(0, 1)
        if show_ylabel:
            ax.set_ylabel("Fraction", fontsize=7)
        if show_legend:
            ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1),
                      loc="upper left", ncol=1, fontsize=7, frameon=False)
        primary = ax

    if title:
        bb = ax.get_position()          # main stacked-bar axis (fig fraction)
        ss_bb = ss.get_position(fig)    # full panel bounding box
        fig.text(0.5 * (bb.x0 + bb.x1), ss_bb.y1 + 0.012, title,
                 ha="center", va="bottom", fontsize=7, fontweight="bold",
                 color="black", transform=fig.transFigure)

    return primary


def add_panel_letters(fig, specs):
    """Add 8 pt bold panel letters at the upper-left of each row panel."""
    for letter, ss in zip(PANEL_LETTERS, specs):
        bb = ss.get_position(fig)
        fig.text(max(0.005, bb.x0 - 0.035), bb.y1 + 0.006, letter,
                 ha="left", va="bottom", fontsize=8, fontweight="bold",
                 family="sans-serif", color="black",
                 transform=fig.transFigure)

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Combined 8-panel MED4/HM2 figure from cache matrices")
    p.add_argument("--med4_cache", required=True)
    p.add_argument("--hm2_cache",  required=True)
    p.add_argument("--out_dir", default="combined_8panel")
    p.add_argument("--nonref_alpha", type=float, default=0.25)
    p.add_argument("--nonref_lw",    type=float, default=0.8)
    p.add_argument("--ref_lw",       type=float, default=1.0)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load caches ───────────────────────────────────────────────────────────
    med4_rmsd = load_line(args.med4_cache, "rmsd")
    hm2_rmsd  = load_line(args.hm2_cache,  "rmsd")
    med4_rg   = load_line(args.med4_cache, "rg")
    hm2_rg    = load_line(args.hm2_cache,  "rg")
    med4_rmsf_arr, med4_rmsf_labs = load_rmsf(args.med4_cache)
    hm2_rmsf_arr,  hm2_rmsf_labs  = load_rmsf(args.hm2_cache)
    med4_ss_sys, med4_ss_fracs = load_secstruct(args.med4_cache)
    hm2_ss_sys,  hm2_ss_fracs  = load_secstruct(args.hm2_cache)

    # ── Shared axis limits ────────────────────────────────────────────────────
    def _ylim(da, db, pad=0.05):
        all_y = [v for d in (da, db)
                 for _, y in d.values() for v in y[np.isfinite(y)]]
        lo, hi = min(all_y), max(all_y)
        m = (hi - lo) * pad
        return (lo - m, hi + m)

    rmsd_ylim = _ylim(med4_rmsd, hm2_rmsd)
    rmsd_ylim = (0.0, rmsd_ylim[1])  # clamp RMSD axis to start at 0
    rg_ylim   = _ylim(med4_rg,   hm2_rg)
    rmsf_vmax = max(1.0, float(max(
        np.nanpercentile(med4_rmsf_arr, 99),
        np.nanpercentile(hm2_rmsf_arr,  99),
    )))

    print(f"RMSD ylim: {rmsd_ylim}")
    print(f"Rg   ylim: {rg_ylim}")
    print(f"RMSF vmax: {rmsf_vmax:.2f}")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(6.8, 8.5))
    outer = gridspec.GridSpec(
        4, 2, figure=fig,
        height_ratios=[1, 1.8, 1, 1.8],
        hspace=0.38, wspace=0.22,
        left=0.09, right=0.91,
        bottom=0.05, top=0.94,
    )

    lkw = dict(nonref_alpha=args.nonref_alpha,
               nonref_lw=args.nonref_lw,
               ref_lw=args.ref_lw)

    # Row 0 – RMSD
    ax_A = fig.add_subplot(outer[0, 0])
    draw_line(ax_A, med4_rmsd, STRAIN_COLORS["MED4"],
              ylim=rmsd_ylim, ylabel="RMSD (Å)", xlabel="Time (ns)", **lkw)

    ax_B = fig.add_subplot(outer[0, 1])
    draw_line(ax_B, hm2_rmsd, STRAIN_COLORS["HM2"],
              ylim=rmsd_ylim, xlabel="Time (ns)", **lkw)

    # Row 1 – RMSF
    ax_C = draw_rmsf(fig, outer[1, 0],
                     med4_rmsf_arr, med4_rmsf_labs, 0.0, rmsf_vmax,
                     show_ylabel=False, show_xlabel=True, show_colorbar=False)
    ax_D = draw_rmsf(fig, outer[1, 1],
                     hm2_rmsf_arr, hm2_rmsf_labs, 0.0, rmsf_vmax,
                     show_ylabel=True, show_xlabel=True)

    # Row 2 – Rg
    ax_E = fig.add_subplot(outer[2, 0])
    draw_line(ax_E, med4_rg, STRAIN_COLORS["MED4"],
              ylim=rg_ylim, ylabel="Rg (Å)", xlabel="Time (ns)", **lkw)

    ax_F = fig.add_subplot(outer[2, 1])
    draw_line(ax_F, hm2_rg, STRAIN_COLORS["HM2"],
              ylim=rg_ylim, xlabel="Time (ns)", **lkw)

    # Row 3 – Secondary structure
    ax_G = draw_secstruct(fig, outer[3, 0], med4_ss_sys, med4_ss_fracs,
                          show_ylabel=True, show_legend=False)
    ax_H = draw_secstruct(fig, outer[3, 1], hm2_ss_sys, hm2_ss_fracs,
                          show_ylabel=False, show_legend=True)

    # ── Column headers ────────────────────────────────────────────────────────
    ax_A.set_title("MED4", fontsize=7, fontweight="bold", color="black", pad=4)
    ax_B.set_title("P-HM2", fontsize=7, fontweight="bold", color="black", pad=4)

    add_panel_letters(fig, [outer[i, 0] for i in range(4)])

    # ── Save ──────────────────────────────────────────────────────────────────
    stem = os.path.join(args.out_dir, "combined_8panel")
    for ext in ("pdf", "tiff"):
        path = f"{stem}.{ext}"
        fig.savefig(path, dpi=600, bbox_inches="tight")
        print(f"-> {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
