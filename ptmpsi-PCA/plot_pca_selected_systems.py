#!/usr/bin/env python
# coding: utf-8
# plot_pca_selected_systems.py
"""
Select top/bottom systems and produce:
  - selected_systems.csv       (the 40 systems by global-PCA)
  1) Selected_globalPCA_systems.pdf
  2) Selected_globalPCA-combined_systems.pdf
  3) Selected_combined_systems.pdf
"""
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─── Global styling ──────────────────────────────────────────────────────────
mpl.rcParams["font.size"]   = 16
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["legend.fontsize"] = 16
mpl.rcParams["axes.titlesize"]  = 16
mpl.rcParams["axes.labelsize"]  = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 12
mpl.rcParams["xtick.major.pad"] = 6
mpl.rcParams["ytick.major.pad"] = 6
# ─────────────────────────────────────────────────────────────────────────────

def plot_simple_diverging(x, y, labels, ylabel, out_pdf):
    colors = ["blue" if v>=0 else "orange" for v in y]
    fig, ax = plt.subplots(figsize=(12,6))
    ax.bar(x, y, color=colors, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl.lstrip("0") for lbl in labels], rotation=45, fontsize=14)
    ax.set_xlabel("Ref #", fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.axhline(0, color='black', lw=0.8)
    blue_patch = mpatches.Patch(color='blue',label='+ve PC Score')
    orange_patch = mpatches.Patch(color='orange', label='-ve PC Score')
    ax.legend(handles=[blue_patch, orange_patch], loc='upper right', frameon=False, prop={"family": "monospace"})
    ax.margins(y=0.1) 
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")

def plot_combined_simple_diverging(x, y, labels, ylabel, out_pdf):
    colors = ["indigo" if v>=0 else "gold" for v in y]
    fig, ax = plt.subplots(figsize=(12,6))
    ax.bar(x, y, color=colors, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl.lstrip("0") for lbl in labels], rotation=45, fontsize=14)
    ax.set_xlabel("Ref #", fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.axhline(0, color='black', lw=0.8)
    indigo_patch = mpatches.Patch(color='indigo', label='+ve Combined Score')
    gold_patch   = mpatches.Patch(color='gold',   label='-ve Combined Score')
    ax.legend(handles=[indigo_patch, gold_patch], loc='upper right', frameon=False, prop={"family": "monospace"})
    ax.margins(y=0.1) 
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")

def plot_stacked(x, del_pc1, del_mfpt, del_tau, labels, ylabel, out_pdf):
    pos_bot = np.zeros_like(x, dtype=float)
    neg_bot = np.zeros_like(x, dtype=float)
    fig, ax = plt.subplots(figsize=(12,6))

    for vals, color, label in [
        (del_pc1,  None,        "PC Score"),
        (del_mfpt, "blue",      "MFPT (Mean First-Passage Time)"),
        (del_tau,  "orange",    "τ (Residence Time)")
    ]:
        pos = np.clip(vals,  0, None)
        neg = np.clip(vals, None, 0)

        if label == "PC Score":
            ax.bar(x, pos, bottom=pos_bot, color="blue", label="+ve PC Score")
            ax.bar(x, neg, bottom=neg_bot, color="orange", label="-ve PC Score")
        else:
            ax.bar(x, pos, bottom=pos_bot, color=color, label=label)
            ax.bar(x, neg, bottom=neg_bot, color=color)

        pos_bot += pos
        neg_bot += neg

    ax.axhline(0, color='black', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl.lstrip("0") for lbl in labels], rotation=45, fontsize=14)
    ax.set_xlabel("Ref #", fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend(loc="upper right", frameon=False, ncol =3, prop={"family": "monospace"})
    ax.margins(y=0.1) 
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)
    print(f"Saved {out_pdf}")

def main():
    p = argparse.ArgumentParser(
        description="Select top/bottom 20 systems and compare scores"
    )
    p.add_argument(
        "--base_dir", required=True,
        help="Directory containing PC_score_ranking.csv and combined_score_ranking.csv"
    )
    args = p.parse_args()
    base = args.base_dir

    # load rankings
    gp = pd.read_csv(os.path.join(base, "PC_score_ranking.csv"))
    cs = pd.read_csv(os.path.join(base, "combined_score_ranking.csv"))
    # drop reference
    ref_names = {'REF', 'no_PTMs'}
    gp = gp[~gp['system'].isin(ref_names)].copy()
    cs = cs[~cs['system'].isin(ref_names)].copy()

    # --- choose bottom20 & top20 by global‐PCA ---
    gp_sorted = gp.sort_values("diff_to_ref", ascending=True)
    bottom20  = gp_sorted.head(20)
    top20     = gp_sorted.tail(20)
    sel_gp    = pd.concat([bottom20, top20], ignore_index=True)
    ids_gp    = sel_gp["system"].tolist()

    # write selected_systems.csv
    sel_csv = os.path.join(base, "selected_systems.csv")
    sel_gp.to_csv(sel_csv, index=False)
    # append a final line of just the IDs, comma-separated
    with open(sel_csv, "a") as f:
        f.write(",".join(ids_gp) + "\n")
    print(f"Wrote selection to {sel_csv}")

    # --- Plot #1: simple green/red ΔPC1 ---
    plot_simple_diverging(
        x      = np.arange(len(ids_gp)),
        y      = sel_gp["diff_to_ref"].values,
        labels = ids_gp,
        ylabel = "Δ <PC> Score",
        out_pdf= os.path.join(base, "PC_score_ranking_selected_systems.pdf")
    )

    # --- Plot #1.5: Mapped indigo/gold Δ Combined Score for the same 40 ---
    combined_idx = cs.set_index("system")
    plot_combined_simple_diverging(
        x       = np.arange(len(ids_gp)),
        y       = combined_idx.loc[ids_gp, "diff_to_ref"].values,
        labels  = ids_gp,
        ylabel  = "Δ Combined Score",
        out_pdf = os.path.join(base, "Mapped_Combined_score_selected_systems.pdf")
    )
    # --- Plot #2: stacked ΔPC1, ΔMFPT, Δτ for those same 40 ---
    comp = cs.set_index("system")
    del_pc1  = comp.loc[ids_gp, "delta_PC1"].values
    del_mfpt = comp.loc[ids_gp, "delta_MFPT"].values
    del_tau  = comp.loc[ids_gp, "delta_tau"].values

    plot_stacked(
        x      = np.arange(len(ids_gp)),
        del_pc1  = del_pc1,
        del_mfpt = del_mfpt,
        del_tau  = del_tau,
        labels = ids_gp,
        ylabel = "Δ Combined Score",
        out_pdf= os.path.join(base, "Detailed_Mapped_Combined_score_selected_systems.pdf")
    )

    # --- Plot #3: stacked ΔPC1, ΔMFPT, Δτ for bottom/top 40 by combined ---
    cs_sorted = cs.sort_values("diff_to_ref", ascending=True)
    bot20_cs  = cs_sorted.head(20)
    top20_cs  = cs_sorted.tail(20)
    sel_cs    = pd.concat([bot20_cs, top20_cs], ignore_index=True)
    ids_cs    = sel_cs["system"].tolist()

    plot_stacked(
        x        = np.arange(len(ids_cs)),
        del_pc1  = sel_cs["delta_PC1"].values,
        del_mfpt = sel_cs["delta_MFPT"].values,
        del_tau  = sel_cs["delta_tau"].values,
        labels   = ids_cs,
        ylabel   = "Δ Combined Score",
        out_pdf  = os.path.join(base, "General_Combined_score_ranking.pdf")
    )

if __name__ == "__main__":
    main()
