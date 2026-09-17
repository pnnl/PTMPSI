#!/usr/bin/env python3
"""
Plot MI coupling matrices and network-difference panels for the 4 conditions.
Reads contact_coupling/contact_coupling_mi.npz.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_ORDER = ["med4_ox", "med4_red", "hm2_ox", "hm2_red"]


def make_mi_panel(npz: dict, partner: str, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.3), sharex=True, sharey=True)
    mats = [npz[f"mi_{c}_{partner}"] for c in CONDITION_ORDER]
    vmax = float(max(m.max() for m in mats))
    for ax, cond, M in zip(axes, CONDITION_ORDER, mats):
        im = ax.imshow(M, cmap="magma", vmin=0, vmax=vmax,
                       extent=(0.5, M.shape[1] + 0.5, M.shape[0] + 0.5, 0.5),
                       interpolation="none")
        ax.set_title(f"{cond}  ({(M > 0.05).sum() // 2} edges ≥ 0.05 bits)",
                     fontsize=10)
        ax.set_xlabel("CP12 HMM col")
        ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70])
    axes[0].set_ylabel("CP12 HMM col")
    fig.suptitle(f"Pairwise mutual information of CP12 HMM-column "
                 f"engagement signals — CP12 vs {partner.upper()}",
                 fontsize=11, y=1.02)
    fig.colorbar(im, ax=axes, fraction=0.015, pad=0.02, label="MI (bits)")
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"fig_mi_matrix_{partner}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"fig_mi_matrix_{partner}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def make_host_diff_panel(npz: dict, partner: str, outdir: Path) -> None:
    """Heatmap of MI(MED4_redox) − MI(HM2_redox) at each redox."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    for ax, redox in zip(axes, ["ox", "red"]):
        m_med4 = npz[f"mi_med4_{redox}_{partner}"]
        m_hm2  = npz[f"mi_hm2_{redox}_{partner}"]
        diff = m_med4 - m_hm2
        v = float(np.max(np.abs(diff)))
        im = ax.imshow(diff, cmap="RdBu_r", vmin=-v, vmax=v,
                       extent=(0.5, diff.shape[1] + 0.5, diff.shape[0] + 0.5, 0.5),
                       interpolation="none")
        ax.set_title(f"{redox}: MI(MED4) − MI(HM2)  ({partner.upper()})",
                     fontsize=10)
        ax.set_xlabel("CP12 HMM col")
        ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70])
    axes[0].set_ylabel("CP12 HMM col")
    fig.colorbar(im, ax=axes, fraction=0.04, pad=0.02, label="ΔMI (bits)")
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"fig_mi_host_diff_{partner}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"fig_mi_host_diff_{partner}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def make_per_hmm_strength(npz: dict, outdir: Path) -> None:
    """For each (condition, partner, HMM col), sum its MI to all other columns -
    a "coupling strength" per residue. Plot side-by-side comparison."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 4.5), sharex=True)
    for ax, partner in zip(axes, ["prk", "gap2"]):
        for cond, color in zip(CONDITION_ORDER, ["#1f77b4", "#aec7e8", "#d62728", "#ff9896"]):
            M = npz[f"mi_{cond}_{partner}"]
            strength = M.sum(axis=1)
            ax.plot(range(1, M.shape[0] + 1), strength,
                    "-o", ms=2.5, lw=1, color=color, label=cond, alpha=0.85)
        ax.set_ylabel("Σ MI to other cols (bits)")
        ax.set_title(f"CP12 HMM-column coupling strength — CP12 vs {partner.upper()}",
                     fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="best")
        ax.grid(lw=0.3, alpha=0.4)
    axes[-1].set_xlabel("CP12 HMM alignment position")
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "fig_mi_strength_per_hmm.pdf")
    fig.savefig(outdir / "fig_mi_strength_per_hmm.png", dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="contact_coupling/contact_coupling_mi.npz")
    ap.add_argument("--outdir", default="contact_coupling_plots")
    args = ap.parse_args()

    arr = np.load(args.npz)
    npz = {k: arr[k] for k in arr.files}

    outdir = Path(args.outdir).resolve()
    for partner in ["prk", "gap2"]:
        make_mi_panel(npz, partner, outdir)
        make_host_diff_panel(npz, partner, outdir)
    make_per_hmm_strength(npz, outdir)
    print(f"wrote coupling plots to {outdir}")


if __name__ == "__main__":
    main()
