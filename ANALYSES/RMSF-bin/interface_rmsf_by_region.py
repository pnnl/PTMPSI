#!/usr/bin/env python3
"""
Compute and plot CP12 RMSF by interface region.

This is the RMSF-only extraction from analyze_section31_strengthening.py.
It writes:
  - section31_strengthening/interface_rmsf.csv
  - FINAL_FIGURES/interface_rmsf_by_region.{pdf,png,tiff}
"""

import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy import stats


BASE_DIR = Path("/path/to")
OUT_DIR = BASE_DIR / "CP12_focused" / "section31_strengthening"
FIG_DIR = BASE_DIR / "CP12_focused" / "FINAL_FIGURES"

REPS = ["rep1", "rep2", "rep3"]
STRIDE = 10

MED4_COLOR = "#1f77b4"
HM2_COLOR = "#ff7f0e"

CP12_CHAINS = ["C", "H", "M", "P"]

# Region definitions from the original Section 3.1 RMSF analysis.
PRK_BINDING_RESIDS_MED4 = list(range(1, 46))
GAP2_BINDING_RESIDS_MED4 = list(range(46, 75))
PRK_BINDING_RESIDS_HM2 = list(range(1, 42))
GAP2_BINDING_RESIDS_HM2 = list(range(42, 71))

CONDITIONS = [
    ("MED4", "0003", "ox", "MED4_NAD"),
    ("MED4", "0000", "red", "MED4_NAD"),
    ("HM2", "0003", "ox", "HM2_NAD"),
    ("HM2", "0000", "red", "HM2_NAD"),
]

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#3a3a3a",
        "axes.linewidth": 0.8,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    }
)


def get_traj_paths(strain_dir, sys_id, rep):
    topo = str(BASE_DIR / strain_dir / sys_id / rep / "solute_fixed.pdb")
    traj = str(BASE_DIR / strain_dir / sys_id / rep / "solute_fit.xtc")
    return topo, traj


def save_fig(fig, name):
    FIG_DIR.mkdir(exist_ok=True)
    for ext in ("pdf", "png", "tiff"):
        outpath = FIG_DIR / f"{name}.{ext}"
        fig.savefig(
            outpath,
            dpi=600 if ext == "tiff" else 300,
            format=ext if ext == "tiff" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"Saved: {name}.pdf/.png/.tiff", flush=True)


def compute_interface_rmsf():
    """Compute mean RMSF in PRK-binding, GAP2-binding, and non-interface bins."""
    print("=" * 60, flush=True)
    print("Interface-restricted CP12 RMSF", flush=True)
    print("=" * 60, flush=True)

    results = []

    for strain, sys_id, redox, strain_dir in CONDITIONS:
        prk_resids = PRK_BINDING_RESIDS_MED4 if strain == "MED4" else PRK_BINDING_RESIDS_HM2
        gap2_resids = GAP2_BINDING_RESIDS_MED4 if strain == "MED4" else GAP2_BINDING_RESIDS_HM2

        for rep in REPS:
            topo, traj = get_traj_paths(strain_dir, sys_id, rep)
            if not os.path.exists(topo) or not os.path.exists(traj):
                print(f"  {strain} {redox} {rep}: missing topology or trajectory", flush=True)
                continue

            u = mda.Universe(topo, traj)

            for chain in CP12_CHAINS:
                cp12 = u.select_atoms(f"segid {chain} and protein and name CA")
                if cp12.n_atoms == 0:
                    continue

                coords = []
                for i, ts in enumerate(u.trajectory):
                    if i % STRIDE != 0:
                        continue
                    coords.append(cp12.positions.copy())

                coords = np.array(coords)
                mean_pos = coords.mean(axis=0)
                rmsf = np.sqrt(((coords - mean_pos) ** 2).sum(axis=2).mean(axis=0))
                resids = [r.resid for r in cp12.residues]

                for resid, rmsf_val in zip(resids, rmsf):
                    if resid in prk_resids:
                        region = "PRK-binding"
                    elif resid in gap2_resids:
                        region = "GAP2-binding"
                    else:
                        region = "Non-interface"

                    results.append(
                        {
                            "strain": strain,
                            "redox": redox,
                            "rep": rep,
                            "chain": chain,
                            "resid": resid,
                            "region": region,
                            "rmsf": rmsf_val,
                        }
                    )

            print(f"  {strain} {redox} {rep}: done", flush=True)

    if not results:
        raise RuntimeError("No RMSF results were generated; check input paths.")

    df = pd.DataFrame(results)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DIR / "interface_rmsf.csv", index=False)

    summary = (
        df.groupby(["strain", "redox", "region"])
        .agg(mean_rmsf=("rmsf", "mean"), std_rmsf=("rmsf", "std"), n=("rmsf", "count"))
        .reset_index()
    )
    summary["sem_rmsf"] = summary["std_rmsf"] / np.sqrt(summary["n"])

    print("\nInterface RMSF summary:", flush=True)
    for _, row in summary.iterrows():
        print(
            f"  {row.strain} {row.redox} {row.region}: "
            f"{row.mean_rmsf:.2f} +/- {row.sem_rmsf:.2f} A (n={row.n})",
            flush=True,
        )

    print("\nWelch t-tests (MED4 vs HM2):", flush=True)
    for redox in ["ox", "red"]:
        for region in ["PRK-binding", "GAP2-binding", "Non-interface"]:
            med4 = df[
                (df["strain"] == "MED4")
                & (df["redox"] == redox)
                & (df["region"] == region)
            ]["rmsf"]
            hm2 = df[
                (df["strain"] == "HM2")
                & (df["redox"] == redox)
                & (df["region"] == region)
            ]["rmsf"]
            if len(med4) > 0 and len(hm2) > 0:
                _, p_value = stats.ttest_ind(med4, hm2, equal_var=False)
                print(
                    f"  {redox} {region}: "
                    f"MED4={med4.mean():.2f}, HM2={hm2.mean():.2f}, p={p_value:.4f}",
                    flush=True,
                )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    regions = ["PRK-binding", "GAP2-binding", "Non-interface"]
    x = np.arange(len(regions))
    width = 0.35

    for ax_idx, redox in enumerate(["ox", "red"]):
        ax = axes[ax_idx]
        for i, strain in enumerate(["MED4", "HM2"]):
            vals = []
            errs = []
            for region in regions:
                sub = df[
                    (df["strain"] == strain)
                    & (df["redox"] == redox)
                    & (df["region"] == region)
                ]
                vals.append(sub["rmsf"].mean())
                errs.append(sub["rmsf"].std() / np.sqrt(len(sub)))

            color = MED4_COLOR if strain == "MED4" else HM2_COLOR
            ax.bar(
                x + (i - 0.5) * width,
                vals,
                width,
                yerr=errs,
                color=color,
                alpha=0.85,
                edgecolor="#333",
                linewidth=0.4,
                label=strain,
                error_kw={"lw": 0.8, "capsize": 3},
            )

        ax.set_xticks(x)
        ax.set_xticklabels(regions, fontsize=9)
        ax.set_ylabel("Mean RMSF (A)", fontsize=9)
        ax.set_title("Oxidized" if redox == "ox" else "Reduced", fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("CP12 RMSF by interface region", fontsize=12, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "interface_rmsf_by_region")

    return df


def main():
    compute_interface_rmsf()


if __name__ == "__main__":
    main()
