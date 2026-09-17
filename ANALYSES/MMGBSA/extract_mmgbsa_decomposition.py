#!/usr/bin/env python3
"""Parse `FINAL_DECOMP_MMPBSA.dat` DELTAS block for every top-5% system
and write a long-format per-residue ΔTOTAL table per ensemble.

The decomposition file is produced by `gmx_MMPBSA` with `idecomp = 2`
(per-residue, including 1-4 interactions in EEL and VDW). The DELTAS
block at the end of the file holds Δ = Complex − Receptor − Ligand
broken down per chain-symmetry-collapsed residue. The residue index
in the DELTAS block restarts when the chain changes (GAPDH → PRK,
PRK → CP12 ligand, then NAD cofactor block), so we detect chain
transitions by index decrease.

Chain → subunit mapping
-----------------------
- Block 1 (R::, indices ≤ ~340)   → GAPDH (chain-symmetric average over 8 chains)
- Block 2 (R::, indices ~2-300)   → PRK   (chain-symmetric average over 4 chains)
- Block 3 (L::, indices 1-74)     → CP12  (chain-symmetric average over 4 chains; the ligand)
- Block 4 (R::NAD:1-8)            → NAD   (8 cofactor instances, one per GAPDH chain)

Outputs
-------
PRK_SSG_RERUN/ensemble/mmgbsa_summary_structural_only_90pct/decomposition/
   decomp_per_residue_<ENS>.csv     -- long format: system, subunit, resi, resname, dTotal_*, dVDW_*, dEEL_*, dPS_*, dNPS_*
   decomp_per_subunit_<ENS>.csv     -- per-system per-subunit summed dTOTAL
   decomp_combined.csv              -- both ensembles concatenated
   extraction_summary.txt
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

SCRATCH = Path("/scratch/path/to")
ENS_BASES = {
    "MED4": SCRATCH / "MED4_NAD" / "postprocess_final",
    "HM2":  SCRATCH / "HM2_NAD"  / "postprocess_final",
}
ENS_TOP5_SOURCE = {
    "MED4": SCRATCH / "MED4_NAD" / "postprocess_final" / "ptm_enrichment" / "ptm_occupancy_binary.csv",
    "HM2":  SCRATCH / "HM2_NAD"  / "postprocess_final" / "ptm_enrichment" / "ptm_occupancy_binary.csv",
}
OUT_DIR = SCRATCH / "PRK_SSG_RERUN" / "ensemble" / "mmgbsa_summary_structural_only_90pct" / "decomposition"

# DELTAS row pattern: 'X::RESNAME:idx,avg1,sd1,sem1,avg2,sd2,sem2,...'
# 6 component groups × 3 stats = 18 numeric fields after the residue label
RES_PATTERN = re.compile(
    r"^([RL])::([A-Z0-9]+):(\d+)," +
    r",?".join([r"(-?\d+\.\d+|-?\d+\.\d+e[+-]\d+|\d+\.\d+|0.0)" for _ in range(18)])
)


def parse_deltas(path: Path) -> list[dict]:
    """Yield per-residue dicts from the DELTAS block of a decomp file."""
    with path.open() as f:
        lines = f.readlines()
    # Find "DELTAS:" marker line
    start = None
    for i, line in enumerate(lines):
        if line.strip() == "DELTAS:":
            start = i
            break
    if start is None:
        return []

    rows: list[dict] = []
    prev_idx = -1
    block_idx = 0  # 0 = GAPDH, 1 = PRK, 2 = CP12 (L::), 3 = NAD
    block_label_for_R = ["GAPDH", "PRK"]  # consumed in order for R:: blocks before NAD

    for line in lines[start + 1:]:
        line = line.rstrip("\n")
        if not line or not (line.startswith("R::") or line.startswith("L::")):
            continue
        # Split on commas to extract residue label and 18 floats
        parts = line.split(",")
        head = parts[0]
        try:
            cls, rest = head.split("::", 1)
            resname, resi_str = rest.split(":")
            resi = int(resi_str)
        except ValueError:
            continue
        try:
            vals = [float(x) for x in parts[1:19]]
        except (IndexError, ValueError):
            continue
        if len(vals) < 18:
            continue

        # Detect block transitions
        if cls == "L":
            subunit = "CP12"
        elif resname == "NAD":
            subunit = "NAD"
            block_idx = 3
        else:
            # R:: non-NAD: GAPDH (block 0) until index decreases, then PRK (block 1)
            if prev_idx > 0 and resi < prev_idx and block_idx < 1:
                block_idx = 1
            subunit = block_label_for_R[min(block_idx, 1)]
        prev_idx = resi if cls == "R" and resname != "NAD" else prev_idx

        # Components in DELTAS row order:
        #   Internal, vdW, EEL, Polar Solv, Non-Polar Solv, TOTAL
        # Each has 3 stats: avg, sd, sem
        avgs = vals[0::3]   # 6 values: int, vdw, eel, ps, nps, total
        sds  = vals[1::3]
        sems = vals[2::3]
        rows.append({
            "subunit": subunit,
            "residue": resi,
            "resname": resname,
            "class": cls,
            "dInternal_mean": avgs[0],
            "dVDW_mean":      avgs[1],
            "dEEL_mean":      avgs[2],
            "dPolarSolv_mean": avgs[3],
            "dNonPolarSolv_mean": avgs[4],
            "dTOTAL_mean":    avgs[5],
            "dTOTAL_sd":      sds[5],
            "dTOTAL_sem":     sems[5],
        })
    return rows


def extract_ensemble(ens: str) -> pd.DataFrame:
    occ = pd.read_csv(ENS_TOP5_SOURCE[ens], dtype={"system": str})
    top5 = set(occ.loc[occ["in_top5pct"] == True, "system"])
    base = ENS_BASES[ens]

    all_rows = []
    n_systems = 0
    skipped = []
    for s in sorted(top5):
        path = base / s / "mmgbsa" / "FINAL_DECOMP_MMPBSA.dat"
        if not path.exists():
            skipped.append(s)
            continue
        rows = parse_deltas(path)
        if not rows:
            skipped.append(s)
            continue
        n_systems += 1
        for row in rows:
            row["system"] = s
            row["ensemble"] = ens
            all_rows.append(row)
    if skipped:
        print(f"[{ens}] WARNING: {len(skipped)} top-5% systems with unparseable decomp: {skipped[:5]}")
    df = pd.DataFrame(all_rows)
    print(f"[{ens}] parsed {n_systems} systems x {len(df) // max(n_systems,1)} residues each")
    return df


def per_subunit(df: pd.DataFrame) -> pd.DataFrame:
    """Sum ΔTOTAL across residues in each subunit per system.

    Note: GAPDH and PRK rows in the DELTAS block are already
    chain-symmetric averages (a single residue position represents the
    average over all symmetric chains). The sum across residue positions
    gives the per-chain-average summed contribution; multiplying by the
    chain count would give the absolute total contribution from all
    chains. We report both.
    """
    chain_count = {"GAPDH": 8, "PRK": 4, "CP12": 4, "NAD": 1}
    g = df.groupby(["ensemble", "system", "subunit"]).agg(
        n_residues=("residue", "count"),
        sum_dTOTAL_per_chain=("dTOTAL_mean", "sum"),
    ).reset_index()
    g["chain_count"] = g["subunit"].map(chain_count)
    g["sum_dTOTAL_all_chains"] = g["sum_dTOTAL_per_chain"] * g["chain_count"]
    return g


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frames = []
    sub_frames = []
    summary_lines = ["MM/GBSA decomposition extraction (top-5% systems only)"]
    for ens in ("MED4", "HM2"):
        df = extract_ensemble(ens)
        if df.empty:
            print(f"[{ens}] NO data; skipping")
            continue
        df.to_csv(OUT_DIR / f"decomp_per_residue_{ens}.csv", index=False)
        sub = per_subunit(df)
        sub.to_csv(OUT_DIR / f"decomp_per_subunit_{ens}.csv", index=False)
        frames.append(df)
        sub_frames.append(sub)
        # Summary
        n_sys = df["system"].nunique()
        summary_lines += [
            f"--- {ens} ---",
            f"  systems parsed       : {n_sys}",
            f"  rows in long table   : {len(df)} ({len(df)//max(n_sys,1)} per system)",
            f"  subunit residue counts (chain-collapsed): " +
                ", ".join(f"{s}={(df['subunit']==s).sum()//max(n_sys,1)}" for s in ("GAPDH","PRK","CP12","NAD")),
        ]
        for sub_name in ("GAPDH", "PRK", "CP12", "NAD"):
            vals = sub.loc[sub["subunit"] == sub_name, "sum_dTOTAL_all_chains"]
            if len(vals):
                summary_lines.append(
                    f"  Σ ΔTOTAL ({sub_name}, all chains, kcal/mol): "
                    f"median={vals.median():.1f}, IQR={vals.quantile(0.25):.1f}..{vals.quantile(0.75):.1f}"
                )
        summary_lines.append("")

    if frames:
        combined = pd.concat(frames, ignore_index=True)
        combined.to_csv(OUT_DIR / "decomp_combined.csv", index=False)
        sub_combined = pd.concat(sub_frames, ignore_index=True)
        sub_combined.to_csv(OUT_DIR / "decomp_per_subunit_combined.csv", index=False)

    text = "\n".join(summary_lines) + "\n"
    (OUT_DIR / "extraction_summary.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
