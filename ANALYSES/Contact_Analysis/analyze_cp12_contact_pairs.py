#!/usr/bin/env python3
"""
Differential contact pairs (rewiring map) analysis.

Computes (CP12_residue, partner_residue) pairwise contact frequencies from
MD trajectories, aggregates across replicates, collapses across CP12 copies
via the alignment coordinate system, and produces differential tables comparing
conditions (host × redox).

Reuses infrastructure from interface_contacts_freq.py for trajectory loading,
chain assignments, and atom selection.

Outputs:
  cp12_contact_pairs/<condition>/rep{N}/contact_pairs_{partner}.csv
  cp12_contact_pairs/<condition>/<condition>_contact_pairs_agg_{partner}.csv
  cp12_contact_pairs/<condition>/<condition>_contact_pairs_collapsed_{partner}.csv
  cp12_contact_pairs/differential_pairs_{partner}.csv
  cp12_contact_pairs/top_rewired_pairs_{partner}.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from collections import Counter

# Ensure the directory containing this script is on sys.path so that
# interface_contacts_freq.py can be found regardless of where Python is invoked from.
_SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in dir() else Path("/Path/to/user/root")
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import numpy as np
import pandas as pd

import interface_contacts_freq as icf
from MDAnalysis.lib.nsgrid import FastNS


SYSTEMS = ["med4_ox", "med4_red", "hm2_ox", "hm2_red"]
SYSTEM_DIRS = {
    "med4_ox": Path("/path/to/user/MED4_NAD/0003"),
    "med4_red": Path("/path/to/user/MED4_NAD/0000"),
    "hm2_ox": Path("/path/to/user/HM2_NAD/0003"),
    "hm2_red": Path("/path/to/user/HM2_NAD/0000"),
}
REPS = ["rep1", "rep2", "rep3"]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def contact_freq_pairs(u, A, B, cutoff, stride):
    """
    Compute per-(CP12_residue, partner_residue) contact frequency.

    A contact is counted when ANY heavy atom pair between two residues
    falls within cutoff in a given frame. Each residue pair is counted
    at most once per frame (binary contact per frame).
    """
    A_atom2res, A_resid, A_resname, A_chain, _ = icf.residues_from_atoms(A)
    B_atom2res, B_resid, B_resname, B_chain, _ = icf.residues_from_atoms(B)

    pair_counts = Counter()
    n_frames = 0

    for ts in u.trajectory[::stride]:
        n_frames += 1
        box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None

        ns = FastNS(cutoff, B.positions, box=box)
        pairs = ns.search(A.positions).get_pairs()
        if pairs.size == 0:
            continue

        res_pairs = np.column_stack([A_atom2res[pairs[:, 0]], B_atom2res[pairs[:, 1]]])
        unique_pairs = np.unique(res_pairs, axis=0)
        for ai, bi in unique_pairs:
            pair_counts[(int(ai), int(bi))] += 1

    rows = []
    for (ai, bi), count in pair_counts.items():
        rows.append({
            "cp12_chain": str(A_chain[ai]),
            "cp12_resid": int(A_resid[ai]),
            "cp12_resname": str(A_resname[ai]),
            "partner_chain": str(B_chain[bi]),
            "partner_resid": int(B_resid[bi]),
            "partner_resname": str(B_resname[bi]),
            "n_frames_contact": count,
            "n_frames_total": n_frames,
            "freq": count / max(n_frames, 1),
        })
    return pd.DataFrame(rows), n_frames


def load_alignment(alignment_csv):
    """Load alignment CSV and build resid -> alignment position maps."""
    aln = pd.read_csv(alignment_csv)
    med4_map = {}
    hm2_map = {}
    for _, row in aln.iterrows():
        aln_pos = int(row["aln_col_1based"])
        if pd.notna(row["med4_pos"]):
            med4_map[int(row["med4_pos"])] = aln_pos
        if pd.notna(row["hm2_pos"]):
            hm2_map[int(row["hm2_pos"])] = aln_pos
    return aln, med4_map, hm2_map


def aggregate_pairs(dfs):
    """Mean + SEM across replicates for pairwise contact frequencies."""
    all_df = pd.concat(dfs, ignore_index=True)
    key_cols = ["cp12_chain", "cp12_resid", "cp12_resname",
                "partner_chain", "partner_resid", "partner_resname"]
    grp = all_df.groupby(key_cols, as_index=False)
    return grp.agg(
        mean_freq=("freq", "mean"),
        sem_freq=("freq", lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0),
        n_reps=("freq", "count"),
    )


def collapse_pairs(agg_df, resid_to_aln, partner_chains):
    """
    Collapse CP12 copies (C/H/M/P) to alignment positions and partner copies
    to canonical residue numbers. Takes max freq across copies at each
    (aln_pos, partner_resid) pair.
    """
    work = agg_df.copy()
    work["cp12_aln_pos"] = work["cp12_resid"].map(resid_to_aln)
    work = work.dropna(subset=["cp12_aln_pos"])
    work["cp12_aln_pos"] = work["cp12_aln_pos"].astype(int)

    grp = work.groupby(["cp12_aln_pos", "cp12_resname", "partner_resid", "partner_resname"],
                        as_index=False)
    collapsed = grp.agg(max_freq=("mean_freq", "max"))
    return collapsed


def compute_differential(collapsed_dfs, systems):
    """
    Join collapsed tables from all conditions and compute all pairwise deltas.
    Delta columns are generated dynamically from whatever systems are provided.
    """
    merged = None
    for sys_name in systems:
        df = collapsed_dfs.get(sys_name)
        if df is None:
            continue
        sub = df[["cp12_aln_pos", "partner_resid", "partner_resname",
                   "cp12_resname", "max_freq"]].copy()
        sub = sub.rename(columns={"max_freq": f"{sys_name}_freq"})
        if merged is None:
            merged = sub
        else:
            merged = merged.merge(sub, on=["cp12_aln_pos", "partner_resid",
                                           "partner_resname", "cp12_resname"],
                                  how="outer")

    if merged is None:
        return pd.DataFrame()

    freq_cols = [f"{s}_freq" for s in systems]
    for col in freq_cols:
        if col not in merged.columns:
            merged[col] = 0.0
        merged[col] = merged[col].fillna(0.0)

    # Compute all pairwise deltas between systems
    for i, s1 in enumerate(systems):
        for s2 in systems[i + 1:]:
            col_name = f"delta_{s1}_vs_{s2}"
            merged[col_name] = merged[f"{s1}_freq"] - merged[f"{s2}_freq"]

    return merged.sort_values("cp12_aln_pos")


def run_condition(system, rep_dir, rep_name, args):
    """Run contact pair analysis for one replicate."""
    tpr = rep_dir / "solute.tpr"
    xtc = rep_dir / "solute_fit.xtc"
    pdb = rep_dir / "solute_fixed.pdb"
    for p in (tpr, xtc, pdb):
        if not p.exists():
            raise FileNotFoundError(f"Missing: {p}")

    u = icf.load_universe_with_chains(pdb, tpr, xtc)
    cp12 = icf.select_cp12(u)

    results = {}
    for partner in ["prk", "gap2"]:
        partner_atoms = icf.select_partner(u, partner, include_nad=False)
        df, n_frames = contact_freq_pairs(u, cp12, partner_atoms, args.cutoff, args.stride)
        df.insert(0, "strain", system)
        df.insert(1, "replicate", rep_name)
        df.insert(2, "partner", partner)
        results[partner] = df
        print(f"  [{system}/{rep_name}] {partner}: {len(df)} pairs from {n_frames} frames")
        u.trajectory.rewind()

    return results


def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute CP12↔partner residue-pair contact frequencies and differential rewiring maps."
    )
    ap.add_argument("--alignment_csv", default="CP12_MED4_vs_HM2_alignment_annotation.csv",
                    help="CP12 MED4↔HM2 alignment annotation CSV.")
    ap.add_argument("--outdir", default="cp12_contact_pairs",
                    help="Output directory.")
    ap.add_argument("--cutoff", type=float, default=6.0,
                    help="Heavy-atom contact cutoff in angstroms.")
    ap.add_argument("--stride", type=int, default=10,
                    help="Frame stride (subsample).")
    ap.add_argument("--peri_thresh", type=float, default=0.10,
                    help="Minimum freq in any condition to retain a pair.")
    ap.add_argument("--strain_dir", action="append", default=[],
                    help="Override system dir: e.g. med4_ox=/path/to/dir")
    ap.add_argument("--reps", nargs="+", default=REPS,
                    help="Replicate directory names.")
    ap.add_argument("--systems", nargs="+", default=SYSTEMS,
                    help="System labels to process.")
    return ap.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse any overrides
    dir_overrides = icf.parse_mapping(args.strain_dir, "--strain_dir")
    system_dirs = {k: Path(v) for k, v in dir_overrides.items()}
    for s in args.systems:
        if s not in system_dirs:
            system_dirs[s] = SYSTEM_DIRS[s]

    # Load alignment
    aln, med4_aln_map, hm2_aln_map = load_alignment(args.alignment_csv)
    print(f"Alignment: {len(med4_aln_map)} MED4 positions, {len(hm2_aln_map)} HM2 positions")

    # Collect per-rep results
    per_system = {s: {"prk": [], "gap2": []} for s in args.systems}

    for system in args.systems:
        base_dir = system_dirs[system]
        print(f"\nProcessing {system} from {base_dir}")
        for rep in args.reps:
            rep_dir = base_dir / rep
            if not rep_dir.exists():
                eprint(f"  Skipping missing: {rep_dir}")
                continue
            results = run_condition(system, rep_dir, rep, args)
            for partner, df in results.items():
                # Save per-rep
                rep_out = outdir / system / rep
                rep_out.mkdir(parents=True, exist_ok=True)
                df.to_csv(rep_out / f"contact_pairs_{partner}.csv", index=False)
                per_system[system][partner].append(df)

    # Aggregate and collapse per system
    collapsed_by_system = {s: {} for s in args.systems}

    for system in args.systems:
        resid_to_aln = med4_aln_map if system.startswith("med4") else hm2_aln_map
        sys_out = outdir / system
        sys_out.mkdir(parents=True, exist_ok=True)

        for partner in ["prk", "gap2"]:
            dfs = per_system[system][partner]
            if not dfs:
                continue

            agg = aggregate_pairs(dfs)
            agg.to_csv(sys_out / f"{system}_contact_pairs_agg_{partner}.csv", index=False)

            partner_chains = icf.PRK_CHAINS if partner == "prk" else icf.GAP2_CHAINS
            collapsed = collapse_pairs(agg, resid_to_aln, partner_chains)
            collapsed.to_csv(sys_out / f"{system}_contact_pairs_collapsed_{partner}.csv", index=False)
            collapsed_by_system[system][partner] = collapsed

            n_sig = int((agg["mean_freq"] >= args.peri_thresh).sum())
            print(f"  [{system}] {partner}: {len(agg)} unique pairs, {n_sig} with freq >= {args.peri_thresh}")

    # Differential analysis
    print("\nComputing differential contact maps...")
    for partner in ["prk", "gap2"]:
        coll = {s: collapsed_by_system[s].get(partner) for s in args.systems}
        diff = compute_differential(coll, args.systems)
        if diff.empty:
            eprint(f"  No differential data for {partner}")
            continue

        diff.to_csv(outdir / f"differential_pairs_{partner}.csv", index=False)

        # Filter to significant pairs (any system above threshold)
        freq_cols = [f"{s}_freq" for s in args.systems if f"{s}_freq" in diff.columns]
        any_sig = (diff[freq_cols] >= args.peri_thresh).any(axis=1)
        diff_sig = diff[any_sig].copy()

        # Top rewired pairs by max |delta| across all pairwise comparisons
        delta_cols = [c for c in diff_sig.columns if c.startswith("delta_")]
        diff_sig["max_abs_delta"] = diff_sig[delta_cols].abs().max(axis=1)
        top = diff_sig.nlargest(30, "max_abs_delta")
        top.to_csv(outdir / f"top_rewired_pairs_{partner}.csv", index=False)

        print(f"  {partner}: {len(diff)} total pairs, {len(diff_sig)} above threshold, top delta = {top['max_abs_delta'].max():.3f}" if len(top) else f"  {partner}: no significant pairs")

    print(f"\nDone. Results in: {outdir}")


if __name__ == "__main__":
    main()
