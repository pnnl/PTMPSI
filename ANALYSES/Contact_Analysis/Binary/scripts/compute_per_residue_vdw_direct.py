#!/usr/bin/env python3
"""
Direct per-residue VDW computation between CP12 and PRK (and CP12 ↔ GAP2 for
validation), reading only the existing trajectories and AMBER topology — no
gmx_MMPBSA run required.

For each (condition, rep, CP12_chain, CP12_resid):
  - per frame, find all CP12-atom × partner-atom pairs within 10 Å
  - evaluate L-J(σ_ij, ε_ij) with Lorentz-Berthelot mixing
  - accumulate to the CP12 residue total

σ/ε per atom type are read from `mmgbsa_cp12_prk/COM.prmtop` (AMBER topology
gmx_MMPBSA was built against, same amber99sb.ff). Atom types in the full system
trajectory (`solute.tpr` + `solute_fit.xtc`) are looked up by name. CP12 and
PRK/GAP2 chains are selected by chain ID via interface_contacts_freq's
load_universe_with_chains.

Output: per_residue_vdw_direct_{partner}.csv
  columns: condition, rep, chain, resid, resname, vdw_mean_kcal_mol, n_frames

For validation: the GAP2 partner output should rank-correlate strongly with the
existing gmx_MMPBSA idecomp=2 DELTAS VDW from `mmgbsa_cp12_gap2_nad`.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import parmed
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS

import interface_contacts_freq as icf


DEFAULT_SYSTEM_DIRS = {
    "med4_ox":  "/home/scott.widmann/DOE/Perlmutter/MED4_NAD/0003",
    "med4_red": "/home/scott.widmann/DOE/Perlmutter/MED4_NAD/0000",
    "hm2_ox":   "/home/scott.widmann/DOE/Perlmutter/HM2_NAD/0003",
    "hm2_red":  "/home/scott.widmann/DOE/Perlmutter/HM2_NAD/0000",
}


def build_type_lj_table(prmtop_path: Path) -> dict[str, tuple[float, float]]:
    """Return {atom_type: (sigma_Å, epsilon_kcal/mol)} from an AMBER prmtop."""
    top = parmed.load_file(str(prmtop_path))
    table: dict[str, tuple[float, float]] = {}
    for a in top.atoms:
        if a.type not in table:
            table[a.type] = (float(a.sigma), float(a.epsilon))
    return table


def get_lj_arrays(group: mda.core.groups.AtomGroup,
                  type_lj: dict[str, tuple[float, float]]) -> tuple[np.ndarray, np.ndarray]:
    """Map atomgroup types to (sigma, epsilon) arrays. Raise if any type is unknown."""
    types = list(group.types)
    missing = sorted({t for t in types if t not in type_lj})
    if missing:
        raise SystemExit(f"types missing from prmtop: {missing}")
    sig = np.array([type_lj[t][0] for t in types], dtype=np.float64)
    eps = np.array([type_lj[t][1] for t in types], dtype=np.float64)
    return sig, eps


def atom_to_residue_index(group: mda.core.groups.AtomGroup) -> tuple[np.ndarray, list[str], list[int], list[str]]:
    """Map each atom in group to a contiguous residue index within the group.

    Returns:
        atom2res — int array of length n_atoms, values in [0, n_residues)
        chain    — list of length n_residues with chain (segid) per residue
        resid    — list of length n_residues with resid per residue
        resname  — list of length n_residues with resname per residue
    """
    chains: list[str] = []
    resids: list[int] = []
    resnames: list[str] = []
    atom2res = np.empty(group.n_atoms, dtype=np.int64)
    res_key_to_idx: dict[tuple, int] = {}
    for i, atom in enumerate(group.atoms):
        r = atom.residue
        key = (str(r.segid), int(r.resid))
        idx = res_key_to_idx.get(key)
        if idx is None:
            idx = len(chains)
            res_key_to_idx[key] = idx
            chains.append(str(r.segid))
            resids.append(int(r.resid))
            resnames.append(str(r.resname))
        atom2res[i] = idx
    return atom2res, chains, resids, resnames


def compute_replicate(rep_dir: Path, partner: str,
                      type_lj: dict[str, tuple[float, float]],
                      cutoff_A: float,
                      start_frame: int = 0,
                      stop_frame: int | None = None,
                      stride: int = 1) -> pd.DataFrame:
    """Compute per-CP12-residue VDW interaction with the chosen partner across the trajectory."""
    tpr = rep_dir / "solute.tpr"
    xtc = rep_dir / "solute_fit.xtc"
    pdb = rep_dir / "solute_fixed.pdb"
    for p in (tpr, xtc, pdb):
        if not p.exists():
            raise SystemExit(f"missing input: {p}")

    u = icf.load_universe_with_chains(pdb, tpr, xtc)
    cp12 = icf.select_cp12(u)  # heavy + hydrogen? icf returns no-H atoms only
    # icf.select_cp12 returns heavy atoms only (uses "not name H*"). For full L-J we
    # want all atoms. Override the selection.
    cp12 = u.select_atoms("protein and (segid " + " ".join(icf.CP12_CHAINS) + ")")
    p_atoms = u.select_atoms(
        "protein and (segid " + " ".join(
            icf.PRK_CHAINS if partner == "prk" else icf.GAP2_CHAINS
        ) + ")"
    )
    if cp12.n_atoms == 0 or p_atoms.n_atoms == 0:
        raise SystemExit(f"empty selection in {rep_dir}")

    sig_cp, eps_cp = get_lj_arrays(cp12, type_lj)
    sig_pt, eps_pt = get_lj_arrays(p_atoms, type_lj)
    cp_atom2res, cp_chain, cp_resid, cp_resname = atom_to_residue_index(cp12)
    n_cp_res = len(cp_chain)

    per_res_sum = np.zeros(n_cp_res, dtype=np.float64)
    n_frames = 0
    t0 = time.time()
    if stop_frame is None:
        stop_frame = u.trajectory.n_frames
    sliced = u.trajectory[start_frame:stop_frame:stride]
    for ts in sliced:
        n_frames += 1
        box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None
        # neighbor search at 10 Å
        ns = FastNS(cutoff_A, p_atoms.positions, box=box)
        pairs = ns.search(cp12.positions).get_pairs()
        if pairs.size == 0:
            continue
        cp_idx = pairs[:, 0]
        pt_idx = pairs[:, 1]
        sig_ij = 0.5 * (sig_cp[cp_idx] + sig_pt[pt_idx])
        eps_ij = np.sqrt(eps_cp[cp_idx] * eps_pt[pt_idx])
        # distances (apply PBC if box present)
        diff = cp12.positions[cp_idx] - p_atoms.positions[pt_idx]
        if box is not None:
            box3 = box[:3].astype(np.float64)
            diff -= np.round(diff / box3) * box3
        r2 = (diff * diff).sum(axis=1)
        # guard against r=0 (shouldn't happen between different atoms)
        r2 = np.where(r2 < 1e-8, 1e-8, r2)
        # L-J energy
        s2_over_r2 = (sig_ij * sig_ij) / r2
        s6 = s2_over_r2 ** 3
        s12 = s6 * s6
        e_lj = 4.0 * eps_ij * (s12 - s6)
        # accumulate per CP12 residue
        np.add.at(per_res_sum, cp_atom2res[cp_idx], e_lj)

    elapsed = time.time() - t0
    print(f"    {rep_dir.name}/{partner}: {n_frames} frames, {elapsed:.1f}s",
          file=sys.stderr, flush=True)

    per_res_mean = per_res_sum / max(n_frames, 1)
    return pd.DataFrame({
        "chain": cp_chain,
        "resid": cp_resid,
        "resname": cp_resname,
        "vdw_mean_kcal_mol": per_res_mean,
        "n_frames": n_frames,
    })


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="contact_kinetic_thermo_bridge")
    ap.add_argument("--partners", nargs="+", default=["prk", "gap2"],
                    help="Which partners to compute. 'gap2' is computed for validation.")
    ap.add_argument("--reps", nargs="+", default=["rep1", "rep2", "rep3"])
    ap.add_argument("--cutoff", type=float, default=10.0,
                    help="Pair distance cutoff in Å for L-J neighbor list.")
    ap.add_argument("--prmtop",
                    default="/home/scott.widmann/DOE/Perlmutter/MED4_NAD/0003/rep1/mmgbsa_cp12_prk/COM.prmtop",
                    help="Source of atom-type → (σ, ε) lookup table.")
    ap.add_argument("--start", type=int, default=200,
                    help="First frame index to include (0-based; default 200, matches gmx_MMPBSA's startframe=201).")
    ap.add_argument("--stop", type=int, default=None,
                    help="Stop frame index (exclusive). Default = full trajectory.")
    ap.add_argument("--stride", type=int, default=4,
                    help="Frame stride (default 4, matches gmx_MMPBSA's interval=4 → 200 frames).")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    type_lj = build_type_lj_table(Path(args.prmtop))
    print(f"[lj] loaded {len(type_lj)} atom types from {args.prmtop}", file=sys.stderr)

    for partner in args.partners:
        rows = []
        for cond, base in DEFAULT_SYSTEM_DIRS.items():
            for rep in args.reps:
                rep_dir = Path(base) / rep
                if not rep_dir.exists():
                    continue
                df = compute_replicate(rep_dir, partner, type_lj, args.cutoff,
                                        start_frame=args.start, stop_frame=args.stop,
                                        stride=args.stride)
                df["condition"] = cond
                df["rep"] = rep
                df["partner"] = partner
                rows.append(df)
        if not rows:
            continue
        out = pd.concat(rows, ignore_index=True)
        out_path = outdir / f"per_residue_vdw_direct_{partner}.csv"
        out.to_csv(out_path, index=False)
        print(f"[ok] wrote {out_path}  ({len(out)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
