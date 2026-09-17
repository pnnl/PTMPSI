#!/usr/bin/env python3
"""
HMM-column-level contact coupling networks via mutual information.

For each (condition, rep, partner):
  - load the trajectory
  - per frame, mark CP12 residues that make ANY heavy-atom contact (5.0 A)
    with the partner — this gives a per-(CP12_chain, CP12_resid) binary signal
  - OR across the 4 CP12 chains -> per-(CP12_resid) signal per chain copy
  - map CP12_resid -> HMM_col via the alignment annotation
  - OR all chain copies into one HMM-column-level binary signal per frame

Result per (condition, partner): a binary matrix B[HMM_col, frame] of shape
(70, n_frames_total = 3 * 1001) where the three replicates are concatenated.

Then for each (condition, partner) compute pairwise mutual information between
all HMM-column signals, and a binary "coupled" edge set thresholded at MI >= tau.

Comparison across conditions:
  - Jaccard similarity on coupled-edge sets between every pair of conditions
  - Frobenius distance on the MI matrices
  - Spectral distance on the eigenvalue spectra

Outputs:
  contact_coupling/per_condition_mi_<partner>.npz       # MI matrices
  contact_coupling/coupling_similarity_summary.csv      # pairwise stats
  contact_coupling/coupling_text_summary.md
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS

import interface_contacts_freq as icf


DEFAULT_SYSTEM_DIRS = {
    "med4_ox":  "/home/scott.widmann/DOE/Perlmutter/MED4_NAD/0003",
    "med4_red": "/home/scott.widmann/DOE/Perlmutter/MED4_NAD/0000",
    "hm2_ox":   "/home/scott.widmann/DOE/Perlmutter/HM2_NAD/0003",
    "hm2_red":  "/home/scott.widmann/DOE/Perlmutter/HM2_NAD/0000",
}


def host_alignment(aln_csv: Path) -> tuple[dict[int, int], dict[int, int]]:
    """med4_pos -> hmm_col, hm2_pos -> hmm_col."""
    df = pd.read_csv(aln_csv)
    med4 = {int(p): int(c)
            for p, c in zip(df["med4_pos"], df["aln_col_1based"])
            if pd.notna(p)}
    hm2  = {int(p): int(c)
            for p, c in zip(df["hm2_pos"], df["aln_col_1based"])
            if pd.notna(p)}
    return med4, hm2


def replicate_hmm_signal(rep_dir: Path, partner: str, cutoff: float,
                         min_run: int, host_aln: dict[int, int],
                         n_hmm_cols: int) -> np.ndarray | None:
    """Return B[hmm_col, frame] binary array of shape (n_hmm_cols, n_frames)."""
    tpr = rep_dir / "solute.tpr"
    xtc = rep_dir / "solute_fit.xtc"
    pdb = rep_dir / "solute_fixed.pdb"
    if not (tpr.exists() and xtc.exists() and pdb.exists()):
        return None
    u = icf.load_universe_with_chains(pdb, tpr, xtc)
    cp12 = icf.select_cp12(u)
    p_atoms = icf.select_partner(u, partner, include_nad=False)
    n_frames = u.trajectory.n_frames

    cp_atom2res, cp_resid, _, cp_chain, _ = icf.residues_from_atoms(cp12)
    n_cp = len(cp_resid)

    # CP12_residue -> HMM column (or None if no map / gap)
    res_to_hmm = np.full(n_cp, -1, dtype=int)
    for i, r in enumerate(cp_resid):
        h = host_aln.get(int(r))
        if h is not None:
            res_to_hmm[i] = int(h) - 1  # 0-based

    B = np.zeros((n_hmm_cols, n_frames), dtype=bool)
    for f, ts in enumerate(u.trajectory):
        if f >= n_frames:
            break
        box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None
        ns = FastNS(cutoff, p_atoms.positions, box=box)
        pairs = ns.search(cp12.positions).get_pairs()
        if pairs.size == 0:
            continue
        cp_res_idx = np.unique(cp_atom2res[pairs[:, 0]])
        hmm_idx = res_to_hmm[cp_res_idx]
        hmm_idx = hmm_idx[hmm_idx >= 0]
        if hmm_idx.size:
            B[hmm_idx, f] = True

    # min-run filter per row (HMM col)
    if min_run > 1:
        for k in range(B.shape[0]):
            if B[k].any():
                B[k] = _min_run_filter(B[k], min_run)
    return B


def _min_run_filter(b: np.ndarray, min_run: int) -> np.ndarray:
    out = b.copy()
    # open
    out = _flip_short(out, target=True, min_run=min_run)
    # close
    out = _flip_short(out, target=False, min_run=min_run)
    return out


def _flip_short(b: np.ndarray, target: bool, min_run: int) -> np.ndarray:
    out = b.copy()
    n = out.size
    diff = np.diff(out.astype(np.int8), prepend=out[0])
    starts = np.r_[0, np.where(diff != 0)[0]]
    ends = np.r_[starts[1:], n]
    lengths = ends - starts
    states = out[starts]
    for i in range(1, len(starts) - 1):
        if states[i] == target and lengths[i] < min_run:
            out[starts[i]:ends[i]] = not target
    return out


def mutual_information_matrix(B: np.ndarray) -> np.ndarray:
    """B is (K, T) binary. Returns (K, K) MI matrix in bits."""
    K, T = B.shape
    p = B.mean(axis=1)  # P(X=1)
    H_individual = np.zeros(K)
    for k in range(K):
        px = p[k]
        if 0.0 < px < 1.0:
            H_individual[k] = -(px * np.log2(px) + (1 - px) * np.log2(1 - px))
    # joint via einsum on int counts
    Bi = B.astype(np.int32)
    n11 = Bi @ Bi.T
    n10 = Bi @ (1 - Bi).T
    n01 = (1 - Bi) @ Bi.T
    n00 = (1 - Bi) @ (1 - Bi).T
    total = T
    p11 = n11 / total
    p10 = n10 / total
    p01 = n01 / total
    p00 = n00 / total
    with np.errstate(divide="ignore", invalid="ignore"):
        # P(X)*P(Y)
        pxy_indep = np.outer(p, p)
        # MI = sum p(x,y) log p(x,y) / (p(x)p(y))
        terms = []
        for pj, pind in [(p11, np.outer(p, p)),
                         (p10, np.outer(p, 1 - p)),
                         (p01, np.outer(1 - p, p)),
                         (p00, np.outer(1 - p, 1 - p))]:
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where((pj > 0) & (pind > 0), pj / pind, 1.0)
                t = pj * np.log2(ratio)
                terms.append(np.where((pj > 0) & np.isfinite(t), t, 0.0))
        mi = sum(terms)
    np.fill_diagonal(mi, 0.0)
    return mi


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / max(len(a | b), 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="contact_coupling")
    ap.add_argument("--alignment_csv", default="CP12_MED4_vs_HM2_alignment_annotation.csv")
    ap.add_argument("--conditions", nargs="+",
                    default=["med4_ox", "med4_red", "hm2_ox", "hm2_red"])
    ap.add_argument("--reps", nargs="+", default=["rep1", "rep2", "rep3"])
    ap.add_argument("--partners", nargs="+", default=["prk", "gap2"])
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--min_run", type=int, default=3)
    ap.add_argument("--mi_threshold", type=float, default=0.05,
                    help="MI (bits) above which to count an edge as 'coupled'.")
    args = ap.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    med4_aln, hm2_aln = host_alignment(Path(args.alignment_csv).resolve())
    n_hmm_cols = 70  # by design

    # extract per (condition, partner) HMM-col signals across reps
    signals: dict[tuple[str, str], np.ndarray] = {}
    for cond in args.conditions:
        host_aln = med4_aln if cond.startswith("med4") else hm2_aln
        base = Path(DEFAULT_SYSTEM_DIRS[cond])
        for partner in args.partners:
            rep_arrays = []
            for rep in args.reps:
                print(f"[mi] extracting {cond}/{rep}/{partner}", file=sys.stderr, flush=True)
                B = replicate_hmm_signal(base / rep, partner,
                                         cutoff=args.cutoff,
                                         min_run=args.min_run,
                                         host_aln=host_aln,
                                         n_hmm_cols=n_hmm_cols)
                if B is None:
                    print(f"    skipping (missing files)", file=sys.stderr)
                    continue
                rep_arrays.append(B)
            if not rep_arrays:
                continue
            stacked = np.concatenate(rep_arrays, axis=1)  # (70, 3 * n_frames)
            signals[(cond, partner)] = stacked
            print(f"    stacked shape: {stacked.shape}, fraction-active: "
                  f"{stacked.mean():.3f}", file=sys.stderr)

    # MI matrices
    mi_matrices = {key: mutual_information_matrix(B) for key, B in signals.items()}

    # save numpy bundle
    save = {}
    for (cond, partner), mi in mi_matrices.items():
        save[f"mi_{cond}_{partner}"] = mi
    for (cond, partner), B in signals.items():
        save[f"occupancy_{cond}_{partner}"] = B.mean(axis=1)
    np.savez_compressed(outdir / "contact_coupling_mi.npz", **save)
    print(f"wrote {outdir / 'contact_coupling_mi.npz'}", file=sys.stderr)

    # pairwise comparisons
    rows = []
    keys = list(mi_matrices.keys())
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            (c1, p1) = keys[i]
            (c2, p2) = keys[j]
            if p1 != p2:
                continue
            M1 = mi_matrices[keys[i]]
            M2 = mi_matrices[keys[j]]
            # edge sets at threshold
            e1 = set((a, b) for a in range(M1.shape[0]) for b in range(a + 1, M1.shape[0])
                     if M1[a, b] >= args.mi_threshold)
            e2 = set((a, b) for a in range(M2.shape[0]) for b in range(a + 1, M2.shape[0])
                     if M2[a, b] >= args.mi_threshold)
            jac = jaccard(e1, e2)
            frob = float(np.linalg.norm(M1 - M2))
            ev1 = np.sort(np.linalg.eigvalsh(M1))[::-1]
            ev2 = np.sort(np.linalg.eigvalsh(M2))[::-1]
            spec = float(np.linalg.norm(ev1 - ev2))
            rows.append({
                "partner": p1,
                "cond_a": c1, "cond_b": c2,
                "n_edges_a": len(e1), "n_edges_b": len(e2),
                "n_edges_shared": len(e1 & e2),
                "jaccard": jac,
                "frobenius_mi_diff": frob,
                "spectral_distance": spec,
            })
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "coupling_similarity_summary.csv", index=False)

    # text summary
    lines = ["# Contact-coupling network similarity (MI-thresholded HMM-column graphs)",
             "",
             f"- MI threshold for edge inclusion: {args.mi_threshold} bits",
             ""]
    for partner in args.partners:
        d = df[df["partner"] == partner]
        if d.empty:
            continue
        lines.append(f"## {partner.upper()} coupling graph similarity (Jaccard on edges)")
        lines.append("")
        lines.append("| condition A | condition B | edges A | edges B | shared | Jaccard | Frobenius | Spectral |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            lines.append(f"| {r['cond_a']} | {r['cond_b']} | {r['n_edges_a']} | {r['n_edges_b']} | "
                         f"{r['n_edges_shared']} | {r['jaccard']:.3f} | {r['frobenius_mi_diff']:.3f} | "
                         f"{r['spectral_distance']:.3f} |")
        lines.append("")
    (outdir / "coupling_text_summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {outdir / 'coupling_text_summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
