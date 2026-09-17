#!/usr/bin/env python3
"""
Interface residue contact frequency via MDAnalysis (heavy-atom any-atom-within-cutoff).

Per replicate directory expects:
	solute.tpr
	solute_fit.xtc
	solute_fixed.pdb  (contains chain IDs)

Runs from project root with structure:
	<root>/med4/rep1 ...
	<root>/hm2/rep1 ...

Outputs under: <root>/<outdir>/<strain>/<rep>/

By default computes three partner definitions:
	- receptor: GAP2+PRK (+ optional NAD)
	- gap2: GAP2 only
	- prk: PRK only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS


CP12_CHAINS = ["C", "H", "M", "P"]
PRK_CHAINS  = ["B", "G", "J", "L"]
GAP2_CHAINS = ["A", "D", "E", "F", "I", "K", "N", "O"]
NAD_CHAINS  = ["Q"]


def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


def parse_mapping(entries: list[str], flag_name: str) -> dict[str, str]:
	mapping: dict[str, str] = {}
	for entry in entries:
		if "=" not in entry:
			raise SystemExit(f"ERROR: {flag_name} entries must look like name=/path, got: {entry}")
		key, value = entry.split("=", 1)
		key = key.strip()
		value = value.strip()
		if not key or not value:
			raise SystemExit(f"ERROR: {flag_name} entries must look like name=/path, got: {entry}")
		mapping[key] = value
	return mapping


def _ref_chain_labels(u_ref: mda.Universe) -> np.ndarray:
	"""Return per-atom chain labels from reference PDB (prefer chainIDs; fall back to segids)."""
	labels = None
	try:
		labels = np.array(u_ref.atoms.chainIDs, dtype=object)
		if labels.size and np.any(labels != ""):
			return labels
	except Exception:
		pass
	labels = np.array(u_ref.atoms.segids, dtype=object)
	return labels


def normalize_segids_from_chainids(u: mda.Universe) -> None:
	"""
	Standardize chain labeling by ensuring segids carry chain IDs when available.
	This enables robust selection via 'segid X Y Z'.
	"""
	try:
		chainIDs = np.array(u.atoms.chainIDs, dtype=object)
	except Exception:
		chainIDs = np.array([""] * len(u.atoms), dtype=object)

	segids = np.array(u.atoms.segids, dtype=object)

	# If segids are empty/constant but chainIDs are informative, copy chainIDs -> segids
	if (np.all(segids == segids[0]) or np.all(segids == "")) and np.any(chainIDs != ""):
		u.atoms.segids = chainIDs


def build_chain_map_by_key(u_ref: mda.Universe, u: mda.Universe) -> dict[int, str]:
	"""
	Map atom index in u -> chain label from u_ref using key (resid, resname, name).
	If duplicates exist, require identical ordering and fall back to index-based mapping.
	"""
	if len(u.atoms) != len(u_ref.atoms):
		raise ValueError(f"Atom count mismatch: traj/topology={len(u.atoms)} vs ref_pdb={len(u_ref.atoms)}")

	ref_keys = [(a.resid, a.resname, a.name) for a in u_ref.atoms]
	trg_keys = [(a.resid, a.resname, a.name) for a in u.atoms]

	ref_counts = Counter(ref_keys)
	trg_counts = Counter(trg_keys)

	if ref_counts != trg_counts:
		mism = []
		for k, c in ref_counts.items():
			if trg_counts.get(k, 0) != c:
				mism.append((k, c, trg_counts.get(k, 0)))
				if len(mism) >= 5:
					break
		raise ValueError(
			"Atom key multiset mismatch between solute_fixed.pdb and solute.tpr/xtc.\n"
			f"Example mismatches (key, ref_count, traj_count): {mism}\n"
			"Suggestion: ensure solute_fixed.pdb was generated from the SAME solute.tpr/trajectory and atom order."
		)

	ref_chain = _ref_chain_labels(u_ref)

	# If duplicated keys exist, require identical ordering to be safe.
	if any(v > 1 for v in ref_counts.values()):
		for i, (rk, tk) in enumerate(zip(ref_keys, trg_keys)):
			if rk != tk:
				raise ValueError(
					"Duplicate atom keys exist and atom ordering differs; cannot safely transfer chain IDs.\n"
					f"First mismatch at atom index {i}: ref={rk} traj={tk}\n"
					"Suggestion: regenerate solute_fixed.pdb directly from this trajectory (same atom order)."
				)
		return {i: str(ref_chain[i]) for i in range(len(u.atoms))}

	# Unique keys
	key_to_chain = {ref_keys[i]: str(ref_chain[i]) for i in range(len(u_ref.atoms))}
	return {i: key_to_chain[trg_keys[i]] for i in range(len(u.atoms))}


def load_universe_with_chains(pdb_path: Path, tpr_path: Path, xtc_path: Path) -> mda.Universe:
	"""
	Prefer Universe(pdb, xtc) if it works (fast path).
	Otherwise Universe(tpr, xtc) then transfer chain labels from pdb via atom mapping.
	"""
	# Fast path
	try:
		u = mda.Universe(str(pdb_path), str(xtc_path))
		normalize_segids_from_chainids(u)
		# Require that CP12 chain IDs exist as segids
		if len(u.select_atoms("segid " + " ".join(CP12_CHAINS))) > 0:
			return u
	except Exception:
		pass

	# Fallback: TPR+XTC + mapping
	u = mda.Universe(str(tpr_path), str(xtc_path))
	u_ref = mda.Universe(str(pdb_path))
	chain_map = build_chain_map_by_key(u_ref, u)

	segids = np.array([chain_map[i] for i in range(len(u.atoms))], dtype=object)
	u.atoms.segids = segids
	try:
		u.atoms.chainIDs = segids
	except Exception:
		pass

	# Validate
	if len(u.select_atoms("segid " + " ".join(CP12_CHAINS))) == 0:
		raise ValueError("Chain transfer completed but CP12 segids not found; check mapping and chain labels.")
	return u


def residues_from_atoms(atoms: mda.core.groups.AtomGroup):
	"""
	Return:
		atom_local_residx: per-atom index into local residue table
		resid,resname,chain: per-local-residue arrays
		res_uindex: universe residue indices for local residues (optional keying)
	"""
	res = atoms.residues
	resid = res.resids.astype(int)
	resname = np.array(res.resnames, dtype=object)
	chain = np.array(res.segids, dtype=object)
	res_uindex = np.array(res.resindices, dtype=int)

	uni_to_local = {ur: i for i, ur in enumerate(res_uindex)}
	atom_local_residx = np.array([uni_to_local[ri] for ri in atoms.resindices], dtype=int)
	return atom_local_residx, resid, resname, chain, res_uindex


def contact_freq_any_heavy(
	u: mda.Universe,
	A: mda.core.groups.AtomGroup,
	B: mda.core.groups.AtomGroup,
	cutoff: float,
	stride: int,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
	"""
	Per-residue contact frequency for A side and B side.
	Contact in a frame if ANY heavy atom in residue is within cutoff of partner.
	Uses FastNS neighbor search (periodic if box present).
	"""
	A_atom2res, A_resid, A_resname, A_chain, _ = residues_from_atoms(A)
	B_atom2res, B_resid, B_resname, B_chain, _ = residues_from_atoms(B)

	A_counts = np.zeros(len(A_resid), dtype=int)
	B_counts = np.zeros(len(B_resid), dtype=int)

	n_frames = 0
	for ts in u.trajectory[::stride]:
		n_frames += 1
		box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None

		ns = FastNS(cutoff, B.positions, box=box)
		pairs = ns.search(A.positions).get_pairs()
		if pairs.size == 0:
			continue

		A_hit = np.unique(A_atom2res[pairs[:, 0]])
		B_hit = np.unique(B_atom2res[pairs[:, 1]])
		A_counts[A_hit] += 1
		B_counts[B_hit] += 1

	A_df = pd.DataFrame({
		"chain": A_chain,
		"resid": A_resid,
		"resname": A_resname,
		"n_frames_contact": A_counts,
		"n_frames_total": n_frames,
		"freq": A_counts / max(n_frames, 1),
	})
	B_df = pd.DataFrame({
		"chain": B_chain,
		"resid": B_resid,
		"resname": B_resname,
		"n_frames_contact": B_counts,
		"n_frames_total": n_frames,
		"freq": B_counts / max(n_frames, 1),
	})
	return A_df, B_df, n_frames


def classify_core_peripheral(df: pd.DataFrame, core_thresh: float, peri_thresh: float,
								core_quantile: float | None) -> pd.DataFrame:
	"""
	Add 'class': core | peripheral | non-interface.
	If core_quantile is set, core = top (core_quantile) fraction among residues with freq>=peri_thresh.
	Peripheral = freq>=peri_thresh but not core.
	"""
	out = df.copy()
	if core_quantile is not None:
		eligible = out["freq"] >= peri_thresh
		if eligible.any():
			q = out.loc[eligible, "freq"].quantile(1.0 - core_quantile)
			out["class"] = np.where(out["freq"] >= q, "core",
									np.where(out["freq"] >= peri_thresh, "peripheral", "non-interface"))
		else:
			out["class"] = "non-interface"
		return out

	out["class"] = np.where(out["freq"] >= core_thresh, "core",
							np.where(out["freq"] >= peri_thresh, "peripheral", "non-interface"))
	return out


def cp12_union_consensus(cp12_df: pd.DataFrame, peri_thresh: float, consensus_n: int) -> tuple[pd.DataFrame, pd.DataFrame]:
	"""
	CP12-side union/consensus across copies (C/H/M/P) by (resid,resname).
	Inclusion uses peri_thresh (interface membership), independent of core definition.
	"""
	work = cp12_df[["chain", "resid", "resname", "freq"]].copy()
	work["meets"] = work["freq"] >= peri_thresh

	grp = work.groupby(["resid", "resname"], as_index=False)
	agg = grp.agg(
		freq_max=("freq", "max"),
		freq_mean=("freq", "mean"),
		n_copies_meeting_threshold=("meets", "sum"),
		copies_meeting_threshold=("chain", lambda s: ",".join(sorted(set(s[work.loc[s.index, "meets"]])))),
	)
	union = agg[agg["n_copies_meeting_threshold"] >= 1].copy()
	consensus = agg[agg["n_copies_meeting_threshold"] >= consensus_n].copy()
	return union, consensus


def receptor_union_consensus_by_cp12_copy_singlepass(
	u: mda.Universe,
	partner_atoms: mda.core.groups.AtomGroup,
	cutoff: float,
	stride: int,
	peri_thresh: float,
	consensus_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
	"""
	Receptor-side union/consensus across CP12 copies computed in ONE trajectory pass.
	For each CP12 chain separately, compute partner-residue contact frequency (partner residues that contact that CP12 chain).
	Then summarize per partner residue:
		- freq_max, freq_mean over CP12 chains
		- n_copies_meeting_threshold and list of CP12 chains meeting peri_thresh
	"""
	# Partner residue metadata
	P_atom2res, P_resid, P_resname, P_chain, _ = residues_from_atoms(partner_atoms)
	nP = len(P_resid)

	# CP12 chain atom groups (heavy atoms)
	cp_chain_groups = {}
	for c in CP12_CHAINS:
		g = u.select_atoms(f"protein and segid {c} and not name H*")
		if len(g) > 0:
			cp_chain_groups[c] = g

	chains = [c for c in CP12_CHAINS if c in cp_chain_groups]
	nC = len(chains)
	counts = np.zeros((nC, nP), dtype=int)

	n_frames = 0
	for ts in u.trajectory[::stride]:
		n_frames += 1
		box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None

		ns = FastNS(cutoff, partner_atoms.positions, box=box)
		for ci, c in enumerate(chains):
			pairs = ns.search(cp_chain_groups[c].positions).get_pairs()
			if pairs.size == 0:
				continue
			# pairs[:,1] indexes partner atoms
			hit = np.unique(P_atom2res[pairs[:, 1]])
			counts[ci, hit] += 1

	freqs = counts / max(n_frames, 1)  # shape (nC, nP)

	rows = []
	for pi in range(nP):
		vals = freqs[:, pi] if nC else np.array([0.0])
		meets_chains = [chains[ci] for ci in range(nC) if vals[ci] >= peri_thresh]
		rows.append({
			"chain": str(P_chain[pi]),
			"resid": int(P_resid[pi]),
			"resname": str(P_resname[pi]),
			"freq_max": float(np.max(vals)),
			"freq_mean": float(np.mean(vals)),
			"n_copies_meeting_threshold": int(len(meets_chains)),
			"copies_meeting_threshold": ",".join(meets_chains),
		})

	agg = pd.DataFrame(rows)
	union = agg[agg["n_copies_meeting_threshold"] >= 1].copy()
	consensus = agg[agg["n_copies_meeting_threshold"] >= consensus_n].copy()
	return union, consensus


def aggregate_resfreq_across_reps(dfs: list[pd.DataFrame], key_cols: list[str]) -> pd.DataFrame:
	"""Mean+SEM across replicates for per-residue frequencies."""
	all_df = pd.concat(dfs, ignore_index=True)
	grp = all_df.groupby(key_cols, as_index=False)
	return grp.agg(
		mean_freq=("freq", "mean"),
		sem_freq=("freq", lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0),
		n_reps=("freq", "count"),
	)


def aggregate_union_consensus_across_reps(df: pd.DataFrame, key_cols: list[str]) -> pd.DataFrame:
	"""
	Aggregate union/consensus tables across replicates.
	Computes mean+SEM for freq_mean and freq_max, plus mean n_copies_meeting_threshold.
	"""
	grp = df.groupby(key_cols, as_index=False)
	return grp.agg(
		mean_freq_mean=("freq_mean", "mean"),
		sem_freq_mean=("freq_mean", lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0),
		mean_freq_max=("freq_max", "mean"),
		sem_freq_max=("freq_max", lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0),
		mean_n_copies=("n_copies_meeting_threshold", "mean"),
		n_reps=("freq_mean", "count"),
	)


def top10(df: pd.DataFrame) -> pd.DataFrame:
	return df.sort_values("freq", ascending=False).head(10)[["chain", "resid", "resname", "freq"]]


def select_partner(u: mda.Universe, partner: str, include_nad: bool) -> mda.core.groups.AtomGroup:
	"""
	partner in {"receptor","gap2","prk"}.
	Uses segid as chain label (ensured by mapping/normalization).
	Heavy atoms only. NAD (chain Q) is not protein; include explicitly if requested for "receptor".
	"""
	if partner == "gap2":
		sel = "protein and (segid " + " ".join(GAP2_CHAINS) + ") and not name H*"
	elif partner == "prk":
		sel = "protein and (segid " + " ".join(PRK_CHAINS) + ") and not name H*"
	elif partner == "receptor":
		prot_sel = "(protein and (segid " + " ".join(GAP2_CHAINS + PRK_CHAINS) + ") and not name H*)"
		if include_nad:
			# NAD may not be 'protein'; include heavy atoms from segid Q
			nad_sel = "(segid " + " ".join(NAD_CHAINS) + " and not name H*)"
			sel = prot_sel + " or " + nad_sel
		else:
			sel = prot_sel
	else:
		raise ValueError(f"Unknown partner: {partner}")

	ag = u.select_atoms(sel)
	if len(ag) == 0:
		raise ValueError(f"Empty partner selection for '{partner}'. Selection was: {sel}")
	return ag


def select_cp12(u: mda.Universe) -> mda.core.groups.AtomGroup:
	sel = "protein and (segid " + " ".join(CP12_CHAINS) + ") and not name H*"
	ag = u.select_atoms(sel)
	if len(ag) == 0:
		raise ValueError(f"Empty CP12 selection. Check chain IDs/segids.")
	return ag


def run_replicate(strain: str, rep: str, rep_dir: Path, args, outroot: Path) -> dict:
	tpr = rep_dir / "solute.tpr"
	xtc = rep_dir / "solute_fit.xtc"
	pdb = rep_dir / "solute_fixed.pdb"
	for p in (tpr, xtc, pdb):
		if not p.exists():
			raise FileNotFoundError(f"Missing required file: {p}")

	u = load_universe_with_chains(pdb, tpr, xtc)
	cp12 = select_cp12(u)

	partners = ["receptor", "gap2", "prk"] if args.split_subunits else ["receptor"]

	rep_out = outroot / strain / rep
	rep_out.mkdir(parents=True, exist_ok=True)

	out = {"per_partner": {}}

	for partner in partners:
		partner_atoms = select_partner(u, partner, include_nad=args.include_nad)

		cp_df, p_df, n_frames = contact_freq_any_heavy(
			u=u, A=cp12, B=partner_atoms, cutoff=args.cutoff, stride=args.stride
		)

		# add identifiers
		cp_df.insert(0, "strain", strain); cp_df.insert(1, "replicate", rep); cp_df.insert(2, "partner", partner); cp_df.insert(3, "side", "cp12")
		p_df.insert(0, "strain", strain);  p_df.insert(1, "replicate", rep);  p_df.insert(2, "partner", partner);  p_df.insert(3, "side", "receptor")

		cp_class = classify_core_peripheral(cp_df, args.core_thresh, args.peri_thresh, args.core_quantile)
		p_class  = classify_core_peripheral(p_df,  args.core_thresh, args.peri_thresh, args.core_quantile)

		cp_union, cp_cons = cp12_union_consensus(cp_df, args.peri_thresh, args.consensus_n)
		# receptor union/consensus across CP12 copies (single pass)
		rc_union, rc_cons = receptor_union_consensus_by_cp12_copy_singlepass(
			u=u, partner_atoms=partner_atoms, cutoff=args.cutoff, stride=args.stride,
			peri_thresh=args.peri_thresh, consensus_n=args.consensus_n
		)

		# write
		cp_df.to_csv(rep_out / f"interface_resfreq_cp12_side_{partner}.csv", index=False)
		p_df.to_csv(rep_out / f"interface_resfreq_receptor_side_{partner}.csv", index=False)
		cp_class.to_csv(rep_out / f"interface_core_peripheral_cp12_{partner}.csv", index=False)
		p_class.to_csv(rep_out / f"interface_core_peripheral_receptor_{partner}.csv", index=False)

		cp_union.to_csv(rep_out / f"interface_union_cp12_{partner}.csv", index=False)
		cp_cons.to_csv(rep_out / f"interface_consensus_cp12_{partner}.csv", index=False)
		rc_union.to_csv(rep_out / f"interface_union_receptor_{partner}.csv", index=False)
		rc_cons.to_csv(rep_out / f"interface_consensus_receptor_{partner}.csv", index=False)

		# summary
		cp_core_n = int((cp_class["class"] == "core").sum())
		cp_peri_n = int((cp_class["class"] == "peripheral").sum())
		rc_core_n = int((p_class["class"] == "core").sum())
		rc_peri_n = int((p_class["class"] == "peripheral").sum())

		print(f"\n[{strain}/{rep}] partner={partner} frames={n_frames} stride={args.stride} cutoff={args.cutoff:.2f}Å include_nad={args.include_nad if partner=='receptor' else False}")
		print(f"  CP12: core={cp_core_n} peripheral={cp_peri_n}  | top10:\n{top10(cp_df).to_string(index=False)}")
		print(f"  {partner.upper():<6}: core={rc_core_n} peripheral={rc_peri_n}  | top10:\n{top10(p_df).to_string(index=False)}")
		print(f"  CP12 union={len(cp_union)} consensus(n>={args.consensus_n})={len(cp_cons)}")
		print(f"  {partner.upper():<6} union={len(rc_union)} consensus(n>={args.consensus_n})={len(rc_cons)}")

		out["per_partner"][partner] = {
			"cp_df": cp_df, "rc_df": p_df,
			"cp_union": cp_union, "cp_cons": cp_cons,
			"rc_union": rc_union, "rc_cons": rc_cons,
		}

		# rewind trajectory for next partner
		u.trajectory.rewind()

	return out


def parse_args():
	ap = argparse.ArgumentParser(
		description="Compute CP12 interfaces vs GAP2/PRK from MD trajectories using heavy-atom contact frequency."
	)
	ap.add_argument("--root", default=".", help="Project root containing strain/replicate directories.")
	ap.add_argument("--strains", nargs="+", default=["med4", "hm2"], help="Strain directory names.")
	ap.add_argument(
		"--strain_dir",
		action="append",
		default=[],
		help="Map a strain name to an explicit base directory, e.g. med4_ox=/path/to/MED4_NAD/0003.",
	)
	ap.add_argument("--reps", nargs="+", default=["rep1", "rep2", "rep3"], help="Replicate directory names.")
	ap.add_argument("--cutoff", type=float, default=6.0, help="Heavy-atom contact cutoff in Å.")
	ap.add_argument("--core_thresh", type=float, default=0.50, help="Core interface freq threshold.")
	ap.add_argument("--peri_thresh", type=float, default=0.10, help="Peripheral interface freq threshold.")
	ap.add_argument("--core_quantile", type=float, default=None,
					help="Optional: core = top fraction (e.g., 0.20=top 20%%) among residues with freq>=peri_thresh.")
	ap.add_argument("--consensus_n", type=int, default=3, help="Consensus requires meeting peri_thresh in >=N CP12 copies.")
	ap.add_argument("--stride", type=int, default=1, help="Frame stride (subsample frames).")
	ap.add_argument("--include_nad", action="store_true", help="Include NAD (segid Q) in 'receptor' partner selection.")
	ap.add_argument("--split_subunits", action="store_true", default=True,
					help="Compute separate interfaces for gap2 and prk in addition to combined receptor (default: on).")
	ap.add_argument("--no_split_subunits", dest="split_subunits", action="store_false",
					help="Only compute combined receptor (GAP2+PRK).")
	ap.add_argument("--outdir", default="interface_results", help="Output directory under root.")
	return ap.parse_args()


def main():
	args = parse_args()
	root = Path(args.root).resolve()
	outroot = root / args.outdir
	outroot.mkdir(parents=True, exist_ok=True)
	strain_dir_map = {k: Path(v).expanduser().resolve() for k, v in parse_mapping(args.strain_dir, "--strain_dir").items()}

	# Collect per strain per partner for aggregation
	per_strain = {s: {} for s in args.strains}

	for strain in args.strains:
		per_strain[strain] = {}
		strain_base = strain_dir_map.get(strain, root / strain)
		for rep in args.reps:
			rep_dir = strain_base / rep
			if not rep_dir.exists():
				eprint(f"Skipping missing directory: {rep_dir}")
				continue
			res = run_replicate(strain, rep, rep_dir, args, outroot)

			for partner, d in res["per_partner"].items():
				per_strain[strain].setdefault(partner, {
					"cp": [], "rc": [],
					"cp_union": [], "cp_cons": [],
					"rc_union": [], "rc_cons": [],
				})
				per_strain[strain][partner]["cp"].append(d["cp_df"])
				per_strain[strain][partner]["rc"].append(d["rc_df"])
				per_strain[strain][partner]["cp_union"].append(d["cp_union"].assign(strain=strain, replicate=rep, partner=partner))
				per_strain[strain][partner]["cp_cons"].append(d["cp_cons"].assign(strain=strain, replicate=rep, partner=partner))
				per_strain[strain][partner]["rc_union"].append(d["rc_union"].assign(strain=strain, replicate=rep, partner=partner))
				per_strain[strain][partner]["rc_cons"].append(d["rc_cons"].assign(strain=strain, replicate=rep, partner=partner))

	# Strain-level aggregation
	for strain, partners in per_strain.items():
		if not partners:
			continue
		(outroot / strain).mkdir(parents=True, exist_ok=True)

		for partner, d in partners.items():
			if not d["cp"]:
				continue

			cp_agg = aggregate_resfreq_across_reps(d["cp"], key_cols=["chain", "resid", "resname"])
			rc_agg = aggregate_resfreq_across_reps(d["rc"], key_cols=["chain", "resid", "resname"])

			cp_agg.to_csv(outroot / strain / f"{strain}_interface_resfreq_cp12_agg_{partner}.csv", index=False)
			rc_agg.to_csv(outroot / strain / f"{strain}_interface_resfreq_receptor_agg_{partner}.csv", index=False)

			# union/consensus aggregates
			cp_union_all = pd.concat(d["cp_union"], ignore_index=True) if d["cp_union"] else None
			cp_cons_all  = pd.concat(d["cp_cons"],  ignore_index=True) if d["cp_cons"] else None
			rc_union_all = pd.concat(d["rc_union"], ignore_index=True) if d["rc_union"] else None
			rc_cons_all  = pd.concat(d["rc_cons"],  ignore_index=True) if d["rc_cons"] else None

			if cp_union_all is not None and len(cp_union_all):
				agg = aggregate_union_consensus_across_reps(cp_union_all, key_cols=["resid", "resname"])
				agg.to_csv(outroot / strain / f"{strain}_interface_union_cp12_agg_{partner}.csv", index=False)

			if cp_cons_all is not None and len(cp_cons_all):
				agg = aggregate_union_consensus_across_reps(cp_cons_all, key_cols=["resid", "resname"])
				agg.to_csv(outroot / strain / f"{strain}_interface_consensus_cp12_agg_{partner}.csv", index=False)

			if rc_union_all is not None and len(rc_union_all):
				agg = aggregate_union_consensus_across_reps(rc_union_all, key_cols=["chain", "resid", "resname"])
				agg.to_csv(outroot / strain / f"{strain}_interface_union_receptor_agg_{partner}.csv", index=False)

			if rc_cons_all is not None and len(rc_cons_all):
				agg = aggregate_union_consensus_across_reps(rc_cons_all, key_cols=["chain", "resid", "resname"])
				agg.to_csv(outroot / strain / f"{strain}_interface_consensus_receptor_agg_{partner}.csv", index=False)

			print(f"\n[{strain}] aggregated outputs written for partner='{partner}' -> {outroot/strain}")

	print(f"\nDone. Results in: {outroot}")


if __name__ == "__main__":
	main()
