#!/usr/bin/env python3
"""
Binary contact dynamics for CP12 against PRK and GAP2.

For each (condition, rep, CP12_chain, CP12_resid, partner_chain, partner_resid) pair:
	- Build per-frame binary contact signal B(t) in {0,1} at heavy-atom 5.0 A cutoff.
	- Apply minimum-run filter (>= --min_run frames in either state) to suppress
	near-cutoff thermal flicker. This is the cheap equivalent of 4.0/5.5 A
	hysteresis on a 1001-frame trajectory at 100 ps stride.
	- Compute persistence statistics (n_runs, mean/median/p90/max run length).
	- Compute turnover (transitions per ns), Markov-style k_on, k_off.
	- Compute survival fractions S(1 ns), S(5 ns), S(20 ns).
	- Compute two minimum heavy-atom distance summaries:
		* pair-distance summaries across all sampled trajectory frames for each
		retained contact-capable residue pair.
		* contact-state distance summaries using only frames where that pair is
		observed within the contact cutoff.
	- Tag each pair by realized contact chemistry:
		hydrophobic / polar_hbond / salt_bridge / cation_pi /
		same_charge / charged_neutral / polar_hydrophobic / other_mixed / unknown.
	The original broad five-class label is also retained as contact_class_broad
	for backward comparison.

Inputs per replicate dir:
	solute.tpr
	solute_fit.xtc
	solute_fixed.pdb

Outputs (per condition aggregate + per pair detail):
	contact_dynamics/per_pair_stats.csv.gz
	contact_dynamics/per_pair_hmm_aggregated.csv
	contact_dynamics/per_class_summary.csv
	contact_dynamics/condition_summary.csv
"""

from __future__ import annotations

import argparse
import sys
import gzip
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.lib.nsgrid import FastNS

import interface_contacts_freq as icf


DEFAULT_SYSTEM_DIRS = {
	"med4_ox": "source_data/md_trajectories/0000_0003_analysis/med4/0003",
	"med4_red": "source_data/md_trajectories/0000_0003_analysis/med4/0000",
	"hm2_ox": "source_data/md_trajectories/0000_0003_analysis/hm2/0003",
	"hm2_red": "source_data/md_trajectories/0000_0003_analysis/hm2/0000",
}

HYDROPHOBIC = set("AVLIMFWYP")
POLAR_HBOND = set("STNQHC")
POS_CHARGE  = set("KR")
NEG_CHARGE  = set("DE")
AROMATIC    = set("FWY")
CHARGED     = POS_CHARGE | NEG_CHARGE

THREE_TO_ONE = {
	"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","CYX":"C","CYM":"C",
	"GLN":"Q","GLU":"E","GLY":"G","HIS":"H","HID":"H","HIE":"H","HIP":"H",
	"ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
	"THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def aa_one(resname: str) -> str:
	return THREE_TO_ONE.get(resname.upper(), "X")


def aa_seq_class(aa: str) -> str:
	"""Single-residue chemical class used for sequence/contact comparison."""
	if aa in ("-", "", "X") or pd.isna(aa):
		return "gap_or_unknown"
	if aa in AROMATIC:
		return "aromatic"
	if aa in POS_CHARGE:
		return "positive"
	if aa in NEG_CHARGE:
		return "negative"
	if aa in HYDROPHOBIC:
		return "hydrophobic"
	if aa in POLAR_HBOND:
		return "polar"
	return "other"


def aa_shift_class(aa: str) -> str:
	"""Coarser single-residue class used to classify MED4/HM2 substitutions."""
	if aa in ("-", "", "X") or pd.isna(aa):
		return "gap"
	if aa in AROMATIC:
		return "aromatic"
	if aa in POS_CHARGE:
		return "pos"
	if aa in NEG_CHARGE:
		return "neg"
	if aa in HYDROPHOBIC:
		return "hydrophobic"
	if aa in POLAR_HBOND:
		return "polar"
	return "other"


def substitution_shift_category(med4_aa: str, hm2_aa: str) -> str:
	"""Classify the sequence-level MED4/HM2 chemical shift at one HMM column."""
	if med4_aa == hm2_aa:
		return "conserved"
	cm = aa_shift_class(med4_aa)
	ch = aa_shift_class(hm2_aa)
	if cm == "gap" or ch == "gap":
		return "gap"
	if cm == ch:
		return "within-class"
	# Put aromatic-specific comparisons before generic hydrophobic/charged rules.
	if {cm, ch} == {"aromatic", "pos"}:
		return "aromatic↔charged"
	if {cm, ch} == {"aromatic", "neg"}:
		return "aromatic↔charged"
	if "aromatic" in {cm, ch} and "polar" in {cm, ch}:
		return "hydrophobic↔polar"
	if "aromatic" in {cm, ch} and "hydrophobic" in {cm, ch}:
		return "within-class"
	s = {cm, ch}
	if s == {"hydrophobic", "polar"}:
		return "hydrophobic↔polar"
	if s == {"hydrophobic", "pos"} or s == {"hydrophobic", "neg"}:
		return "hydrophobic↔charged"
	if s == {"polar", "pos"} or s == {"polar", "neg"}:
		return "polar↔charged"
	if s == {"pos", "neg"}:
		return "charge-flip"
	return "other"


def classify_pair_legacy(a: str, b: str) -> str:
	"""Original broad five-class contact chemistry label."""
	if a == "X" or b == "X":
		return "mixed"
	if (a in POS_CHARGE and b in NEG_CHARGE) or (a in NEG_CHARGE and b in POS_CHARGE):
		return "salt_bridge"
	if (a in AROMATIC and b in POS_CHARGE) or (a in POS_CHARGE and b in AROMATIC):
		return "cation_pi"
	if a in HYDROPHOBIC and b in HYDROPHOBIC:
		return "hydrophobic"
	if a in POLAR_HBOND or b in POLAR_HBOND:
		return "polar_hbond"
	return "mixed"


def classify_pair(a: str, b: str) -> str:
	"""Refined realized chemistry of a CP12-partner residue pair.

	This keeps specific interpretable classes first, then splits the old
	catch-all `mixed` class into categories that can be compared to sequence
	chemical shifts.
	"""
	if a == "X" or b == "X":
		return "unknown"
	if (a in POS_CHARGE and b in NEG_CHARGE) or (a in NEG_CHARGE and b in POS_CHARGE):
		return "salt_bridge"
	if (a in AROMATIC and b in POS_CHARGE) or (a in POS_CHARGE and b in AROMATIC):
		return "cation_pi"
	if a in HYDROPHOBIC and b in HYDROPHOBIC:
		return "hydrophobic"
	if (a in POS_CHARGE and b in POS_CHARGE) or (a in NEG_CHARGE and b in NEG_CHARGE):
		return "same_charge"
	if (a in POLAR_HBOND and b in POLAR_HBOND) or (
		(a in POLAR_HBOND and b in CHARGED) or (b in POLAR_HBOND and a in CHARGED)
	):
		return "polar_hbond"
	if (a in CHARGED and b not in CHARGED) or (b in CHARGED and a not in CHARGED):
		return "charged_neutral"
	if (a in POLAR_HBOND and b in HYDROPHOBIC) or (b in POLAR_HBOND and a in HYDROPHOBIC):
		return "polar_hydrophobic"
	return "other_mixed"


def contact_class_broad(refined_class: str) -> str:
	"""Map refined contact classes back to the original broad class family."""
	if refined_class in {"hydrophobic", "polar_hbond", "salt_bridge", "cation_pi"}:
		return refined_class
	return "mixed"


def _flip_short_runs(b: np.ndarray, target_state: bool, min_run: int) -> np.ndarray:
	"""Flip every interior run of `target_state` with length < min_run to `not target_state`.

	Edge runs (at the very start or end of the array) are kept as-is — they may be
	truncated by the trajectory boundary rather than spurious.
	"""
	if b.size == 0 or min_run <= 1:
		return b.astype(bool, copy=True)
	out = b.astype(bool, copy=True)
	n = out.size
	diff = np.diff(out.astype(np.int8), prepend=out[0])
	starts = np.r_[0, np.where(diff != 0)[0]]
	ends = np.r_[starts[1:], n]
	lengths = ends - starts
	states = out[starts]
	for i in range(1, len(starts) - 1):  # skip the two edge runs
		if states[i] == target_state and lengths[i] < min_run:
			out[starts[i]:ends[i]] = not target_state
	return out


def min_run_filter(b: np.ndarray, min_run: int) -> np.ndarray:
	"""Suppress thermal flicker by morphological opening then closing of B(t).

	Opening (flip short 1-runs to 0)  -> kills brief spurious "on" events.
	Closing (flip short 0-runs to 1) -> bridges brief disconnects in a real contact.

	A contact is considered "real" iff it persists for >= min_run frames; brief
	disconnects shorter than min_run do not break it.
	"""
	if b.size == 0 or min_run <= 1:
		return b.astype(bool, copy=True)
	b = _flip_short_runs(b, target_state=True,  min_run=min_run)
	b = _flip_short_runs(b, target_state=False, min_run=min_run)
	return b


def run_length_stats(b: np.ndarray) -> dict:
	"""Run-length statistics for the 1-runs and 0-runs of a binary array."""
	if b.size == 0:
		return dict(n_runs_on=0, n_runs_off=0,
					mean_on=0.0, median_on=0.0, p90_on=0.0, max_on=0,
					mean_off=0.0, median_off=0.0, p90_off=0.0, max_off=0,
					transitions_01=0, transitions_10=0)
	diff = np.diff(b.astype(np.int8), prepend=b[0])
	starts = np.r_[0, np.where(diff != 0)[0]]
	ends = np.r_[starts[1:], b.size]
	lengths = ends - starts
	states = b[starts]
	on_runs = lengths[states == True]
	off_runs = lengths[states == False]

	def stats(arr):
		if arr.size == 0:
			return 0.0, 0.0, 0.0, 0
		return (float(np.mean(arr)),
				float(np.median(arr)),
				float(np.percentile(arr, 90)),
				int(np.max(arr)))
	m_on, med_on, p90_on, max_on = stats(on_runs)
	m_off, med_off, p90_off, max_off = stats(off_runs)

	# transitions
	transitions_01 = int(np.sum((np.diff(b.astype(np.int8))) == 1))
	transitions_10 = int(np.sum((np.diff(b.astype(np.int8))) == -1))

	return dict(n_runs_on=int(on_runs.size), n_runs_off=int(off_runs.size),
				mean_on=m_on, median_on=med_on, p90_on=p90_on, max_on=max_on,
				mean_off=m_off, median_off=med_off, p90_off=p90_off, max_off=max_off,
				transitions_01=transitions_01, transitions_10=transitions_10)


def survival_at(b: np.ndarray, frame_indices: list[int]) -> dict:
	"""Empirical survival S(t) of 1-runs measured at requested lags (frames).

	S(t) := P(a contact still on at time t after first formed | contact formed).
	Computed as: take all 1-runs of length L; for each lag k, fraction of runs
	with L >= k. Returns NaN if no 1-runs.
	"""
	diff = np.diff(b.astype(np.int8), prepend=b[0])
	starts = np.r_[0, np.where(diff != 0)[0]]
	ends = np.r_[starts[1:], b.size]
	lengths = ends - starts
	states = b[starts]
	on_runs = lengths[states == True]
	out = {}
	if on_runs.size == 0:
		for k in frame_indices:
			out[k] = float("nan")
		return out
	for k in frame_indices:
		out[k] = float(np.mean(on_runs >= k))
	return out


def extract_pair_signals(
	u: mda.Universe,
	cp12_atoms: mda.core.groups.AtomGroup,
	partner_atoms: mda.core.groups.AtomGroup,
	cutoff: float,
	n_frames_total: int,
) -> tuple[dict, dict, dict]:
	"""One pass over the trajectory.

	Returns:
		contact_signals: dict[(cp12_chain, cp12_resid, p_chain, p_resid)] -> np.array(bool, n_frames_total)
		contact_distances: same keys -> np.array(float, n_frames_total), minimum heavy-atom distance in A
		cp12_resname  : dict[(cp12_chain, cp12_resid)] -> resname
		p_resname     : dict[(p_chain,    p_resid)]   -> resname
	"""
	cp12_atom2res, cp12_resid, cp12_resname, cp12_chain, _ = icf.residues_from_atoms(cp12_atoms)
	p_atom2res, p_resid, p_resname, p_chain, _           = icf.residues_from_atoms(partner_atoms)

	contact_signals: dict[tuple, np.ndarray] = {}
	contact_distances: dict[tuple, np.ndarray] = {}
	cp12_resname_map: dict[tuple, str] = {}
	p_resname_map: dict[tuple, str] = {}

	for ci, cr, crn in zip(cp12_chain, cp12_resid, cp12_resname):
		cp12_resname_map[(str(ci), int(cr))] = str(crn)
	for pi, pr, prn in zip(p_chain, p_resid, p_resname):
		p_resname_map[(str(pi), int(pr))] = str(prn)

	frame_idx = 0
	for ts in u.trajectory:
		if frame_idx >= n_frames_total:
			break
		box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None
		ns = FastNS(cutoff, partner_atoms.positions, box=box)
		results = ns.search(cp12_atoms.positions)
		pairs = results.get_pairs()
		if pairs.size:
			distances = results.get_pair_distances()
			c_res_idx = cp12_atom2res[pairs[:, 0]]
			p_res_idx = p_atom2res[pairs[:, 1]]
			# de-dup to unique residue-pair contacts this frame
			stacked = (c_res_idx.astype(np.int64) << 32) | p_res_idx.astype(np.int64)
			uniq, inverse = np.unique(stacked, return_inverse=True)
			min_dists = np.full(uniq.size, np.inf, dtype=float)
			np.minimum.at(min_dists, inverse, distances)
			c_unique = (uniq >> 32).astype(np.int64)
			p_unique = (uniq & 0xFFFFFFFF).astype(np.int64)
			for cri, pri, d_min in zip(c_unique, p_unique, min_dists):
				key = (str(cp12_chain[cri]), int(cp12_resid[cri]),
						str(p_chain[pri]),    int(p_resid[pri]))
				arr = contact_signals.get(key)
				if arr is None:
					arr = np.zeros(n_frames_total, dtype=bool)
					contact_signals[key] = arr
					contact_distances[key] = np.full(n_frames_total, np.nan, dtype=float)
				arr[frame_idx] = True
				contact_distances[key][frame_idx] = float(d_min)
		frame_idx += 1

	return contact_signals, contact_distances, cp12_resname_map, p_resname_map


def distance_stats_for_pair(b: np.ndarray, distances: np.ndarray) -> dict:
	"""Distance summaries for raw observed contact frames retained by the clean signal.

	The minimum-run filter can bridge very short off gaps. Those bridged frames
	do not have observed atom-pair distances, so distance summaries use only
	frames where the pair was actually found within the cutoff.
	"""
	if distances is None or distances.size == 0:
		vals = np.array([], dtype=float)
	else:
		vals = distances[b & np.isfinite(distances)]
	if vals.size == 0:
		return dict(
			n_frames_with_distance=0,
			mean_contact_dist_A=np.nan,
			sd_contact_dist_A=np.nan,
			median_contact_dist_A=np.nan,
			min_contact_dist_A=np.nan,
			p90_contact_dist_A=np.nan,
		)
	return dict(
		n_frames_with_distance=int(vals.size),
		mean_contact_dist_A=float(np.mean(vals)),
		sd_contact_dist_A=float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
		median_contact_dist_A=float(np.median(vals)),
		min_contact_dist_A=float(np.min(vals)),
		p90_contact_dist_A=float(np.percentile(vals, 90)),
	)


def residue_atomgroups(atoms: mda.core.groups.AtomGroup) -> dict[tuple[str, int], mda.core.groups.AtomGroup]:
	"""Map (chain, resid) to the heavy-atom AtomGroup used by contact detection."""
	atom2res, resid, _resname, chain, _ = icf.residues_from_atoms(atoms)
	groups: dict[tuple[str, int], mda.core.groups.AtomGroup] = {}
	for local_idx, (chain_id, resid_i) in enumerate(zip(chain, resid)):
		atom_idx = np.where(atom2res == local_idx)[0]
		groups[(str(chain_id), int(resid_i))] = atoms[atom_idx]
	return groups


def _summarise_distance_values(values: list[float], prefix: str) -> dict:
	vals = np.asarray(values, dtype=float)
	vals = vals[np.isfinite(vals)]
	if vals.size == 0:
		return {
			f"n_frames_{prefix}_distance": 0,
			f"mean_{prefix}_dist_A": np.nan,
			f"sd_{prefix}_dist_A": np.nan,
			f"median_{prefix}_dist_A": np.nan,
			f"min_{prefix}_dist_A": np.nan,
			f"p90_{prefix}_dist_A": np.nan,
			f"max_{prefix}_dist_A": np.nan,
		}
	return {
		f"n_frames_{prefix}_distance": int(vals.size),
		f"mean_{prefix}_dist_A": float(np.mean(vals)),
		f"sd_{prefix}_dist_A": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
		f"median_{prefix}_dist_A": float(np.median(vals)),
		f"min_{prefix}_dist_A": float(np.min(vals)),
		f"p90_{prefix}_dist_A": float(np.percentile(vals, 90)),
		f"max_{prefix}_dist_A": float(np.max(vals)),
	}


def full_trajectory_distance_stats(
	u: mda.Universe,
	cp12_atoms: mda.core.groups.AtomGroup,
	partner_atoms: mda.core.groups.AtomGroup,
	keys: Iterable[tuple],
) -> dict[tuple, dict]:
	"""Distance summaries over all trajectory frames for retained contact-capable pairs.

	This reproduces the original Layer 1 distance meaning: for each residue pair
	that passed the binary contact filter, measure the minimum heavy-atom
	distance in every sampled frame, not only in frames where the pair is within
	the contact cutoff. The resulting `mean_pair_dist_A` is therefore a
	trajectory-averaged pair-proximity metric.
	"""
	keys = list(keys)
	if not keys:
		return {}

	cp_groups = residue_atomgroups(cp12_atoms)
	partner_groups = residue_atomgroups(partner_atoms)
	values_by_key: dict[tuple, list[float]] = {key: [] for key in keys}

	u.trajectory.rewind()
	for ts in u.trajectory:
		box = ts.dimensions if ts.dimensions is not None and np.all(ts.dimensions[:3] > 0) else None
		for key in keys:
			cp_chain, cp_resid, p_chain, p_resid = key
			cp_ag = cp_groups.get((str(cp_chain), int(cp_resid)))
			p_ag = partner_groups.get((str(p_chain), int(p_resid)))
			if cp_ag is None or p_ag is None or len(cp_ag) == 0 or len(p_ag) == 0:
				values_by_key[key].append(np.nan)
				continue
			values_by_key[key].append(float(np.nanmin(distance_array(cp_ag.positions, p_ag.positions, box=box))))

	return {key: _summarise_distance_values(vals, "pair") for key, vals in values_by_key.items()}


def stats_for_pair(
	b: np.ndarray,
	distances: np.ndarray,
	dt_ps: float,
	surv_frames: list[int],
) -> dict:
	n = int(b.size)
	occ = float(np.mean(b)) if n else float("nan")
	rls = run_length_stats(b)
	# convert frame counts to ps
	on_to_ps = lambda x: float(x * dt_ps)
	# k_on, k_off in 1/ns: 1 / mean dwell in ns
	k_off = (1.0 / (rls["mean_on"]  * dt_ps / 1000.0)) if rls["mean_on"]  > 0 else float("nan")
	k_on  = (1.0 / (rls["mean_off"] * dt_ps / 1000.0)) if rls["mean_off"] > 0 else float("nan")
	sv = survival_at(b, surv_frames)

	return dict(
		n_frames=n,
		occupancy=occ,
		n_runs_on=rls["n_runs_on"],
		mean_on_ps=on_to_ps(rls["mean_on"]),
		median_on_ps=on_to_ps(rls["median_on"]),
		p90_on_ps=on_to_ps(rls["p90_on"]),
		max_on_ps=on_to_ps(rls["max_on"]),
		mean_off_ps=on_to_ps(rls["mean_off"]),
		median_off_ps=on_to_ps(rls["median_off"]),
		p90_off_ps=on_to_ps(rls["p90_off"]),
		transitions_01=rls["transitions_01"],
		transitions_10=rls["transitions_10"],
		turnover_per_ns=(rls["transitions_01"] + rls["transitions_10"]) / (n * dt_ps / 1000.0)
							if n > 0 else float("nan"),
		k_off_per_ns=k_off,
		k_on_per_ns=k_on,
		survival_1ns=sv.get(int(round(1000.0 / dt_ps)), float("nan")),
		survival_5ns=sv.get(int(round(5000.0 / dt_ps)), float("nan")),
		survival_20ns=sv.get(int(round(20000.0 / dt_ps)), float("nan")),
		**distance_stats_for_pair(b, distances),
	)


def process_replicate(
	condition: str,
	rep: str,
	rep_dir: Path,
	partner: str,
	cutoff: float,
	min_run: int,
	surv_frames: list[int],
	dt_ps: float,
) -> pd.DataFrame:
	"""Returns a per-(pair) DataFrame for this (condition, rep, partner)."""
	tpr = rep_dir / "solute.tpr"
	xtc = rep_dir / "solute_fit.xtc"
	pdb = rep_dir / "solute_fixed.pdb"
	if not (tpr.exists() and xtc.exists() and pdb.exists()):
		return pd.DataFrame()

	u = icf.load_universe_with_chains(pdb, tpr, xtc)
	cp12 = icf.select_cp12(u)
	p_atoms = icf.select_partner(u, partner, include_nad=False)
	n_frames_total = u.trajectory.n_frames

	signals, distances, cp12_resname_map, p_resname_map = extract_pair_signals(
		u, cp12, p_atoms, cutoff, n_frames_total
	)

	retained_signals: dict[tuple, np.ndarray] = {}
	for key, raw in signals.items():
		# noise cleanup
		b_clean = min_run_filter(raw, min_run)
		# if pair never actually contacts after cleanup, skip
		if not b_clean.any():
			continue
		retained_signals[key] = b_clean

	pair_distance_stats = full_trajectory_distance_stats(
		u, cp12, p_atoms, retained_signals.keys()
	)

	rows = []
	for key, b_clean in retained_signals.items():
		cp_chain, cp_resid, p_chain_, p_resid_ = key
		raw = signals[key]
		stats = stats_for_pair(b_clean, distances.get(key), dt_ps, surv_frames)
		cp_resn = cp12_resname_map.get((cp_chain, cp_resid), "UNK")
		p_resn  = p_resname_map.get((p_chain_,    p_resid_),   "UNK")
		cp_aa = aa_one(cp_resn)
		p_aa  = aa_one(p_resn)
		refined_class = classify_pair(cp_aa, p_aa)
		rows.append({
			"condition": condition,
			"rep": rep,
			"partner": partner,
			"cp12_chain": cp_chain,
			"cp12_resid": cp_resid,
			"cp12_resname": cp_resn,
			"cp12_aa": cp_aa,
			"partner_chain": p_chain_,
			"partner_resid": p_resid_,
			"partner_resname": p_resn,
			"partner_aa": p_aa,
			"cp12_seq_class": aa_seq_class(cp_aa),
			"partner_seq_class": aa_seq_class(p_aa),
			"contact_class": refined_class,
			"contact_class_broad": contact_class_broad(refined_class),
			"contact_class_legacy": classify_pair_legacy(cp_aa, p_aa),
			**pair_distance_stats.get(key, _summarise_distance_values([], "pair")),
			**stats,
			"raw_occupancy": float(np.mean(raw)),
			"n_frames_raw_contact": int(raw.sum()),
		})
	return pd.DataFrame(rows)


def annotate_hmm(df: pd.DataFrame, aln_csv: Path) -> pd.DataFrame:
	"""Add aln_col_1based, med4_aa, hm2_aa from the HMM alignment."""
	if df.empty:
		return df
	aln = pd.read_csv(aln_csv)
	aln_med4 = aln[["aln_col_1based", "med4_aa", "hm2_aa", "med4_pos"]].dropna(subset=["med4_pos"]).copy()
	aln_med4["med4_pos"] = aln_med4["med4_pos"].astype(int)
	aln_hm2  = aln[["aln_col_1based", "med4_aa", "hm2_aa", "hm2_pos"]].dropna(subset=["hm2_pos"]).copy()
	aln_hm2["hm2_pos"]   = aln_hm2["hm2_pos"].astype(int)

	out = df.copy()
	out["hmm_col"] = np.nan
	out["host_aa_at_hmm"] = ""
	out["partner_host_aa_at_hmm"] = ""
	out["med4_aa_at_hmm"] = ""
	out["hm2_aa_at_hmm"] = ""

	# Strain split
	med_mask = out["condition"].str.startswith("med4")
	hm2_mask = out["condition"].str.startswith("hm2")

	if med_mask.any():
		merged = out.loc[med_mask].merge(
			aln_med4.rename(columns={"med4_pos": "cp12_resid",
										"aln_col_1based": "hmm_col_m",
										"med4_aa": "host_aa_m",
										"hm2_aa": "alt_aa_m"}),
			on="cp12_resid", how="left",
		)
		out.loc[med_mask, "hmm_col"] = merged["hmm_col_m"].values
		out.loc[med_mask, "host_aa_at_hmm"] = merged["host_aa_m"].values
		out.loc[med_mask, "partner_host_aa_at_hmm"] = merged["alt_aa_m"].values
		out.loc[med_mask, "med4_aa_at_hmm"] = merged["host_aa_m"].values
		out.loc[med_mask, "hm2_aa_at_hmm"] = merged["alt_aa_m"].values

	if hm2_mask.any():
		merged = out.loc[hm2_mask].merge(
			aln_hm2.rename(columns={"hm2_pos": "cp12_resid",
									"aln_col_1based": "hmm_col_h",
									"hm2_aa": "host_aa_h",
									"med4_aa": "alt_aa_h"}),
			on="cp12_resid", how="left",
		)
		out.loc[hm2_mask, "hmm_col"] = merged["hmm_col_h"].values
		out.loc[hm2_mask, "host_aa_at_hmm"] = merged["host_aa_h"].values
		out.loc[hm2_mask, "partner_host_aa_at_hmm"] = merged["alt_aa_h"].values
		out.loc[hm2_mask, "hm2_aa_at_hmm"] = merged["host_aa_h"].values
		out.loc[hm2_mask, "med4_aa_at_hmm"] = merged["alt_aa_h"].values

	out["med4_seq_class_at_hmm"] = out["med4_aa_at_hmm"].map(aa_seq_class)
	out["hm2_seq_class_at_hmm"] = out["hm2_aa_at_hmm"].map(aa_seq_class)
	out["substitution_shift_category"] = [
		substitution_shift_category(m, h)
		for m, h in zip(out["med4_aa_at_hmm"], out["hm2_aa_at_hmm"])
	]
	out["is_substituted_hmm"] = (
		out["med4_aa_at_hmm"].ne(out["hm2_aa_at_hmm"])
		& out["med4_aa_at_hmm"].ne("")
		& out["hm2_aa_at_hmm"].ne("")
	)

	return out


def aggregate_per_hmm(per_pair: pd.DataFrame) -> pd.DataFrame:
	"""Aggregate per-pair stats up to (condition, hmm_col, partner, contact_class).

	Aggregates over CP12 chains (C/H/M/P), partner residue+chain, and replicates.
	Reports mean + SEM (treating each (chain,rep,partner_residue) row as an
	independent observation).
	"""
	if per_pair.empty:
		return pd.DataFrame()
	df = per_pair.copy()
	df = df.dropna(subset=["hmm_col"])
	df["hmm_col"] = df["hmm_col"].astype(int)
	keys = ["condition", "partner", "hmm_col", "contact_class"]

	def agg_metric(s):
		s = s.dropna()
		if s.empty:
			return pd.Series(dict(mean=np.nan, sem=np.nan, n=0))
		return pd.Series(dict(
			mean=float(s.mean()),
			sem=float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0,
			n=int(len(s)),
		))

	metrics = ["occupancy", "mean_on_ps", "median_on_ps", "p90_on_ps",
				"max_on_ps", "turnover_per_ns", "k_off_per_ns", "k_on_per_ns",
				"survival_1ns", "survival_5ns", "survival_20ns",
				"mean_pair_dist_A", "sd_pair_dist_A",
				"median_pair_dist_A", "min_pair_dist_A",
				"p90_pair_dist_A", "max_pair_dist_A",
				"n_frames_pair_distance",
				"mean_contact_dist_A", "sd_contact_dist_A",
				"median_contact_dist_A", "min_contact_dist_A",
				"p90_contact_dist_A", "n_frames_with_distance"]

	out = []
	for keys_vals, g in df.groupby(keys):
		row = dict(zip(keys, keys_vals))
		row["host_aa_at_hmm"] = g["host_aa_at_hmm"].mode().iat[0] if not g["host_aa_at_hmm"].mode().empty else ""
		row["partner_host_aa_at_hmm"] = g["partner_host_aa_at_hmm"].mode().iat[0] if not g["partner_host_aa_at_hmm"].mode().empty else ""
		for col in ["med4_aa_at_hmm", "hm2_aa_at_hmm", "med4_seq_class_at_hmm",
					"hm2_seq_class_at_hmm", "substitution_shift_category"]:
			row[col] = g[col].mode().iat[0] if col in g and not g[col].mode().empty else ""
		row["n_pairs"] = int(len(g))
		for m in metrics:
			ag = agg_metric(g[m])
			row[f"{m}_mean"] = ag["mean"]
			row[f"{m}_sem"]  = ag["sem"]
		out.append(row)
	return pd.DataFrame(out)


def aggregate_per_class(per_pair: pd.DataFrame) -> pd.DataFrame:
	"""Class-level summary: per (condition, partner, contact_class)."""
	if per_pair.empty:
		return pd.DataFrame()
	keys = ["condition", "partner", "contact_class"]
	metrics = ["occupancy", "mean_on_ps", "median_on_ps", "p90_on_ps",
				"max_on_ps", "turnover_per_ns", "k_off_per_ns",
				"survival_1ns", "survival_5ns", "survival_20ns",
				"mean_pair_dist_A", "sd_pair_dist_A",
				"median_pair_dist_A", "min_pair_dist_A",
				"p90_pair_dist_A", "max_pair_dist_A",
				"n_frames_pair_distance",
				"mean_contact_dist_A", "sd_contact_dist_A",
				"median_contact_dist_A", "min_contact_dist_A",
				"p90_contact_dist_A", "n_frames_with_distance"]
	out = []
	for keys_vals, g in per_pair.groupby(keys):
		row = dict(zip(keys, keys_vals))
		row["n_pairs"] = int(len(g))
		for m in metrics:
			v = g[m].dropna()
			row[f"{m}_mean"]   = float(v.mean()) if len(v) else np.nan
			row[f"{m}_median"] = float(v.median()) if len(v) else np.nan
			row[f"{m}_sem"]    = float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
		out.append(row)
	return pd.DataFrame(out)


def aggregate_condition(per_pair: pd.DataFrame) -> pd.DataFrame:
	"""High-level (condition, partner) summary."""
	if per_pair.empty:
		return pd.DataFrame()
	keys = ["condition", "partner"]
	metrics = ["occupancy", "mean_on_ps", "p90_on_ps", "turnover_per_ns",
				"k_off_per_ns", "survival_5ns", "survival_20ns",
				"mean_pair_dist_A", "sd_pair_dist_A",
				"mean_contact_dist_A", "sd_contact_dist_A"]
	out = []
	for keys_vals, g in per_pair.groupby(keys):
		row = dict(zip(keys, keys_vals))
		row["n_pairs"] = int(len(g))
		for m in metrics:
			v = g[m].dropna()
			row[f"{m}_mean"]   = float(v.mean()) if len(v) else np.nan
			row[f"{m}_median"] = float(v.median()) if len(v) else np.nan
		out.append(row)
	return pd.DataFrame(out)


def parse_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("--outdir", default="contact_dynamics")
	ap.add_argument("--conditions", nargs="+",
					default=["med4_ox", "med4_red", "hm2_ox", "hm2_red"])
	ap.add_argument("--reps", nargs="+", default=["rep1", "rep2", "rep3"])
	ap.add_argument("--partners", nargs="+", default=["prk", "gap2"])
	ap.add_argument("--system_dir", action="append", default=[],
					help="Override condition->dir mapping, e.g. med4_ox=/path/.../0003")
	ap.add_argument(
		"--alignment_csv",
		default="source_data/family_alignment/CP12_MED4_vs_HM2_alignment_annotation.csv",
	)
	ap.add_argument("--cutoff", type=float, default=5.0,
					help="Heavy-atom contact cutoff in A (default 5.0).")
	ap.add_argument("--min_run", type=int, default=3,
					help="Minimum-run filter length, in frames. 3 frames at 100 ps/frame "
							"= 300 ps minimum contact event.")
	ap.add_argument("--dt_ps", type=float, default=100.0,
					help="Trajectory frame stride in ps.")
	ap.add_argument("--save_per_pair", action="store_true", default=True,
					help="Write per-pair CSV (gzipped, can be large).")
	return ap.parse_args()


def main() -> None:
	args = parse_args()
	outdir = Path(args.outdir).resolve()
	outdir.mkdir(parents=True, exist_ok=True)

	# condition -> base dir
	cond_dirs = {k: Path(v) for k, v in DEFAULT_SYSTEM_DIRS.items()}
	for entry in args.system_dir:
		if "=" not in entry:
			raise SystemExit(f"--system_dir entries must be name=path, got: {entry}")
		k, v = entry.split("=", 1)
		cond_dirs[k.strip()] = Path(v.strip())

	# survival frame lags
	surv_frames = [int(round(1000.0  / args.dt_ps)),
					int(round(5000.0  / args.dt_ps)),
					int(round(20000.0 / args.dt_ps))]

	all_pairs: list[pd.DataFrame] = []
	for cond in args.conditions:
		base = cond_dirs.get(cond)
		if base is None or not base.exists():
			print(f"[skip] {cond}: directory missing ({base})", file=sys.stderr)
			continue
		for rep in args.reps:
			rep_dir = base / rep
			if not rep_dir.exists():
				print(f"[skip] {cond}/{rep}: missing", file=sys.stderr)
				continue
			for partner in args.partners:
				print(f"[run ] {cond}/{rep}/{partner}", file=sys.stderr, flush=True)
				df = process_replicate(
					cond, rep, rep_dir, partner,
					cutoff=args.cutoff,
					min_run=args.min_run,
					surv_frames=surv_frames,
					dt_ps=args.dt_ps,
				)
				if not df.empty:
					all_pairs.append(df)
					print(f"        {len(df):d} pairs after min-run filter", file=sys.stderr)

	if not all_pairs:
		raise SystemExit("No pair data produced. Check inputs.")

	per_pair = pd.concat(all_pairs, ignore_index=True)
	per_pair = annotate_hmm(per_pair, Path(args.alignment_csv).resolve())

	if args.save_per_pair:
		out_path = outdir / "per_pair_stats.csv.gz"
		per_pair.to_csv(out_path, index=False, compression="gzip")
		print(f"wrote {out_path}  ({len(per_pair)} rows)", file=sys.stderr)

	# aggregations
	per_hmm = aggregate_per_hmm(per_pair)
	if not per_hmm.empty:
		out_path = outdir / "per_pair_hmm_aggregated.csv"
		per_hmm.to_csv(out_path, index=False)
		print(f"wrote {out_path}  ({len(per_hmm)} rows)", file=sys.stderr)

	per_class = aggregate_per_class(per_pair)
	if not per_class.empty:
		out_path = outdir / "per_class_summary.csv"
		per_class.to_csv(out_path, index=False)
		print(f"wrote {out_path}  ({len(per_class)} rows)", file=sys.stderr)

	cond_summary = aggregate_condition(per_pair)
	if not cond_summary.empty:
		out_path = outdir / "condition_summary.csv"
		cond_summary.to_csv(out_path, index=False)
		print(f"wrote {out_path}  ({len(cond_summary)} rows)", file=sys.stderr)

	print("done.", file=sys.stderr)


if __name__ == "__main__":
	main()
