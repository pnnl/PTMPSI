#!/usr/bin/env python3
"""
Compute complex-level SASA and CP12-partner buried surface area with GROMACS.

This script is intended for the existing Dark Complex trajectory layout:

  {analysis_dir}/{med4,hm2}/{0000,0003}/rep*/solute_fixed.pdb
  {analysis_dir}/{med4,hm2}/{0000,0003}/rep*/solute.tpr
  {analysis_dir}/{med4,hm2}/{0000,0003}/rep*/solute_fit.xtc

It builds a chain-based GROMACS index file from solute_fixed.pdb and runs
`gmx sasa` for:

  - whole protein complex SASA
  - CP12, PRK, and GAP2/GAPDH SASA in the assembled complex from atom-level `-oa`
  - per-residue protein SASA summaries by summing atom-level `-oa`
  - isolated CP12, PRK, and GAP2/GAPDH SASA
  - CP12-PRK, CP12-GAP2, and CP12-(PRK+GAP2) pair-complex SASA

Buried surface area is then calculated as:

  BSA(A,B) = SASA(A isolated) + SASA(B isolated) - SASA(A+B complex)
"""

from __future__ import annotations

import argparse
import math
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


GAP2_CHAINS = set("ADEFI KNO".replace(" ", ""))
PRK_CHAINS = set("BGJL")
CP12_CHAINS = set("CHMP")
NAD_CHAINS = set("Q")
CONDITIONS = {"0000": "reduced", "0003": "oxidized"}


@dataclass(frozen=True)
class AtomRecord:
	index1: int
	name: str
	resname: str
	chain: str
	resid: int


@dataclass(frozen=True)
class ReplicateInput:
	strain: str
	system: str
	condition: str
	replicate: str
	rep_dir: Path
	pdb: Path
	tpr: Path
	xtc: Path


@dataclass(frozen=True)
class SasaJob:
	metric: str
	surface_group: str
	output_group: str


SASA_JOBS: Tuple[SasaJob, ...] = (
	SasaJob("protein_complex_sasa_A2", "ProteinComplex", "ProteinComplex"),
	SasaJob("full_complex_with_nad_sasa_A2", "FullComplex", "FullComplex"),
	SasaJob("cp12_isolated_sasa_A2", "CP12", "CP12"),
	SasaJob("prk_isolated_sasa_A2", "PRK", "PRK"),
	SasaJob("gap2_isolated_sasa_A2", "GAP2", "GAP2"),
	SasaJob("prk_gap2_pair_sasa_A2", "PRK_GAP2", "PRK_GAP2"),
	SasaJob("cp12_prk_pair_sasa_A2", "CP12_PRK", "CP12_PRK"),
	SasaJob("cp12_gap2_pair_sasa_A2", "CP12_GAP2", "CP12_GAP2"),
    SasaJob("cp12_prk_gap2_pair_sasa_A2", "CP12_PRK_GAP2", "CP12_PRK_GAP2"),
)

ATOM_COMPONENT_JOBS: Tuple[SasaJob, ...] = (
	SasaJob("cp12_in_complex_sasa_A2", "ProteinComplex", "CP12"),
	SasaJob("prk_in_complex_sasa_A2", "ProteinComplex", "PRK"),
	SasaJob("gap2_in_complex_sasa_A2", "ProteinComplex", "GAP2"),
)

BSA_FORMULAS: Dict[str, Tuple[str, str, str]] = {
	"bsa_cp12_prk_A2": (
		"cp12_isolated_sasa_A2",
		"prk_isolated_sasa_A2",
		"cp12_prk_pair_sasa_A2",
	),
	"bsa_cp12_gap2_A2": (
		"cp12_isolated_sasa_A2",
		"gap2_isolated_sasa_A2",
		"cp12_gap2_pair_sasa_A2",
	),
	"bsa_cp12_prk_gap2_A2": (
		"cp12_isolated_sasa_A2",
		"prk_gap2_pair_sasa_A2",
		"cp12_prk_gap2_pair_sasa_A2",
	),
}


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(
		description="Run GROMACS SASA and calculate CP12-partner buried surface area."
	)
	p.add_argument("--analysis-dir", default="data_dir/0000_0003_analysis")
	p.add_argument(
		"--out-dir",
		default="data_dir/0000_0003_analysis/complex_sasa_bsa_results",
	)
	p.add_argument("--strains", nargs="+", default=["med4", "hm2"])
	p.add_argument("--systems", nargs="+", default=["0000", "0003"])
	p.add_argument("--replicates", nargs="+", default=["rep1", "rep2", "rep3"])
	p.add_argument("--pdb-name", default="solute_fixed.pdb")
	p.add_argument("--tpr-name", default="solute.tpr")
	p.add_argument("--xtc-name", default="solute_fit.xtc")
	p.add_argument(
		"--mode",
		choices=["static", "dynamic", "both"],
		default="dynamic",
		help="static uses solute_fixed.pdb; dynamic uses solute_fit.xtc.",
	)
	p.add_argument(
		"--static-source",
		choices=["pdb", "trajectory"],
		default="pdb",
		help=(
			"Coordinate source for static mode. 'pdb' uses solute_fixed.pdb. "
			"'trajectory' uses one frame from solute_fit.xtc with solute.tpr, "
			"which can avoid GROMACS PDB radius/mass guessing warnings."
		),
	)
	p.add_argument(
		"--static-time-ps",
		type=float,
		default=0.0,
		help="Frame time used when --static-source trajectory.",
	)
	p.add_argument(
		"--begin-ps",
		type=float,
		default=None,
		help="Optional start time for dynamic analysis, e.g. 20000 to omit first 20 ns.",
	)
	p.add_argument("--end-ps", type=float, default=None)
	p.add_argument(
		"--dt-ps",
		type=float,
		default=None,
		help="Optional GROMACS -dt sampling interval in ps. Use 2000 for every ~20th frame if frames are 100 ps apart.",
	)
	p.add_argument(
		"--probe-radius",
		type=float,
		default=0.14,
		help="Probe radius in nm for gmx sasa. 0.14 nm equals 1.4 A.",
	)
	p.add_argument(
		"--gmx-bin",
		default="gmx",
		help="GROMACS executable, e.g. gmx, gmx_mpi, or /path/to/gmx.",
	)
	p.add_argument(
		"--skip-existing",
		action="store_true",
		help="Reuse existing XVG outputs instead of rerunning gmx sasa.",
	)
	p.add_argument(
		"--dry-run",
		action="store_true",
		help="Write index files and command list, but do not run gmx sasa.",
	)
	p.add_argument(
		"--skip-residue-sasa",
		action="store_true",
		help="Skip atom-level per-residue SASA output.",
	)
	p.add_argument(
		"--med4-cp12-hmm-map",
		default=None,
		help="Optional MED4 CP12 residue-to-HMM table with CP12_Chain, ResID, HMM_col columns.",
	)
	p.add_argument(
		"--hm2-cp12-hmm-map",
		default=None,
		help="Optional HM2 CP12 residue-to-HMM table with CP12_Chain, ResID, HMM_col columns.",
	)
	p.add_argument(
		"--cp12-combined-hmm-map",
		default=None,
		help=(
			"Optional family-wide CP12 residue-to-HMM table with Accession, "
			"seq_residue_index, and HMM_col columns. Preferred for MED4/HM2 "
			"cross-system residue pairing."
		),
	)
	p.add_argument("--med4-cp12-accession", default="WP_011131859.1")
	p.add_argument("--hm2-cp12-accession", default="YP_004323596.1")
	p.add_argument(
		"--skip-empty-xtc",
		action="store_true",
		default=True,
		help="Skip missing or 0-byte trajectories in dynamic mode.",
	)
	return p.parse_args()


def parse_pdb_atoms(pdb: Path) -> List[AtomRecord]:
	atoms: List[AtomRecord] = []
	with pdb.open() as handle:
		for line in handle:
			if not line.startswith(("ATOM", "HETATM")):
				continue
			atoms.append(
				AtomRecord(
					index1=len(atoms) + 1,
					name=line[12:16].strip(),
					resname=line[17:20].strip(),
					chain=line[21].strip() or "_",
					resid=int(line[22:26]),
				)
			)
	return atoms


def write_ndx(atoms: Sequence[AtomRecord], ndx: Path) -> Dict[str, List[int]]:
	groups: Dict[str, List[int]] = {
		"CP12": [],
		"PRK": [],
		"GAP2": [],
		"NAD": [],
	}
	for atom in atoms:
		if atom.chain in CP12_CHAINS:
			groups["CP12"].append(atom.index1)
		elif atom.chain in PRK_CHAINS:
			groups["PRK"].append(atom.index1)
		elif atom.chain in GAP2_CHAINS:
			groups["GAP2"].append(atom.index1)
		elif atom.chain in NAD_CHAINS or atom.resname.upper() in {"NAD", "NAP"}:
			groups["NAD"].append(atom.index1)

	groups["ProteinComplex"] = sorted(groups["CP12"] + groups["PRK"] + groups["GAP2"])
	groups["FullComplex"] = sorted(groups["ProteinComplex"] + groups["NAD"])
	groups["CP12_PRK"] = sorted(groups["CP12"] + groups["PRK"])
	groups["CP12_GAP2"] = sorted(groups["CP12"] + groups["GAP2"])
	groups["PRK_GAP2"] = sorted(groups["PRK"] + groups["GAP2"])
	groups["CP12_PRK_GAP2"] = sorted(groups["CP12"] + groups["PRK"] + groups["GAP2"])

	missing = [name for name in ("CP12", "PRK", "GAP2") if not groups[name]]
	if missing:
		raise ValueError(f"Missing required chain-derived groups in {ndx}: {', '.join(missing)}")

	with ndx.open("w") as handle:
		for name, indices in groups.items():
			if not indices:
				continue
			handle.write(f"[ {name} ]\n")
			for start in range(0, len(indices), 15):
				handle.write(" ".join(f"{idx:8d}" for idx in indices[start : start + 15]) + "\n")
			handle.write("\n")
	return groups


def discover_inputs(args: argparse.Namespace) -> List[ReplicateInput]:
	root = Path(args.analysis_dir)
	reps: List[ReplicateInput] = []
	for strain in args.strains:
		for system in args.systems:
			system = str(system).zfill(4) if str(system).isdigit() else str(system)
			for rep in args.replicates:
				rep_dir = root / strain / system / rep
				pdb = rep_dir / args.pdb_name
				tpr = rep_dir / args.tpr_name
				xtc = rep_dir / args.xtc_name
				if not pdb.exists():
					continue
				reps.append(
					ReplicateInput(
						strain=strain,
						system=system,
						condition=CONDITIONS.get(system, system),
						replicate=rep,
						rep_dir=rep_dir,
						pdb=pdb,
						tpr=tpr,
						xtc=xtc,
					)
				)
	return reps


def selection(group: str) -> str:
	return f'group "{group}"'


def gmx_sasa_command(
	args: argparse.Namespace,
	rep: ReplicateInput,
	input_file: Path,
	ndx: Path,
	job: SasaJob,
	xvg: Path,
	analysis: str,
) -> List[str]:
	cmd = [
		args.gmx_bin,
		"sasa",
		"-s",
		str(rep.tpr),
		"-f",
		str(input_file),
		"-n",
		str(ndx),
		"-surface",
		selection(job.surface_group),
		"-output",
		selection(job.output_group),
		"-o",
		str(xvg),
		"-probe",
		str(args.probe_radius),
	]
	if analysis == "static" and args.static_source == "trajectory":
		cmd.extend(["-b", str(args.static_time_ps), "-e", str(args.static_time_ps)])
	elif analysis == "dynamic":
		if args.begin_ps is not None:
			cmd.extend(["-b", str(args.begin_ps)])
		if args.end_ps is not None:
			cmd.extend(["-e", str(args.end_ps)])
		if args.dt_ps is not None:
			cmd.extend(["-dt", str(args.dt_ps)])
	return cmd


def gmx_atom_sasa_command(
	args: argparse.Namespace,
	rep: ReplicateInput,
	input_file: Path,
	ndx: Path,
	surface_group: str,
	output_group: str,
	xvg: Path,
	analysis: str,
) -> List[str]:
	cmd = [
		args.gmx_bin,
		"sasa",
		"-s",
		str(rep.tpr),
		"-f",
		str(input_file),
		"-n",
		str(ndx),
		"-surface",
		selection(surface_group),
		"-output",
		selection(output_group),
		"-oa",
		str(xvg),
		"-probe",
		str(args.probe_radius),
	]
	if analysis == "static" and args.static_source == "trajectory":
		cmd.extend(["-b", str(args.static_time_ps), "-e", str(args.static_time_ps)])
	elif analysis == "dynamic":
		if args.begin_ps is not None:
			cmd.extend(["-b", str(args.begin_ps)])
		if args.end_ps is not None:
			cmd.extend(["-e", str(args.end_ps)])
		if args.dt_ps is not None:
			cmd.extend(["-dt", str(args.dt_ps)])
	return cmd


def command_to_text(cmd: Sequence[str]) -> str:
	return " ".join(shlex.quote(part) for part in cmd)


def parse_xvg(xvg: Path, metric: str) -> pd.DataFrame:
	rows: List[Tuple[float, float]] = []
	with xvg.open() as handle:
		for line in handle:
			if not line.strip() or line.startswith(("#", "@")):
				continue
			parts = line.split()
			if len(parts) < 2:
				continue
			rows.append((float(parts[0]), float(parts[1]) * 100.0))
	return pd.DataFrame(rows, columns=["time_ps", metric])


def parse_atom_xvg(xvg: Path, output_indices: Sequence[int]) -> pd.DataFrame:
	raw_rows: List[Tuple[int, float, float]] = []
	with xvg.open() as handle:
		for line in handle:
			if not line.strip() or line.startswith(("#", "@")):
				continue
			parts = line.split()
			if len(parts) < 2:
				continue
			atom_token = int(float(parts[0]))
			mean_sasa_A2 = float(parts[1]) * 100.0
			sd_sasa_A2 = float(parts[2]) * 100.0 if len(parts) >= 3 else float("nan")
			raw_rows.append((atom_token, mean_sasa_A2, sd_sasa_A2))
	if not raw_rows:
		return pd.DataFrame(columns=["atom_index1", "mean_atom_sasa_A2", "sd_atom_sasa_A2"])

	output_indices = list(output_indices)
	tokens = [row[0] for row in raw_rows]
	output_index_set = set(output_indices)
	token_set = set(tokens)

	if len(raw_rows) == len(output_indices) and tokens == list(range(1, len(output_indices) + 1)):
		atom_indices = [output_indices[token - 1] for token in tokens]
	elif token_set.issubset(output_index_set):
		atom_indices = tokens
	elif token_set.issubset(set(range(1, len(output_indices) + 1))):
		atom_indices = [output_indices[token - 1] for token in tokens]
	elif len(raw_rows) == len(output_indices):
		atom_indices = output_indices
	else:
		raise ValueError(
			f"Cannot map atom SASA rows from {xvg} to output group atoms "
			f"(rows={len(raw_rows)}, output_atoms={len(output_indices)})"
		)

	return pd.DataFrame(
		{
			"atom_index1": atom_indices,
			"mean_atom_sasa_A2": [row[1] for row in raw_rows],
			"sd_atom_sasa_A2": [row[2] for row in raw_rows],
		}
	)


def protein_for_chain(chain: str) -> str:
	if chain in CP12_CHAINS:
		return "CP12"
	if chain in PRK_CHAINS:
		return "PRK"
	if chain in GAP2_CHAINS:
		return "GAP2"
	if chain in NAD_CHAINS:
		return "NAD"
	return "Other"


def load_cp12_hmm_maps(args: argparse.Namespace) -> Dict[Tuple[str, str, int], object]:
	maps: Dict[Tuple[str, str, int], object] = {}
	if args.cp12_combined_hmm_map:
		path = Path(args.cp12_combined_hmm_map)
		if not path.exists():
			raise FileNotFoundError(f"Combined CP12 HMM map not found: {path}")
		df = pd.read_csv(path, sep="\t", low_memory=False)
		required = {"Accession", "seq_residue_index", "HMM_col"}
		missing = required - set(df.columns)
		if missing:
			raise ValueError(f"Combined CP12 HMM map {path} is missing columns: {', '.join(sorted(missing))}")
		accessions = {
			"med4": args.med4_cp12_accession,
			"hm2": args.hm2_cp12_accession,
		}
		for strain, accession in accessions.items():
			sub = df[df["Accession"].astype(str) == accession]
			if sub.empty:
				raise ValueError(f"Accession {accession} for {strain} was not found in {path}")
			for _, row in sub.iterrows():
				hmm_col = row["HMM_col"]
				if pd.isna(hmm_col) or str(hmm_col).upper() == "NA":
					continue
				resid = int(row["seq_residue_index"])
				for chain in CP12_CHAINS:
					maps[(strain, chain, resid)] = int(float(hmm_col))
	for strain, map_path in (("med4", args.med4_cp12_hmm_map), ("hm2", args.hm2_cp12_hmm_map)):
		if not map_path:
			continue
		path = Path(map_path)
		if not path.exists():
			raise FileNotFoundError(f"CP12 HMM map not found for {strain}: {path}")
		df = pd.read_csv(path, sep="\t")
		required = {"CP12_Chain", "ResID", "HMM_col"}
		missing = required - set(df.columns)
		if missing:
			raise ValueError(f"CP12 HMM map {path} is missing columns: {', '.join(sorted(missing))}")
		for _, row in df.iterrows():
			hmm_col = row["HMM_col"]
			if pd.isna(hmm_col) or str(hmm_col).upper() == "NA":
				continue
			maps.setdefault((strain, str(row["CP12_Chain"]), int(row["ResID"])), int(float(hmm_col)))
	return maps


def atom_lookup(
	atoms: Sequence[AtomRecord],
	strain: str,
	cp12_hmm_maps: Dict[Tuple[str, str, int], object],
) -> pd.DataFrame:
	rows = []
	residue_index: Dict[Tuple[str, int, str], int] = {}
	for atom in atoms:
		if atom.chain in NAD_CHAINS:
			continue
		protein = protein_for_chain(atom.chain)
		if protein not in {"CP12", "PRK", "GAP2"}:
			continue
		key = (atom.chain, atom.resid, atom.resname)
		if key not in residue_index:
			residue_index[key] = len(residue_index) + 1
		hmm_col = cp12_hmm_maps.get((strain, atom.chain, atom.resid)) if protein == "CP12" else None
		residue_match_key = (
			f"CP12:HMM{hmm_col}" if hmm_col is not None else f"{protein}:RES{atom.resid}"
		)
		rows.append(
			{
				"atom_index1": atom.index1,
				"residue_index1": residue_index[key],
				"atom_name": atom.name,
				"chain": atom.chain,
				"resid": atom.resid,
				"resname": atom.resname,
				"protein": protein,
				"hmm_col": hmm_col,
				"site": f"{protein}:{atom.chain}:{atom.resname}{atom.resid}",
				"site_no_chain": f"{protein}:{atom.resname}{atom.resid}",
				"residue_match_key": residue_match_key,
			}
		)
	return pd.DataFrame(rows)


def residue_lookup(
	atoms: Sequence[AtomRecord],
	strain: str,
	cp12_hmm_maps: Dict[Tuple[str, str, int], object],
) -> pd.DataFrame:
	lookup = atom_lookup(atoms, strain, cp12_hmm_maps)
	if lookup.empty:
		return lookup
	residue_cols = [
		"residue_index1",
		"chain",
		"resid",
		"resname",
		"protein",
		"hmm_col",
		"site",
		"site_no_chain",
		"residue_match_key",
	]
	return lookup[residue_cols].drop_duplicates().reset_index(drop=True)


def aggregate_atom_sasa_to_residues(
	atom_sasa: pd.DataFrame,
	lookup: pd.DataFrame,
	rep: ReplicateInput,
	analysis: str,
	pdb: Path,
	input_file: Path,
) -> pd.DataFrame:
	if atom_sasa.empty or lookup.empty:
		return pd.DataFrame()
	merged = atom_sasa.merge(lookup, on="atom_index1", how="inner")
	if merged.empty:
		return pd.DataFrame()
	id_cols = [
		"residue_index1",
		"chain",
		"resid",
		"resname",
		"protein",
		"hmm_col",
		"site",
		"site_no_chain",
		"residue_match_key",
	]
	out = (
		merged.groupby(id_cols, dropna=False)
		.agg(
			n_residue_atoms=("atom_index1", "count"),
			mean_residue_sasa_A2=("mean_atom_sasa_A2", "sum"),
			sd_residue_sasa_A2=("sd_atom_sasa_A2", lambda x: float(np.sqrt(np.nansum(np.square(x))))),
		)
		.reset_index()
	)
	out.insert(0, "replicate", rep.replicate)
	out.insert(0, "condition", rep.condition)
	out.insert(0, "system", rep.system)
	out.insert(0, "strain", rep.strain)
	out.insert(0, "analysis", analysis)
	out["pdb"] = str(pdb)
	out["trajectory_or_structure"] = str(input_file)
	return out


def run_jobs_for_analysis(
	args: argparse.Namespace,
	rep: ReplicateInput,
	analysis: str,
	input_file: Path,
	out_rep_dir: Path,
	ndx: Path,
	groups: Dict[str, List[int]],
	commands: List[str],
) -> pd.DataFrame:
	xvg_dir = out_rep_dir / analysis / "xvg"
	xvg_dir.mkdir(parents=True, exist_ok=True)
	metric_tables = []
	for job in SASA_JOBS:
		xvg = xvg_dir / f"{job.metric}.xvg"
		cmd = gmx_sasa_command(args, rep, input_file, ndx, job, xvg, analysis)
		commands.append(command_to_text(cmd))
		if args.dry_run:
			continue
		if not (args.skip_existing and xvg.exists() and xvg.stat().st_size > 0):
			subprocess.run(cmd, check=True)
		metric_tables.append(parse_xvg(xvg, job.metric))

	if args.dry_run:
		for job in ATOM_COMPONENT_JOBS:
			atom_xvg = xvg_dir.parent / "atom_xvg" / f"{job.metric}_atoms.xvg"
			cmd = gmx_atom_sasa_command(
				args, rep, input_file, ndx, job.surface_group, job.output_group, atom_xvg, analysis
			)
			commands.append(command_to_text(cmd))
		return pd.DataFrame()

	if not metric_tables:
		return pd.DataFrame()

	merged = metric_tables[0]
	for table in metric_tables[1:]:
		merged = merged.merge(table, on="time_ps", how="outer")

	atom_xvg_dir = out_rep_dir / analysis / "atom_xvg"
	atom_xvg_dir.mkdir(parents=True, exist_ok=True)
	for job in ATOM_COMPONENT_JOBS:
		atom_xvg = atom_xvg_dir / f"{job.metric}_atoms.xvg"
		cmd = gmx_atom_sasa_command(
			args, rep, input_file, ndx, job.surface_group, job.output_group, atom_xvg, analysis
		)
		commands.append(command_to_text(cmd))
		if not (args.skip_existing and atom_xvg.exists() and atom_xvg.stat().st_size > 0):
			subprocess.run(cmd, check=True)
		atom_sasa = parse_atom_xvg(atom_xvg, groups[job.surface_group])
		output_atoms = set(groups[job.output_group])
		atom_sasa = atom_sasa[atom_sasa["atom_index1"].isin(output_atoms)]
		merged[job.metric] = float(atom_sasa["mean_atom_sasa_A2"].sum()) if not atom_sasa.empty else np.nan

	for metric, (a_metric, b_metric, pair_metric) in BSA_FORMULAS.items():
		merged[metric] = merged[a_metric] + merged[b_metric] - merged[pair_metric]
	merged.insert(0, "analysis", analysis)
	merged.insert(1, "strain", rep.strain)
	merged.insert(2, "system", rep.system)
	merged.insert(3, "condition", rep.condition)
	merged.insert(4, "replicate", rep.replicate)
	merged["pdb"] = str(rep.pdb)
	merged["trajectory_or_structure"] = str(input_file)
	return merged


def run_residue_job_for_analysis(
	args: argparse.Namespace,
	rep: ReplicateInput,
	analysis: str,
	input_file: Path,
	out_rep_dir: Path,
	ndx: Path,
	lookup: pd.DataFrame,
	groups: Dict[str, List[int]],
	commands: List[str],
) -> pd.DataFrame:
	if args.skip_residue_sasa:
		return pd.DataFrame()
	residue_dir = out_rep_dir / analysis / "atom_xvg"
	residue_dir.mkdir(parents=True, exist_ok=True)
	xvg = residue_dir / "protein_complex_residue_atom_sasa_A2.xvg"
	cmd = gmx_atom_sasa_command(
		args, rep, input_file, ndx, "ProteinComplex", "ProteinComplex", xvg, analysis
	)
	commands.append(command_to_text(cmd))
	if args.dry_run:
		return pd.DataFrame()
	if not (args.skip_existing and xvg.exists() and xvg.stat().st_size > 0):
		subprocess.run(cmd, check=True)
	atom_sasa = parse_atom_xvg(xvg, groups["ProteinComplex"])
	residue = aggregate_atom_sasa_to_residues(atom_sasa, lookup, rep, analysis, rep.pdb, input_file)
	if not residue.empty:
		residue.to_csv(
			out_rep_dir / analysis / "protein_complex_residue_sasa.tsv",
			sep="\t",
			index=False,
		)
	return residue


def process_replicate(
	args: argparse.Namespace,
	rep: ReplicateInput,
	out_dir: Path,
	cp12_hmm_maps: Dict[Tuple[str, str, int], object],
) -> pd.DataFrame:
	if not rep.tpr.exists():
		raise FileNotFoundError(f"Missing TPR required by gmx sasa: {rep.tpr}")
	if args.mode in {"static", "both"} and args.static_source == "trajectory":
		if not rep.xtc.exists() or (args.skip_empty_xtc and rep.xtc.stat().st_size == 0):
			raise FileNotFoundError(
				f"Static trajectory source requested but trajectory is missing/empty: {rep.xtc}"
			)
	if args.mode in {"dynamic", "both"}:
		if not rep.xtc.exists() or (args.skip_empty_xtc and rep.xtc.stat().st_size == 0):
			print(f"SKIP dynamic empty/missing trajectory: {rep.xtc}", flush=True)
			if args.mode == "dynamic":
				return pd.DataFrame()

	out_rep_dir = out_dir / rep.strain / rep.system / rep.replicate
	out_rep_dir.mkdir(parents=True, exist_ok=True)
	atoms = parse_pdb_atoms(rep.pdb)
	lookup = atom_lookup(atoms, rep.strain, cp12_hmm_maps)
	ndx = out_rep_dir / "complex_sasa_groups.ndx"
	groups = write_ndx(atoms, ndx)
	group_sizes = pd.DataFrame(
		[{"group": name, "n_atoms": len(indices)} for name, indices in groups.items()]
	)
	group_sizes.to_csv(out_rep_dir / "complex_sasa_group_sizes.tsv", sep="\t", index=False)

	commands: List[str] = []
	frames: List[pd.DataFrame] = []
	residue_frames: List[pd.DataFrame] = []
	if args.mode in {"static", "both"}:
		static_input = rep.xtc if args.static_source == "trajectory" else rep.pdb
		frames.append(run_jobs_for_analysis(args, rep, "static", static_input, out_rep_dir, ndx, groups, commands))
		residue_frames.append(
			run_residue_job_for_analysis(args, rep, "static", static_input, out_rep_dir, ndx, lookup, groups, commands)
		)
	if args.mode in {"dynamic", "both"} and rep.xtc.exists():
		frames.append(run_jobs_for_analysis(args, rep, "dynamic", rep.xtc, out_rep_dir, ndx, groups, commands))
		residue_frames.append(
			run_residue_job_for_analysis(args, rep, "dynamic", rep.xtc, out_rep_dir, ndx, lookup, groups, commands)
		)

	(out_rep_dir / "complex_sasa_commands.sh").write_text(
		"#!/usr/bin/env bash\nset -euo pipefail\n\n" + "\n".join(commands) + "\n"
	)
	residue_frames = [frame for frame in residue_frames if not frame.empty]
	if residue_frames:
		pd.concat(residue_frames, ignore_index=True).to_csv(
			out_rep_dir / "protein_complex_residue_sasa.tsv",
			sep="\t",
			index=False,
		)
	frames = [frame for frame in frames if not frame.empty]
	if not frames:
		return pd.DataFrame()
	wide = pd.concat(frames, ignore_index=True)
	wide.to_csv(out_rep_dir / "complex_sasa_bsa_raw_wide.tsv", sep="\t", index=False)
	return wide


def collect_residue_tables(out_dir: Path) -> pd.DataFrame:
	tables = []
	for path in out_dir.glob("*/*/*/protein_complex_residue_sasa.tsv"):
		tables.append(pd.read_csv(path, sep="\t", dtype={"system": str}))
	if not tables:
		return pd.DataFrame()
	out = pd.concat(tables, ignore_index=True)
	if "system" in out.columns:
		out["system"] = out["system"].astype(str).str.zfill(4)
	return out


def wide_to_long(raw_wide: pd.DataFrame) -> pd.DataFrame:
	if raw_wide.empty:
		return pd.DataFrame()
	id_cols = [
		"analysis",
		"strain",
		"system",
		"condition",
		"replicate",
		"time_ps",
		"pdb",
		"trajectory_or_structure",
	]
	metric_cols = [col for col in raw_wide.columns if col.endswith("_A2")]
	return raw_wide.melt(
		id_vars=id_cols,
		value_vars=metric_cols,
		var_name="metric",
		value_name="value_A2",
	)


def sem(values: pd.Series) -> float:
	vals = pd.to_numeric(values, errors="coerce").dropna()
	if len(vals) <= 1:
		return 0.0 if len(vals) == 1 else float("nan")
	return float(vals.std(ddof=1) / math.sqrt(len(vals)))


def summarize(raw_long: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
	if raw_long.empty:
		return pd.DataFrame(), pd.DataFrame()
	rep_cols = ["analysis", "strain", "system", "condition", "replicate", "metric"]
	rep_summary = (
		raw_long.groupby(rep_cols, dropna=False)["value_A2"]
		.agg(["count", "mean", "std", "median", "min", "max"])
		.reset_index()
		.rename(columns={"count": "n_frames", "mean": "mean_value_A2", "std": "sd_value_A2"})
	)
	group_cols = ["analysis", "strain", "system", "condition", "metric"]
	group_summary = (
		rep_summary.groupby(group_cols, dropna=False)["mean_value_A2"]
		.agg(["count", "mean", "std", "median", "min", "max"])
		.reset_index()
		.rename(
			columns={
				"count": "n_replicates",
				"mean": "mean_value_A2",
				"std": "replicate_sd_value_A2",
			}
		)
	)
	group_summary["replicate_sem_value_A2"] = group_summary[
		"replicate_sd_value_A2"
	].fillna(0.0) / np.sqrt(group_summary["n_replicates"].clip(lower=1))
	return rep_summary, group_summary


def make_delta_tables(group_summary: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
	if group_summary.empty:
		return pd.DataFrame(), pd.DataFrame()
	strain_rows = []
	key_cols = ["analysis", "system", "condition", "metric"]
	for key, sub in group_summary.groupby(key_cols, dropna=False):
		med4 = sub[sub["strain"] == "med4"]
		hm2 = sub[sub["strain"] == "hm2"]
		if med4.empty or hm2.empty:
			continue
		strain_rows.append(
			{
				**dict(zip(key_cols, key)),
				"med4_mean_value_A2": float(med4.iloc[0]["mean_value_A2"]),
				"hm2_mean_value_A2": float(hm2.iloc[0]["mean_value_A2"]),
				"hm2_minus_med4_A2": float(hm2.iloc[0]["mean_value_A2"])
				- float(med4.iloc[0]["mean_value_A2"]),
			}
		)

	redox_rows = []
	key_cols = ["analysis", "strain", "metric"]
	for key, sub in group_summary.groupby(key_cols, dropna=False):
		reduced = sub[sub["condition"] == "reduced"]
		oxidized = sub[sub["condition"] == "oxidized"]
		if reduced.empty or oxidized.empty:
			continue
		redox_rows.append(
			{
				**dict(zip(key_cols, key)),
				"reduced_mean_value_A2": float(reduced.iloc[0]["mean_value_A2"]),
				"oxidized_mean_value_A2": float(oxidized.iloc[0]["mean_value_A2"]),
				"reduced_minus_oxidized_A2": float(reduced.iloc[0]["mean_value_A2"])
				- float(oxidized.iloc[0]["mean_value_A2"]),
			}
		)
	return pd.DataFrame(strain_rows), pd.DataFrame(redox_rows)


def summarize_residue_sasa(residue: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
	if residue.empty:
		return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
	id_cols = [
		"analysis",
		"strain",
		"system",
		"condition",
		"protein",
		"chain",
		"resid",
		"resname",
		"site",
		"site_no_chain",
		"residue_match_key",
		"hmm_col",
	]
	group = (
		residue.groupby(id_cols, dropna=False)["mean_residue_sasa_A2"]
		.agg(["count", "mean", "std", "median", "min", "max"])
		.reset_index()
		.rename(
			columns={
				"count": "n_replicate_chain_observations",
				"mean": "mean_residue_sasa_A2",
				"std": "replicate_chain_sd_residue_sasa_A2",
			}
		)
	)
	group["replicate_chain_sem_residue_sasa_A2"] = group[
		"replicate_chain_sd_residue_sasa_A2"
	].fillna(0.0) / np.sqrt(group["n_replicate_chain_observations"].clip(lower=1))

	strain_rows = []
	key_cols = ["analysis", "system", "condition", "protein", "residue_match_key"]
	for key, sub in group.groupby(key_cols, dropna=False):
		med4 = sub[sub["strain"] == "med4"]
		hm2 = sub[sub["strain"] == "hm2"]
		if med4.empty or hm2.empty:
			continue
		med4_row = med4.iloc[0]
		hm2_row = hm2.iloc[0]
		strain_rows.append(
			{
				**dict(zip(key_cols, key)),
				"resid": hm2_row.get("resid", med4_row.get("resid", "")),
				"hmm_col": hm2_row.get("hmm_col", med4_row.get("hmm_col", "")),
				"med4_resname": med4_row.get("resname", ""),
				"hm2_resname": hm2_row.get("resname", ""),
				"med4_resid": med4_row.get("resid", ""),
				"hm2_resid": hm2_row.get("resid", ""),
				"med4_site_no_chain": med4_row.get("site_no_chain", ""),
				"hm2_site_no_chain": hm2_row.get("site_no_chain", ""),
				"med4_mean_residue_sasa_A2": float(med4_row["mean_residue_sasa_A2"]),
				"hm2_mean_residue_sasa_A2": float(hm2_row["mean_residue_sasa_A2"]),
				"hm2_minus_med4_residue_sasa_A2": float(hm2_row["mean_residue_sasa_A2"])
				- float(med4_row["mean_residue_sasa_A2"]),
			}
		)

	redox_rows = []
	key_cols = ["analysis", "strain", "protein", "residue_match_key"]
	for key, sub in group.groupby(key_cols, dropna=False):
		reduced = sub[sub["condition"] == "reduced"]
		oxidized = sub[sub["condition"] == "oxidized"]
		if reduced.empty or oxidized.empty:
			continue
		red_row = reduced.iloc[0]
		ox_row = oxidized.iloc[0]
		redox_rows.append(
			{
				**dict(zip(key_cols, key)),
				"resid": ox_row.get("resid", red_row.get("resid", "")),
				"hmm_col": ox_row.get("hmm_col", red_row.get("hmm_col", "")),
				"reduced_resname": red_row.get("resname", ""),
				"oxidized_resname": ox_row.get("resname", ""),
				"reduced_resid": red_row.get("resid", ""),
				"oxidized_resid": ox_row.get("resid", ""),
				"reduced_site_no_chain": red_row.get("site_no_chain", ""),
				"oxidized_site_no_chain": ox_row.get("site_no_chain", ""),
				"reduced_mean_residue_sasa_A2": float(red_row["mean_residue_sasa_A2"]),
				"oxidized_mean_residue_sasa_A2": float(ox_row["mean_residue_sasa_A2"]),
				"reduced_minus_oxidized_residue_sasa_A2": float(red_row["mean_residue_sasa_A2"])
				- float(ox_row["mean_residue_sasa_A2"]),
			}
		)
	return group, pd.DataFrame(strain_rows), pd.DataFrame(redox_rows)


def write_notes(out_dir: Path, args: argparse.Namespace, elapsed: float) -> None:
	(out_dir / "complex_sasa_bsa_notes.txt").write_text(
		"Complex SASA and BSA notes\n"
		"==========================\n\n"
		f"mode: {args.mode}\n"
		f"static_source: {args.static_source}\n"
		f"static_time_ps: {args.static_time_ps}\n"
		f"begin_ps: {args.begin_ps}\n"
		f"end_ps: {args.end_ps}\n"
		f"dt_ps: {args.dt_ps}\n"
		f"probe_radius_nm: {args.probe_radius}\n"
		f"gmx_bin: {args.gmx_bin}\n"
		f"dry_run: {args.dry_run}\n"
		f"skip_residue_sasa: {args.skip_residue_sasa}\n"
		f"cp12_combined_hmm_map: {args.cp12_combined_hmm_map}\n"
		f"med4_cp12_accession: {args.med4_cp12_accession}\n"
		f"hm2_cp12_accession: {args.hm2_cp12_accession}\n"
		f"med4_cp12_hmm_map: {args.med4_cp12_hmm_map}\n"
		f"hm2_cp12_hmm_map: {args.hm2_cp12_hmm_map}\n"
		f"elapsed_seconds: {elapsed:.1f}\n\n"
		"Chain groups: CP12=C,H,M,P; PRK=B,G,J,L; GAP2=A,D,E,F,I,K,N,O; NAD=Q.\n"
		"Component-in-complex and residue SASA are summed from atom-level gmx sasa -oa output.\n"
		"BSA is calculated from isolated component SASA and pair-complex SASA.\n"
		"For viral-substitution tests, compare HM2 minus MED4 within the same redox state,\n"
		"especially the oxidized 0003 reference condition.\n"
	)


def main() -> int:
	args = parse_args()
	start = time.time()
	out_dir = Path(args.out_dir)
	tables_dir = out_dir / "tables"
	tables_dir.mkdir(parents=True, exist_ok=True)

	reps = discover_inputs(args)
	if not reps:
		raise SystemExit("No replicate inputs found.")
	cp12_hmm_maps = load_cp12_hmm_maps(args)
	if cp12_hmm_maps:
		print(f"Loaded {len(cp12_hmm_maps)} CP12 residue-to-HMM mappings", flush=True)

	raw_frames = []
	for rep in reps:
		print(f"Processing {rep.strain} {rep.system} {rep.replicate}", flush=True)
		raw_frames.append(process_replicate(args, rep, out_dir, cp12_hmm_maps))

	if args.dry_run:
		write_notes(out_dir, args, time.time() - start)
		print(f"Dry run complete. Command files are under {out_dir}", flush=True)
		return 0

	raw_wide = (
		pd.concat([frame for frame in raw_frames if not frame.empty], ignore_index=True)
		if raw_frames
		else pd.DataFrame()
	)
	raw_wide.to_csv(tables_dir / "complex_sasa_bsa_raw_wide.tsv", sep="\t", index=False)
	raw_long = wide_to_long(raw_wide)
	raw_long.to_csv(tables_dir / "complex_sasa_bsa_raw_long.tsv", sep="\t", index=False)
	rep_summary, group_summary = summarize(raw_long)
	rep_summary.to_csv(tables_dir / "complex_sasa_bsa_replicate_summary.tsv", sep="\t", index=False)
	group_summary.to_csv(tables_dir / "complex_sasa_bsa_group_summary.tsv", sep="\t", index=False)
	strain_delta, redox_delta = make_delta_tables(group_summary)
	strain_delta.to_csv(tables_dir / "complex_sasa_bsa_hm2_minus_med4.tsv", sep="\t", index=False)
	redox_delta.to_csv(tables_dir / "complex_sasa_bsa_reduced_minus_oxidized.tsv", sep="\t", index=False)
	residue = collect_residue_tables(out_dir)
	residue.to_csv(tables_dir / "protein_complex_residue_sasa_raw.tsv", sep="\t", index=False)
	residue_summary, residue_strain_delta, residue_redox_delta = summarize_residue_sasa(residue)
	residue_summary.to_csv(
		tables_dir / "protein_complex_residue_sasa_group_summary.tsv",
		sep="\t",
		index=False,
	)
	residue_strain_delta.to_csv(
		tables_dir / "protein_complex_residue_sasa_hm2_minus_med4.tsv",
		sep="\t",
		index=False,
	)
	residue_redox_delta.to_csv(
		tables_dir / "protein_complex_residue_sasa_reduced_minus_oxidized.tsv",
		sep="\t",
		index=False,
	)
	write_notes(out_dir, args, time.time() - start)
	print(f"Wrote outputs to {out_dir}", flush=True)
	print(f"Elapsed seconds: {time.time() - start:.1f}", flush=True)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
