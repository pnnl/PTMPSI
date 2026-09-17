#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
from collections import defaultdict, Counter

import pandas as pd


AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA20_SET = set(AA20)

# gap-like characters in Stockholm alignments
GAP_SET = set(["-", ".", "~"])

# ChemClass mapping for Panel D style features
AA_TO_CHEM = {
	"C": "Cysteine",
	"K": "Positive", "R": "Positive", "H": "Positive",
	"D": "Negative", "E": "Negative",
	"S": "Polar", "T": "Polar", "N": "Polar", "Q": "Polar",
	"A": "Hydrophobic", "V": "Hydrophobic", "I": "Hydrophobic", "L": "Hydrophobic",
	"M": "Hydrophobic", "F": "Hydrophobic", "W": "Hydrophobic", "Y": "Hydrophobic",
	"G": "Gly_Pro", "P": "Gly_Pro"
}


def safe_str(x):
	"""
	Convert values (including pandas objects) into safe strings.
	Handles:
	- scalar
	- Series (take first element)
	- DataFrame (take [0,0])
	"""
	if x is None:
		return ""
	if isinstance(x, pd.Series):
		if x.size == 0:
			return ""
		return safe_str(x.iloc[0])
	if isinstance(x, pd.DataFrame):
		if x.shape[0] == 0 or x.shape[1] == 0:
			return ""
		return safe_str(x.iloc[0, 0])
	if pd.isna(x):
		return ""
	return str(x)


def parse_stockholm(sto_path: Path):
	"""
	Parse Stockholm alignment into {aln_name: aligned_sequence} and RF string.
	Concatenates wrapped blocks.
	"""
	seq_chunks = defaultdict(list)
	rf_chunks = []

	with sto_path.open("r", encoding="utf-8", errors="replace") as f:
		for line in f:
			line = line.rstrip("\n")
			if not line:
				continue

			# RF line (match-state annotation)
			if line.startswith("#=GC RF"):
				parts = line.split()
				if len(parts) >= 3:
					rf_chunks.append(parts[2])
				continue

			# other annotations
			if line.startswith("#"):
				continue

			if line.strip() == "//":
				break

			parts = line.split()
			if len(parts) < 2:
				continue

			name = parts[0]
			seq = parts[1]
			seq_chunks[name].append(seq)

	seqs = {k: "".join(v) for k, v in seq_chunks.items()}
	if not seqs:
		raise ValueError(f"No sequences parsed from Stockholm: {sto_path}")

	aln_lens = {len(s) for s in seqs.values()}
	if len(aln_lens) != 1:
		raise ValueError(f"Inconsistent alignment lengths: {sorted(aln_lens)}")
	aln_len = next(iter(aln_lens))

	rf = "".join(rf_chunks)
	if rf == "":
		raise ValueError("No #=GC RF line found in Stockholm file, cannot map match columns to HMM positions.")
	if len(rf) != aln_len:
		raise ValueError(f"RF length ({len(rf)}) does not match alignment length ({aln_len}).")

	# Match-state columns, RF uses '.' for inserts in many HMMER outputs.
	match_col_idx = [i for i, ch in enumerate(rf) if ch not in [".", "-"]]

	if len(match_col_idx) == 0:
		raise ValueError("RF parsing produced 0 match columns, check RF encoding in the Stockholm file.")

	return seqs, aln_len, rf, match_col_idx


def resolve_accession(aln_name: str, meta_acc_set: set):
	raw = safe_str(aln_name)

	candidates = []
	candidates.append(raw)

	if "|" in raw:
		candidates.append(raw.split("|")[-1])
	if "/" in raw:
		candidates.append(raw.split("/")[0])
	if "." in raw:
		candidates.append(raw.split(".")[0])

	for c in candidates:
		if c in meta_acc_set:
			return c

	return None


def shannon_entropy_from_counts(counts: Counter, total: int):
	if total <= 0:
		return float("nan")

	h = 0.0
	for aa, c in counts.items():
		if c <= 0:
			continue
		p = c / total
		h -= p * math.log(p)
	return h


def add_primary_groups(meta: pd.DataFrame):
	def _primary_group(row):
		if safe_str(row["SourceType"]) == "Reference":
			return "reference"
		ss = safe_str(row["SequenceStatus"]).lower()
		if ss == "full":
			return "metagenome_full"
		if ss == "partial":
			return "metagenome_partial"
		return "metagenome_unknown"

	m = meta.copy()
	m["PrimaryGroup"] = m.apply(_primary_group, axis=1)
	m["PrimaryGroupCombined"] = m["SourceType"].apply(lambda x: "reference" if safe_str(x) == "Reference" else "metagenome_combined")
	m["LineageCleaned"] = m["Lineage_consensus_cleaned"].astype(str)
	return m


def build_accession_long(seqs_by_acc: dict, match_col_idx: list, meta_by_acc: pd.DataFrame):
	"""
	Return per-accession per HMM match-position long table:
	Accession, HMM_col, AA_raw, AA, IsGap, ChemClass, plus metadata labels.
	"""
	rows = []
	for acc, seq in seqs_by_acc.items():
		seq_list = list(seq)

		for i_match, col0 in enumerate(match_col_idx):
			hmm_col = i_match + 1

			aa_raw = seq_list[col0]
			is_gap = aa_raw in GAP_SET
			aa = aa_raw.upper()
			chem = AA_TO_CHEM.get(aa, "")

			row = {
				"Accession": acc,
				"HMM_col": hmm_col,
				"AA_raw": aa_raw,
				"AA": aa,
				"IsGap": int(is_gap),
				"ChemClass": chem
			}

			if acc in meta_by_acc.index:
				md = meta_by_acc.loc[acc]
				if isinstance(md, pd.DataFrame):
					md = md.iloc[0]

				row["SourceType"] = safe_str(md.get("SourceType", ""))
				row["SequenceStatus"] = safe_str(md.get("SequenceStatus", ""))
				row["PrimaryGroup"] = safe_str(md.get("PrimaryGroup", ""))
				row["PrimaryGroupCombined"] = safe_str(md.get("PrimaryGroupCombined", ""))
				row["LineageCleaned"] = safe_str(md.get("LineageCleaned", ""))
				row["TaxID"] = safe_str(md.get("TaxID", ""))
			else:
				row["SourceType"] = ""
				row["SequenceStatus"] = ""
				row["PrimaryGroup"] = ""
				row["PrimaryGroupCombined"] = ""
				row["LineageCleaned"] = ""
				row["TaxID"] = ""

			rows.append(row)

	return pd.DataFrame(rows)


def compute_group_tables(accession_long: pd.DataFrame, out_prefix: str, out_dir: Path):
	group_defs = [
		("GLOBAL", lambda df: pd.Series(["GLOBAL"] * len(df), index=df.index)),
		("PrimaryGroup", lambda df: df["PrimaryGroup"]),
		("PrimaryGroupCombined", lambda df: df["PrimaryGroupCombined"]),
		("LineageCleaned", lambda df: df["LineageCleaned"])
	]

	aa_counts_rows = []
	stats_rows = []

	for group_type, group_func in group_defs:
		tmp = accession_long.copy()
		tmp["GroupType"] = group_type
		tmp["GroupValue"] = group_func(tmp).astype(str)

		if group_type != "GLOBAL":
			tmp = tmp[tmp["GroupValue"].fillna("").astype(str) != ""]

		for (gv, hmm_col), g in tmp.groupby(["GroupValue", "HMM_col"], sort=False):
			accs = g["Accession"].unique().tolist()
			numseq = len(accs)

			numnongap = int((g["IsGap"] == 0).sum())
			gap_count = numseq - numnongap

			aa_ct = Counter()
			for aa in g.loc[g["IsGap"] == 0, "AA"].tolist():
				if aa in AA20_SET:
					aa_ct[aa] += 1

			aa_row = {
				"HMM_col": int(hmm_col),
				"GroupType": group_type,
				"GroupValue": gv,
				"NumSeq": numseq,
				"NumNonGap": numnongap,
				"Gap_count": gap_count
			}
			for aa in AA20:
				aa_row[f"{aa}_count"] = int(aa_ct.get(aa, 0))
			aa_counts_rows.append(aa_row)

			occupancy = (numnongap / numseq) if numseq > 0 else float("nan")
			c_frac = (aa_ct.get("C", 0) / numnongap) if numnongap > 0 else 0.0
			shannon = shannon_entropy_from_counts(aa_ct, numnongap)

			stats_base = {
				"HMM_col": int(hmm_col),
				"GroupType": group_type,
				"GroupValue": gv,
				"NumSeq": numseq,
				"NumNonGap": numnongap,
				"MED4_0003_Interface": pd.NA,
				"MED4_0003_Type": pd.NA,
				"HM2_0003_Interface": pd.NA,
				"HM2_0003_Type": pd.NA,
				"ChangeClass": pd.NA
			}
			stats_rows.append({**stats_base, "Metric": "Shannon", "Value": shannon})
			stats_rows.append({**stats_base, "Metric": "Occupancy", "Value": occupancy})
			stats_rows.append({**stats_base, "Metric": "C_fraction", "Value": c_frac})

	aa_counts_df = pd.DataFrame(aa_counts_rows)
	stats_df = pd.DataFrame(stats_rows)

	aa_out = out_dir / f"{out_prefix}_aa_counts_step2v6.txt"
	st_out = out_dir / f"{out_prefix}_stats_long_step2v6.txt"

	aa_counts_df.to_csv(aa_out, sep="\t", index=False)
	stats_df.to_csv(st_out, sep="\t", index=False)

	print(f"Wrote: {aa_out}")
	print(f"Wrote: {st_out}")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--sto", required=True)
	ap.add_argument("--meta", required=True)
	ap.add_argument("--out_prefix", default="figure2")
	ap.add_argument("--out_dir", default=None)
	args = ap.parse_args()

	sto_path = Path(args.sto).resolve()
	meta_path = Path(args.meta).resolve()

	out_prefix_path = Path(args.out_prefix)
	if args.out_dir is None:
		if out_prefix_path.parent == Path("."):
			out_dir = Path.cwd()
			out_prefix = out_prefix_path.name
		else:
			out_dir = out_prefix_path.parent.resolve()
			out_prefix = out_prefix_path.name
	else:
		out_dir = Path(args.out_dir).resolve()
		out_prefix = out_prefix_path.name
	out_dir.mkdir(parents=True, exist_ok=True)

	meta = pd.read_csv(meta_path, sep="\t", dtype=str)

	required = ["Accession", "SourceType", "SequenceStatus", "Lineage_consensus_cleaned"]
	for col in required:
		if col not in meta.columns:
			raise ValueError(f"Metadata missing required column: {col}")

	meta["Accession"] = meta["Accession"].astype(str)

	dup_n = int(meta["Accession"].duplicated().sum())
	if dup_n > 0:
		print(f"Warning: metadata contains {dup_n} duplicate Accession rows, keeping the first occurrence for join.")
	meta = meta.drop_duplicates(subset=["Accession"], keep="first").reset_index(drop=True)

	meta = add_primary_groups(meta)
	meta_by_acc = meta.set_index("Accession", drop=False)

	seqs_by_name, aln_len, rf, match_col_idx = parse_stockholm(sto_path)

	meta_acc_set = set(meta["Accession"])
	name_to_acc = {}
	unmatched = 0

	for aln_name in seqs_by_name.keys():
		acc = resolve_accession(aln_name, meta_acc_set)
		if acc is None:
			unmatched += 1
			continue
		name_to_acc[aln_name] = acc

	if len(name_to_acc) == 0:
		raise ValueError("No alignment names could be matched to metadata Accession.")

	if unmatched > 0:
		print(f"Warning: {unmatched} alignment names could not be matched to metadata, skipped.")

	seqs_by_acc = {}
	dup = 0
	for aln_name, acc in name_to_acc.items():
		if acc in seqs_by_acc:
			dup += 1
			continue
		seqs_by_acc[acc] = seqs_by_name[aln_name]

	if dup > 0:
		print(f"Warning: {dup} duplicate accessions encountered in alignment names, kept first occurrence.")

	accession_long = build_accession_long(seqs_by_acc, match_col_idx, meta_by_acc)

	acc_out = out_dir / f"{out_prefix}_accession_per_column_step2v6.txt"
	accession_long.to_csv(acc_out, sep="\t", index=False)
	print(f"Wrote: {acc_out}")

	compute_group_tables(accession_long, out_prefix, out_dir)

	print("Unique accessions:", accession_long["Accession"].nunique())
	print("Raw alignment length:", aln_len)
	print("HMM match columns (RF-based):", len(match_col_idx))


if __name__ == "__main__":
	main()