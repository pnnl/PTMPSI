#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import AlignIO
import pandas as pd
import sys


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def parse_args():
	p = argparse.ArgumentParser(
		description="Derive per-HMM-column amino acid counts from Stockholm alignment"
	)
	p.add_argument("--sto", required=True, help="Stockholm alignment file")
	p.add_argument("--meta", required=True, help="Combined CP12 metadata table")
	p.add_argument("--out", required=True, help="Output AA count table")
	return p.parse_args()


def extract_rf_string(sto_file):
	rf_lines = []
	with open(sto_file) as f:
		for line in f:
			if line.startswith("#=GC RF"):
				rf_lines.append(line.rstrip().split()[-1])
	if not rf_lines:
		raise ValueError("Stockholm file missing #=GC RF annotation")

	rf = "".join(rf_lines)
	return rf


def main():
	args = parse_args()

	aln = AlignIO.read(args.sto, "stockholm")
	nseq = len(aln)

	rf = extract_rf_string(args.sto)
	aln_len = len(rf)

	hmm_cols = [i for i, c in enumerate(rf) if c != "."]
	n_hmm_cols = len(hmm_cols)

	print(f"[INFO] RF length: {aln_len}")
	print(f"[INFO] HMM match columns: {n_hmm_cols}")
	print(f"[INFO] Aligned sequences: {nseq}", file=sys.stderr)

	# Map alignment index → HMM column number (1-based)
	aln_to_hmm = {}
	hmm_index = 0
	for i in range(aln_len):
		if rf[i] != ".":
			hmm_index += 1
			aln_to_hmm[i] = hmm_index

	meta = pd.read_csv(args.meta, sep="\t", dtype=str)

	meta["PrimaryGroup"] = meta["PrimaryGroup"].fillna("Unknown")
	meta["PrimaryGroupCombined"] = meta["PrimaryGroupCombined"].fillna("Unknown")
	meta["Lineage_consensus_cleaned"] = meta["Lineage_consensus_cleaned"].fillna("Unknown")

	meta_map = meta.set_index("Accession")

	aa_counts = defaultdict(int)
	gap_counts = defaultdict(int)
	numseq = defaultdict(int)
	numnongap = defaultdict(int)

	for record in aln:
		acc = record.id
		seq = str(record.seq)

		if acc not in meta_map.index:
			continue

		groups = [
			("GLOBAL", "GLOBAL"),
			("PrimaryGroup", meta_map.at[acc, "PrimaryGroup"]),
			("PrimaryGroupCombined", meta_map.at[acc, "PrimaryGroupCombined"]),
			("Lineage", meta_map.at[acc, "Lineage_consensus_cleaned"]),
		]

		for aln_pos, aa in enumerate(seq):
			if aln_pos not in aln_to_hmm:
				continue

			hmm_col = aln_to_hmm[aln_pos]

			for gtype, gval in groups:
				key = (hmm_col, gtype, gval)
				numseq[key] += 1

				if aa == "-" or aa == ".":
					gap_counts[key] += 1
				else:
					numnongap[key] += 1
					if aa in AMINO_ACIDS:
						aa_counts[(hmm_col, gtype, gval, aa)] += 1

	rows = []

	for (hmm_col, gtype, gval), nseq_val in numseq.items():
		nongap = numnongap.get((hmm_col, gtype, gval), 0)
		gap = gap_counts.get((hmm_col, gtype, gval), 0)

		row = {
			"HMM_col": hmm_col,
			"GroupType": gtype,
			"GroupValue": gval,
			"NumSeq": nseq_val,
			"NumNonGap": nongap,
			"Gap_count": gap,
		}

		for aa in AMINO_ACIDS:
			row[f"{aa}_count"] = aa_counts.get((hmm_col, gtype, gval, aa), 0)

		rows.append(row)

	df_out = pd.DataFrame(rows)
	df_out.sort_values(["HMM_col", "GroupType", "GroupValue"], inplace=True)

	df_out.to_csv(args.out, sep="\t", index=False)

	print(f"[INFO] Written per-column AA counts to {args.out}", file=sys.stderr)


if __name__ == "__main__":
	main()