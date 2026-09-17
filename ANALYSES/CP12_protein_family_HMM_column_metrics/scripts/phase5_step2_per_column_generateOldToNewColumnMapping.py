#!/usr/bin/env python3

import argparse
from Bio import AlignIO
import pandas as pd


def parse_args():
	p = argparse.ArgumentParser(
		description="Export per-sequence residue index to HMM column mapping from a Stockholm alignment"
	)
	p.add_argument("--sto", required=True, help="Stockholm alignment file")
	p.add_argument("--meta", required=True, help="Metadata table with Accession column")
	p.add_argument("--out", required=True, help="Output mapping TSV")
	return p.parse_args()


def extract_rf_string(sto_file):
	rf_lines = []
	with open(sto_file) as f:
		for line in f:
			if line.startswith("#=GC RF"):
				rf_lines.append(line.rstrip().split()[-1])

	if not rf_lines:
		raise ValueError("Stockholm file missing #=GC RF annotation")

	return "".join(rf_lines)


def main():
	args = parse_args()

	aln = AlignIO.read(args.sto, "stockholm")
	rf = extract_rf_string(args.sto)

	if len(rf) != aln.get_alignment_length():
		raise ValueError(
			f"RF length {len(rf)} does not match alignment length {aln.get_alignment_length()}"
		)

	# -------------------------------
	# Build alignment position -> HMM column mapping
	# Only RF positions not equal to "." are HMM match columns
	# HMM columns are 1-based
	# -------------------------------
	aln_to_hmm = {}
	hmm_col = 0
	for i, c in enumerate(rf):
		if c != ".":
			hmm_col += 1
			aln_to_hmm[i] = hmm_col

	print(f"[INFO] RF length: {len(rf)}")
	print(f"[INFO] HMM match columns: {hmm_col}")
	print(f"[INFO] Aligned sequences: {len(aln)}")


	meta = pd.read_csv(args.meta, sep="\t", dtype=str)

	if "Accession" not in meta.columns:
		raise ValueError(
			f"Metadata file missing required column 'Accession'. Found columns: {list(meta.columns)}"
		)

	meta = meta.fillna("")

	if meta["Accession"].duplicated().any():
		n_dup_accessions = meta.loc[meta["Accession"].duplicated(keep=False), "Accession"].nunique()
		print(f"[WARN] Duplicate Accessions found in metadata, keeping first row for {n_dup_accessions} accessions")

	meta = meta.drop_duplicates(subset="Accession", keep="first").copy()
	meta_map = meta.set_index("Accession").to_dict(orient="index")

	# -------------------------------
	# Generate mapping rows
	# Output only non-gap residues at HMM match columns
	# seq_residue_index is 1-based in the ungapped sequence
	# -------------------------------
	rows = []

	for record in aln:
		acc = record.id
		seq = str(record.seq)

		seq_residue_index = 0
		meta_row = meta_map.get(acc, {})

		for aln_pos, aa in enumerate(seq):
			if aa != "-" and aa != ".":
				seq_residue_index += 1

			if aln_pos not in aln_to_hmm:
				continue

			if aa == "-" or aa == ".":
				continue

			row = {
				"Accession": acc,
				"aln_pos_1based": aln_pos + 1,
				"seq_residue_index": seq_residue_index,
				"AA": aa,
				"HMM_col": aln_to_hmm[aln_pos],
			}

			row.update(meta_row)
			rows.append(row)

	df = pd.DataFrame(rows)

	if df.empty:
		print("[WARN] No mapping rows generated")
	else:
		df.sort_values(
			["Accession", "HMM_col", "seq_residue_index", "aln_pos_1based"],
			inplace=True
		)

	df.to_csv(args.out, sep="\t", index=False)

	print(f"[INFO] Wrote mapping table to {args.out}")
	print(f"[INFO] Rows: {len(df)}")
	if not df.empty:
		print(f"[INFO] Unique sequences: {df['Accession'].nunique()}")
		print(f"[INFO] Unique HMM columns: {df['HMM_col'].nunique()}")


if __name__ == "__main__":
	main()