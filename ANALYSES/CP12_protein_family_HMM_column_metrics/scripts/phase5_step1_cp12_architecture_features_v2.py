#!/usr/bin/env python

import argparse
import pandas as pd
import ast


def parse_list_safe(x):
	if pd.isna(x):
		return []
	if isinstance(x, list):
		return x
	if isinstance(x, str):
		x = x.strip()
		if x in ["", "NA"]:
			return []
		try:
			return list(ast.literal_eval(x))
		except Exception:
			return []
	return []

def define_primary_group(row):
	if row["SourceType"] == "Reference":
		return "Reference"
	if row["SourceType"] == "Metagenome":
		if row["SequenceStatus"] == "full":
			return "Metagenome_full"
		if row["SequenceStatus"] == "partial":
			return "Metagenome_partial"
	return "Other"


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--meta", required=True, help="Combined CP12 metadata TSV")
	parser.add_argument("--out", required=True, help="Output architecture feature TSV")
	args = parser.parse_args()

	df = pd.read_csv(args.meta, sep="\t")

	required_cols = [
		"Accession",
		"SourceType",
		"SequenceStatus",
		"NumCys",
		"CysPositions",
		"CysSpacings",
		"AcidicFraction"
	]

	missing = [c for c in required_cols if c not in df.columns]
	if missing:
		raise ValueError(f"Missing required columns: {missing}")

	df["CysPositions_parsed"] = df["CysPositions"].apply(parse_list_safe)
	df["CysSpacings_parsed"] = df["CysSpacings"].apply(parse_list_safe)

	df["NumCysPairs"] = df["CysSpacings_parsed"].apply(
		lambda x: len(x) if isinstance(x, list) else 0
	)

	df["MeanCysSpacing"] = df["CysSpacings_parsed"].apply(
		lambda x: sum(x) / len(x) if isinstance(x, list) and len(x) > 0 else float("nan")
	)

	df["PrimaryGroup"] = df.apply(define_primary_group, axis=1)

	print("\nPrimaryGroup counts:")
	print(df["PrimaryGroup"].value_counts())

	out_cols = [
		"Accession",
		"PrimaryGroup",
		"SourceType",
		"SequenceStatus",
		"NumCys",
		"NumCysPairs",
		"MeanCysSpacing",
		"AcidicFraction"
	]

	df[out_cols].to_csv(args.out, sep="\t", index=False)

	print(f"\nWrote architecture table to: {args.out}")

if __name__ == "__main__":
	main()