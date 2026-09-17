#!/usr/bin/env python

import pandas as pd
import argparse

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--meta", required=True)
	ap.add_argument("--out", required=True)
	args = ap.parse_args()

	df = pd.read_csv(args.meta, sep="\t")

	# PrimaryGroup
	def assign_primary(row):
		if row["SourceType"] == "Reference":
			return "Reference"
		if row["SourceType"] == "Metagenome":
			if row["SequenceStatus"] == "full":
				return "Metagenome_full"
			else:
				return "Metagenome_partial"
		return "Other"

	df["PrimaryGroup"] = df.apply(assign_primary, axis=1)

	# PrimaryGroupCombined
	df["PrimaryGroupCombined"] = df["SourceType"].map(
		lambda x: "Reference" if x == "Reference" else "Metagenome_combined"
	)

	out = df[[
		"Accession",
		"PrimaryGroup",
		"PrimaryGroupCombined",
		"SourceType",
		"SequenceStatus",
		"Lineage_consensus_cleaned"
	]].drop_duplicates()

	out.to_csv(args.out, sep="\t", index=False)

	print("PrimaryGroup counts:")
	print(out["PrimaryGroup"].value_counts())
	print("\nPrimaryGroupCombined counts:")
	print(out["PrimaryGroupCombined"].value_counts())

if __name__ == "__main__":
	main()