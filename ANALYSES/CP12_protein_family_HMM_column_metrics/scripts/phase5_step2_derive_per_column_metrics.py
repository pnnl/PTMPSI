#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from math import log2

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def shannon_entropy(counts):
	total = sum(counts)
	if total == 0:
		return np.nan
	probs = [c / total for c in counts if c > 0]
	return -sum(p * log2(p) for p in probs)

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--counts", required=True,
					help="figure2_per_column_aa_counts.tsv")
	ap.add_argument("--interfaces", required=True,
					help="MED4_vs_HM2_cp12_interface_hmm.tsv")
	ap.add_argument("--out", required=True)
	ap.add_argument("--cys_thresh", type=float, default=0.7)
	args = ap.parse_args()

	df = pd.read_csv(args.counts, sep="\t")

	# Load interface annotation
	iface = pd.read_csv(args.interfaces, sep="\t")

	keep_iface_cols = [
		"HMM_col",
		"MED4_0003_Interface",
		"MED4_0003_Type",
		"HM2_0003_Interface",
		"HM2_0003_Type",
		"ChangeClass"
	]

	iface = iface[keep_iface_cols].drop_duplicates()
	df = df.merge(iface, on="HMM_col", how="left")

	records = []

	for _, r in df.iterrows():
		num_non_gap = r["NumNonGap"]
		num_seq = r["NumSeq"]

		aa_counts = [r[f"{aa}_count"] for aa in AA_LIST]

		# Core metrics
		c_frac = np.nan
		acidic_frac = np.nan
		entropy = np.nan
		num_unique = 0

		if num_non_gap > 0:
			c_frac = r["C_count"] / num_non_gap
			acidic_frac = (r["D_count"] + r["E_count"]) / num_non_gap
			entropy = shannon_entropy(aa_counts)
			num_unique = sum(c > 0 for c in aa_counts)

		occupancy = num_non_gap / num_seq if num_seq > 0 else np.nan

		metrics = {
			"C_fraction": c_frac,
			"IsConservedCys": c_frac >= args.cys_thresh if not np.isnan(c_frac) else False,
			"NumUnique": num_unique,
			"Shannon": entropy,
			"AcidicFraction": acidic_frac,
			"Occupancy": occupancy
		}

		for metric, value in metrics.items():
			records.append({
				"HMM_col": r["HMM_col"],
				"GroupType": r["GroupType"],
				"GroupValue": r["GroupValue"],
				"Metric": metric,
				"Value": value,
				"NumSeq": num_seq,
				"NumNonGap": num_non_gap,
				"MED4_0003_Interface": r["MED4_0003_Interface"],
				"MED4_0003_Type": r["MED4_0003_Type"],
				"HM2_0003_Interface": r["HM2_0003_Interface"],
				"HM2_0003_Type": r["HM2_0003_Type"],
				"ChangeClass": r["ChangeClass"]
			})

	out = pd.DataFrame(records)
	out.to_csv(args.out, sep="\t", index=False)

	print("Wrote:", args.out)
	print("Metrics:", out["Metric"].unique())

if __name__ == "__main__":
	main()