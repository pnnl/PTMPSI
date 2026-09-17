#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict

def parse_args():
	p = argparse.ArgumentParser(
		description="Phase 4.3 + 4.4: Project CP12 interfaces to HMM columns and collapse across chains"
	)
	p.add_argument("--interfaces", required=True, help="Phase 3 CP12 interface TSV")
	p.add_argument("--resid-map", required=True, help="Residue-to-HMM column mapping TSV")
	p.add_argument("--system", required=True, help="System name (e.g. MED4_0003)")
	p.add_argument("--out-prefix", required=True, help="Output prefix")
	return p.parse_args()

def main():
	args = parse_args()

	# Load inputs
	iface = pd.read_csv(args.interfaces, sep="\t")
	resmap = pd.read_csv(args.resid_map, sep="\t")

	# Phase 3 output is already CP12-only, no filtering needed
	iface = iface.copy()

	# Merge interface residues with HMM mapping
	merged = iface.merge(
		resmap,
		left_on=["System", "CP12_Chain", "CP12_ResID"],
		right_on=["System", "CP12_Chain", "ResID"],
		how="left"
	)

	merged = merged[merged["HMM_col"].notna()].copy()

	# =========================
	# Phase 4.3: chain-resolved
	# =========================
	out_chain = f"{args.out_prefix}_cp12_interface_hmm_chain.tsv"
	chain_cols = [
		"System",
		"HMM_col",
		"CP12_Chain",
		"CP12_ResID",
		"AA",
		"InterfaceType"
	]
	merged[chain_cols].to_csv(out_chain, sep="\t", index=False)

	print(f"Wrote chain-resolved interface table: {out_chain}")

	# =========================
	# Phase 4.4: collapse chains
	# =========================
	rows = []

	for (system, hmm_col), sub in merged.groupby(["System", "HMM_col"]):
		types = set(sub["InterfaceType"])
		chains = sorted(set(sub["CP12_Chain"]))

		if len(types) == 1:
			iface_type = list(types)[0]
		else:
			iface_type = "Mixed"

		rows.append({
			"System": system,
			"HMM_col": hmm_col,
			"Interface_present": 1,
			"InterfaceType_collapsed": iface_type,
			"ChainsSupporting": ",".join(chains),
			"NumChainsSupporting": len(chains)
		})

	collapsed = pd.DataFrame(rows)

	out_collapsed = f"{args.out_prefix}_cp12_interface_hmm_collapsed.tsv"
	collapsed.to_csv(out_collapsed, sep="\t", index=False)

	print(f"Wrote collapsed interface table: {out_collapsed}")
	print("Collapsed interface type counts:")
	print(collapsed["InterfaceType_collapsed"].value_counts().to_string())

if __name__ == "__main__":
	main()