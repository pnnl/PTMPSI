#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import AlignIO

def parse_args():
	p = argparse.ArgumentParser()
	p.add_argument("--sto", required=True, help="HMM-aligned Stockholm file")
	p.add_argument("--interfaces", required=True, help="Interface annotation TSV")
	p.add_argument("--out", required=True, help="Output TSV")
	return p.parse_args()

def main():
	args = parse_args()

	aln = AlignIO.read(args.sto, "stockholm")
	aln_len = aln.get_alignment_length()
	seq_ids = [rec.id for rec in aln]

	rows = []
	for i in range(aln_len):
		col = aln[:, i]
		residues = [aa for aa in col if aa not in "-."]
		uniq = set(residues)
		rows.append({
			"HMM_col": i + 1,
			"NumSeq": len(residues),
			"NumUnique": len(uniq),
			"Diversity": len(uniq) / len(residues) if residues else 0
		})

	div = pd.DataFrame(rows)

	iface = pd.read_csv(args.interfaces, sep="\t")

	# ---- MINIMAL FIX STARTS HERE ----
	med4_iface = [c for c in iface.columns if c.startswith("MED4") and c.endswith("_Interface")]
	med4_type  = [c for c in iface.columns if c.startswith("MED4") and c.endswith("_Type")]
	hm2_iface  = [c for c in iface.columns if c.startswith("HM2")  and c.endswith("_Interface")]
	hm2_type   = [c for c in iface.columns if c.startswith("HM2")  and c.endswith("_Type")]

	keep_cols = ["HMM_col", "ChangeClass"] + med4_iface + med4_type + hm2_iface + hm2_type

	missing = [c for c in keep_cols if c not in iface.columns]
	if missing:
		raise ValueError(f"Missing expected interface columns: {missing}")

	iface = iface[keep_cols]

	out = div.merge(iface, on="HMM_col", how="left")
	out.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
	main()