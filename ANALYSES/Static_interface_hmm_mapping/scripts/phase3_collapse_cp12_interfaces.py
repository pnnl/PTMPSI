#!/usr/bin/env python3
import argparse
import pandas as pd

def parse_args():
	p = argparse.ArgumentParser(
		description="Phase 3: collapse CP12 residue contacts into interface annotations"
	)
	p.add_argument(
		"--contacts",
		required=True,
		help="Phase 2 contact table (TSV)"
	)
	p.add_argument(
		"--system",
		required=True,
		help="System name (e.g. MED4_0003 or HM2_0003)"
	)
	p.add_argument(
		"--focal-protein",
		default="CP12",
		help="Focal protein name (default: CP12)"
	)
	p.add_argument(
		"--out-prefix",
		default="phase3_cp12_interface",
		help="Output prefix"
	)
	return p.parse_args()

def main():
	args = parse_args()

	df = pd.read_csv(args.contacts, sep="\t")

	# Filter focal protein rows only
	df = df[df["FocalProtein"] == args.focal_protein].copy()

	if df.empty:
		raise RuntimeError("No focal protein residues found in contact table")

	# Group by CP12 residue
	grouped = df.groupby(
		["FocalChain", "FocalResID", "FocalResName"],
		as_index=False
	)

	records = []

	for _, g in grouped:
		partners = sorted(set(g["PartnerProtein"]))
		min_dist = g["MinDistance"].min()
		n_contacts = g.shape[0]

		if len(partners) == 1:
			interface_type = f"{partners[0]}_only"
		else:
			interface_type = "+".join(partners)

		records.append({
			"System": args.system,
			"CP12_Chain": g["FocalChain"].iloc[0],
			"CP12_ResID": g["FocalResID"].iloc[0],
			"CP12_ResName": g["FocalResName"].iloc[0],
			"InterfaceType": interface_type,
			"PartnerProteins": ",".join(partners),
			"MinDistance": round(min_dist, 3),
			"NumContacts": n_contacts,
			"IsInterface": True
		})

	interface_df = pd.DataFrame(records)

	# Optional: add non-interface residues later via sequence map
	interface_df.sort_values("CP12_ResID", inplace=True)

	out_tsv = f"{args.out_prefix}.tsv"
	interface_df.to_csv(out_tsv, sep="\t", index=False)

	print(f"Wrote interface table: {out_tsv}")
	print(f"Residues with contacts: {interface_df.shape[0]}")
	print("Interface type counts:")
	print(interface_df["InterfaceType"].value_counts().to_string())

if __name__ == "__main__":
	main()