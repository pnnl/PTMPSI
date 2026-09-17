#!/usr/bin/env python3
import argparse
from collections import OrderedDict
from Bio.PDB import PDBParser

def parse_args():
	p = argparse.ArgumentParser(
		description="Phase 1: flexible chain-to-protein mapping for PTM-PSI complexes"
	)
	p.add_argument("--pdb", required=True, help="Input PDB, e.g. solute_0_chains.pdb")
	p.add_argument("--system", required=True, help="System name, e.g. MED4_0003")
	p.add_argument(
		"--protein-map",
		required=True,
		nargs="+",
		help="Protein to chain mapping, e.g. CP12:C,H,M,P PRK:A,B GAPDH:D,E,F"
	)
	p.add_argument("--out", default="system_chain_map.tsv", help="Output TSV")
	return p.parse_args()

def parse_protein_map(entries):
	chain_to_protein = {}
	for entry in entries:
		if ":" not in entry:
			raise ValueError(f"Invalid protein-map entry: {entry}")
		protein, chains = entry.split(":", 1)
		for c in chains.split(","):
			c = c.strip()
			if not c:
				continue
			if c in chain_to_protein:
				raise ValueError(f"Chain {c} assigned twice")
			chain_to_protein[c] = protein
	return chain_to_protein

def count_chain_residues(structure):
	chain_to_res = OrderedDict()
	for model in structure:
		for chain in model:
			seen = set()
			for res in chain.get_residues():
				if res.id[0] != " ":
					continue
				key = (res.id[1], res.id[2])
				seen.add(key)
			chain_to_res[chain.id] = len(seen)
		break
	return chain_to_res

def main():
	args = parse_args()
	chain_to_protein = parse_protein_map(args.protein_map)

	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("complex", args.pdb)
	chain_lengths = count_chain_residues(structure)

	with open(args.out, "w") as f:
		f.write("System\tChain\tProtein\tResidueCount\n")
		for chain_id, nres in chain_lengths.items():
			protein = chain_to_protein.get(chain_id, "UNASSIGNED")
			f.write(f"{args.system}\t{chain_id}\t{protein}\t{nres}\n")

	print(f"Wrote {args.out}")
	print("Protein assignments:")
	for c, p in sorted(chain_to_protein.items()):
		print(f"  {c} -> {p}")

if __name__ == "__main__":
	main()