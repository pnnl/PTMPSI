#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

import numpy as np

import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS

def parse_args():
	p = argparse.ArgumentParser(
		description="Phase 2: extract residue contacts between proteins in PTM-PSI complexes"
	)
	p.add_argument("--pdb", required=True, help="Input PDB with chain IDs")
	p.add_argument("--chain-map", required=True, help="Phase 1 system_chain_map.tsv")
	p.add_argument("--system", required=True, help="System name")
	p.add_argument("--focal-protein", required=True, help="Protein of interest, e.g. CP12")
	p.add_argument(
		"--partner-proteins",
		required=True,
		nargs="+",
		help="One or more partner proteins, e.g. PRK GAPDH"
	)
	p.add_argument(
		"--distance-cutoff",
		type=float,
		default=5.0,
		help="Distance cutoff in Å (default: 5.0)"
	)
	p.add_argument("--out", default="phase2_contacts.tsv")
	return p.parse_args()

def load_chain_map(tsv):
	chain_to_protein = {}
	with open(tsv) as f:
		reader = csv.DictReader(f, delimiter="\t")
		for row in reader:
			chain_to_protein[row["Chain"]] = row["Protein"]
	return chain_to_protein

def _atom_is_heavy(atom):
	element = str(getattr(atom, "element", "") or "").upper().strip()
	name = str(getattr(atom, "name", "") or "").upper().strip()
	type_name = str(getattr(atom, "type", "") or "").upper().strip()

	if element in ("H", "D"):
		return False
	if name.startswith("H"):
		return False
	if type_name.startswith("H"):
		return False
	return True

def _chain_id(atom):
	chain = str(getattr(atom, "chainID", "") or "").strip()
	if chain:
		return chain
	return str(getattr(atom, "segid", "") or "").strip()

def main():
	args = parse_args()
	chain_map = load_chain_map(args.chain_map)
	partner_set = set(args.partner_proteins)

	u = mda.Universe(args.pdb)
	residue_atoms = defaultdict(list)

	for atom in u.atoms:
		chain = _chain_id(atom)
		if not chain:
			continue
		protein = chain_map.get(chain, "UNASSIGNED")
		if protein != args.focal_protein and protein not in partner_set:
			continue
		if not _atom_is_heavy(atom):
			continue

		res = atom.residue
		key = (chain, protein, int(res.resid), str(res.resname).strip())
		residue_atoms[key].append(atom.index)

	focal_residues = []
	partner_residues = []

	for key, idxs in residue_atoms.items():
		chain, protein, resid, resname = key
		coords = u.atoms[idxs].positions.astype(np.float32)
		if coords.shape[0] == 0:
			continue
		rec = {
			"chain": chain,
			"protein": protein,
			"resid": resid,
			"resname": resname,
			"coords": coords,
		}
		if protein == args.focal_protein:
			focal_residues.append(rec)
		else:
			partner_residues.append(rec)

	if not focal_residues:
		raise RuntimeError(f"No heavy-atom residues found for focal protein: {args.focal_protein}")
	if not partner_residues:
		raise RuntimeError("No heavy-atom residues found for partner proteins")

	partner_coord_blocks = []
	partner_key_by_atom = []

	for p in partner_residues:
		partner_coord_blocks.append(p["coords"])
		partner_key_by_atom.extend([p] * p["coords"].shape[0])

	partner_coords = np.vstack(partner_coord_blocks).astype(np.float32)

	all_coords = np.vstack([u.atoms.positions.astype(np.float32), partner_coords])
	mins = all_coords.min(axis=0) - (args.distance_cutoff + 1.0)
	maxs = all_coords.max(axis=0) + (args.distance_cutoff + 1.0)
	box_lengths = np.maximum(maxs - mins, np.array([1.0, 1.0, 1.0], dtype=np.float32))
	box = np.array([box_lengths[0], box_lengths[1], box_lengths[2], 90.0, 90.0, 90.0], dtype=np.float32)

	partner_coords_shift = partner_coords - mins
	ns = FastNS(args.distance_cutoff, partner_coords_shift, box, pbc=False)

	with open(args.out, "w") as out:
		out.write(
			"System\tFocalProtein\tFocalChain\tFocalResID\tFocalResName\t"
			"PartnerProtein\tPartnerChain\tPartnerResID\tPartnerResName\t"
			"MinDistance\n"
		)

		for f in focal_residues:
			f_coords_shift = (f["coords"] - mins).astype(np.float32)
			results = ns.search(f_coords_shift)
			pairs = results.get_pairs()
			dists = results.get_pair_distances()

			if pairs.size == 0:
				continue

			min_by_partner = {}
			for pair, dist in zip(pairs, dists):
				p_atom_idx = int(pair[1])
				p = partner_key_by_atom[p_atom_idx]
				p_key = (p["chain"], p["protein"], p["resid"], p["resname"])
				dist = float(dist)
				if p_key not in min_by_partner or dist < min_by_partner[p_key]:
					min_by_partner[p_key] = dist

			for (p_chain, p_prot, p_resid, p_resname), min_dist in sorted(min_by_partner.items()):
				out.write(
					f"{args.system}\t{f['protein']}\t{f['chain']}\t{f['resid']}\t{f['resname']}\t"
					f"{p_prot}\t{p_chain}\t{p_resid}\t{p_resname}\t"
					f"{min_dist:.3f}\n"
				)

	print(f"Wrote {args.out}")
	print(f"Distance cutoff: {args.distance_cutoff} Å")
	print(f"Focal protein: {args.focal_protein}")
	print(f"Partners: {', '.join(args.partner_proteins)}")

if __name__ == "__main__":
	main()