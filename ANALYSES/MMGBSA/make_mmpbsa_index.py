#!/usr/bin/env python
"""
Create GROMACS index file from PDB with specific chain groups
Usage: python make_index_groups.py reference.pdb output.ndx
"""

import sys

def read_pdb_atoms(pdb_file):
    """Read PDB and return list of (atom_index, chain_id, atom_name, resname)"""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # PDB format: columns 7-11 for atom serial number, column 22 for chain
                atom_num = int(line[6:11].strip())
                chain_id = line[21].strip()
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                atoms.append((atom_num, chain_id, atom_name, resname))
    return atoms

def write_index_group(f, group_name, atom_indices):
    """Write index group in GROMACS .ndx format"""
    f.write(f"[ {group_name} ]\n")
    for i, idx in enumerate(atom_indices):
        f.write(f"{idx:6d} ")
        if (i + 1) % 15 == 0:  # 15 indices per line
            f.write("\n")
    if len(atom_indices) % 15 != 0:
        f.write("\n")
    f.write("\n")

def make_index_file(pdb_file, ndx_file):
    """Create index file with CP12, PRK, GAP2, NAD, LIGAND, and RECEPTOR groups"""
    
    # Define chain assignments
    CP12_CHAINS = ["C", "H", "M", "P"]
    PRK_CHAINS = ["B", "G", "J", "L"]
    GAP2_CHAINS = ["A", "D", "E", "F", "I", "K", "N", "O"]
    NAD_CHAINS = ["Q"]
    
    # Read PDB
    atoms = read_pdb_atoms(pdb_file)
    print(f"Read {len(atoms)} atoms from {pdb_file}")
    
    # Collect atom indices for each group
    cp12_atoms = []
    prk_atoms = []
    gap2_atoms = []
    nad_atoms = []
    
    for atom_num, chain_id, atom_name, resname in atoms:
        if chain_id in CP12_CHAINS:
            cp12_atoms.append(atom_num)
        elif chain_id in PRK_CHAINS:
            prk_atoms.append(atom_num)
        elif chain_id in GAP2_CHAINS:
            gap2_atoms.append(atom_num)
        elif chain_id in NAD_CHAINS:
            nad_atoms.append(atom_num)
    
    # RECEPTOR = PRK + GAP2 + NAD
    receptor_atoms = sorted(prk_atoms + gap2_atoms + nad_atoms)
    
    # Write index file
    with open(ndx_file, 'w') as f:
        write_index_group(f, "CP12", cp12_atoms)
        write_index_group(f, "PRK", prk_atoms)
        write_index_group(f, "GAP2", gap2_atoms)
        write_index_group(f, "NAD", nad_atoms)
        write_index_group(f, "LIGAND", cp12_atoms)
        write_index_group(f, "RECEPTOR", receptor_atoms)
    
    # Print summary
    print(f"\nIndex groups created in {ndx_file}:")
    print(f"  CP12:     {len(cp12_atoms)} atoms")
    print(f"  PRK:      {len(prk_atoms)} atoms")
    print(f"  GAP2:     {len(gap2_atoms)} atoms")
    print(f"  NAD:      {len(nad_atoms)} atoms")
    print(f"  LIGAND:   {len(cp12_atoms)} atoms (= CP12)")
    print(f"  RECEPTOR: {len(receptor_atoms)} atoms (= PRK + GAP2 + NAD)")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python make_index_groups.py <reference.pdb> <output.ndx>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    ndx_file = sys.argv[2]
    
    make_index_file(pdb_file, ndx_file)
    print(f"\nDone. Use with: gmx_mmpbsa -n {ndx_file}")
