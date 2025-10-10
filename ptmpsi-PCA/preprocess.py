#!/usr/bin/env python
# coding: utf-8
"""
Module: preprocess.py

This module converts .gro files to PDB, fixes chain IDs and periodic boundaries (PBC).
"""

import os
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.rms import rmsd as mda_rmsd
from MDAnalysis import transformations

def gro_to_pdb(gro_file, pdb_file):
    """
    Converts a GRO file to a PDB file using MDAnalysis.
    """
    u = mda.Universe(gro_file)
    from MDAnalysis.coordinates.PDB import PDBWriter
    with PDBWriter(pdb_file, n_atoms=u.atoms.n_atoms) as writer:
        writer.write(u.atoms)
    print(f"Converted {gro_file} → {pdb_file}")

def add_chain_ids(input_pdb, output_pdb, chain_ids="ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    Reads a PDB file, assigns chain IDs at residue number resets,
    and replaces the SEGID field with the chain ID.
    """
    chain_idx = 0
    na_chain = None
    cl_chain = None
    prev_resid = None
    out_lines = []
    with open(input_pdb, 'r') as fin:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                resid_str = line[22:26].strip()
                try:
                    resid = int(resid_str)
                except ValueError:
                    out_lines.append(line)
                    continue
                if resname == "NA":
                    if na_chain is None:
                        chain_idx += 1
                        na_chain = chain_ids[chain_idx]
                    chain_id = na_chain
                elif resname == "CL":
                    if cl_chain is None:
                        chain_idx += 1
                        cl_chain = chain_ids[chain_idx]
                    chain_id = cl_chain
                else:
                    if prev_resid is not None and resid < prev_resid:
                        chain_idx += 1
                    chain_id = chain_ids[chain_idx]
                    prev_resid = resid
                line = line.rstrip("\n").ljust(80)
                new_line = line[:21] + chain_id + line[22:72] + chain_id.ljust(4) + line[76:]
                out_lines.append(new_line + "\n")
            else:
                out_lines.append(line)
    with open(output_pdb, 'w') as fout:
        fout.writelines(out_lines)
    print(f"Chain IDs assigned. Output written to {output_pdb}")

def fix_pbc(in_pdb, out_pdb):
    """
    Fully fixes periodic boundaries so that each chain is made whole.
    Steps:
      1. Parse CRYST1 (or use a default box if missing).
      2. Compute the box matrix.
      3. Convert atomic coordinates to fractional coordinates.
      4. For each chain, adjust fractional coordinates so that the chain is contiguous.
      5. Convert back to Cartesian coordinates and write the new PDB.
    """
    def parse_cryst1(line):
        a = float(line[6:15])
        b = float(line[15:24])
        c = float(line[24:33])
        alpha = float(line[33:40])
        beta = float(line[40:47])
        gamma = float(line[47:54])
        return a, b, c, alpha, beta, gamma

    def box_matrix_from_abc(a, b, c, alpha, beta, gamma):
        alpha_r = np.radians(alpha)
        beta_r  = np.radians(beta)
        gamma_r = np.radians(gamma)
        ax = a; ay = 0.0; az = 0.0
        bx = b * np.cos(gamma_r); by = b * np.sin(gamma_r); bz = 0.0
        cx = c * np.cos(beta_r)
        cy = c * ((np.cos(alpha_r) - np.cos(gamma_r)*np.cos(beta_r)) / np.sin(gamma_r))
        cz = c * np.sqrt(1 - np.cos(alpha_r)**2 - np.cos(beta_r)**2 - np.cos(gamma_r)**2 + 2*np.cos(alpha_r)*np.cos(beta_r)*np.cos(gamma_r)) / np.sin(gamma_r)
        return np.array([[ax,ay,az],[bx,by,bz],[cx,cy,cz]], dtype=float)

    with open(in_pdb, 'r') as fin:
        lines = fin.readlines()
    a = b = c = alpha = beta = gamma = None
    for line in lines:
        if line.startswith("CRYST1"):
            a, b, c, alpha, beta, gamma = parse_cryst1(line)
            break
    if a is None:
        print(f"Warning: No CRYST1 record in {in_pdb}. Using default box.")
        a, b, c, alpha, beta, gamma = (100.0, 100.0, 100.0, 90.0, 90.0, 90.0)
    box = box_matrix_from_abc(a, b, c, alpha, beta, gamma)

    coords = []
    chainIDs = []
    resnames = []
    atom_line_indices = []
    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
            chainIDs.append(line[21])
            resnames.append(line[17:20].strip())
            atom_line_indices.append(i)
    coords = np.array(coords, dtype=float)

    def fix_pbc_by_chain(coords, chainIDs, resnames, box):
        inv_box = np.linalg.inv(box)
        frac = np.dot(coords, inv_box)
        new_frac = frac.copy()
        ion_residues = {"NA", "CL"}
        chain_map = {}
        for i, (ch, resn) in enumerate(zip(chainIDs, resnames)):
            key = (ch, i) if resn in ion_residues else ch
            chain_map.setdefault(key, []).append(i)
        for key, indices in chain_map.items():
            if len(indices) == 1:
                new_frac[indices[0]] = frac[indices[0]] % 1.0
            else:
                for i_local in range(1, len(indices)):
                    i_prev = indices[i_local - 1]
                    i_curr = indices[i_local]
                    ref_frac = new_frac[i_prev]
                    cur_frac = frac[i_curr]
                    best_shift = np.array([0, 0, 0], dtype=float)
                    best_dist = 1e9
                    for sx in (-1, 0, 1):
                        for sy in (-1, 0, 1):
                            for sz in (-1, 0, 1):
                                trial = cur_frac + np.array([sx, sy, sz])
                                dist = np.linalg.norm(np.dot(trial - ref_frac, box))
                                if dist < best_dist:
                                    best_dist = dist
                                    best_shift = np.array([sx, sy, sz])
                    new_frac[i_curr] = cur_frac + best_shift
        return np.dot(new_frac, box)

    fixed_coords = fix_pbc_by_chain(coords, chainIDs, resnames, box)
    out_lines = list(lines)
    for i_atom, idx in enumerate(atom_line_indices):
        old_line = out_lines[idx]
        x, y, z = fixed_coords[i_atom]
        new_line = old_line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + old_line[54:]
        out_lines[idx] = new_line
    with open(out_pdb, 'w') as fout:
        fout.writelines(out_lines)
    print(f"PBC fixed for {in_pdb}. Output written to {out_pdb}")

########################################################################
# Atom Selection for Analysis
########################################################################

def get_unique(ag):
    # Some versions of MDAnalysis support callable unique()
    return ag.unique() if callable(ag.unique) else ag.unique

def select_atoms_for_analysis(u, ptm_resnames=None, standard_resnames={
    "ABU", "ACE", "AIB", "ALA", "ARG", "ARGN", "ASN", "ASN1", 
    "ASP", "ASP1", "ASPH", "ASPP", "ASH", "CT3", "CYS", "CYS1", 
    "CYS2", "CYSH", "DALA", "GLN", "GLU", "GLUH", "GLUP", "GLH", 
    "GLY", "HIS", "HIS1", "HISA", "HISB", "HISH", "HISD", "HISE", 
    "HISP", "HSD", "HSE", "HSP", "HYP", "ILE", "LEU", "LSN", "LYS", 
    "LYSH", "MELEU", "MET", "MEVAL", "NAC", "NME", "NHE", "NH2", "PHE", 
    "PHEH", "PHEU", "PHL", "PRO", "SER", "THR", "TRP", "TRPH", "TRPU", 
    "TYR", "TYRH", "TYRU", "VAL", "PGLU", "HID", "HIE", "HIP", "LYP", 
    "LYN", "CYN", "CYM", "CYX", "DAB", "ORN", "HYP", "NALA", "NGLY", 
    "NSER", "NTHR", "NLEU", "NILE", "NVAL", "NASN", "NGLN", "NARG", 
    "NHID", "NHIE", "NHIP", "NHISD", "NHISE", "NHISH", "NTRP", "NPHE", 
    "NTYR", "NGLU", "NASP", "NLYS", "NORN", "NDAB", "NLYSN", "NPRO", 
    "NHYP", "NCYS", "NCYS2", "NMET", "NASPH", "NGLUH", "CALA", "CGLY", 
    "CSER", "CTHR", "CLEU", "CILE", "CVAL", "CASN", "CGLN", "CARG", 
    "CHID", "CHIE", "CHIP", "CHISD", "CHISE", "CHISH", "CTRP", "CPHE", 
    "CTYR", "CGLU", "CASP", "CLYS", "CORN", "CDAB", "CLYSN", "CPRO", 
    "CHYP", "CCYS", "CCYS2", "CMET", "CASPH", "CGLUH",
    "HOH", "WAT", "TIP3", "UNK", "SOL",
    "NAD", "NADH", "NADPH", "NAI",
    "MG", "ZN", "NA", "CL", "K", "CA"
}, gap_segids=["A","D","E","F","I","K","N","O"],
                              prk_segids=["B","G","J","L"],
                              cp12_segids=["C","H","M","P"],
                              nad_resnames=["NAD","NADH","NADPH","NAI"],
                              radius_ptm=3.0, radius_nad=3.0, radius_interfaces=2.0):
    if isinstance(ptm_resnames, str) and ptm_resnames.upper() == "AUTO":
        auto_list = []
        for res in u.residues:
            if res.resname.upper() not in standard_resnames:
                auto_list.append(res.resname.upper())
        ptm_resnames = list(set(auto_list))
        print("Auto-detected non-standard PTM resnames:", ptm_resnames)
    ptm_group = u.atoms[[]]
    for rname in ptm_resnames:
        ptm_group += u.select_atoms(f"resname {rname}")
    ptm_neighbors = u.select_atoms(f"byres (around {radius_ptm} group ptm_group)", ptm_group=ptm_group)
    ptm_total = ptm_group + ptm_neighbors
    nad_group = u.atoms[[]]
    for rname in nad_resnames:
        nad_group += u.select_atoms(f"resname {rname}")
    nad_neighbors = u.select_atoms(f"byres (around {radius_nad} group nad_group)", nad_group=nad_group)
    nad_total = nad_group + nad_neighbors
    def segid_to_query(segids):
        return " or ".join([f"segid {s}" for s in segids])
    gap_sel  = u.select_atoms(segid_to_query(gap_segids))
    prk_sel  = u.select_atoms(segid_to_query(prk_segids))
    cp12_sel = u.select_atoms(segid_to_query(cp12_segids))
    gap_around_prk = u.select_atoms(f"byres (around {radius_interfaces} group gap_sel)", gap_sel=gap_sel) & prk_sel
    prk_around_gap = u.select_atoms(f"byres (around {radius_interfaces} group prk_sel)", prk_sel=prk_sel) & gap_sel
    gap_prk_if = gap_around_prk + prk_around_gap
    prk_around_cp12 = u.select_atoms(f"byres (around {radius_interfaces} group prk_sel)", prk_sel=prk_sel) & cp12_sel
    cp12_around_prk = u.select_atoms(f"byres (around {radius_interfaces} group cp12_sel)", cp12_sel=cp12_sel) & prk_sel
    prk_cp12_if = prk_around_cp12 + cp12_around_prk
    gap_around_cp12 = u.select_atoms(f"byres (around {radius_interfaces} group gap_sel)", gap_sel=gap_sel) & cp12_sel
    cp12_around_gap = u.select_atoms(f"byres (around {radius_interfaces} group cp12_sel)", cp12_sel=cp12_sel) & gap_sel
    gap_cp12_if = gap_around_cp12 + cp12_around_gap
    interfaces = {
        "gap_prk": get_unique(gap_prk_if),
        "prk_cp12": get_unique(prk_cp12_if),
        "gap_cp12": get_unique(gap_cp12_if)
    }
    final_selection = ptm_total + nad_total + prk_cp12_if + gap_cp12_if
    try:
        final_unique = final_selection.unique()
    except TypeError:
        final_unique = final_selection.unique
    return {"final": final_unique, "interfaces": interfaces, "ptm": ptm_total, "nad": "nad_total"}

########################################################################


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: preprocess.py <structure.(gro)>")
        sys.exit(1)

    u   = mda.Universe(sys.argv[1])

    # grab the AtomGroup you’ll feed to PCA / analyses
    sel = select_atoms_for_analysis(u, ptm_resnames="AUTO")["final"]

    print("Selected atoms:", sel.n_atoms)
