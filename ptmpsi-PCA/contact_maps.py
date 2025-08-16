#!/usr/bin/env python
# coding: utf-8
"""
contact_maps.py

Computes residue–residue contact‐probability heatmaps for selected interface segment pairs
and for PTM residues vs their neighbors, arranging them in a single row of subplots:

  CP12 vs GAP    CP12 vs PRK    [PTM vs Neighbors]

Usage:
  python contact_maps.py \
    --top     system.pdb \
    --xtc     traj1.xtc [traj2.xtc …] \
    --name    LABEL \
    [--cutoff 3.0] [--skip 1] [--rep CA|centroid] [--outdir figures]
"""
import os
import numpy as np
import MDAnalysis as mda
import matplotlib as mpl
import matplotlib
mpl.use("Agg")
import matplotlib.pyplot as plt
import argparse

# ─── Global styling ──────────────────────────────────────────────────────────
mpl.rcParams["font.size"]   = 14
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["legend.fontsize"] = 14
mpl.rcParams["axes.titlesize"]  = 14
mpl.rcParams["axes.labelsize"]  = 14
mpl.rcParams["xtick.labelsize"] = 14
mpl.rcParams["ytick.labelsize"] = 14
mpl.rcParams["axes.labelpad"] = 12
mpl.rcParams["xtick.major.pad"] = 6
mpl.rcParams["ytick.major.pad"] = 6
# ─────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
GAP_SEGIDS   = ["A","D","E","F","I","K","N","O"]
PRK_SEGIDS   = ["B","G","J","L"]
CP12_SEGIDS  = ["C","H","M","P"]
STANDARD_RES = {
    "ABU","ACE","AIB","ALA","ARG","ARGN","ASN","ASN1","ASP","ASP1","ASPH","ASPP","ASH",
    "CT3","CYS","CYS1","CYS2","CYSH","DALA","GLN","GLU","GLUH","GLUP","GLH","GLY",
    "HIS","HIS1","HISA","HISB","HISH","HISD","HISE","HISP","HSD","HSE","HSP","HYP",
    "ILE","LEU","LSN","LYS","LYSH","MELEU","MET","MEVAL","NAC","NME","NHE","NH2","PHE",
    "PHEH","PHEU","PHL","PRO","SER","THR","TRP","TRPH","TRPU","TYR","TYRH","TYRU","VAL",
    "PGLU","HID","HIE","HIP","LYP","LYN","CYN","CYM","CYX","DAB","ORN","NALA","NGLY",
    "NSER","NTHR","NLEU","NILE","NVAL","NASN","NGLN","NARG","NHID","NHIE","NHIP",
    "NHISD","NHISE","NHISH","NTRP","NPHE","NTYR","NGLU","NASP","NLYS","NORN","NDAB",
    "NLYSN","NPRO","NHYP","NCYS","NCYS2","NMET","NASPH","NGLUH","CALA","CGLY","CSER",
    "CTHR","CLEU","CILE","CVAL","CASN","CGLN","CARG","CHID","CHIE","CHIP","CHISD",
    "CHISE","CHISH","CTRP","CPHE","CTYR","CGLU","CASP","CLYS","CORN","CDAB","CLYSN",
    "CPRO","CHYP","CCYS","CCYS2","CMET","CASPH","CGLUH",
    "HOH","WAT","TIP3","UNK","SOL",
    "NAD","NADH","NADPH","NAI","MG","ZN","NA","CL","K","CA"
}
# ──────────────────────────────────────────────────────────────────────────────

def get_segment_residues(u, segids):
    sel = u.select_atoms(" or ".join(f"segid {s}" for s in segids))
    return list(sel.residues)

def get_ptm_residues(u):
    return [r for r in u.residues if r.resname.upper() not in STANDARD_RES]

def preload_coordinates(u, residues, frames, rep):
    """
    Load representative coordinates per frame for each residue:
      rep=="CA" uses Cα if present, else centroid.
    """
    F, N = len(frames), len(residues)
    coords = np.zeros((F, N, 3), float)
    for i, res in enumerate(residues):
        if rep=="CA":
            ca = res.atoms.select_atoms("name CA")
            if len(ca)==1:
                idx = ca.indices[0]
                for fi, frm in enumerate(frames):
                    u.trajectory[frm]
                    coords[fi,i] = u.atoms[idx].position
                continue
        for fi, frm in enumerate(frames):
            u.trajectory[frm]
            coords[fi,i] = res.atoms.center_of_mass()
    return coords

def compute_pairwise_prob_chunked(c1, c2, cutoff, block_size=200):
    """
    Compute contact probabilities between coords1 (F×N1×3) and coords2 (F×N2×3)
    by processing frames in blocks to limit memory.
    """
    F, N1, _ = c1.shape
    N2 = c2.shape[1]
    counts = np.zeros((N1, N2), int)
    for start in range(0, F, block_size):
        end = min(start+block_size, F)
        diff = c1[start:end,:,None,:] - c2[start:end,None,:,:]
        d = np.linalg.norm(diff, axis=-1)
        counts += (d < cutoff).sum(axis=0)
    return counts / float(F)

def plot_maps(maps, groups, titles, name, outpath):
    """
    maps: list of (Ni×Nj) arrays
    groups: list of (group1_residues, group2_residues)
    """
    n = len(maps)
    fig, axes = plt.subplots(1, n, figsize=(6*n,6), constrained_layout=True)
    if n == 1:
        axes = [axes]
    cmap = plt.get_cmap("Reds")

    for idx, ax in enumerate(axes):
        mat, (grp1, grp2) = maps[idx], groups[idx]

        im = ax.imshow(mat, aspect='auto', cmap=cmap, vmin=0, vmax=1)

        # residue tick labels
        N1, N2 = mat.shape
        ax.set_yticks(np.arange(N1))
        ax.set_yticklabels([f"{r.resname}-{r.resid}" for r in grp1], fontsize=6)
        ax.set_xticks(np.arange(N2))
        ax.set_xticklabels([f"{r.resname}-{r.resid}" for r in grp2],
                           rotation=90, fontsize=6)

        # segment boundaries & centred region labels (segname–segid)
        # Y axis
        segids1 = [r.segid for r in grp1]
        starts1, curr = [0], segids1[0]
        for i, s in enumerate(segids1):
            if s != curr:
                starts1.append(i)
                curr = s
        mids1, labs1 = [], []
        for i, st in enumerate(starts1):
            end = starts1[i+1]-1 if i+1<len(starts1) else N1-1
            mids1.append((st+end)/2)
            seg = segids1[st]
            name = ("GAP" if seg in GAP_SEGIDS else
                    "PRK" if seg in PRK_SEGIDS else
                    "CP12" if seg in CP12_SEGIDS else seg)
            labs1.append(f"{name}-{seg}")
        for b in starts1[1:]:
            ax.hlines(b-0.5, xmin=-0.5, xmax=N2-0.5,
                      color='black', linewidth=1, linestyle=':')
        for m, l in zip(mids1, labs1):
            ax.text(1.01, m, l, ha='left', va='center', rotation=90, 
                    fontsize=14, transform=ax.get_yaxis_transform())

        # X axis
        segids2 = [r.segid for r in grp2]
        starts2, curr = [0], segids2[0]
        for i, s in enumerate(segids2):
            if s != curr:
                starts2.append(i)
                curr = s
        mids2, labs2 = [], []
        for i, st in enumerate(starts2):
            end = starts2[i+1]-1 if i+1<len(starts2) else N2-1
            mids2.append((st+end)/2)
            seg = segids2[st]
            name = ("GAP" if seg in GAP_SEGIDS else
                    "PRK" if seg in PRK_SEGIDS else
                    "CP12" if seg in CP12_SEGIDS else seg)
            labs2.append(f"{name}-{seg}")
        for b in starts2[1:]:
            ax.vlines(b-0.5, ymin=-0.5, ymax=N1-0.5,
                      color='black', linewidth=1, linestyle=':')
        for m, l in zip(mids2, labs2):
            ax.text(m, 1.02, l, ha='center', va='bottom',
                    rotation=0, fontsize=6,
                    transform=ax.get_xaxis_transform())

        ax.set_title(titles[idx], pad=14)

    cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                        fraction=0.05, pad=0.02)
    cbar.set_label("Contact probability")
    fig.suptitle(f"Contact Maps", fontsize=14)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(f"Saved contact maps → {outpath}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top",    required=True, help="System topology (PDB/GRO)")
    parser.add_argument("--xtc",    required=True, nargs="+", help="Trajectory XTC(s)")
    parser.add_argument("--name",   required=True, help="Output label")
    parser.add_argument("--cutoff", type=float, default=3.0, help="Contact cutoff (Å)")
    parser.add_argument("--skip",   type=int,   default=1,   help="Frame skip")
    parser.add_argument("--rep",    choices=["CA","centroid"], default="CA",
                        help="Representative atom per residue")
    parser.add_argument("--outdir", default="figures_contactmaps", help="Output directory")
    args = parser.parse_args()
    iface_cut = args.cutoff - 1.0

    os.makedirs(args.outdir, exist_ok=True)
    u      = mda.Universe(args.top, *args.xtc, skip=args.skip)
    frames = list(range(len(u.trajectory)))

    # full segment AtomGroups & PTMs
    gap_all  = u.select_atoms(" or ".join(f"segid {s}" for s in GAP_SEGIDS))
    prk_all  = u.select_atoms(" or ".join(f"segid {s}" for s in PRK_SEGIDS))
    cp12_all = u.select_atoms(" or ".join(f"segid {s}" for s in CP12_SEGIDS))
    ptm_all  = get_ptm_residues(u)
    ptm_atoms = u.atoms[[]]
    for r in ptm_all: ptm_atoms += r.atoms
    neigh_sel = u.select_atoms(f"byres (around {args.cutoff} group PTM)", PTM=ptm_atoms)
    ptm_neigh = [r for r in neigh_sel.residues if r not in ptm_all]

    # always: CP12–GAP, CP12–PRK
    groups = []
    groups.append((
        list(u.select_atoms(
            f"(segid {' or segid '.join(CP12_SEGIDS)}) and byres (around {iface_cut} group gap_all)",
            gap_all=gap_all).residues),
        list(u.select_atoms(
            f"(segid {' or segid '.join(GAP_SEGIDS)}) and byres (around {iface_cut} group cp12_all)",
            cp12_all=cp12_all).residues)
    ))
    groups.append((
        list(u.select_atoms(
            f"(segid {' or segid '.join(CP12_SEGIDS)}) and byres (around {iface_cut} group prk_all)",
            prk_all=prk_all).residues),
        list(u.select_atoms(
            f"(segid {' or segid '.join(PRK_SEGIDS)}) and byres (around {iface_cut} group cp12_all)",
            cp12_all=cp12_all).residues)
    ))
    titles = ["CP12 vs GAP", "CP12 vs PRK"]

    # optional PTM vs Neighbors
    if ptm_all:
        groups.append((ptm_all, ptm_neigh))
        titles.append("PTM vs Neighbors")

    # preload coords
    coords = [preload_coordinates(u, g, frames, args.rep) for g,_ in groups]
    coords += [preload_coordinates(u, g, frames, args.rep) for _,g in groups]

    # compute maps
    maps = [compute_pairwise_prob_chunked(coords[i], coords[i+len(groups)], args.cutoff)
            for i in range(len(groups))]

    # plot and save
    outpath = os.path.join(args.outdir, f"{args.name}_contact_maps.pdf")
    plot_maps(maps, groups, titles, args.name, outpath)

if __name__ == "__main__":
    main()
