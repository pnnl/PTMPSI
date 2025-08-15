#!/usr/bin/env python
# coding: utf-8
"""
rmsf_secstruct.py

Compare one system vs reference:
  • Top panel: stacked bar of secondary‐structure fractions
  • Bottom panel: discrete‐green heatmap of per-residue RMSF (0–20 Å), grouped by segment
Produces:
  • PDF figure in --outdir
  • <name>_secstruct.dat
  • <name>_rmsf.dat

Usage:
  python rmsf_secstruct.py \
    --top      system.pdb \
    --xtc      sys.xtc [...] \
    --ref_top  ref.pdb \
    --ref_xtc  ref.xtc [...] \
    (--all | --mode ptm|nad) \
    [--cutoff 3.0] [--skip 1] \
    [--outdir figures_rmsf_secstruct] [--name LABEL]
"""
import os
import argparse
import tempfile

import numpy as np
import MDAnalysis as mda
import mdtraj as md
import matplotlib as mpl
import matplotlib
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

# ─── Global styling ──────────────────────────────────────────────────────────
mpl.rcParams["font.size"]   = 16
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["legend.fontsize"] = 16
mpl.rcParams["axes.titlesize"]  = 16
mpl.rcParams["axes.labelsize"]  = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 12
mpl.rcParams["xtick.major.pad"] = 6
mpl.rcParams["ytick.major.pad"] = 6
# ─────────────────────────────────────────────────────────────────────────────

# DSSP codes plus “O” for other/coil
DSSP_CODES = ['H','B','E','G','I','P','S','T']
ALL_CODES  = DSSP_CODES + ['O']

# Standard residues (including water/ions) to exclude from PTM detection
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

# segment definitions
GAP_SEGIDS  = set(["A","D","E","F","I","K","N","O"])
PRK_SEGIDS  = set(["B","G","J","L"])
CP12_SEGIDS = set(["C","H","M","P"])


def compute_secfrac(u, sel, frames):
    """
    For selection 'sel' in Universe u over frames, compute
    the mean fraction of each DSSP code in ALL_CODES.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
    with mda.Writer(tmp, sel.n_atoms) as W:
        u.trajectory[frames[0]]; W.write(sel)
    topo = md.load_pdb(tmp).topology
    os.unlink(tmp)

    F = len(frames)
    counts = {c: 0 for c in ALL_CODES}

    for frm in frames:
        u.trajectory[frm]
        coords = sel.positions / 10.0
        dssp   = md.compute_dssp(md.Trajectory(coords[None,:,:], topo))[0]
        N      = dssp.size
        cdict  = {c: (dssp==c).sum() for c in DSSP_CODES}
        cdict['O'] = N - sum(cdict.values())
        for c in ALL_CODES:
            counts[c] += cdict[c]

    base = F * max(1, len(sel.residues))
    return {c: counts[c] / base for c in ALL_CODES}


def compute_rmsf(u, sel, frames):
    """
    Compute per-residue RMSF (COM) for selection sel over frames.
    Returns (labels, rmsf_array).
    """
    residues = list(sel.residues)
    labels   = [f"{r.resname}-{r.segid}-{r.resid}" for r in residues]
    F, A     = len(frames), sel.n_atoms

    coords = np.zeros((F, A, 3), float)
    for i, frm in enumerate(frames):
        u.trajectory[frm]
        coords[i] = sel.positions

    idxmap = {a: i for i, a in enumerate(sel.indices)}
    R = len(residues)
    coms = np.zeros((F, R, 3), float)
    for ri, r in enumerate(residues):
        atom_inds = [idxmap[a] for a in r.atoms.indices if a in idxmap]
        coms[:,ri,:] = coords[:, atom_inds, :].mean(axis=1)

    mean_com = coms.mean(axis=0)
    sqr      = ((coms - mean_com)**2).sum(axis=2)
    rmsf     = np.sqrt(sqr.mean(axis=0))

    return labels, rmsf


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--top",     required=True, help="System topology")
    p.add_argument("--xtc",     required=True, nargs="+", help="System XTC(s)")
    p.add_argument("--ref_top", required=True, help="Reference topology")
    p.add_argument("--ref_xtc", required=True, nargs="+", help="Reference XTC(s)")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--all",   action="store_true", help="All non-water/ions")
    grp.add_argument("--mode",  choices=["ptm","nad"], help="Shell mode")
    p.add_argument("--cutoff",  type=float, default=3.0, help="Shell cutoff (Å)")
    p.add_argument("--skip",    type=int,   default=1,   help="Frame skip")
    p.add_argument("--outdir",  default="figures_rmsf_secstruct", help="Output directory")
    p.add_argument("--name",    default="system",  help="Label for system")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Reference selection
    u_ref = mda.Universe(args.ref_top, *args.ref_xtc, skip=args.skip)
    fr_ref = list(range(len(u_ref.trajectory)))
    if args.all:
        sel_ref = u_ref.select_atoms("not resname SOL NA CL NAD")
    else:
        if args.mode=="ptm":
            pg = None
            for r in u_ref.residues:
                if r.resname.upper() not in STANDARD_RES:
                    pg = r.atoms if pg is None else pg + r.atoms
            if not pg:
                sel_ref = u_ref.select_atoms("not resname SOL NA CL NAD")
            else:
                sel_ref = pg + u_ref.select_atoms(
                    f"byres (around {args.cutoff} group pg)", pg=pg)
        else:  # nad
            ng = None
            for r in u_ref.residues:
                if r.resname.upper() in {"NAD","NADH","NADPH","NAI"}:
                    ng = r.atoms if ng is None else ng + r.atoms
            if not ng:
                sel_ref = u_ref.select_atoms("not resname SOL NA CL NAD")
            else:
                sel_ref = ng + u_ref.select_atoms(
                    f"byres (around {args.cutoff} group ng)", ng=ng)

    # ── System selection
    u_sys = mda.Universe(args.top, *args.xtc, skip=args.skip)
    fr_sys = list(range(len(u_sys.trajectory)))
    if args.all:
        sel_sys = u_sys.select_atoms("not resname SOL NA CL NAD")
    else:
        if args.mode=="ptm":
            pg = None
            for r in u_sys.residues:
                if r.resname.upper() not in STANDARD_RES:
                    pg = r.atoms if pg is None else pg + r.atoms
            if not pg:
                sel_sys = u_sys.select_atoms("not resname SOL NA CL NAD")
            else:
                sel_sys = pg + u_sys.select_atoms(
                    f"byres (around {args.cutoff} group pg)", pg=pg)
        else:
            ng = None
            for r in u_sys.residues:
                if r.resname.upper() in {"NAD","NADH","NADPH","NAI"}:
                    ng = r.atoms if ng is None else ng + r.atoms
            if not ng:
                sel_sys = u_sys.select_atoms("not resname SOL NA CL NAD")
            else:
                sel_sys = ng + u_sys.select_atoms(
                    f"byres (around {args.cutoff} group ng)", ng=ng)

    # ── Compute secondary‐structure fractions
    sec_ref = compute_secfrac(u_ref, sel_ref, fr_ref)
    sec_sys = compute_secfrac(u_sys, sel_sys, fr_sys)

    # ── Compute RMSF per residue
    labels, rmsf_ref = compute_rmsf(u_ref, sel_ref, fr_ref)
    _,      rmsf_sys = compute_rmsf(u_sys, sel_sys, fr_sys)
    N = len(labels)

    # ── Write .dat files
    out1 = os.path.join(args.outdir, f"{args.name}_secstruct.dat")
    with open(out1, "w") as f:
        f.write("# code\tfrac_REF\tfrac_SYS\n")
        for c in ALL_CODES:
            f.write(f"{c}\t{sec_ref[c]:.6f}\t{sec_sys[c]:.6f}\n")
    print("Wrote →", out1)

    out2 = os.path.join(args.outdir, f"{args.name}_rmsf.dat")
    with open(out2, "w") as f:
        f.write("# residue\tRMSF_REF\tRMSF_SYS\n")
        for lbl, r1, r2 in zip(labels, rmsf_ref, rmsf_sys):
            f.write(f"{lbl}\t{r1:.6f}\t{r2:.6f}\n")
    print("Wrote →", out2)

#    # ── Plot
#    greens = ["#edf8e9","#bae4b3","#74c476","#31a354","#006d2c"]
#    bounds = np.linspace(0,20,6)
#    norm = BoundaryNorm(bounds, len(greens))
#    cmap_r = ListedColormap(greens)


    M = np.vstack([rmsf_ref, rmsf_sys])
    vmax  = np.max(M) # maximum RMSF in your data
    n_shades = 5 # how many discrete green levels

    # define a light→dark green palette (pick as many as n_shades)
    greens = ["#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c"][:n_shades]

    # build colormap and normalization from 0→vmax in n_shades bins
    cmap_r = ListedColormap(greens)
    bounds = np.linspace(0, vmax, n_shades+1)
    norm   = BoundaryNorm(bounds, n_shades)

    fig = plt.figure(figsize=(14,10))
    gs  = fig.add_gridspec(2,1, height_ratios=[1.5,2], hspace=0.3)

    # Top: stacked bar
    ax0 = fig.add_subplot(gs[0])
    x0 = np.array([0,1])
    cmap_ss = {c: plt.get_cmap("tab10")(i) for i,c in enumerate(ALL_CODES)}
    bottom = np.zeros(2)
    barw = 0.15
    for code in ALL_CODES:
        vals = [sec_ref[code], sec_sys[code]]
        ax0.bar(x0, vals, bottom=bottom, color=cmap_ss[code], width=barw)
        bottom += vals
    ax0.set_xticks(x0);
    ax0.set_xticklabels(["REF", args.name])
    ax0.set_ylabel("Fraction of Secondary Structure")
    #ax0.set_title("Secondary‐Structure Composition")
    ax0.legend(ALL_CODES, bbox_to_anchor=(1,1))

    # Physically shrink top subplot to half width
    pos0 = ax0.get_position()
    ax0.set_position([pos0.x0, pos0.y0, pos0.width * 0.2, pos0.height])

    # Bottom: RMSF heatmap
    ax1 = fig.add_subplot(gs[1])
    M = np.vstack([rmsf_ref, rmsf_sys])
    im = ax1.imshow(M, aspect='auto', cmap=cmap_r, norm=norm)
    ax1.axhline(0.5, color='white', linewidth=4)

    segids = [lbl.split('-')[1] for lbl in labels]
    groups = []
    current = None
    for i, s in enumerate(segids):
        if s != current:
            current = s
            groups.append([s, [i]])
        else:
            groups[-1][1].append(i)

    mids, xtlabs = [], []
    for seg, idxs in groups:
        start, end = idxs[0], idxs[-1]
        mids.append((start+end)//2)
        if seg in GAP_SEGIDS:
            grp="GAP"
        elif seg in PRK_SEGIDS:
            grp="PRK"
        elif seg in CP12_SEGIDS:
            grp="CP12"
        else:
            grp=seg
        xtlabs.append(f"{seg}\n{grp}")
        ax1.axvline(end+0.5, color='black', linestyle=':', linewidth=2)

    ax1.set_xticks(mids);
    ax1.set_xticklabels(xtlabs, fontsize=16)
    ax1.set_yticks([0,1]);
    ax1.set_yticklabels(["REF", args.name])
    ax1.set_ylabel("Ref #")
    ax1.set_title("Per‐residue RMSF (Å)")
    # ─── add numeric residue indices along the TOP of the same plot ───
    ax1_top = ax1.twiny()               # create a second x‐axis
    N = len(labels)                     # total residues
    step = max(1, N // 8)               # ~8 ticks
    num_ticks = np.arange(0, N, step)   # 0‐based
    num_labels = (num_ticks + 1).tolist()  # convert to 1‐based

    ax1_top.set_xticks(num_ticks)
    ax1_top.set_xticklabels(num_labels, fontsize=16, rotation=0)
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xlabel("Residue #")
    
    cbar = plt.colorbar(im, ax=ax1, boundaries=bounds, ticks=bounds)
    cbar.set_label("RMSF (Å)")

    plt.tight_layout()
    outfig = os.path.join(args.outdir, f"{args.name}_rmsf_secstruct.pdf")
    plt.savefig(outfig, dpi=300)
    print("Saved →", outfig)


if __name__ == "__main__":
    main()
