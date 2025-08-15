#!/usr/bin/env python
# coding: utf-8
"""
plot_nonpca_selected_systems.py

Plot combined RMSD, Rg, RMSF heatmap & SecStruct for REF + subsets:

  * explicit:     --systems 0000 0100
  * range:        --start 0000 --end 0019
  * default:      blocks of 20 (0–19,20–39,…)

Generates combined RMSD, Rg, RMSF, and secondary structure plots
across blocks of systems (e.g., REF+0000–0019, REF+0020–0039…).

Usage:
    python multi_system_plots.py \
      --base_dir /path/to/postprocess \
      --output_dir /path/to/postprocess/combined_analysis \
      --block-size 20 \
      --skip 1 \
      --mode ptm \
      --cutoff 3.0 \
      --systems 0000 0001 0002 \
      --start 0001 \
      --end 0020 \
"""
import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import MDAnalysis as mda

from rmsd_rg        import compute_rmsd_sel, compute_rg_sel
from rmsf_secstruct import compute_rmsf, compute_secfrac, ALL_CODES, STANDARD_RES
from rmsf_secstruct import GAP_SEGIDS, PRK_SEGIDS, CP12_SEGIDS
from collections import defaultdict

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

def discover_systems(base_dir):
    return sorted(
        [d for d in os.listdir(base_dir) if d == "REF" or d.isdigit()],
        key=lambda x: (x != "REF", int(x) if x.isdigit() else x)
    )

def select_shell(u, mode, cutoff):
    if mode == "all":
        return u.select_atoms("not resname SOL NA CL NAD")
    if mode == "ptm":
        pts = [r.atoms for r in u.residues if r.resname.upper() not in STANDARD_RES]
        if not pts:
            return u.select_atoms("not resname SOL NA CL NAD")
        pg = pts[0]
        for grp in pts[1:]:
            pg = pg.union(grp)
        return pg.union(u.select_atoms(f"byres (around {cutoff} group pg)", pg=pg))
    # NAD mode
    ngs = [r.atoms for r in u.residues if r.resname.upper() in {"NAD","NADH","NADPH","NAI"}]
    if not ngs:
        return u.select_atoms("not resname SOL NA CL NAD")
    ng = ngs[0]
    for grp in ngs[1:]:
        ng = ng.union(grp)
    return ng.union(u.select_atoms(f"byres (around {cutoff} group ng)", ng=ng))

def plot_rmsd(base, systems, skip, mode, cutoff, dt, topo, traj, blk, outdir):
    idxs = [int(s) for s in systems if s != "REF"]
    lo, hi = (min(idxs), max(idxs)) if idxs else (0, 0)
    prefix = f"{blk}-" if blk else ""
    outfile = os.path.join(outdir, f"{prefix}{lo}_{hi}_combined_rmsd.pdf")

    fig, ax = plt.subplots(figsize=(10,5))
    cmap = plt.get_cmap("tab20")
    for s in systems:
        u = mda.Universe(os.path.join(base, s, topo),
                         os.path.join(base, s, traj), skip=skip)
        sel  = select_shell(u, mode, cutoff)
        n    = len(u.trajectory)
        t    = np.arange(n) * dt
        data = compute_rmsd_sel(u, sel, list(range(n)))
        avg  = np.cumsum(data) / (np.arange(n) + 1)

        if s == "REF":
            lbl, color, lw, ls = "REF", "k", 2, "--"
        else:
            lbl   = str(int(s))
            color = cmap(int(s) % cmap.N)
            lw, ls = 1, "-"

        ax.plot(t, avg, linestyle=ls, linewidth=lw,
                color=color, alpha=0.8, label=lbl)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.grid(ls="--", alpha=0.3)
    ticks = ax.get_xticks()
    ax.set_xticklabels((ticks).astype(int))
    ax.legend(bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"→ {outfile}")
    plt.close(fig)

def plot_rg(base, systems, skip, mode, cutoff, dt, topo, traj, blk, outdir):
    idxs = [int(s) for s in systems if s != "REF"]
    lo, hi = (min(idxs), max(idxs)) if idxs else (0, 0)
    prefix = f"{blk}-" if blk else ""
    outfile = os.path.join(outdir, f"{prefix}{lo}_{hi}_combined_rg.pdf")

    fig, ax = plt.subplots(figsize=(10,5))
    cmap = plt.get_cmap("tab20")
    for s in systems:
        u = mda.Universe(os.path.join(base, s, topo),
                         os.path.join(base, s, traj), skip=skip)
        sel  = select_shell(u, mode, cutoff)
        n    = len(u.trajectory)
        t    = np.arange(n) * dt
        data = compute_rg_sel(u, sel, list(range(n)))
        avg  = np.cumsum(data) / (np.arange(n) + 1)

        if s == "REF":
            lbl, color, lw, ls = "REF", "k", 2, "--"
        else:
            lbl   = str(int(s))
            color = cmap(int(s) % cmap.N)
            lw, ls = 1, "-"

        ax.plot(t, avg, linestyle=ls, linewidth=lw,
                color=color, alpha=0.8, label=lbl)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rg (Å)")
    ax.grid(ls="--", alpha=0.3)
    ticks = ax.get_xticks()
    ax.set_xticklabels((ticks).astype(int))
    ax.legend(bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"→ {outfile}")
    plt.close(fig)

def plot_rmsf(base, systems, skip, mode, cutoff, topo, traj, blk, outdir):
    # 0) Output filename
    idxs   = [int(s) for s in systems if s!="REF"]
    lo, hi = (min(idxs), max(idxs)) if idxs else (0,0)
    prefix = f"{blk}-" if blk else ""
    outpath = os.path.join(outdir, f"{prefix}{lo}_{hi}_combined_rmsf_heatmap.pdf")

    # 1) Load REF & compute RMSF
    u0        = mda.Universe(
                   os.path.join(base,"REF",topo),
                   os.path.join(base,"REF",traj),
                   skip=skip
               )
    sel0      = select_shell(u0, mode, cutoff)
    labs0, vals0 = compute_rmsf(u0, sel0, list(range(len(u0.trajectory))))
    sel_idx   = sel0.indices
    L         = len(vals0)

    # 2) Compute/pad/truncate RMSF for variants
    all_vals = [vals0]
    for s in systems[1:]:
        u      = mda.Universe(
                     os.path.join(base,s,topo),
                     os.path.join(base,s,traj),
                     skip=skip
                 )
        selsys = u.atoms[sel_idx]
        _, v   = compute_rmsf(u, selsys, list(range(len(u.trajectory))))
        if   len(v) < L: v = np.pad(v, (0, L-len(v)), constant_values=0.0)
        elif len(v) > L: v = v[:L]
        all_vals.append(v)
    arr = np.vstack(all_vals)  # shape = (1 + n_variants) × n_residues

    # 3) Figure: two rows (REF / variants)
    fig = plt.figure(figsize=(max(8, len(systems)*0.4), 6))
    G   = gs.GridSpec(2, 1, height_ratios=[1, len(systems)-1], hspace=0.05)
    ax0 = fig.add_subplot(G[0])
    ax1 = fig.add_subplot(G[1], sharex=ax0)

    # 4) Plot heatmaps with fixed 0–40 Å limits
    vmin, vmax = 0.0, 40.0
    ax0.imshow(arr[0:1], aspect="auto", origin="lower",
               cmap="YlGn", vmin=vmin, vmax=vmax, interpolation="none")
    im = ax1.imshow(arr[1:], aspect="auto", origin="lower",
                    cmap="YlGn", vmin=vmin, vmax=vmax, interpolation="none")

    # 5) Draw chain‐boundary lines (dotted)
    # extract chain ID from each label :contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}
    chain_ids = [lbl.split("-", 2)[1] for lbl in labs0]
    bounds    = [i for i in range(1, len(chain_ids))
                      if chain_ids[i] != chain_ids[i-1]]
    for b in bounds:
        # position at b - 0.5 to sit between residues
        ax0.axvline(b-0.5, linestyle="--", color="k", linewidth=0.5)
        ax1.axvline(b-0.5, linestyle="--", color="k", linewidth=0.5)

    # 6) Y-axis labels
    ax0.set_yticks([0])
    ax0.set_yticklabels(["REF"], fontsize=16)
    yt = list(range(len(systems)-1))
    yl = [str(int(s)) for s in systems[1:]]
    ax1.set_yticks(yt)
    ax1.set_yticklabels(yl, fontsize=16)
    ax1.invert_yaxis()

    # 7) Suppress x-ticks on REF panel
    ax0.xaxis.set_tick_params(bottom=False, labelbottom=False)

    # 8) Bottom‐axis ticks: one per chain midpoint, label “chainID\nGROUP(count)”
    starts = [0] + bounds
    ends   = bounds + [len(chain_ids)]
    mids   = [(s + e - 1)//2 for s, e in zip(starts, ends)]

    # count how many times each group has appeared so far
    group_counts = defaultdict(int)
    bottom_labels = []
    for i in mids:
        ch = chain_ids[i]
        # map chain → group name
        if   ch in GAP_SEGIDS:
            grp = "GAP"
        elif ch in PRK_SEGIDS:
            grp = "PRK"
        elif ch in CP12_SEGIDS:
            grp = "CP12"
        else:
            grp = ch
        # increment and fetch this group’s occurrence
        group_counts[grp] += 1
        occ = group_counts[grp]
        bottom_labels.append(f"{ch}\n{grp}({occ})")

    ax1.set_xticks(mids)
    ax1.set_xticklabels(bottom_labels, rotation=0, fontsize=16)
    for lbl in ax1.get_xticklabels():
        lbl.set_bbox(dict(facecolor='none', edgecolor='none', pad=3))
    ax1.set_xlim(0, len(labs0)-1)

    # 9) Numeric residue-index ticks *only* atop REF panel
    N    = arr.shape[1]
    step = 500 if N >= 500 else max(1, N//8)
    ticks = np.arange(0, N, step)
    labs  = (ticks + 1).tolist()

    ax0_top = ax0.twiny()
    ax0_top.set_xticks(ticks)
    ax0_top.set_xticklabels(labs, fontsize=16, rotation=0)
    ax0_top.set_xlim(0, len(labs0)-1)
    ax0_top.set_xlabel("Residue index", labelpad=4)

    # 10) Single colorbar fixed 0–40 Å
    cax  = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cax, label="RMSF (Å)")
    cbar.set_ticks(np.linspace(vmin, vmax, 9))

    plt.tight_layout(rect=[0,0,0.9,1])
    plt.savefig(outpath, dpi=300)
    plt.close(fig)
    print(f"→ {outpath}")

def plot_secstruct(base, systems, skip, mode, cutoff, topo, traj, blk, outdir):
    idxs = [int(s) for s in systems if s != "REF"]
    lo, hi = (min(idxs), max(idxs)) if idxs else (0, 0)
    prefix = f"{blk}-" if blk else ""
    outfile = os.path.join(outdir, f"{prefix}{lo}_{hi}_combined_secstruct.pdf")

    fracs = {}
    for s in systems:
        u = mda.Universe(os.path.join(base, s, topo),
                         os.path.join(base, s, traj), skip=skip)
        fracs[s] = compute_secfrac(u, select_shell(u, mode, cutoff),
                                   list(range(len(u.trajectory))))

    codes = ALL_CODES
    x = np.arange(len(systems))
    bottom = np.zeros(len(systems))
    fig, ax = plt.subplots(figsize=(max(8,len(systems)*0.4),5))
    cmap = plt.get_cmap("tab10")
    for i, code in enumerate(codes):
        vals = [fracs[s][code] for s in systems]
        ax.bar(x, vals, bottom=bottom, width=0.8, color=cmap(i), label=code)
        bottom += np.array(vals)

    labels = ["REF"] + [str(int(s)) for s in systems[1:]]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0, fontsize=16)
    ax.set_ylabel("Fraction of Contacts")
    ax.set_title("SecStruct")
    ax.legend(bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"→ {outfile}")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Plot RMSD, Rg, RMSF & SecStruct")
    parser.add_argument("--base_dir",   required=True, help="Root folder with REF/, 0000/, …")
    parser.add_argument("--output_dir", default=None, help="Where to save PDFs")
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("--systems", nargs="+", help="Explicit systems, e.g. REF 0000 0019")
    grp.add_argument("--start", help="First system in a contiguous range")
    parser.add_argument("--end",   help="Last system in a contiguous range (if using --start)")
    parser.add_argument("--block-size", type=int, default=20, help="Systems per block")
    parser.add_argument("--skip",      type=int, default=1, help="Trajectory frame skip")
    parser.add_argument("--mode",      choices=["all","ptm","nad"], default="all", help="Shell mode")
    parser.add_argument("--cutoff",    type=float, default=3.0, help="Cutoff (Å)")
    parser.add_argument("--dt",        type=float, default=1.0, help="Time step (ns)")
    parser.add_argument("--topo",      default="solute_fixed.pdb", help="Topology file")
    parser.add_argument("--traj",      default="solute_fit.xtc",   help="Trajectory file")
    args = parser.parse_args()

    base   = args.base_dir
    outdir = args.output_dir or base
    os.makedirs(outdir, exist_ok=True)

    all_sys = discover_systems(base)
    if args.systems:
        blocks = [["REF"] + args.systems]
    elif args.start:
        if not args.end:
            parser.error("--end is required when using --start")
        i0 = all_sys.index(args.start)
        i1 = all_sys.index(args.end)
        blocks = [["REF"] + all_sys[i0:i1+1]]
    else:
        nums = [s for s in all_sys if s!="REF"]
        blocks = [nums[i:i+args.block_size] for i in range(0, len(nums), args.block_size)]
        blocks = [["REF"] + b for b in blocks]

    for idx, systems in enumerate(blocks, start=1):
        plot_rmsd(     base, systems, args.skip, args.mode,
                      args.cutoff, args.dt, args.topo,
                      args.traj, idx, outdir)
        plot_rg(       base, systems, args.skip, args.mode,
                      args.cutoff, args.dt, args.topo,
                      args.traj, idx, outdir)
        plot_rmsf(     base, systems, args.skip, args.mode,
                      args.cutoff, args.topo, args.traj,
                      idx, outdir)
        plot_secstruct(base, systems, args.skip, args.mode,
                      args.cutoff, args.topo, args.traj,
                      idx, outdir)

if __name__ == "__main__":
    main()
