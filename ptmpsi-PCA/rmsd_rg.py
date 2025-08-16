#!/usr/bin/env python
# coding: utf-8
"""
Module: rmsd_rg.py

Computes RMSD and radius of gyration (Rg) over time for a system and reference,
plots cumulative (window=1) or running‐window (window>1) averages versus time (ns),
and writes a .dat file of the resulting time series—all with no dropped points.
"""
import os
import numpy as np
import MDAnalysis as mda
import matplotlib as mpl
import matplotlib
mpl.use("Agg")
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import rmsd as mda_rmsd
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

def moving_average(data, window=10):
    """
    Compute a running‐window average of `data` using at most `window` points
    up to each index.  For i < window, averages data[0:i+1]; thereafter
    averages data[i-window+1:i+1].  Returns an array of same length.
    """
    N = len(data)
    out = np.zeros(N)
    for i in range(N):
        start = max(0, i - window + 1)
        out[i] = data[start:i+1].mean()
    return out

def compute_rmsd_sel(u, sel_atoms, frame_indices, reference_frame=0):
    """
    Compute RMSD of sel_atoms against the first reference frame
    for each frame in frame_indices.
    """
    rmsd_vals = []
    u.trajectory[frame_indices[reference_frame]]
    ref_coords = sel_atoms.positions.copy()
    for fi in frame_indices:
        u.trajectory[fi]
        coords = sel_atoms.positions.copy()
        rmsd_vals.append(
            mda_rmsd(coords, ref_coords, center=True, superposition=True)
        )
    return np.array(rmsd_vals)

def compute_rg_sel(u, sel_atoms, frame_indices):
    """
    Compute radius of gyration of sel_atoms for each frame in frame_indices.
    """
    rg_vals = []
    for fi in frame_indices:
        u.trajectory[fi]
        coords = sel_atoms.positions.copy()
        com = np.mean(coords, axis=0)
        rg_vals.append(np.sqrt(np.mean(np.sum((coords - com)**2, axis=1))))
    return np.array(rg_vals)

def plot_rmsd_rg(topology, traj_list,
                 ref_topology, ref_traj_list,
                 sel_str="not resname SOL NA CL NAD",
                 replicates=True,
                 skip=1, window=10,
                 dt=0.1,
                 outdir="figures_rmsd_rg",
                 system_name="system"):
    """
    Generate RMSD and Rg plots (cumulative or running‐window average) versus
    time in ns for both the system and a reference, and write a .dat of the
    time series.  No points are dropped.
    """
    os.makedirs(outdir, exist_ok=True)

    # ── Load reference and compute its RMSD & Rg ─────────────────────────
    u_ref = mda.Universe(ref_topology, *ref_traj_list, skip=skip)
    sel_ref = u_ref.select_atoms(sel_str)
    idx_ref = list(range(len(u_ref.trajectory)))
    rmsd_ref_raw = compute_rmsd_sel(u_ref, sel_ref, idx_ref)
    rg_ref_raw   = compute_rg_sel(u_ref, sel_ref, idx_ref)

    # cumulative vs running-window for REF
    if window == 1:
        cr_ref = np.cumsum(rmsd_ref_raw) / (np.arange(len(rmsd_ref_raw)) + 1)
        cg_ref = np.cumsum(rg_ref_raw)   / (np.arange(len(rg_ref_raw))   + 1)
    else:
        cr_ref = moving_average(rmsd_ref_raw, window)
        cg_ref = moving_average(rg_ref_raw,   window)

    frames_ref = np.arange(len(cr_ref))
    time_ref   = frames_ref * skip * dt
    dt_point   = skip * dt

    # ── System branch: replicates vs single trajectory ───────────────────
    if replicates and len(traj_list) > 1:
        # gather per-trajectory arrays
        rmsd_reps = []
        rg_reps   = []
        lengths   = []
        for xtc in traj_list:
            u = mda.Universe(topology, xtc, skip=skip)
            sel = u.select_atoms(sel_str)
            idx = list(range(len(u.trajectory)))
            rmsd_reps.append(compute_rmsd_sel(u, sel, idx))
            rg_reps.append(compute_rg_sel(u, sel, idx))
            lengths.append(len(idx))

        # truncate to the shortest
        Lmin = min(lengths)
        arr_r = np.vstack([r[:Lmin] for r in rmsd_reps])
        arr_g = np.vstack([g[:Lmin] for g in rg_reps])

        rmsd_mean = arr_r.mean(axis=0)
        rmsd_std  = arr_r.std(axis=0)
        rg_mean   = arr_g.mean(axis=0)
        rg_std    = arr_g.std(axis=0)

        # cumulative vs running-window for system
        if window == 1:
            cr_sys = np.cumsum(rmsd_mean) / (np.arange(Lmin) + 1)
            cg_sys = np.cumsum(rg_mean)   / (np.arange(Lmin) + 1)
        else:
            cr_sys = moving_average(rmsd_mean, window)
            cg_sys = moving_average(rg_mean,   window)

        frames_sys = np.arange(len(cr_sys))
        time_sys   = frames_sys * dt_point

        # ── Plot RMSD (mean±std) + REF ────────────────────────────────
        plt.figure(figsize=(10,5))
        plt.plot(time_sys, cr_sys, label=system_name, color='blue')
        plt.fill_between(time_sys,
                         cr_sys - rmsd_std[:len(cr_sys)],
                         cr_sys + rmsd_std[:len(cr_sys)],
                         color='blue', alpha=0.2)
        plt.plot(time_ref, cr_ref, 'k--', label='REF')
        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        plt.title(f"{system_name} — Root Mean Square Deviation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{system_name}_rmsd.pdf"), dpi=300)
        plt.close()

        # ── Plot Rg (mean±std) + REF ─────────────────────────────────
        plt.figure(figsize=(10,5))
        plt.plot(time_sys, cg_sys, label=system_name, color='red')
        plt.fill_between(time_sys,
                         cg_sys - rg_std[:len(cg_sys)],
                         cg_sys + rg_std[:len(cg_sys)],
                         color='red', alpha=0.2)
        plt.plot(time_ref, cg_ref, 'k--', label='REF')
        plt.xlabel("Time (ns)")
        plt.ylabel("Rg (Å)")
        plt.title(f"{system_name} — Radius of Gyration")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{system_name}_rg.pdf"), dpi=300)
        plt.close()

        # ── Write combined time-series .dat ───────────────────────────
        datfn = os.path.join(outdir, f"{system_name}_rmsd_rg_timeseries.dat")
        with open(datfn, 'w') as df:
            df.write("# time_ns\tRMSD_sys\tRMSD_REF\tRg_sys\tRg_REF\n")
            for i, t in enumerate(time_sys):
                df.write(f"{t:.3f}\t"
                         f"{cr_sys[i]:.6f}\t{cr_ref[i]:.6f}\t"
                         f"{cg_sys[i]:.6f}\t{cg_ref[i]:.6f}\n")
        print("Wrote time series →", datfn)

    else:
        # ── Single trajectory cumulative avg or running-window ────────
        u = mda.Universe(topology, *traj_list, skip=skip)
        sel = u.select_atoms(sel_str)
        N = len(u.trajectory)
        idx = list(range(N))
        raw_rmsd = compute_rmsd_sel(u, sel, idx)
        raw_rg   = compute_rg_sel(u, sel, idx)

        if window == 1:
            cr_sys = np.cumsum(raw_rmsd) / (np.arange(N) + 1)
            cg_sys = np.cumsum(raw_rg)   / (np.arange(N) + 1)
        else:
            cr_sys = moving_average(raw_rmsd, window)
            cg_sys = moving_average(raw_rg,   window)

        L = len(cr_sys)
        time_sys = np.arange(L) * dt_point

        # ── Plot RMSD + REF ──────────────────────────────────────────
        plt.figure(figsize=(10,5))
        plt.plot(time_sys, cr_sys, label=system_name, color='blue')
        plt.plot(time_ref, cr_ref, 'k--', label='REF')
        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        plt.title(f"{system_name} — Root Mean Square Deviation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{system_name}_rmsd.pdf"), dpi=300)
        plt.close()

        # ── Plot Rg + REF ───────────────────────────────────────────
        plt.figure(figsize=(10,5))
        plt.plot(time_sys, cg_sys, label=system_name, color='red')
        plt.plot(time_ref, cg_ref, 'k--', label='REF')
        plt.xlabel("Time (ns)")
        plt.ylabel("Rg (Å)")
        plt.title(f"{system_name} — Radius of Gyration")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{system_name}_rg.pdf"), dpi=300)
        plt.close()

        # ── Write combined time-series .dat ─────────────────────────
        datfn = os.path.join(outdir, f"{system_name}_rmsd_rg_timeseries.dat")
        with open(datfn, 'w') as df:
            df.write("# time_ns\tRMSD_sys\tRMSD_REF\tRg_sys\tRg_REF\n")
            for i, t in enumerate(time_sys):
                df.write(f"{t:.3f}\t"
                         f"{cr_sys[i]:.6f}\t{cr_ref[i]:.6f}\t"
                         f"{cg_sys[i]:.6f}\t{cg_ref[i]:.6f}\n")
        print("Wrote time series →", datfn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute & plot RMSD & Rg for system and REF"
    )
    parser.add_argument("--top",     required=True,
                        help="System topology (TPR/GRO)")
    parser.add_argument("--xtc",     required=True, nargs="+",
                        help="System trajectory XTC(s)")
    parser.add_argument("--ref_top", required=True,
                        help="Reference topology (TPR/GRO)")
    parser.add_argument("--ref_xtc", required=True, nargs="+",
                        help="Reference trajectory XTC(s)")
    parser.add_argument("--sel",     default="not resname SOL NA CL NAD",
                        help="Atom selection string")
    parser.add_argument("--no-rep",  action="store_false", dest="replicates",
                        help="Treat all XTCs as one concatenated trajectory")
    parser.add_argument("--skip",    type=int, default=1,
                        help="Frame skip interval")
    parser.add_argument("--window",  type=int, default=10,
                        help="Window length (1=cumulative; >1=running‐window)")
    parser.add_argument("--dt",      type=float, default=0.1,
                        help="Time per raw frame in ns")
    parser.add_argument("--outdir",  default="figures_rmsd_rg",
                        help="Directory for output plots and .dat")
    parser.add_argument("--name",    default="system",
                        help="System name tag for titles and filenames")
    args = parser.parse_args()

    plot_rmsd_rg(
        topology=args.top,
        traj_list=args.xtc,
        ref_topology=args.ref_top,
        ref_traj_list=args.ref_xtc,
        sel_str=args.sel,
        replicates=args.replicates,
        skip=args.skip,
        window=args.window,
        dt=args.dt,
        outdir=args.outdir,
        system_name=args.name
    )
