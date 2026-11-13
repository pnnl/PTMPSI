"""Command-line interface for PTM-Psi trajectory preprocessing and analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from .trajectory import check_pbc_jumps, load_universe, save_standard_report, unwrap_trajectory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m AnalysisTools",
        description="Unwrap GROMACS trajectories and run PTM-aware MDAnalysis structural analyses.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_unwrap = sub.add_parser("unwrap", help="Unwrap PBC jumps via gmx/gmx_mpi trjconv (-pbc whole -> -pbc nojump).")
    p_unwrap.add_argument("-s", "--topology", required=True, help="Input GROMACS .tpr (preferred) or structure file.")
    p_unwrap.add_argument("-f", "--trajectory", required=True, help="Input raw .xtc trajectory.")
    p_unwrap.add_argument("-o", "--output", required=True, help="Output unwrapped .xtc path.")
    p_unwrap.add_argument("--gmx", default=None, help="Explicit gmx or gmx_mpi binary.")
    p_unwrap.add_argument("--fit-group", default=None, help="Optional trjconv -fit rot+trans group (e.g. Backbone).")

    p_run = sub.add_parser("analyze", help="Run RMSD, mass-weighted RoG, RMSF, contact map, and mass-weighted PCA.")
    p_run.add_argument("-s", "--topology", required=True, help="GROMACS .tpr topology file.")
    p_run.add_argument("-f", "--trajectory", required=True, help="Unwrapped .xtc trajectory file.")
    p_run.add_argument("-o", "--output-dir", required=True, help="Directory for CSV, NPY, and PNG outputs.")
    p_run.add_argument("--prefix", default="analysis", help="Filename prefix inside --output-dir.")
    p_run.add_argument("--contact-rep", choices=("ca", "heavy"), default="ca")
    p_run.add_argument("--n-components", type=int, default=3)
    p_run.add_argument("--no-plots", action="store_true", help="Skip PNG summary plot generation.")

    args = parser.parse_args(argv)
    if args.command == "unwrap":
        out = unwrap_trajectory(
            args.topology,
            args.trajectory,
            args.output,
            gmx=args.gmx,
            fit_group=args.fit_group,
        )
        print(f"Wrote unwrapped trajectory: {out}")
        return 0

    u = load_universe(args.topology, args.trajectory)
    pbc_report = check_pbc_jumps(u)
    if pbc_report["has_pbc_jumps"]:
        print(
            "WARNING: Possible un-unwrapped PBC jumps detected "
            f"(max step {pbc_report['max_step_angstrom']:.1f} Å; "
            f"split-chain frames: {len(pbc_report['split_chain_frames'])}). "
            "Run 'python -m AnalysisTools unwrap' first for continuous coordinates."
        )
    written = save_standard_report(
        args.topology,
        args.trajectory,
        Path(args.output_dir),
        prefix=args.prefix,
        contact_rep=args.contact_rep,
        n_pca_components=args.n_components,
        make_plots=not args.no_plots,
    )
    for name, path in written.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
