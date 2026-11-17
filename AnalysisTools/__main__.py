"""Command-line interface for PTM-Psi trajectory preprocessing and analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from .trajectory import check_pbc_jumps, load_universe, save_standard_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m AnalysisTools",
        description="Check, unwrap and analyze PTM-Psi GROMACS trajectories.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_unwrap = sub.add_parser(
        "unwrap",
        help="Make whole, remove jumps and fit via ptmpsi.gromacs.analysis.post_process (gmx_mpi or gmx).",
    )
    p_unwrap.add_argument("-s", "--topology", required=True, help="Input GROMACS .tpr.")
    p_unwrap.add_argument("-f", "--trajectory", required=True, help="Input raw .xtc trajectory.")
    p_unwrap.add_argument("-o", "--output-dir", required=True, help="Directory for the processed trajectory and .gro/.tpr.")
    p_unwrap.add_argument("--group", default="Protein", help="gmx group written out (default Protein).")
    p_unwrap.add_argument("--fit-group", type=int, default=4, help="gmx group index the fit uses (3 C-alpha, 4 Backbone; default 4).")
    p_unwrap.add_argument("--gmx", default=None, help="Explicit gmx or gmx_mpi binary.")

    p_run = sub.add_parser("analyze", help="Run the PTM-aware contact map and mass-weighted Cartesian PCA.")
    p_run.add_argument("-s", "--topology", required=True, help="GROMACS .tpr topology file.")
    p_run.add_argument("-f", "--trajectory", required=True, help="Unwrapped .xtc trajectory file.")
    p_run.add_argument("--contacts-trajectory", help="Non-fitted .xtc for minimum-image contacts (default: --trajectory).")
    p_run.add_argument("-o", "--output-dir", required=True, help="Directory for CSV, NPY, and PNG outputs.")
    p_run.add_argument("--prefix", default="analysis", help="Filename prefix inside --output-dir.")
    p_run.add_argument("--contact-rep", choices=("ca", "heavy"), default="ca")
    p_run.add_argument("--n-components", type=int, default=3)
    p_run.add_argument("--no-plots", action="store_true", help="Skip PNG summary plot generation.")

    args = parser.parse_args(argv)
    if args.command == "unwrap":
        from ptmpsi.gromacs.analysis import post_process

        try:
            outputs = post_process(
                args.trajectory,
                args.topology,
                prefix=args.output_dir,
                group=args.group,
                fit_group=args.fit_group,
                gmx=args.gmx,
                interactive=True,
            )
        except (RuntimeError, OSError) as error:
            parser.exit(1, f"Unwrap failed: {error}\n")
        for out in outputs:
            name = str(out).removesuffix("_whole_nojump_fit.xtc")
            print(f"fitted_trajectory: {out}")
            print(f"topology: {name}_{args.group}.tpr")
            print(f"contacts_trajectory: {name}_whole_nojump.xtc")
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
        contacts_trajectory=args.contacts_trajectory,
    )
    for name, path in written.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
