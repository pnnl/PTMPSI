#!/usr/bin/env python3
"""Turn a generated cluster script (<NNNN>_slurm.sbatch) into a plain shell script
that runs the same GROMACS steps on a workstation (macOS or Linux, no Slurm, no MPI
launcher, no container, no GPU flags).

  python localize_sbatch.py <uid>/1tuples/0000/0000_slurm.sbatch            > run_local.sh
  python localize_sbatch.py <uid>/1tuples/0000/0000_slurm.sbatch --smoke 10 > smoke.sh

Run the output from inside the combination folder (the one that holds 0000.pdb).

--smoke N adds "-nsteps N" to every mdrun, so the whole chain can be checked in minutes.
Without it, the heating, NPT and free-NPT stages run for their full length.

The script stops after "grompp ... md.tpr". The production run is in md.sbatch on the
cluster; locally run it yourself:  gmx mdrun -deffnm md   (add -cpi md.cpt to resume).
"""
import argparse
import re
import sys

LAUNCHER = re.compile(r"\b(?:srun|mpirun|mpiexec|flux run)\b[^|;]*?(?=\bgmx(?:_mpi)?\b)")
GPU_FLAGS = [
    r"-gpu_id\s+\S+", r"-gputasks\s+\S+", r"-nb\s+gpu", r"-pme\s+gpu",
    r"-npme\s+\d+", r"-bonded\s+gpu", r"-update\s+gpu",
]


SED_INSERT = re.compile(r"^sed -i '/; Include water topology/i #ifdef P(\d+)")


def perl_insert(force):
    # BSD sed (macOS) cannot do this edit; perl behaves the same everywhere.
    return ('perl -0pi -e \'s/; Include water topology/#ifdef P%s\\n  #include "%sp.itp"\\n#endif\\n'
            '; Include water topology/\' topol.top' % (force, force))


def localize(line, smoke):
    m = SED_INSERT.match(line)
    if m:
        return perl_insert(m.group(1))
    line = LAUNCHER.sub("", line)
    line = re.sub(r"\bgmx_mpi\b", "gmx", line)
    for flag in GPU_FLAGS:
        line = re.sub(flag, "", line)
    if smoke and " mdrun " in line:
        line = line.rstrip() + f" -nsteps {smoke} -ntmpi 1"
    return re.sub(r"[ \t]{2,}", " ", line).rstrip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sbatch")
    ap.add_argument("--smoke", type=int, default=0, metavar="N")
    args = ap.parse_args()

    out = ["#!/bin/bash", "set -euo pipefail", "export GMX_MAXBACKUP=-1"]
    for raw in open(args.sbatch).read().splitlines():
        raw = raw.strip()
        gromacs = re.match(
            r"^(echo\b.*\|\s*)?(?:srun|mpirun|mpiexec|flux run|gmx(?:_mpi)?)\b", raw
        ) and re.search(r"\bgmx(?:_mpi)?\s+\w+\b", raw)
        keep = raw.startswith(("perl -0pi", "sed -i '/; Include water topology")) or gromacs
        if keep:
            out.append(localize(raw, args.smoke))
        if gromacs and re.search(r"\bgrompp\b.*(?:^|\s)-o\s+md\.tpr(?:\s|$)", raw):
            break
    else:
        sys.exit("did not find the 'grompp ... -o md.tpr' line; is this a <NNNN>_slurm.sbatch file?")
    print("\n".join(out))


if __name__ == "__main__":
    main()
