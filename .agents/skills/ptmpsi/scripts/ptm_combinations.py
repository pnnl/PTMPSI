#!/usr/bin/env python3
"""Tutorial step 9: PTM combinations for PDB 6gve and a GROMACS run directory.

Run from a scratch directory with the ptmpsi venv active (pdb2pqr30 on PATH).
Needs internet for the PDB download.

  python ptm_combinations.py            # list combinations only
  python ptm_combinations.py --generate # also write the run directory
"""
import sys

import ptmpsi
from ptmpsi.protein import Protein


def print_tuple_table(combos):
    for k, level in enumerate(combos, start=1):
        print(f"{k}-tuple ({len(level)} combinations)")
        for i, combo in enumerate(level, start=1):
            print(f"  {i:3d}  " + "   ".join(f"{res} {ptm}" for res, ptm in combo))


def main():
    protein = Protein(pdbid="6gve")
    protein.findresidue("CYS")

    both = [("CYS", "nitrosylation"), ("CYS", "sulfhydration")]
    c1 = protein.get_ptm_combinations(both, ntuple=1)
    exclude = [f"{ch}:{r}" for ch in "ADEFIKNO" for r in ("CYS19", "CYS77")]
    c2 = protein.get_ptm_combinations(both, exclude=exclude, ntuple=2)

    protein.findresidue("CYX")
    c3 = protein.get_ptm_combinations([("M:CYX", "reduction")], ntuple=4)

    print_tuple_table(c3)
    print("sizes:", [len(x) for x in c1], [len(x) for x in c2], [len(x) for x in c3])
    assert [len(x) for x in c3] == [4, 6, 4, 1], "unexpected combination counts for M:CYX"

    if "--generate" in sys.argv:
        opts = ptmpsi.slurm.aqe_ldrd.partitions["premium"].options["gromacs"]
        uid = protein.gen_ptm_files([c3[3]], machine=ptmpsi.slurm.aqe_ldrd, **opts)
        print("run directory:", uid)


if __name__ == "__main__":
    main()
