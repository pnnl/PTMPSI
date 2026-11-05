#!/usr/bin/env python3
"""Tutorial steps 7 and 8: dock a nylon-6 fragment into the ICCG cutinase active site.

Run from a scratch directory with the ptmpsi venv active. Needs internet (PDB 6tht),
`vina` and `pdb2pqr30` on PATH, and `mk_prepare_receptor.py` (installed with Meeko).

The receptor is always protonated first. Without hydrogens Meeko stops with
"Updated 1 H positions but deleted 5", no receptor .pdbqt is written, and
dock() then fails with FileNotFoundError for the Vina output.

  python dock_cutinase.py            # step 7: docking only
  python dock_cutinase.py --qmmm     # also step 8: write NWChem QM/MM inputs
"""
import sys

import ptmpsi
from ptmpsi.protein import Protein


def main():
    qmmm = "--qmmm" in sys.argv

    cutinase = Protein(pdbid="6tht")
    cutinase.mutate("A:ALA130", "SER")
    cutinase.protonate()             # adds hydrogens; also sets cutinase.charge
    cutinase.write_pdb("cutinase.pdb")

    nylon6 = Protein()
    nylon6.prepend("A", "EAC")
    nylon6.prepend("A", "ACE")
    nylon6.append("A", "EAC")
    nylon6.append("A", "NME")
    nylon6.write_xyz("nylon6.xyz")

    cutinase.dock(ligand="nylon6.xyz", receptor="cutinase.pdb",
                  boxcenter="A:SER130", boxsize=20.0)

    if qmmm:
        ptmpsi.docking.write_pdb(cutinase.docking, filename="test2.pdb", pose=0, template=nylon6)
        ptmpsi.nwchem.qmmm_optimize(filename="test2.pdb",
                                    qmres=[130, 207, 259, 260, 261, 262],
                                    counter=True, totcharge=cutinase.charge,
                                    xcfun="pbe96", grid="medium", aobasis="def2-svp")


if __name__ == "__main__":
    main()
