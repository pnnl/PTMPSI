import subprocess, os
import numpy as np
from shutil import which, copyfileobj
from ptmpsi.exceptions import MyDockingError
from ptmpsi.constants import nwchem_input
from ptmpsi.protein.tools import get_residue
from xyz2mol import xyz2mol, read_xyz_file
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTMolecule

def dock_ligand(protein,ligand,receptor,boxcenter,boxsize,output,flexible=None,charge=0,mgltools=None):

    # Check if prepare_receptor and vina are in the path
    if which("prepare_receptor") is None and mgltools is None:
        raise MyDockingError("Cannot find prepare_receptor")
    elif mgltools is not None:
        prepare_receptor = f"{mgltools}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        if not os.path.isFile(prepare_receptor):
            raise MyDockingError("Cannot find prepare_receptor4.py")
    else:
        prepare_receptor = "prepare_receptor"
        if mgltools is None: prepare_receptor = "prepare_receptor"

    if which("vina") is None:
        raise MyDockingError("Cannot find vina")

    # Prepare ligand
    ligandpdbqt = ligand[:-4]+".pdbqt"
    atoms, charge, xyz_coords = read_xyz_file(ligand)
    mol = xyz2mol(atoms, xyz_coords, charge=charge,
            use_graph=True, allow_charged_fragments=True,
            embed_chiral=True, use_huckel=True)
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol[0])
    preparator = MoleculePreparation(hydrate=False,remove_smiles=True)
    preparator.prepare(mol[0])
    preparator.write_pdbqt_file(ligandpdbqt)

    # Prepare receptor
    receptorpdbqt = receptor[:-4]+".pdbqt"
    subprocess.run([prepare_receptor,"-r",receptor,"-o",receptor[:-4]+".pdbqt"])

    if flexible is not None:
        receptorrigid = receptor[:-4]+"_rigid.pdbqt"
        receptorflex = receptor[:-4]+"_flex.pdbqt"
        subprocess.run(["pythonsh","prepare_flexreceptor4.py","-r",receptorpdbqt,"-s",flexible])

    # Vina configuration file
    if boxcenter is None:
        xcenter = 0; ycenter = 0; zcenter = 0
        for chain in protein.chains:
            for residue in chain.residues:
                for xyz in residue.coordinates:
                    xcenter += xyz[0]
                    ycenter += xyz[1]
                    zcenter += xyz[2]
        xcenter /= protein.natoms
        ycenter /= protein.natoms
        zcenter /= protein.natoms
    elif isinstance(boxcenter,str):
        residue = get_residue(protein,boxcenter)
        xcenter, ycenter, zcenter = np.mean(residue.coordinates,axis=0)
    else:
        [xcenter,ycenter,zcenter] = boxcenter


    if isinstance(boxsize,float) or isinstance(boxsize,int):
        xsize = ysize = zsize = float(boxsize)
    elif len(boxsize) == 3:
        [xsize,ysize,zsize] = boxsize
    else:
        raise MyDockingError("Boxsize must be a 3-vector")

    with open("config.txt","w") as fh:
        fh.write("""
size_x = {}
size_y = {}
size_z = {}
center_x = {}
center_y = {}
center_z = {}
""".format(xsize,ysize,zsize,xcenter,ycenter,zcenter))


    # Dock structure 
    output = "_".join([receptor[:-4],ligand[:-4],"vina.pdbqt"])
    if flexible is None:
        subprocess.run(["vina",
        "--receptor",receptorpdbqt,
        "--ligand",ligandpdbqt,
        "--config","config.txt",
        "--exhaustiveness=32",
        "--cpu=8",
        "--num_modes=9",
        "--out",output])
    else:
        subprocess.run(["vina",
        "--receptor",receptorrigid,
        "--flex",receptorflex,
        "--ligand",ligandpdbqt,
        "--config","config.txt",
        "--exhaustiveness=32",
        "--cpu=8",
        "--num_modes=9",
        "--out",output])


    return
