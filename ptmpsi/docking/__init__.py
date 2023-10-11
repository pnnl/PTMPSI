import subprocess, os
import re
import numpy as np
from shutil import which, copyfileobj
from ptmpsi.exceptions import MyDockingError
from ptmpsi.constants import nwchem_input
from ptmpsi.protein.tools import get_residue
from ptmpsi.math import norm
from xyz2mol import xyz2mol, read_xyz_file
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate

class Dock:
    def __init__(self):
        self.docked = False
        self.ligand = None
        self.receptor = None
        self.flexible = None
        self.output = None
        self.engine = "vina"
        self.boxcenter = None
        self.boxsize = None
        self.exhaustiveness = 32

def dock_ligand(cls):
    receptorpdbqt = cls.receptor[:-4]+".pdbqt"

    if which(cls.engine) is None:
        raise MyDockingError(f"Cannot find {cls.engine}")

    # Prepare ligand
    ligandpdbqt = cls.ligand[:-4]+".pdbqt"
    atoms, charge, xyz_coords = read_xyz_file(cls.ligand)
    mol = xyz2mol(atoms, xyz_coords, charge=charge,
            use_graph=True, allow_charged_fragments=True,
            embed_chiral=True, use_huckel=True)
    mol = mol[0]
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
    preparator = MoleculePreparation(hydrate=False)
    mol_setups = preparator.prepare(mol)
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
        if is_ok:
            with open(ligandpdbqt,"w") as fh: fh.write(pdbqt_string)
        else:
            raise MyDockingError(error_msg)

    # Prepare receptor
    command = ["mk_prepare_receptor.py","--pdb",cls.receptor,"-o",f"{cls.receptor[:-4]}","--skip_gpf"]
    receptorrigid = cls.receptor[:-4]+"_rigid.pdbqt"
    receptorflex = cls.receptor[:-4]+"_flex.pdbqt"
    if cls.flexible is not None:
        for flexres in cls.flexible:
            command.extend(["-f",flexres])
    subprocess.run(command)

    # Vina configuration file
    [xcenter,ycenter,zcenter] = cls.boxcenter

    if isinstance(cls.boxsize,float) or isinstance(cls.boxsize,int):
        xsize = ysize = zsize = float(cls.boxsize)
    elif len(cls.boxsize) == 3:
        [xsize,ysize,zsize] = cls.boxsize
    else:
        raise MyDockingError("Boxsize must be a 3-vector")

    with open("config.txt","w") as fh:
        fh.write(f"""
size_x = {xsize}
size_y = {ysize}
size_z = {zsize}
center_x = {xcenter}
center_y = {ycenter}
center_z = {zcenter}""")


    # Dock structure 
    if cls.flexible is None:
        subprocess.run(["vina",
        "--receptor",receptorpdbqt,
        "--ligand",ligandpdbqt,
        "--config","config.txt",
        f"--exhaustiveness={cls.exhaustiveness}",
        "--cpu=8",
        "--num_modes=9",
        "--out",cls.output])
    else:
        subprocess.run(["vina",
        "--receptor",receptorrigid,
        "--flex",receptorflex,
        "--ligand",ligandpdbqt,
        "--config","config.txt",
        f"--exhaustiveness={cls.exhaustiveness}",
        "--cpu=8",
        "--num_modes=9",
        "--out",cls.output])

    # We got a docked structure
    cls.docked = True

    # Extract ligand results into XYZ file
    with open(cls.output,"r") as fh: results_pdbqt = fh.read()
    pdbqt_mol = PDBQTMolecule.from_file(cls.output, skip_typing=True)
    docked_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    for i in range(docked_mols[0].GetNumConformers()):
        Chem.MolToXYZFile(docked_mols[0], 
                filename=f"{cls.receptor[:-4]}_{cls.ligand[:-4]}_pose{i+1:02d}_vina.xyz", confId=i)

    return

def write_pdb(cls, filename=None, pose=0, template=None):
    if filename is None:
        filename = f"{cls.receptor[:-4]}_{cls.ligand[:-4]}_complex.pdb"

    with open(filename,"w") as outf:
        if not cls.flexible:
            with open(cls.receptor,"r") as recf: reclines = recf.readlines()
            natom = 0
            ires  = 0
            for line in reclines:
                if line[:4] == "ATOM":
                    natom += 1
                    newline = line[:56] + "0.00  0.00"
                    ires = int(line[23:26])
                    if line[76:78] == "  ":
                        if re.search(r"C[A,B,G,D,E,Z,H]| C$",line[12:16]):
                            newline += "           C\n"
                        elif re.search(r"O[G,D,E,H,X]| O$",line[12:16]):
                            newline += "           O\n"
                        elif re.search(r"N[H,Z,E,D]| N$",line[12:16]):
                            newline += "           N\n"
                        elif re.search(r"S[G,D,H]| S$",line[12:16]):
                            newline += "           S\n"
                        elif re.search(r"H[H,A,B,G,Z,E,D,O]| H$",line[12:16]):
                            newline += "           H\n"
                    else:
                        newline += line[66:]
                    outf.write(newline)
            with open(f"{cls.receptor[:-4]}_{cls.ligand[:-4]}_pose{pose+1:02d}_vina.xyz","r") as ligf:
                liglines = ligf.readlines()
            iatom = 0
            if template is None:
                ires += 1
                for line in liglines[2:]:
                    natom += 1
                    iatom += 1
                    element = line.split()[0] + str(iatom)
                    coords = [float(x) for x in line.split()[1:]]
                    newline = f"ATOM  {natom:5d} {element:>4s} LIG   {ires:3d}    {coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  0.00  0.00          {line.split()[0]:>2s}\n"
                    outf.write(newline)
            else:
                iline = 2
                for chain in template.chains:
                    for residue in chain.residues:
                        ires += 1
                        for element,name,coord in zip(residue.elements,residue.names,residue.coordinates):
                            if element == "H":
                                if liglines[iline].split()[0] == "H": iline += 1
                                continue
                            else:
                                natom += 1
                                coords = np.array([float(x) for x in liglines[iline].split()[1:]])
                                outf.write(f"ATOM  {natom:5d} {name:>4s} {residue.name:3s}   {ires:3d}    {coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  0.00  0.00          {element:>2s}\n")
                                iline += 1
                                if name in ["C","O"]: continue
                                ih = 0
                                for jline in range(2,len(liglines)):
                                    if liglines[jline].split()[0] == "H":
                                        jcoord = np.array([float(x) for x in liglines[jline].split()[1:]])
                                        if norm(coords-jcoord) < 1.2:
                                            natom += 1
                                            ih += 1
                                            if name in ["NH3","CH3"]:
                                                hname = "HH3" + str(ih)
                                            elif name in ["N"]:
                                                hname = "H"
                                            else:
                                                hname = "H" + name[1] + str(ih+1)
                                            outf.write(f"ATOM  {natom:5d} {hname:>4s} {residue.name}   {ires:3d}    {jcoord[0]:8.3f}{jcoord[1]:8.3f}{jcoord[2]:8.3f}  0.00  0.00           H\n")
                                            if hname == "H" and ih == 1: break
                                            elif "HH3" in hname and ih == 3: break
                                            elif "HH3" not in hname and ih == 2: break




            


