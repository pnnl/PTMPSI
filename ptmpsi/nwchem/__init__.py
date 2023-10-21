from ptmpsi.protein import Protein, Chain
from ptmpsi.residues import Residue
from ptmpsi.nwchem.templates import * 
from ptmpsi.nwchem.bonded import bondedks
from ptmpsi.nwchem.reader import readoptim
from ptmpsi.nwchem.qmmm import *
from ptmpsi.math import get_torsion, norm
from ptmpsi.constants import ang2bohr, covrad
import copy
import numpy as np
import os
import subprocess

ACE = {"CH3" : -0.366200,
          "C":  0.597200,
          "O": -0.567900,
       "HH31":  0.112300,
       "HH32":  0.112300,
       "HH33":  0.112300}
NME = {"N"   : -0.415700,
       "H"   :  0.271900,
       "CH3" : -0.149000,
       "HH31":  0.097600,
       "HH32":  0.097600,
       "HH33":  0.097600}
amide = {"N": -0.415700,
         "H":  0.271900,
         "C":  0.597200,
         "O": -0.567900}

class Slurm:
    def __init__(self,**kwargs):
        self.time      = kwargs.get("time","12:00:00")
        self.partition = kwargs.get("partition","normal")
        self.account   = kwargs.get("account","emsls60202")
        self.nnodes    = kwargs.get("nnodes",1)
        self.ntasks    = kwargs.get("ntasks",36)
        self.scratch   = kwargs.get("scratch","/big_scratch")
        self.jobname   = kwargs.get("jobname","ptmpsi")
        self.header    = slurm_header.format(
            account=self.account, time=self.time, nodes=self.nnodes,
            ntasks=self.ntasks, jname=self.jobname, partition=self.partition,
            np=self.nnodes*self.ntasks, scratch=self.scratch)
        return

class NWChem:
    def __init__(self,**kwargs):
        self.mult    = kwargs.get("mult",1)
        self.charge  = kwargs.get("charge",0)
        self.memory  = kwargs.get("memory",2000)
        self.aobasis = kwargs.get("aobasis","def2-tzvp")
        self.tdbasis = kwargs.get("tdbasis","def2-svp")
        self.cdbasis = kwargs.get("cdbasis","def2-universal-jfit")
        self.xcfun   = kwargs.get("xcfun","r2scan")
        self.grid    = kwargs.get("grid","lebedev 120 14")
        self.tdgrid  = kwargs.get("tdgrid","lebedev 100 14")
        self.nscf    = kwargs.get("nscf",100)
        self.nopt    = kwargs.get("nopt",60)
        self.disp    = kwargs.get("disp","disp vdw 4")
        self.delta   = kwargs.get("delta",0.0189)
        self.lshift  = kwargs.get("lshift",0.1)
        self.bqzone  = kwargs.get("bqzone",20.0)
        self.boxsize = kwargs.get("boxsize",10.0)
        return

class TorsionDrive:
    def __init__(self,**kwargs):
        self.spacing  = kwargs.get("spacing",10)
        self.torsions = kwargs.get("torsions",[])
        for torsion in self.torsions:
            assert len(torsion) == 4, f"Invalid torsion {torsion}"
        self.torsions = [[x-1 for x in torsion] for torsion in self.torsions]



def get_qdata(files):
    njobs = len(files)
    with open("qdata.txt","w") as qdata:
        for ifile,filename in enumerate(files):
            coords = ""
            forces = ""
            with open(filename,"r") as logfile:
                found_gradients, energy = False, None
                for line in logfile:
                    line = line.strip()
                    if line.startswith("Total DFT energy"):
                        ls = line.split()
                        energy = float(ls[4])
                    elif line.startswith("DFT ENERGY GRADIENTS"):
                        for i in range(3): line = next(logfile)
                        found_gradients = True
                    elif found_gradients:
                        ls = line.split()
                        if len(ls) == 0: break
                        coords += f" {float(ls[2])/ang2bohr}  {float(ls[3])/ang2bohr}  {float(ls[4])/ang2bohr} "
                        forces += f" {float(ls[5])}  {float(ls[6])}  {float(ls[7])} "
            assert found_gradients, "Could not find Gradients"
            assert energy is not None, "Could not find energy"
            qdata.write(f"JOB {ifile}\n")
            qdata.write(f"COORDS {coords} \n")
            qdata.write(f"ENERGY {energy} \n")
            qdata.write(f"FORCES {forces} \n\n")


def get_qm_data(residue,ligand=False,metal=False,ff="AMBER99",**kwargs):
    """ Get all QM data needed to parameterize a non-standard
    amino acid or a new ligand.
    """

    # Process all keyword arguments
    slurm  = Slurm(**kwargs)
    nwchem = NWChem(**kwargs)
    tdrive = TorsionDrive(**kwargs)

    if ligand and metal:
        raise KeyError("ligand and metal cannot be both True")
    single = not (ligand or metal) 

    # Generate the appropriate conformations
    alpha = copy.deepcopy(residue)
    if ligand:
        if isinstance(alpha, Protein):
            conformers = alpha.chains[0].residues
        elif isinstance(alpha, Chain):
            conformers = alpha.residues
        elif isinstance(alpha, Residue):
            conformers = [alpha]
        else:
            raise KeyError("residue must be a Protein, Chain, or Residue instance")
        if len(conformers) > 1:
            raise KeyError("Ligand parameterization only accepts one residue")
    elif metal:
        _alpha = []
        if isinstance(alpha, Protein):
            for _chain in alpha.chains:
                for _residue in _chain.residues:
                    _alpha.append(_residue)
        elif isinstance(alpha, Chain):
            for _residue in alpha.residues:
                _alpha.append(_residue)
        elif isinstance(alpha, Residue):
            _alpha.append(alpha)
        else:
            raise KeyError("residue must be a Protein, Chain, or Residue instance")

        alpha = copy.deepcopy(_alpha)
        conformers = [alpha]
    else:
        if not isinstance(alpha, Protein):
            raise KeyError("residue must be a Protein instance")

        alpha.prepend(chain="A",residue="ACE",phi=-60.0)
        alpha.append(chain="A",residue="NME",psi=-45.0)
        alpha = alpha.chains[0].residues
        beta = copy.deepcopy(residue)
        beta.prepend(chain="A",residue="ACE",phi=-135.0)
        beta.append(chain="A",residue="NME",psi=135.0)
        beta = beta.chains[0].residues
        conformers = [alpha, beta]

    # Get phi and psi indices
    phi, psi = None, None
    if single:
        offset1 = len(alpha[0].names)
        phi = [alpha[0].find("C"),
               offset1+alpha[1].find("N"),
               offset1+alpha[1].find("CA"),
               offset1+alpha[1].find("C")]
        offset2 = len(alpha[1].names)
        psi = [phi[1], phi[2], phi[3],
               offset1+offset2+alpha[2].find("N")]
    
    # Default torsion
    if single and len(tdrive.torsions) == 0:
        if alpha[1].chi2 is not None:
            tdrive.torsions = [offset1+alpha[1].chi2]
    dotorsions = False if len(tdrive.torsions)==0 else True

    # Update torsion list
    _list = [phi, psi]
    for torsion in tdrive.torsions: _list.append(torsion)
    tdrive.torsions = _list

    # Get coordinates and number of atoms
    names  = [name for residue in alpha for name in residue.names]
    elems  = [element for residue in alpha for element in residue.elements]
    coords = [np.array([coord for residue in alpha for coord in residue.coordinates])]
    if single:
        coords.append(np.array([coord for residue in beta for coord in residue.coordinates]))
    natoms = len(names)

    # Get ESP constraints
    consnw, conspy, printpy, ncons = "", "", "", 0
    if not ligand:
        iatom = 0
        for _residue in alpha:
            for atom in _residue.names:
                iatom += 1
                if _residue.name in ["ACE"]:
                    ncons   += 1
                    consnw  += f" constrain  {ACE[atom]: 10.6f}  {iatom:5d}\n"
                    conspy  += f"[{iatom}, {ACE[atom]}],\n"
                elif _residue.name in ["NME"]:
                    ncons += 1
                    consnw += f" constrain  {NME[atom]: 10.6f}  {iatom:5d}\n"
                    conspy += f"[{iatom}, {NME[atom]}],\n"
                else:
                    if ff == "AMBER99" and atom in amide:
                        ncons += 1
                        consnw += f" constrain  {amide[atom]: 10.6f}  {iatom:5d}\n"
                        conspy += f"[{iatom}, {amide[atom]}],\n"
                printpy += f"""print(f"{iatom: 3d} {atom:4s}: {{q[{iatom-1}]: 10.5f}}")\n"""

    else:
        iatom = 0
        for ires in alpha:
            for atom in ires.names:
                iatom += 1
                printpy += f"""print(f"{iatom: 3d} {atom:4s}: {{q[{iatom-1}]: 10.5f}}")\n"""

    # Write files
    for idx,coord in enumerate(coords):
        geometry = ""
        for ele,atom in zip(elems,coord):
            geometry += f"{ele}   {atom[0]: 14.8f}   {atom[1]: 14.8f}   {atom[2]: 14.8f}\n"
        geometry = geometry[:-1]
        filename = f"conf{str(idx)}.nw"
        zcoord   = ""
        gcons    = ""
        noautoz  = ""
        if single:
            phi_val, psi_val = get_torsion(*coord[phi]), get_torsion(*coord[psi])
            _phi = [x+1 for x in phi]
            _psi = [x+1 for x in psi]
            zcoord  += " zcoord\n"
            zcoord  += "  torsion {} {} {} {} {:8.3f} constant\n".format(*_phi,phi_val)
            zcoord  += "  torsion {} {} {} {} {:8.3f} constant\n".format(*_psi,psi_val)
            if alpha[1].chi2 is not None:
                zcoord  += "  torsion {} {} {} {}\n end".format(*tdrive.torsions[-1]+1,120.0)
            else:
                zcoord  += " end"
        elif metal:
            noautoz = "noautoz"
            gcons += "constraints\n"
            iatom = 0
            for _residue in alpha:
                for atom in _residue.names:
                    iatom +=1
                    if atom in ["CA","ZN","CU","FE","CO"]:
                        gcons += f"  fix atom {iatom}\n"
            gcons += "end"
        with open(filename,"w") as infile:
            infile.write(respnw.format(
                name=f"conf{str(idx)}", zcoord=zcoord, charge=nwchem.charge,
                mult=nwchem.mult, memory=nwchem.memory, xcfun=nwchem.xcfun,
                grid=nwchem.grid, aobasis=nwchem.aobasis, cdbasis=nwchem.cdbasis,
                nscf=nwchem.nscf, nopt=nwchem.nopt, disp=nwchem.disp, gcons=gcons,
                constraint=consnw,geometry=geometry,lshift=nwchem.lshift,
                noautoz=noautoz))
        filename = filename[:-3] + "_hess.nw"
        _geometry = f""" load "conf{str(idx)}.xyz" """
        with open(filename,"w") as infile:
            infile.write(hessnw.format(
                name=f"conf{str(idx)}-hess", zcoord=zcoord, charge=nwchem.charge,
                mult=nwchem.mult, memory=nwchem.memory, xcfun=nwchem.xcfun,
                grid=nwchem.grid, aobasis=nwchem.aobasis, cdbasis=nwchem.cdbasis,
                nscf=nwchem.nscf, disp=nwchem.disp, delta=nwchem.delta,
                geometry=_geometry))
        if dotorsions and single:
            filename = f"conf{str(idx)}_tdrive.nw"
            with open(filename,"w") as infile:
                infile.write(torsnw.format(
                    name=f"conf{str(idx)}-tdrive", charge=nwchem.charge,
                    mult=nwchem.mult, memory=nwchem.memory, xcfun=nwchem.xcfun,
                    grid=nwchem.tdgrid, aobasis=nwchem.tdbasis, cdbasis=nwchem.cdbasis,
                    nscf=nwchem.nscf, disp=nwchem.disp, lshift=nwchem.lshift,
                    geometry="@geometry@"))
            with open(f"dihedrals{str(idx)}.txt","w") as infile:
                infile.write("{} {} {} {}".format(*tdrive.torsions[-1]+1))
            with open(f"extras{str(idx)}.txt","w") as infile:
                infile.write("$set\n dihedral  {}  {}  {}  {}  {}\n".format(*_phi,round(phi_val)))
                infile.write(" dihedral  {}  {}  {}  {}  {}".format(*_psi,round(psi_val)))

    with open(f"get_qm_data.sbatch","w") as infile:
        infile.write(slurm.header)
        infile.write(slurm_copy.format(filename="espfit.py"))
        infile.write(slurm_copy.format(filename="modseminario.py"))
        for idx in range(len(coords)):
            tail = f"conf{str(idx)}"
            infile.write(slurm_copy.format(filename=f"{tail}.nw"))
            infile.write(f"""echo "Running conf{str(idx)} optimization"\n""")
            infile.write(runsingularity.format(scratch=slurm.scratch,name=tail))
           
            infile.write(slurm_copy.format(filename=f"{tail}_hess.nw"))
            infile.write(f"""echo "Running conf{str(idx)} hessian"\n""")
            infile.write(runsingularity.format(scratch=slurm.scratch,name=f"{tail}_hess"))
        infile.write("\npython espfit.py\n")
        infile.write("\npython modseminario.py\n")
        if dotorsions and single: 
            infile.write(slurm_tdrive.format(scratch=slurm.scratch))
        infile.write("\n")
        if dotorsions and single:
            for idx in range(len(coords)):
                tail = f"conf{str(idx)}"
                infile.write(slurm_copy.format(filename=f"dihedrals{str(idx)}.txt"))
                infile.write(slurm_copy.format(filename=f"extras{str(idx)}.txt"))
                infile.write(slurm_copy.format(filename=f"{tail}_tdrive.nw"))
                infile.write(slurm_tdrive_run.format(
                    tail=tail, charge=nwchem.charge, mult=nwchem.mult, aobasis=nwchem.aobasis,
                    cdbasis=nwchem.cdbasis, nscf=nwchem.nscf, xcfun=nwchem.xcfun,
                    grid=nwchem.grid, disp=nwchem.disp, memory=nwchem.memory,
                    spacing=tdrive.spacing, idx=str(idx)))
        infile.write("\n cleanup")
    with open("espfit.py","w") as fitting:
        hcons = []
        for i,iname in enumerate(names):
            if iname[0:2] not in ["HA","HB","HG","HD","HE","HH","H1","H2","H3"]: continue
            if metal and iname[0:2] in ["HH"]: continue

            bonded = None
            for j in range(natoms):
                if i == j: continue
                if norm(coords[0][i]-coords[0][j]) < 1.15:
                    bonded = j
                    break
            
            if bonded is None:
                raise KeyError(f"Could not find bond to hydrogen atom {iname}")

            for j,jname in enumerate(names):
                if j == i: continue
                if jname[0:2] not in ["HA","HB","HG","HD","HE","HH","H1","H2","H3"]: continue
                if norm(coords[0][j]-coords[0][bonded]) < 1.15:
                    found = False
                    for icons in hcons:
                        if icons[0] == i and icons[1] == j: found = True
                        if icons[0] == j and icons[1] == i: found = True
                        if found: break
                    if found: continue
                    hcons.append([i,j])

        fitting.write(fit.format(
            natoms=natoms, ncons=ncons, nconf=len(coords), charge=nwchem.charge,
            printing=printpy, names=names, cons=conspy, hcons=hcons,
            files=[f"conf{str(idx)}" for idx in range(len(coords))]))
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),"bonds_and_angles.py")
    with open(filename,"r") as fh: modseminario = fh.read()
    with open("modseminario.py","w") as seminario:
        seminario.write(modseminario.format(
            elements=elems, names=names, nconf=len(coords)))
    return

def get_bond_graph(coords, elems, factor=1.3):
    natoms = len(coords)
    bonds = list()
    for iatom in range(natoms):
        for jatom in range(iatom+1,natoms):
            rcovsum = covrad[elems[iatom]] + covrad[elems[jatom]]
            if norm(coords[iatom]-coords[jatom]) <= factor*rcovsum:
                bonds.append((iatom,jatom))
    return bonds

def get_angle_graph(bonds):
    angles = list()
    _bonds = copy.copy(bonds)
    if len(_bonds) < 3: return angles

    for i,ibond in enumerate(bonds):
        if len(_bonds) < 3: break
        iatom,jatom = _bonds.pop(0)
        for jbond in _bonds:
            if jatom not in jbond: continue
            katom = jbond[0] if jatom == jbond[1] else jbond[1]
            angles.append((iatom,jatom,katom))
    return angles


def get_torsion_graph(bonds):
    torsions = list()
    _bonds = copy.copy(bonds)
    if len(_bonds) < 4: return torsions

    
    for i,ibond in enumerate(bonds):
        if len(_bonds) < 4: break
        iatom,jatom = _bonds.pop(0)
        for jbond in _bonds:
            if jatom not in jbond: continue
            katom = jbond[0] if jatom == jbond[1] else jbond[1]
            for kbond in _bonds:
                if kbond == jbond: continue
                if katom not in kbond: continue
                latom = kbond[0] if katom == kbond[1] else kbond[1]
                torsions.append((iatom,jatom,katom,latom))
    return torsions


def forcebalance(residue,**kwargs):
    import os
    import json
    import subprocess

    # Process all keyword arguments
    slurm  = Slurm(**kwargs)
    tdrive = TorsionDrive(**kwargs)

    suffix = kwargs.get("suffix","")

    # Generate a copy of the Residue
    alpha = copy.deepcopy(residue)
    alpha.prepend(chain="A",residue="ACE",phi=-60.0)
    alpha.append(chain="A",residue="NME",psi=-45.0)
    from ptmpsi.gromacs.utils import amber_to_gromacs_names
    amber_to_gromacs_names(alpha)
    alpha = alpha.chains[0].residues
    
    # Get phi and psi indices
    phi, psi = None, None
    offset1 = len(alpha[0].names)
    phi = [alpha[0].find("C"),
               offset1+alpha[1].find("N"),
               offset1+alpha[1].find("CA"),
               offset1+alpha[1].find("C")]
    offset2 = len(alpha[1].names)
    psi = [phi[1], phi[2], phi[3],
               offset1+offset2+alpha[2].find("N")]

    # Default torsion
    if len(tdrive.torsions) == 0:
        if alpha[1].chi2 is not None:
            tdrive.torsions = [(offset1+alpha[1].chi2).tolist()]

    # Create metadata dictionary
    metadata = {"dihedrals": tdrive.torsions,
                "spacing"  : [tdrive.spacing],
                "torsion_grid_ids": [[idx] for idx in range(-180+tdrive.spacing,180+tdrive.spacing,tdrive.spacing)]}
    with open("metadata.json","w") as meta:
        metastring = json.dumps(metadata, indent=4)
        print(metastring,file=meta)
    
    # Update torsion list
    _list = [phi, psi]
    for torsion in tdrive.torsions: _list.append(torsion)
    tdrive.torsions = _list

    # Get coordinates and number of atoms
    names  = [name for residue in alpha for name in residue.names]
    natoms = len(names)

    # Generate GRO files for all conformers in torsion scan
    dotopol = True
    confs = ["conf0", "conf1"]
    phi_val = {}
    psi_val = {}
    _phi = [x+1 for x in phi]
    _psi = [x+1 for x in psi]
    _chi = [x+1 for x in tdrive.torsions[-1]]
    for conf in confs:
        with open(f"scan-{conf}{suffix}.xyz", "r") as xyzfile:
            for snapshot in range(len(metadata["torsion_grid_ids"])):
                xyzfile.readline()
                title = xyzfile.readline().strip()
                coords = np.zeros((natoms,3))
                for iatom in range(natoms):
                    xyz = xyzfile.readline().split()[1:]
                    coords[iatom] = np.array([ float(x) for x in xyz ])/10.0
                phi_val[conf] = get_torsion(*coords[phi])
                psi_val[conf] = get_torsion(*coords[psi])
                with open(f"scan-{conf}-{snapshot:02d}.g96","w") as grofile:
                    print("TITLE", file=grofile)
                    print(title, file=grofile)
                    print("END", file=grofile)
                    print("POSITION", file=grofile)
                    iatom = 0
                    iresidue = 0
                    for residue in alpha:
                        for atom in residue.names:
                            print(f"{iresidue+1:5d} {residue.name:3s}   {str(atom):<4s}   {iatom+1:>5d}{coords[iatom,0]:15.9f}{coords[iatom,1]:15.9f}{coords[iatom,2]:15.9f}", file=grofile)
                            iatom += 1
                        iresidue += 1
                if dotopol:
                    dotopol = False
                    p = subprocess.Popen(["gmx","pdb2gmx","-f",f"scan-{conf}-{snapshot:02d}.g96",
                                    "-o","dummy.gro","-water","none","-p","topol.top"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    p.stdin.write("1".encode('ascii'))
                    p.stdin.close()
                    p.wait()
                    os.remove("dummy.gro")
                    with open("topol.top","r") as topfile:
                        topcontents = topfile.readlines()
                    search_forcefield = True
                    search_system = True
                    with open("topol.top","w") as topfile:
                        for line in topcontents:
                            if search_forcefield and "forcefield.itp" in line: 
                                topfile.write("#define _FF_AMBER\n#define _FF_AMBER99SB\n")
                                topfile.write("""#include "ffnonbonded.itp"\n#include "ffbonded.itp"\n\n""")
                                search_forcefield = False
                                continue
                            elif search_system and "[ system ]" in line :
                                topfile.write("#ifdef DIHRES\n")
                                topfile.write("[ dihedral_restraints ]\n")
                                topfile.write("{:5d}{:5d}{:5d}{:5d}   1   @phival@   0   14000\n".format(*_phi))
                                topfile.write("{:5d}{:5d}{:5d}{:5d}   1   @psival@   0   14000\n".format(*_psi))
                                topfile.write("{:5d}{:5d}{:5d}{:5d}   1   @chival@   0  100000\n\n".format(*_chi))
                                topfile.write("#endif\n\n")
                                topfile.write(line)
                                search_system = False
                                continue
                            else:
                                topfile.write(line)
                subprocess.run(["gmx","editconf","-f",f"scan-{conf}-{snapshot:02d}.g96",
                                "-o",f"scan-{conf}-{snapshot:02d}_box.g96",
                                "-c","-box","999.9"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                os.remove(f"scan-{conf}-{snapshot:02d}.g96")
        with open(f"{conf}_all.g96","w") as grofile:
            for snapshot in range(len(metadata["torsion_grid_ids"])):
                with open(f"scan-{conf}-{snapshot:02d}_box.g96","r") as infile:
                    grofile.write(infile.read())
                os.remove(f"scan-{conf}-{snapshot:02d}_box.g96")
    with open("posre.itp","w") as posre:
       print("[ position_restraints ]", file=posre)
       iatom = 0
       for residue in alpha:
           for atom in residue.elements:
               iatom += 1
               if atom == "H": continue
               force = 200
               print(f"{iatom:6d}     1 {force:6d} {force:6d} {force:6d}",file=posre)
    with open("torsionprofile.in","w") as infile:
        infile.write(forcebalance_input.format(name1=confs[0],name2=confs[1]))
    with open("shot.mdp","w") as infile:
        infile.write(shotmdp)
    with open("forcebalance.sbatch","w") as sfile:
        sfile.write(slurm.header)
        sfile.write("module purge\n")
        sfile.write("module load python\n")
        sfile.write("module use /tahoma/emsle60550/modulefiles\n")
        sfile.write("module load gromacs/2023 \n")
        sfile.write("module load gcc/10.2.0 \n\n")
        sfile.write(installptmpsi)
        sfile.write("mkdir -p fbdir/forcefield \n")
        sfile.write("cp ${SLURM_SUBMIT_DIR}/posre.itp fbdir/forcefield \n")
        sfile.write("cp ${SLURM_SUBMIT_DIR}/amber99sb.ff/ffbonded.itp fbdir/forcefield \n")
        sfile.write("cp ${SLURM_SUBMIT_DIR}/amber99sb.ff/ffnonbonded.itp fbdir/forcefield \n")
        sfile.write("cp ${SLURM_SUBMIT_DIR}/torsionprofile.in fbdir \n")
        sfile.write("sed -i '1i [ defaults ]' fbdir/forcefield/ffnonbonded.itp\n")
        sfile.write("sed -i '2i 1               2               yes             0.5     0.8333' fbdir/forcefield/ffnonbonded.itp\n")
        for conf in confs:
            sfile.write(f"mkdir -p fbdir/targets/{conf} \n")
            sfile.write(f"cp ${{SLURM_SUBMIT_DIR}}/{conf}_all.g96 fbdir/targets/{conf} \n")
            sfile.write(f"cp ${{SLURM_SUBMIT_DIR}}/topol.top fbdir/targets/{conf} \n")
            sfile.write(f"cp ${{SLURM_SUBMIT_DIR}}/metadata.json fbdir/targets/{conf} \n")
            sfile.write(f"cp ${{SLURM_SUBMIT_DIR}}/qdata-{conf}{suffix}.txt fbdir/targets/{conf}/qdata.txt \n")
            sfile.write(f"cp ${{SLURM_SUBMIT_DIR}}/shot.mdp fbdir/targets/{conf} \n")
            sfile.write(f"sed -i 's/@phival@/{phi_val[conf]:8.2f}/' fbdir/targets/{conf}/topol.top \n")
            sfile.write(f"sed -i 's/@psival@/{psi_val[conf]:8.2f}/' fbdir/targets/{conf}/topol.top \n")
        sfile.write("\n\n cd fbdir \n")
        sfile.write("ForceBalance -b torsionprofile.in \n")
        sfile.write("cp result/*/ffbonded.itp ${SLURM_SUBMIT_DIR}/ffbonded_new.itp\n") 
    return

def torsionscan(elements,coordinates,torsions,ligand=False,charge=0,mult=1):
    if not ligand:
        assert len(torsions) == 3, "For non-standard residues, 3 torsions must be given"
        phi, psi, scan = torsions
    for torsion in torsions:
        assert len(torsion) == 4, "each torsion should have 4 indices"
    for coord in coordinates:
        assert len(coord) == 3, "Could not find X,Y,Z triplet"

    natoms = len(coordinates)

    geometry = np.array(coordinates,dtype=float)
    phi_value = get_torsion(*geometry[phi])
    psi_value = get_torsion(*geometry[psi])

    zcoord = ""
    for torsion in torsions:
        zcoord +=  "torsion  {}  {}  {}  {}\n".format(torsion)

    coords = ""
    for e,c in zip(elements,geometry):
        coords += f"{e} {c[0]}  {c[1]}  {c[2]} \n"

    text = torsnw.format(memory=4000, name="", torsions=zcoord,
                         geometry=coords, charge=charge, mult=mult)
    with open("dihedrals.txt") as fh: fh.write("{}  {}  {}  {}".format(*scan))
    with open("extras.txt") as fh:
        fh.write("$set\n dihedral  {}  {}  {}  {}  {}\n".format(*phi,phi_value))
        fh.write(" dihedral  {}  {}  {}  {}  {}".format(*psi,psi_value))
    with open("torsiondrive.nw") as fh: fh.write(text)
    with open("torsiondrive.sbatch") as fh:
        fh.write(slurm_torsiondrive.format(memory=memory,mult=mult,charge=charge))
    return



def bonds_and_angles(protein,ligand=False):
    alpha = copy.deepcopy(protein)
    if not ligand:
        alpha.prepend(chain="A", residue="ACE", phi=-60)
        alpha.append(chain="A", residue="NME", psi=-45)

    names = []
    elements = []
    for residue in alpha.chains[0].residues:
        for i in range(len(residue.elements)):
            names.append(residue.names[i])
            elements.append(residue.elements[i])

    natoms,charge,mult,geometry = readoptim("alpha.log")
    _geometry = geometry[:-1]
    _geometry = _geometry.split("\n")
    geometry = np.zeros((natoms,3))
    for i,line in enumerate(_geometry):
        _line = line.split()
        geometry[i] = np.array(_line[1:]).astype(float)
    akab, aktheta, ablist, aalist = bondedks("alpha",elements,geometry)

    natoms,charge,mult,geometry = readoptim("beta.log")
    _geometry = geometry[:-1]
    _geometry = _geometry.split("\n")
    geometry = np.zeros((natoms,3))
    for i,line in enumerate(_geometry):
        _line = line.split()
        geometry[i] = np.array(_line[1:]).astype(float)    
    bkab, bktheta, bblist, balist = bondedks("beta",elements,geometry)

    if not np.array_equal(ablist,bblist):
        print(" Error: bonded lists are not equal")
    if not np.array_equal(aalist,balist):
        print(" Error: angle lists are not equal")

    kab = 0.5*(akab + bkab)
    ktheta = 0.5*(aktheta + bktheta)

    print(" Bond force constants [kJ/mol/nm^2]")
    print(" Atom1     Atom2       kAB")
    print(" ------------------------------")
    for i in range(len(ablist)):
        print(" {:5s}     {:5s}    {:10.1f}".format(names[ablist[i,0]],names[ablist[i,1]],kab[i]))

    print("")
    print(" Angle force constants [kJ/mol/rad^2]")
    print(" Atom1     Atom2     Atom3      kAB")
    print(" --------------------------------------")
    for i in range(len(aalist)):
        print(" {:5s}      {:5s}      {:5s}    {:8.3f}".format(names[aalist[i,0]],names[aalist[i,1]],names[aalist[i,2]],ktheta[i]))
    print("")
    return

def qmmm_optimize(filename, qmres=None, counter=False, center=False, orient=False, totcharge=0, **kwargs):
    
    # Process all keyword arguments
    slurm  = Slurm(**kwargs)
    nwchem = NWChem(**kwargs)

    # Remove extension
    name = filename[:-4]

    # Clean PDB file
    _pdb = f"{name}_clean.pdb"
    _pqr = f"{name}_clean.pqr"
    fh = open("pdb2pqr_clean.log","w")
    subprocess.run(["pdb2pqr30",
            "--ff","AMBER",
            "--ffout","AMBER",
            "--pdb-output",_pdb,
            "--clean",
            filename, _pqr], stdout=fh, stderr=subprocess.STDOUT)
    fh.close()

    # Convert AMBER to NWChem format
    with open("amber2nwchem.sh","w") as fh:
        fh.write(amber2nwchem)

    subprocess.run(["bash","amber2nwchem.sh",_pdb])

    modify = ""
    for res in qmres:
        modify += f" modify segment {res} quantum\n"

    counter_string = ""
    if counter:
        if totcharge > 0:
            counter_string = f"counter {int(round(totcharge))} Cl\n"
        elif totcharge < 0:
            counter_string = f"counter {int(round(totcharge))} Na\n"

    center_string = "center\n" if center else ""
    orient_string = "orient\n" if orient else ""

    # Write NWChem input file for geometry minimization
    with open("prepare.nw","w") as fh:
        fh.write(prepare.format(name=name, complex=_pdb, charge=nwchem.charge,
            aobasis=nwchem.aobasis, cdbasis=nwchem.cdbasis, nscf=nwchem.nscf,
            nopt=nwchem.nopt, bqzone=nwchem.bqzone, size=nwchem.boxsize,
            xcfun=nwchem.xcfun, grid=nwchem.grid, lshift=nwchem.lshift,
            modify=modify, counter=counter_string, memory=nwchem.memory,
            disp=nwchem.disp, center=center_string, orient=orient_string))


    # If needed, write fragment files
    with open(filename,"r") as fh:
        line = fh.readline()
        while line:
            if line[17:20] == "EAC":
                with open("EAC.frg","w") as frg:
                    frg.write(EAC_frg)
                break
            line = fh.readline()





    return
