import requests
import copy
import numpy as np
import subprocess
from shutil import which
from ..exceptions import FeatureError
from ptmpsi.residues import resdict, Residue
from ptmpsi.math import find_clashes, find_clashes_residue, appendc, prependn
from ptmpsi.protein.mutate import point_mutation, post_translational_modification
from ptmpsi.protein.tools import get_residue, ptm_combination
from ptmpsi.io import digestpdb, writepdb, writexyz
from ptmpsi.docking import Dock, dock_ligand
from ptmpsi.gromacs.utils import amber_to_gromacs_names
from ptmpsi.gromacs import generate as generate_gromacs
from ptmpsi.gromacs.templates import qlambdas, vdwlambdas


class Chain:
    def __init__(self,name):
        self.name = name
        self.nresidues = 0
        self.natoms    = 0
        self.residues  = None



class Protein:
    def __init__(self,filename=None,pdbid=None,uniprotid=None,interactive=False,delwat=True,delhet=True):
        self.filename = filename
        self.pdbid = pdbid
        self.uniprotid = uniprotid
        self.pdbfile = None
        self.chains = None
        self.natoms = 0
        self.nresidues = 0
        self.nonstandard = None
        self.missing     = None
        self.nssbonds = 0
        self.ssbonds = None
        self.charge = None
        self.docking = Dock()
        
        # Download file from the PDB
        if self.pdbid is not None:
            print("\t Downloading file from the Protein Databank")
            response = requests.get("https://files.rcsb.org/download/"+self.pdbid.upper()+".pdb")
            response.raise_for_status()
            self.pdbfile = response.text.splitlines()
            with open(self.pdbid+".pdb","w") as fh:
                for line in self.pdbfile:
                    fh.write(line+"\n")
            del(response)
            
        # Download file from AlphaFold Database
        elif self.uniprotid is not None:
            print("\t Downloading file from the AlphaFold Protein Structure Database")
            response = requests.get("https://alphafold.ebi.ac.uk/files/AF-"+self.uniprotid.upper()+"-F1-model_v4.pdb")
            response.raise_for_status()
            self.pdbfile = response.text.splitlines()
            with open(self.uniprotid+".pdb","w") as fh:
                for line in self.pdbfile:
                    fh.write(line+"\n")
        
        # Read local file
        elif self.filename is not None:
            print("\t Reading local PDB file")
            with open(self.filename,'r') as fh:
                self.pdbfile = fh.readlines()

        # Process PDB file
        if self.pdbfile is not None:
            digestpdb(self,interactive,delwat,delhet)

        # Initialize a dummy chain
        else:
            self.chains = [Chain("A")]
            self.chains[0].residues = []
            self.chains[0].natoms = 0
            self.chains[0].nresidues = 0

        return


    def write_pdb(self,pdbfile):
        writepdb(self,pdbfile)

    def write_xyz(self,xyzfile):
        writexyz(self,xyzfile)

    def update(self):
        self.nresidues = 0
        self.natoms = 0
        self.nchains = len(self.chains)
        for chain in self.chains:
            natoms = 0
            for ires,residue in enumerate(chain.residues):
                residue.resid = ires+1
                natoms += residue.natoms
            chain.natoms = natoms
            chain.nresidues = len(chain.residues)
            self.nresidues += chain.nresidues
            self.natoms += natoms
        return


    def delwaters(self):
        for chain in self.chains:
            waters = []
            for ires,residue in enumerate(chain.residues):
                if residue.name in ['HOH','WAT']:
                    waters.append(ires)
            for water in reversed(waters):
                del chain.residues[water]
        return


    @staticmethod
    def addres(residues,residue,natoms,names,elements,coordinates,chain,backbone):
        _residue = Residue(residue,natoms)
        _residue.names = np.array(names)
        _residue.elements = np.array(elements)
        _residue.coordinates = np.array(coordinates)
        _residue.chain = chain
        _residue.backbone = backbone
        _residue.resid = len(residues) + 1
        residues.append(_residue)
        return


    @staticmethod
    def addchain(chains,chain,residues,natoms,resid,nmissing):
        _chain = Chain(chain)
        _chain.residues = copy.deepcopy(residues)
        _chain.natoms = natoms
        _chain.nresidues = resid
        chains.append(_chain)
        return [], 0 ,0 ,0


    def savessbond(self,ssbonds,residue,resid,nmissing,chain):
        if (self.nssbonds > 0) and (residue in ['CYS','CYX']):
            for issbond in range(self.nssbonds):
                if ((resid+nmissing == ssbonds[issbond][0]) and (chain == ssbonds[issbond][1])):
                        residue = 'CYX'
                        self.ssbonds[issbond][0] = resid
                        self.ssbonds[issbond][1] = chain
                elif ((resid+nmissing == ssbonds[issbond][2]) and (chain == ssbonds[issbond][3])):
                        residue = 'CYX'
                        self.ssbonds[issbond][2] = resid
                        self.ssbonds[issbond][3] = chain
        return residue


    def mutate(self,original,new):
        point_mutation(self,original,new)
        return


    def append(self,chain,residue,psi=None):
        try:
            _residue = copy.deepcopy(resdict[residue])
        except:
            raise MyDockingError("There is no residue with name '{}'".format(residue))


        if psi is None:
            psi = -40.0
        elif isinstance(psi,str):
            psi = phi.lower()
            if psi == "alpha":
                psi = -40
            elif psi == "beta":
                psi = 130.0
            else:
                raise MyDockingError()

        for _chain in self.chains:
            if _chain.name == chain:
                newcoords = appendc(_chain,_residue,psi)
                natoms = len(_residue.coordinates)
                newres = Residue(residue,natoms)
                newres.names = _residue.elements[:,1]
                newres.coordinates = newcoords
                newres.elements = _residue.elements[:,0]
                newres.chain = chain
                newres.backbone = _residue.backbone
                newres.resid = len(_chain.residues) + 1
                _chain.residues.append(newres)
                break
        self.update()
        find_clashes_residue(newres,[self])
        return


    def prepend(self,chain,residue,phi=None):
        try:
            _residue = copy.deepcopy(resdict[residue])
        except:
            raise MyDockingError("There is no residue with  name '{}'".format(residue))

        if phi is None:
            phi = -60.0
        elif isinstance(phi,str):
            phi = phi.lower()
            if phi == "alpha":
                phi = -60
            elif phi == "beta":
                phi = -140.0
            else:
                raise MyDockingError()

        for _chain in self.chains:
            if _chain.name == chain:
                newcoords = prependn(_chain,_residue,phi)
                natoms = len(_residue.coordinates)
                newres = Residue(residue,natoms)
                newres.names = _residue.elements[:,1]
                newres.coordinates = newcoords
                newres.elements = _residue.elements[:,0]
                newres.chain = chain
                newres.backbone = _residue.backbone
                newres.resid = len(_chain.residues) + 1
                newres.chi1 = _residue.chi1
                newres.chi2 = _residue.chi2
                newres.cattach = _residue.cattach
                newres.nattach = _residue.nattach
                _chain.residues.insert(0,newres)
                break
        self.update()
        find_clashes_residue(newres,[self])
        return


    def delchain(self,chain):
        newchains = []
        for i, _chain in enumerate(self.chains):
            if _chain.name == chain: continue
            newchains.append(_chain)
        self.chains = newchains
        self.update()
        return


    def findresidue(self,name):
        for _chain in self.chains:
            for residue in _chain.residues:
                if residue.name == name:
                    print("{}:{}{}".format(_chain.name,name,residue.resid))
        return


    def modify(self,original,modification):
        post_translational_modification(self,original,modification)
        return

    def get_ptm_combinations(self, ptms, exclude=None, ntuple=-1):
        return ptm_combination(self, ptms, ntuple, exclude)

    def gen_ptm_files(self, combinations, path=None, prefix=None, ff='amber99sb', **kwargs):
        import os
        import uuid
        import time
        import shutil
        cwd = os.getcwd()
        uid = str(uuid.uuid4())

        lenlambda = kwargs.pop("lenlambda", 2.0)
        timestep  = kwargs.get("timestep",  2.0)
        nsteps    = int(lenlambda*1000000/timestep)
        temp      = kwargs.get("temp", 300.0)

        ff = ff.lower()
        if ff not in ["amber99sb", "amber99zn", "amber14sb"]:
            KeyError(f"Forcefield '{ff}' not available")

        parent_dir = path if path is not None else cwd
        path = os.path.join(parent_dir, uid)
        os.mkdir(path)

        ff_prefix = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../gromacs/forcefield")
        shutil.copytree(f"{ff_prefix}/{ff}.ff", f"{path}/{ff}.ff")
        shutil.copy(f"{ff_prefix}/residuetypes.dat", path)
        shutil.copy(f"{ff_prefix}/specbond.dat", path)

        prefix = "" if prefix is None else f"{prefix}_"

        self.protonate(pdb=f"{path}/{prefix}protonated.pdb", pqr=f"{path}/{prefix}protonated.pqr")


        submit = open(f"{path}/submit.sh", "w")
        submit.write("#!/bin/bash\n")
        with open(f"{path}/README", "w") as fh:
            today = time.ctime()
            fh.write(f"""Files generated on {today}\n""")
            for i in range(len(combinations)):
                ipath = os.path.join(path, f"{i+1}tuples")
                os.mkdir(ipath)
                for j,combi in enumerate(combinations[i]):
                    jpath = os.path.join(ipath,f"{j:04d}")
                    os.mkdir(jpath)
                    _protein = Protein(filename=f"{path}/{prefix}protonated.pdb")
                    string = ""
                    for _ptm in combi:
                        _protein.modify(_ptm[0], _ptm[1])
                        string += f" {_ptm[0]} -> {_ptm[1]}"
                    os.chdir(jpath)
                    os.symlink(os.path.relpath(f"{path}/{ff}.ff", "./"), f"{ff}.ff")
                    os.symlink(os.path.relpath(f"{path}/residuetypes.dat", "./"), f"residuetypes.dat")
                    os.symlink(os.path.relpath(f"{path}/specbond.dat", "./"), f"specbond.dat")

                    # Thermodynamic integration, symlinks will be broken until file generation
                    os.mkdir("./dualti")
                    kpath = os.path.join(jpath, "dualti")
                    os.chdir(kpath)
                    os.symlink(os.path.relpath(f"{jpath}/topol.top", "./"), f"topol.top")
                    for k in range(13):
                        os.mkdir(f"{kpath}/lam-{k:02d}")
                        qpath = os.path.join(kpath, f"lam-{k:02d}/01-q")
                        os.mkdir(qpath)
                        os.chdir(qpath)
                        os.symlink(os.path.relpath(f"{path}/{ff}.ff", "./"), f"{ff}.ff")
                        os.symlink(os.path.relpath(f"{path}/residuetypes.dat", "./"), f"residuetypes.dat")
                        os.symlink(os.path.relpath(f"{jpath}/index.ndx", "./"), f"index.ndx")
                        os.symlink(os.path.relpath(f"{jpath}/md.gro", "./"), f"md.gro")
                        os.symlink(os.path.relpath(f"{jpath}/md.cpt", "./"), f"md.cpt")
                        os.symlink(os.path.relpath(f"{kpath}/TItop.top", "./"), "topol.top")
                        with open("grompp.mdp", "w") as grompp:
                            grompp.write(qlambdas.format(lambda_state=k, nsteps=nsteps, timestep=timestep, temp=temp))
                        os.chdir(jpath)
                        qpath = os.path.join(kpath, f"lam-{k:02d}/02-vdw")
                        os.mkdir(qpath)
                        os.chdir(qpath)
                        os.symlink(os.path.relpath(f"{path}/{ff}.ff", "./"), f"{ff}.ff")
                        os.symlink(os.path.relpath(f"{path}/residuetypes.dat", "./"), f"residuetypes.dat")
                        os.symlink(os.path.relpath(f"{jpath}/index.ndx", "./"), f"index.ndx")
                        os.symlink(os.path.relpath(f"{jpath}/md.gro", "./"), f"md.gro")
                        os.symlink(os.path.relpath(f"{kpath}/TItop.top", "./"), "topol.top")
                        with open("grompp.mdp", "w") as grompp:
                            grompp.write(vdwlambdas.format(lambda_state=k, nsteps=nsteps, timestep=timestep, temp=temp))
                        os.chdir(kpath)
                    os.chdir(jpath)

                    # Generate SLURM script
                    amber_to_gromacs_names(_protein)
                    generate_gromacs(_protein, filename=f"{prefix}{j:04d}.pdb", **kwargs)
                    fh.write(f"{i+1}tuples/{j:04d}/{prefix}{j:04d}.pdb: {string}\n")

                    # Update submission script
                    submit.write(f"cd {os.path.relpath(jpath, path)} \n")
                    submit.write(f"jobid=$(sbatch {prefix}{j:04d}_slurm.sbatch | sed 's/Submitted batch job //') \n")
                    submit.write(f"cd dualti\n")
                    for k in range(13):
                        submit.write(f"cd lam-{k:02d}\n")
                        submit.write(f"sbatch --dependency=afterok:$jobid {prefix}{j:04d}_lam{k:02d}_slurm.sbatch\n")
                        submit.write(f"cd ../ \n")
                        submit.write(f"sleep 1s \n")
                    submit.write(f"cd ../ \n")
                    submit.write(f"cd {os.path.relpath(path, jpath)} \n")
                    submit.write(f"sleep 1s \n")
                    submit.write(f"\n\n")
        submit.close()
        os.chdir(cwd)

        return uid

    def dock(self, ligand=None, receptor=None, boxcenter=None, boxsize=10, output=None, flexible=None, engine=None, exhaustiveness=None):
        if ligand is None:
            raise MyDockingError("No ligand was specified")
        if receptor is None:
            raise MyDockingError("No receptor was specified")

        if engine is None:
            self.docking.engine = "vina"
        else:
            if engine not in ["vina", "adgpu"]:
                raise MyDockingError("Only vina or AutoDockGPU can be used as Docking engines")
            self.docking.engine = engine

        if isinstance(exhaustiveness,int):
            self.docking.exhaustiveness = exhaustiveness
        else:
            raise MyDockingError("Exhaustiveness should be an integer value")

        if output is None:
            if self.docking.engine == 'vina':
                output = f"{receptor[:-4]}_{ligand[:-4]}_{self.docking.engine}.pdbqt"
            elif self.docking.engine == 'adgpu':
                output = f"{receptor[:-4]}_{ligand[:-4]}_{self.docking.engine}.dlg"


        if boxcenter is None:
            xcenter = 0; ycenter = 0; zcenter = 0
            for chain in self.chains:
                for residue in chain.residues:
                    for xyz in residue.coordinates:
                        xcenter += xyz[0]
                        ycenter += xyz[1]
                        zcenter += xyz[2]
            xcenter /= self.natoms
            ycenter /= self.natoms
            zcenter /= self.natoms
        elif isinstance(boxcenter,str):
            residue = get_residue(self, boxcenter)
            xcenter, ycenter, zcenter = np.mean(residue.coordinates,axis=0)
        else:
            [xcenter,ycenter,zcenter] = boxcenter

        self.docking.ligand = ligand
        self.docking.receptor = receptor
        self.docking.flexible = flexible
        self.docking.output = output
        self.docking.boxsize = boxsize
        self.docking.boxcenter = [xcenter, ycenter, zcenter]

        dock_ligand(self.docking)
        return

    def protonate(self,pdbin=None,pdb=None,pqr=None,ph=7):
        # Check if pdb2pqr30 is in the path
        if which("pdb2pqr30") is None:
            raise MyDockingError("Cannot find PDB2PQR")
        if pdbin is None:
            pdbin = ".tmp.pdb"
            self.write_pdb(pdbin)
        if self.pdbid is not None:
            _pdb = f"{self.pdbid}_H.pdb" if pdb is None else pdb
            _pqr = f"{self.pdbid}_H.pqr" if pqr is None else pqr
        elif self.uniprotid is not None:
            _pdb = f"{self.uniprotid}_H.pdb" if pdb is None else pdb
            _pqr = f"{self.uniprotid}_H.pqr" if pqr is None else pqr
        elif self.filename is not None:
            _pdb = f"{self.filename[:-4]}_H.pdb" if pdb is None else pdb
            _pqr = f"{self.filename[:-4]}_H.pqr" if pqr is None else pqr
        else:
            if (pdb is None) or (pqr is None):
                raise MyDockingError("Provide a PDB and PQR filename for protonoation output")
        fh = open("pdb2pqr.log","w")
        subprocess.run(["pdb2pqr30",
            "--ff","AMBER",
            "--ffout","AMBER",
            "--pdb-output",_pdb,
            "--titration-state-method","propka",
            "--with-ph",str(ph),
            "-o",str(ph),
            "--protonate-all",
            pdbin, _pqr], stdout=fh, stderr=subprocess.STDOUT)
        fh.close()
        self.charge = 0
        with open(_pqr,"r") as fh:
            lines = fh.readlines()
        for line in lines:
            try:
                self.charge += float(line.split()[8])
            except:
                pass

        with open(_pdb,"r") as fh:
            self.pdbfile = fh.readlines()
        digestpdb(self,interactive=False,delwat=False,delhet=False)

        return
