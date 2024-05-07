from ptmpsi.gromacs.templates import ions,minim,heating,npt,md
from ptmpsi.gromacs.utils import amber_to_gromacs_names
from ptmpsi.slurm import Slurm
import numpy as np

_navogadro = 6.0221408E23
_cm2nm = 1.0E7
_volumes = {
"ALA":  88.3,
"CYS": 113.4,
"CYM": 113.4,
"CYX": 113.4,
"ASP": 117.1,
"GLU": 139.1,
"PHE": 192.4,
"GLY":  65.2,
"HIS": 159.2,
"HIE": 159.2,
"HIP": 159.2,
"HID": 159.2,
"ILE": 157.9,
"LYS": 163.7,
"LEU": 159.0,
"MET": 164.5,
"ASN": 125.4,
"ASH": 125.4,
"PRO": 121.8,
"GLN": 148.3,
"GLH": 148.3,
"ARG": 187.9,
"SER":  95.9,
"THR": 119.1,
"VAL": 134.6,
"TRP": 227.4,
"TYR": 198.1
}

def generate_mdp(temp=300,posres=[1000.0,500.0,100.0,50.0,10.0,5.0,1.0,0.0],**kwargs):

    lennvt = kwargs.pop("lennvt", 500.0)
    lennpt = kwargs.pop("lennpt", 500.0)
    lenmd  = kwargs.pop("lenmd",  100.0)
    timestep = kwargs.pop("timestep", 2.0)


    # Generate MDP file for ions addition
    with open("ions.mdp","w") as fh:
        fh.write(ions)

    # Generate MDP files for minimization
    for force in posres:
        name = "minim{}.mdp".format(int(force))
        if force > 0.0:
            restraint = "define          = -DP{}".format(int(force))
        else:
            restraint = ""
        with open(name,"w") as fh:
            fh.write(minim.format(restraint=restraint))

    # Generate MDP file NVT equilibration
    with open("heating.mdp","w") as fh:
        fh.write(heating.format(temp=temp, timestep=timestep/1000.0, nsteps=int(lennvt/timestep*1000)))

    # Generate MDP files for NPT equilibration
    with open("npt.mdp", "w") as fh:
        fh.write(npt.format(temp=temp, restraint="", timestep=timestep/1000.0, nsteps=int(lennpt/timestep*1000)))

    # Generate MDP file for MD production
    with open("md.mdp","w") as fh:
        fh.write(md.format(temp=temp, timestep=timestep/1000.0, nsteps=int(lenmd*1000000/timestep)))

    return

def generate_slurm(infile, posres=[1000.0,500.0,100.0,50.0,10.0,5.0,1.0,0.0],
                           results_dir="results",
                           bt="dodecahedron",
                           npos=None,
                           nneg=None,
                           conc=None,
                           d=1.0,
                           center=True,
                           **kwargs):

    c = "-c" if center else ""

    jobname   = kwargs.pop("jobname", f"{infile[:-4]}")
    slurm = Slurm("gromacs", jobname=jobname, **kwargs)
    gmx       = kwargs.pop("gmx", "gmx")
    container = kwargs.pop("container", "")
    gpu_id    = kwargs.pop("gpu_id", "")
    gpu_id    = f"-gpu_id {gpu_id}" if len(gpu_id) > 0 else ""
    mpirun    = kwargs.pop("mpirun", f"mpirun -np {slurm.ncpus}")

    if conc is None and (npos is None or nneg is None):
        raise KeyError("Specify total ion concentration or a number of positive and negative ions to add")

    addion = f"-np {npos} -nn {nneg}" if conc is None else "-neutral -conc {}".format(conc)

    with open(f"{infile[:-4]}_slurm.sbatch","w") as fh:
        fh.write(slurm.header)

        # PDB2GMX
        fh.write(f"""echo -e "1\\n1" | mpirun -np 1 {container} {gmx} pdb2gmx -f {infile} -o step1.gro -merge all \n""")

        # EDITCONF
        fh.write(f"""mpirun -np 1 {container} {gmx} editconf -f step1.gro -o step2.gro -bt {bt} -d {d} {c}\n""")

        # SOLVATE
        fh.write(f"""mpirun -np 1 {container} {gmx} solvate -cp step2.gro -cs spc216.gro -p topol.top -o step3.gro\n""")

        # IONS
        fh.write(f"""mpirun -np 1 {container} {gmx} grompp -f ions.mdp -c step3.gro -p topol.top -o ions.tpr\n""")
        fh.write(f"""echo SOL | mpirun -np 1 {container} {gmx} genion -s ions.tpr -o step4.gro -p topol.top -pname NA -nname CL {addion}\n""")

        # INDEX
        fh.write(f"""echo q | mpirun -np 1 {container} {gmx} make_ndx -f step4.gro -o index.ndx\n""")

        # INCLUDE POSRES DEFINITIONS
        for force in posres:
            fh.write(f"echo Protein | mpirun -np 1 {container} {gmx} genrestr -f step4.gro -n index.ndx -o {int(force)}p.itp -fc {force} {force} {force}\n")
            fh.write("""sed -i '/; Include water topology/i #ifdef P{force}\\n  #include "{force}p.itp"\\n#endif' topol.top\n""".format(force=int(force)))


        # RUN ALL MINIMIZATIONS
        oldforce = -1
        for force in posres:
            if force == posres[0]:
                fh.write(f"mpirun -np 1 {container} {gmx} grompp -f minim{int(force)}.mdp -c step4.gro -r step4.gro -p topol.top -n index.ndx -o minim{int(force)}.tpr\n")
            else:
                fh.write(f"mpirun -np 1 {container} {gmx} grompp -f minim{int(force)}.mdp -c minim{int(posres[oldforce])}.gro -r minim{int(posres[oldforce])}.gro -p topol.top -n index.ndx -o minim{int(force)}.tpr\n")
            if slurm.ngpus >= 1:
                fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -deffnm minim{int(force)} -v \n")
            else:
                fh.write(f"{mpirun} {container} {gmx} mdrun -deffnm minim{int(force)} -v \n")
            oldforce += 1

        # HEATING
        fh.write(f"mpirun -np 1 {container} {gmx} grompp -f heating.mdp -c minim{int(posres[-1])}.gro -r minim{int(posres[-1])}.gro -p topol.top -n index.ndx -o heating.tpr\n")
        if slurm.ngpus > 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm heating \n")
        elif slurm.ngpus == 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -bonded gpu -deffnm heating \n")
        else:
            fh.write(f"{mpirun} {container} {gmx} mdrun -deffnm heating \n")


        # NPT
        fh.write(f"mpirun -np 1 {container} {gmx} grompp -f npt.mdp -c heating.gro -r heating.gro -p topol.top -n index.ndx -o npt.tpr\n")
        if slurm.ngpus > 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm npt \n")
        elif slurm.ngpus == 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -bonded gpu -deffnm npt \n")
        else:
            fh.write(f"{mpirun} {container} {gmx} mdrun -deffnm npt \n")


        # PRODUCTION
        fh.write(f"mpirun -np 1 {container} {gmx} grompp -f md.mdp -c npt.gro -p topol.top -n index.ndx -o md.tpr\n")
        if slurm.ngpus > 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm md \n")
        elif slurm.ngpus == 1:
            fh.write(f"{mpirun} {container} {gmx} mdrun {gpu_id} -nb gpu -bonded gpu -deffnm md \n")
        else:
            fh.write(f"{mpirun} {container} {gmx} mdrun -deffnm md \n")

    return




def sltcap(protein, bt="dodecahedron", size=None, conc=0.154, charge=0, **kwargs):
    volume = 0
    for chain in protein.chains:
        for residue in chain.residues:
            volume += _volumes.get(residue.name,118.9)
    volume /= 1000.0

    total_volume = size**3
    if bt == "dodecahedron": total_volume *= 0.7071
    if bt == "octahedron": total_volume *= 0.7698

    water_volume = total_volume - volume
    nparticles = conc*water_volume/_cm2nm**3 * _navogadro/1000.0
    factor = np.arcsinh(-0.5*charge/nparticles)
    npositive = int(np.exp(factor)*nparticles)
    nnegative = int(np.exp(-factor)*nparticles)

    return npositive,nnegative



def generate(protein, filename=None, size=None, conc=None, charge=0, **kwargs):

    if size is None or conc is None:
        npos = None; nneg = None
    else:
        npos,nneg = sltcap(protein,size=size,conc=conc,charge=charge,**kwargs)

    if conc is None:
        _conc = 0.154 if npos is None else None
    else:
        _conc = conc if npos is None else None

    filename = "protein.pdb" if filename is None else filename
    protein.write_pdb(filename)

    generate_mdp(**kwargs)
    generate_slurm(filename, conc=_conc, npos=npos, nneg=nneg, **kwargs)
    return
