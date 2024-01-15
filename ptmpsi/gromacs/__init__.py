from ptmpsi.gromacs.templates import ions,minim,heating,npt,md,SLURM_header,SLURM_tail
from ptmpsi.gromacs.utils import amber_to_gromacs_names
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

def generate_mdp(temp=300,posres=[1000.0,500.0,100.0,50.0,10.0,5.0,1.0,0.0],posresnpt=[1000.0,0.0],time=100.0,**kwargs):
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
        fh.write(heating.format(temp=temp))

    # Generate MDP files for NPT equilibration
    for force in posresnpt:
        name = "npt{}.mdp".format(int(force))
        nsteps = int(50000/len(posresnpt))
        if force > 0.0:
            restraint = "define          = -DP{}".format(int(force))
        else:
            restraint = ""
        with open(name,"w") as fh:
            fh.write(npt.format(temp=temp,restraint=restraint,nsteps=nsteps))

    # Generate MDP file for MD production
    nsteps = int(1000.0*time/0.002)
    with open("md.mdp","w") as fh:
        fh.write(md.format(temp=temp,nsteps=nsteps))
    return

def generate_slurm(infile,posres=[1000.0,500.0,100.0,50.0,10.0,5.0,1.0,0.0],posresnpt=[1000.0,0.0],results_dir="results",time="3:00:00",account="emsls60202",partition="normal",nodes=1,ntasks=1,nthreads=None,gpu=True,bt="dodecahedron",npos=None,nneg=None,conc=None,ff="amber99sb.ff",**kwargs):
    gpuline = "#SBATCH --gres=gpu:1" if gpu else ""
    ompline = "" if nthreads is None else "export OMP_NUM_THREADS={}".format(nthreads)
    addion = "-np {} -nn {}".format(npos,nneg) if conc is None else "-neutral -conc {}".format(conc)
    _partition = "analysis" if gpu else partition
    with open("slurm.sbatch","w") as fh:
        fh.write(SLURM_header.format(
            results=results_dir,
            time=time,
            account=account,
            partition=_partition,
            nodes=nodes,
            ntasks=ntasks,
            gpu=gpuline,
            omp=ompline,
            infile=infile,
            bt=bt,
            addion=addion,
            jname=infile,
            ff=ff))
        for force in posres:
            fh.write("echo Protein | gmx genrestr -f step4.gro -n index.ndx -o {}p.itp -fc {} {} {}\n".format(int(force),force,force,force))
            fh.write("""sed -i '/; Include water topology/i #ifdef P{force}\\n  #include "{force}p.itp"\\n#endif' topol.top\n""".format(force=int(force)))

        oldforce = -1
        for force in posres:
            if force == posres[0]:
                fh.write("gmx grompp -f minim{force}.mdp -c step4.gro -r step4.gro -p topol.top -n index.ndx -o minim{force}.tpr\n".format(force=int(force)))
            else:
                fh.write("gmx grompp -f minim{force}.mdp -c minim{oldforce}.gro -r minim{oldforce}.gro -p topol.top -n index.ndx -o minim{force}.tpr\n".format(force=int(force),oldforce=int(posres[oldforce])))
            fh.write("gmx_mpi mdrun -nb gpu -deffnm minim{force} \n".format(force=int(force)))
            oldforce += 1

        fh.write("gmx grompp -f heating.mdp -c minim{oldforce}.gro -r minim{oldforce}.gro -p topol.top -n index.ndx -o heating.tpr\n".format(oldforce=int(posres[-1])))
        fh.write("gmx_mpi mdrun -nb gpu -pme gpu -deffnm heating \n")

        oldforce = -1
        for force in posresnpt:
            if force == posresnpt[0]:
                fh.write("gmx grompp -f npt{force}.mdp -c heating.gro -r heating.gro -p topol.top -n index.ndx -o npt{force}.tpr\n".format(force=int(force)))
            else:
                fh.write("gmx grompp -f npt{force}.mdp -c npt{oldforce}.gro -r npt{oldforce}.gro -p topol.top -n index.ndx -o npt{force}.tpr\n".format(force=int(force),oldforce=int(posresnpt[oldforce])))
            fh.write("gmx_mpi mdrun -nb gpu -pme gpu -deffnm npt{force} \n".format(force=int(force)))
            oldforce += 1

        fh.write("gmx grompp -f md.mdp -c npt{force}.gro -p topol.top -n index.ndx -o md.tpr\n".format(force=int(posresnpt[-1])))
        fh.write("gmx_mpi mdrun -nb gpu -pme gpu -deffnm md\n")

        fh.write(SLURM_tail)
    return




def sltcap(protein,bt="dodecahedron",size=None,conc=0.154,charge=0,**kwargs):
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



def generate(protein,size=None,conc=None,charge=0,**kwargs):
    if size is None or conc is None:
        npos = None; nneg = None
    else:
        npos,nneg = sltcap(protein,size=size,charge=charge,**kwargs)

    if conc is None:
        _conc = 0.154 if npos is None else None
    else:
        _conc = conc if npos is None else None

    protein.write_pdb("protein.pdb")

    generate_mdp(**kwargs)
    generate_slurm("protein.pdb",conc=_conc,npos=npos,nneg=nneg,**kwargs)
    return
