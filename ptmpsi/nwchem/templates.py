hfresp = """
unset "dft:cd*"
unset "basis:cd*"

basis "ao basis" spherical
 * library 6-31g*
end

dft
 xc hfexch 1.0
 vectors input atomic
end

task dft
"""

respnw = """start
title "RESP {name}"
set driver:linopt 0
set int:cando_txs f
set int:cando_nw f

memory total {memory} mb noverify


geometry nocenter {noautoz}
{geometry}
{zcoord}
end

{gcons}

charge {charge}

xtb
 acc 0.1
end

driver
 maxiter 100
end

task xtb optimize ignore

basis "ao basis" spherical
 * library {aobasis}
end

basis "cd basis" spherical bse
 * library {cdbasis}
end


dft
  noio
  {disp}
  xc {xcfun}
  grid {grid} nodisk
  mult {mult}
  maxiter {nscf}
  convergence lshift {lshift} energy 1d-7
  noprint "final vectors analysis"
end

driver
 maxiter {nopt}
end

task dft optimize ignore

{dohfresp}

esp
 {constraint}
 restrain hfree
end

task esp
"""


fit = """#!/usr/bin/env python3
import numpy as np
import copy
import math

def norm2(vec):
    return math.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 )

# Parameters
natoms = {natoms}
names  = {names}
nconf  = {nconf}
files  = {files}
ncons  = {ncons}
charge = {charge}

# Constrain hydrogens bonded to the same atom to have same charge
hcons = {hcons}
hncons = len(hcons)

# Constrain NME, ACE, and amide bond charges (for AMBER99)
cons = [
{cons}
]

grids = []
geometries = np.zeros((nconf,natoms,3))
npoints = np.zeros(nconf, dtype=int)
A = np.zeros((natoms+ncons+hncons+1,natoms+ncons+hncons+1))
B = np.zeros(natoms+ncons+hncons+1)

for i,file in enumerate(files):
    filename = file + ".xyz"
    with open(filename,"r") as fh:
        fh.readline(); fh.readline()
        for atom in range(natoms):
            line = fh.readline().split()
            geometries[i,atom] = [float(x)/0.529177 for x in line[1:4]]
    filename = file + ".grid"
    with open(filename,"r") as fh:
        npoints[i] = int(fh.readline().split()[0])
        grids.append(np.zeros((npoints[i],4)))
        for point in range(npoints[i]):
            line = fh.readline().split()
            grids[i][point] = [float(x) for x in line]

for iconf in range(nconf):
    dists = np.zeros((natoms,npoints[iconf]))
    for iatom in range(natoms):
        for k in range(npoints[iconf]):
            dists[iatom,k] = 1.0/norm2(geometries[iconf,iatom]-grids[iconf][k,:3])
        B[iatom] += np.dot(grids[iconf][:,3],dists[iatom])
    for iatom in range(natoms):
        for jatom in range(iatom,natoms):
            A[iatom,jatom] += np.dot(dists[iatom],dists[jatom])

# Symmetrize matrix
for iatom in range(natoms):
    for jatom in range(iatom,natoms):
        A[jatom,iatom] = copy.copy(A[iatom,jatom])

# Total charge constraint
A[:natoms,natoms] = 1.0
A[natoms,:natoms] = 1.0
B[natoms] = charge

# NME, ACE, and amide bond constraints
for icons in range(ncons):
    A[cons[icons][0]-1,natoms+icons+1] = 1.0
    A[natoms+icons+1,cons[icons][0]-1] = 1.0
    B[natoms+icons+1] = cons[icons][1]

# Hydrogen bond constraints
for icons in range(hncons):
    A[hcons[icons][0],natoms+ncons+icons+1] = 1.0
    A[hcons[icons][1],natoms+ncons+icons+1] = -1.0
    A[natoms+ncons+icons+1,hcons[icons][0]] = 1.0
    A[natoms+ncons+icons+1,hcons[icons][1]] = -1.0

# Start from solution without restraints
qold, _, _, _ = np.linalg.lstsq(A,B)

# Hyperbolic restraints are non-linear. Do 50 iterations at most
for iter in range(50):
    Acur = copy.deepcopy(A)

    # Add restraint contribution to matrix
    for i in range(natoms):
        
        # Hydrogens are free of restraints
        if names[i][0] == "H": continue

        # Skip charges already constrained
        skip = False
        for j in range(ncons):
            if i == cons[j][0]-1: 
                skip = True
                break
        if skip: continue
        Acur[i,i] += 0.005 / np.sqrt(qold[i]**2 + 0.01)

    # Solve linear equation system
    q, _, _, _ = np.linalg.lstsq(Acur,B)

    # Check convergence
    delta = np.amax(np.abs(q-qold))
    print("iter {{}}, delta: {{}}".format(iter,delta))
    if delta < 0.000001: break

    # Copy current solution to qold
    qold = copy.deepcopy(q)

# Print charges to STDOUT
print("")
print(" RESP charges")
{printing}
print("")
"""

slurm_header = {}
slurm_header["Tahoma"] = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --time={time}
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --job-name={jname}
#SBATCH --error={jname}-%j.err
#SBATCH --output={jname}-%j.out
#SBATCH --partition={partition}

cleanup()
{{
cp *.xyz $SLURM_SUBMIT_DIR || :
cp *.log $SLURM_SUBMIT_DIR || :
cp *.txt $SLURM_SUBMIT_DIR || :
cp *.json $SLURM_SUBMIT_DIR || :
cp *.grid $SLURM_SUBMIT_DIR || :
cp *.qrs $SLURM_SUBMIT_DIR || :
cp *.pdb $SLURM_SUBMIT_DIR || :
cp *.rst $SLURM_SUBMIT_DIR || :
cp *.top $SLURM_SUBMIT_DIR || :
cp *.trj $SLURM_SUBMIT_DIR || :
cp *.out $SLURM_SUBMIT_DIR || :
}}

trap cleanup SIGINT SIGTERM SIGKILL SIGSEGV SIGCONT
source /etc/profile.d/modules.sh
module purge
module load python
module load gcc/9.3.0
module load openmpi

#export NWCHEM_BASIS_LIBRARY=/cluster/apps/nwchem/nwchem/src/basis/libraries/
#export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export NWC_RANKS_PER_DEVICE=0
export ARMCI_OPENIB_DEVICE=mlx5_0
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export https_proxy="http://proxy.emsl.pnl.gov:3128"
export http_proxy="http://proxy.emsl.pnl.gov:3128"
export NWBIN=/big_scratch/nwchems_`id -u`.img
#export NWCHEM_IMAGE="ghcr.io/edoapra/nwchem-singularity/nwchem-720.ompi41x:latest"
export NWCHEM_IMAGE="ghcr.io/edoapra/nwchem-singularity/nwchem-dev.mpi3.ompi5x:latest"

srun -N $SLURM_NNODES -n $SLURM_NNODES apptainer pull -F --name $NWBIN --disable-cache oras://$NWCHEM_IMAGE
export APPTAINERENV_SCRATCH_DIR={scratch}
export APPTAINERENV_OMP_NUM_THREADS=${OMP_NUM_THREADS}
#export APPTAINERENV_NWCHEM_BASIS_LIBRARY=$NWCHEM_BASIS_LIBRARY


cd {scratch}
"""

slurm_header["Frontier"] = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --time={time}
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --job-name={jname}
#SBATCH --error={jname}-%j.err
#SBATCH --output={jname}-%j.out
#SBATCH --partition={partition}

cleanup()
{{
cp *.xyz $SLURM_SUBMIT_DIR || :
cp *.log $SLURM_SUBMIT_DIR || :
cp *.txt $SLURM_SUBMIT_DIR || :
cp *.json $SLURM_SUBMIT_DIR || :
cp *.grid $SLURM_SUBMIT_DIR || :
cp *.qrs $SLURM_SUBMIT_DIR || :
cp *.pdb $SLURM_SUBMIT_DIR || :
cp *.rst $SLURM_SUBMIT_DIR || :
cp *.top $SLURM_SUBMIT_DIR || :
cp *.trj $SLURM_SUBMIT_DIR || :
cp *.out $SLURM_SUBMIT_DIR || :
}}

export SCRATCH="/lustre/orion/{account}/scratch/${{USER}}"

trap cleanup SIGINT SIGTERM SIGKILL SIGSEGV SIGCONT
module load python
module load rocm/5.7.1
module load cray-mpich-abi

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=hybrid
export APPTAINERENV_LD_LIBRARY_PATH="${{CRAY_MPICH_ROOTDIR}}/lib-abi-mpich:${{CRAY_MPICH_ROOTDIR}}/gtl/lib:${{CRAY_LD_LIBRARY_PATH}}:${{PE_PERFTOOLS_MPICH_LIBDIR}}:/opt/cray/lib64:/opt/cray/pe/lib64:/opt/cray/xpmem/2.8.4-1.0_7.3__ga37cbd9.shasta/lib64:${{ROCM_PATH}}/lib:${{ROCM_PATH}}/lib64:/opt/amdgpu/lib64:/opt/cray/pe/lib64/cce/:/opt/cray/pe/gcc-libs/:${{LD_LIBRARY_PATH}}"
export APPTAINER_CONTAINLIBS="/usr/lib64/libcxi.so.1,/usr/lib64/libjson-c.so.3,/lib64/libtinfo.so.6,/usr/lib64/libnl-3.so.200,/usr/lib64/libgfortran.so.5,/usr/lib64/libjansson.so.4"
export APPTAINERENV_LD_PRELOAD=$CRAY_MPICH_ROOTDIR/gtl/lib/libmpi_gtl_hsa.so.0:
MYFS=$(findmnt -r -T . | tail -1 |cut -d ' ' -f 1)
export BINDS=/usr/share/libdrm,/var/spool/slurmd,/opt/cray,${{MYFS}},${{ROCM_PATH}},/opt/amdgpu

export MKL_NUM_THREADS=1
export http_proxy="http://proxy.ccs.ornl.gov:3128/"
export https_proxy="http://proxy.ccs.ornl.gov:3128/"
export NWBIN=/${{SCRATCH}}/nwchems_`id -u`.img
export NWCHEM_IMAGE="ghcr.io/edoapra/nwchem-singularity/nwchem-dev.mpich3.4.2:latest"

FC=gfortran ARMCI_NETWORK=MPI-PR MPICH=3.4.2 srun -N $SLURM_NNODES -n $SLURM_NNODES apptainer pull -F --name $NWBIN --disable-cache oras://$NWCHEM_IMAGE
export APPTAINERENV_SCRATCH_DIR=/lustre/orion/bip258/scratch/$USER
export APPTAINERENV_OMP_NUM_THREADS=1


cd {scratch}
"""

slurm_header["Polaris"] = """#!/bin/bash -l
#PBS -l select={nnodes}:system=polaris
#PBS -l walltime={time}
#PBS -l filesystems=home:grand
#PBS -q {partition}
#PBS -N {jname}
#PBS -A {account}
#PBS -l place=scatter

cleanup()
{{
cp *.xyz $PBS_O_WORKDIR || :
cp *.log $PBS_O_WORKDIR || :
cp *.txt $PBS_O_WORKDIR || :
cp *.json $PBS_O_WORKDIR || :
cp *.grid $PBS_O_WORKDIR || :
cp *.qrs $PBS_O_WORKDIR || :
cp *.pdb $PBS_O_WORKDIR || :
cp *.rst $PBS_O_WORKDIR || :
cp *.top $PBS_O_WORKDIR || :
cp *.trj $PBS_O_WORKDIR || :
cp *.out $PBS_O_WORKDIR || :
}}

trap cleanup SIGINT SIGTERM SIGKILL SIGSEGV SIGCONT

export SCRATCH="/grand/TwinHostPath/${{USER}}/scratch"
mkdir -p $SCRATCH

# Load the appropriate environment
source "/lus/grand/projects/TwinHostPath/software/nwchem/share/export.sh"
module list
module unload xalt
echo $LD_LIBRARY_PATH

export MPICH_GPU_SUPPORT_ENABLED=0
export COMEX_MAX_NB_OUTSTANDING=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=hybrid
export COMEX_EAGER_THRESHOLD=16384
export FI_CXI_RDZV_THRESHOLD=16384
export FI_CXI_OFLOW_BUF_COUNT=6

export NWBIN=/grand/projects/TwinHostPath/dmejiar/nwchem/bin/LINUX64/nwchem

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS_PER_NODE={ntasks}
NDEPTH=1
NTOTRANKS=$(($NNODES * $NRANKS_PER_NODE))
NTHREADS=1

export OMP_NUM_THREADS=$NTHREADS

echo "NUM_OF_NODES= ${{NNODES}} TOTAL_NUM_RANKS= ${{NTOTRANKS}} RANKS_PER_NODE= ${{NRANKS_PER_NODE}} THREADS_PER_RANK= ${{NTHREADS}}"


cd ${{SCRATCH}}
"""

shotmdp = """integrator     = md
dt              = 0.001
nsteps          = 0             ; Maximum number of (minimization) steps to perform
nstxout         = 0
nstfout         = 1             ; Forces (important)
nstenergy       = 1             ; Write energies to disk every nstenergy steps
nstxtcout       = 0             ; Write coordinates to disk every nstxtcout steps
xtc_grps        = System               ; Which coordinate group(s) to write to disk
energygrps      = System               ; Which energy group(s) to write to disk
constraints     = none
nstlist         = 1
rlist           = 333.3
vdwtype         = cut-off
coulombtype     = cut-off
rcoulomb        = 333.3
rvdw            = 333.3
pbc             = xyz
"""

slurm_tdrive = """
# Create a Virtual Environment
if [ -d "venv" ]; then
  echo "Virtual environment already exists"
else
  python -m venv venv
fi
source venv/bin/activate

# Set proxy server
export https_proxy=http://proxy.emsl.pnl.gov:3128
export http_proxy=http://proxy.emsl.pnl.gov:3128

python -m pip install --upgrade pip
git clone git@code.emsl.pnl.gov:cheung_group/ptmpsi.git
cd ptmpsi
python -m pip install .
cd ..

export NWCHEM_COMMAND="{runsingularity_prefix_tahoma}"
"""

slurm_torsiondrive = """
cp ${{SLURM_SUBMIT_DIR}}/dihedrals.txt .
cp ${{SLURM_SUBMIT_DIR}}/extras.txt .
cp ${{SLURM_SUBMIT_DIR}}/nwchem.nw .

torsiondrive-launch torsiondrive.nw dihedrals.txt -c extras.txt -g 15 -e nwchem --native_opt -v

nlines=$(( $(head -1 scan.txt) + 2 ))
split -d -a2 -l$nlines --additional-suffix=.xyz scan.xyz scan-

for file in `ls scan-*.xyz`; do
name=${{file%.xyz}}.nw
cat <<EOF >$name 
start
memory {memory} mb
set int:cando_txs f
set int:cando_nw f

geometry
load "$file"
end

charge {charge}

basis "ao basis" spherical
 * library def2-tzvp
end
basis "cd basis" spherical bse
 * library def2-universal-jfit
end

dft
  mult {mult}
  grid nodisk fine
  xc r2scan
  disp vdw 4
end

task dft gradient
EOF

$NWCHEM_COMMAND $name > ${{name%.nw}}.log

done

cleanup
"""

slurm_resp = """
cp ${SLURM_SUBMIT_DIR}/alpha.nw .
cp ${SLURM_SUBMIT_DIR}/beta.nw .
cp ${SLURM_SUBMIT_DIR}/fit.py .

echo "Running alpha-helix conformer"
{runsingularity_prefix_tahoma} alpha.nw > alpha.log

echo "Running beta-strand conformer"
{runsingularity_prefix_tahoma} beta.nw > beta.log

# Create a Virtual Environment
if [ -d "venv" ]; then
  echo "Virtual environment already exists"
else
  python -m venv venv
fi
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy

echo "\\n Starting RESP fitting"
python fit.py

cp alpha.log $SLURM_SUBMIT_DIR
cp beta.log $SLURM_SUBMIT_DIR
"""

slurm_hess = """
cp ${SLURM_SUBMIT_DIR}/alpha_hess.nw .
cp ${SLURM_SUBMIT_DIR}/beta_hess.nw .

echo "\\n Running alpha-helix hessian"
{runsingularity_prefix_tahoma} alpha_hess.nw > alpha_hess.log

echo "\\n Running beta-sheet hessian"
{runsingularity_prefix_tahoma} beta_hess.nw > beta_hess.log

cp alpha_hess.log $SLURM_SUBMIT_DIR
cp beta_hess.log $SLURM_SUBMIT_DIR
"""


hessnw = """start
memory total {memory} mb noverify
title "{name} hessian calculation"

set int:cando_txs f
set int:cando_nw f


geometry nocenter
{geometry}
{zcoord}
end

basis "ao basis" spherical
 * library {aobasis}
end

basis "cd basis" spherical bse
 * library {cdbasis}
end

charge {charge}

dft
 noio
 {disp}
 mult {mult}
 xc {xcfun}
 grid {grid} nodisk
 maxiter {nscf}
 convergence energy 1d-7
 noprint "final vectors analysis"
end

freq
 fd_delta {delta}
end

task dft frequencies
"""

torsnw = """start
memory total {memory} mb
title "{name} torsion scan"

set int:cando_txs f
set int:cando_nw f
set driver:linopt 0

geometry
{geometry}
 zcoord
 end
end

basis "ao basis" spherical
 * library {aobasis}
end

basis "cd basis" spherical bse
 * library {cdbasis}
end

charge {charge}

dft
 noio
 grid {grid} nodisk
 mult {mult}
 {disp}
 xc {xcfun}
 maxiter {nscf}
 convergence lshift {lshift} energy 1d-7
 noprint "final vectors analysis"
end

task dft optimize ignore
"""

slurm_tdrive_run = """
nlines=$(( $( cat {tail}.xyz | wc -l)-2 ))
geom=`cat {tail}.xyz | tail -$nlines`
geom="${{geom//$'\\n'/\\\\n}}"
sed -i "s/@geometry@/$geom/" {tail}_tdrive.nw

torsiondrive-launch {tail}_tdrive.nw dihedrals{idx}.txt -c extras{idx}.txt -g {spacing} -e nwchem --native_opt -v

mv scan.xyz scan-{tail}.xyz
mv qdata.txt qdata-{tail}.txt

rm -rf opt_tmp

"""

forcebalance_input = """$options
jobtype newton
forcefield ffbonded.itp ffnonbonded.itp posre.itp
trust0 0.1
penalty_type L1
search_tolerance 1e-2
duplicate_pnames 1
lm_guess 100.0
penalty_additive 1.0
priors
  PDIHMULS1B : 100.0
  PDIHMULS2B : 100.0
  PDIHMULS3B : 100.0
  PDIHMULS4B : 100.0
  PDIHMULS5B : 100.0
  PDIHMULS6B : 100.0
/priors
$end

$target
type         torsionprofile_gmx
name         {name1}
coords       {name1}_all.g96
attenuate
energy_denom 1.5
energy_upper 10.0
writelevel   2
energy_mode  qm_minimum
$end

$target
type         torsionprofile_gmx
name         {name2}
coords       {name2}_all.g96
attenuate
energy_denom 1.5
energy_upper 10.0
writelevel   2
energy_mode  qm_minimum
$end
"""

installptmpsi = """
if [ -d "venv" ]; then
  echo "Virtual environment already exists"
else
  python -m venv _venv
fi
source _venv/bin/activate

python -m pip install --upgrade pip
git clone git@code.emsl.pnl.gov:cheung_group/ptmpsi.git
cd ptmpsi
python -m pip install .
cd ..
"""

qmmm_slurm = {}

qmmm_slurm["Tahoma"] = """
cat <<EOF >nwchemrc
ffield amber
amber_1 /cluster/apps/nwchem/nwchem/src/data/amber_s/
amber_2 /cluster/apps/nwchem/nwchem/src/data/amber_x/
amber_3 /cluster/apps/nwchem/nwchem/src/data/amber_q/
spce /cluster/apps/nwchem/nwchem/src/data/solvents/spce.rst
EOF

cp ${{SLURM_SUBMIT_DIR}}/*.frg /big_scratch
cp ${{SLURM_SUBMIT_DIR}}/{complex} /big_scratch
cp ${{SLURM_SUBMIT_DIR}}/prepare.nw /big_scratch

srun --mpi=pmi2 -N $SLURM_NNODES -n $SLURM_NPROCS apptainer exec --bind /big_scratch,$NWCHEM_BASIS_LIBRARY,/cluster/apps/nwchem/nwchem/src/data/ $NWBIN nwchem prepare.nw > prepare.log

cleanup
"""

qmmm_slurm["Frontier"] = """
cat <<EOF >nwchemrc
ffield amber
amber_1 /opt/nwchem/share/data/amber_s/
amber_2 /opt/nwchem/share/data/amber_q/
amber_3 /opt/nwchem/share/data/amber_x/
amber_4 /opt/nwchem/share/data/amber_u/
spce /opt/nwchem/share/data/solvents/spce.rst
EOF

cp ${{SLURM_SUBMIT_DIR}}/*.frg ${{SCRATCH}}
cp ${{SLURM_SUBMIT_DIR}}/{complex} ${{SCRATCH}}
cp ${{SLURM_SUBMIT_DIR}}/prepare.nw ${{SCRATCH}}

srun -N $SLURM_NNODES -n $SLURM_NPROCS $NWBIN prepare.nw > prepare.log

cleanup
"""

qmmm_slurm["Polaris"] = """
cat <<EOF >nwchemrc
ffield amber
amber_1 /grand/projects/TwinHostPath/software/nwchem/share/data/amber_s/
amber_2 /grand/projects/TwinHostPath/software/nwchem/share/data/amber_q/
amber_3 /grand/projects/TwinHostPath/software/nwchem/share/data/amber_x/
amber_4 /grand/projects/TwinHostPath/software/nwchem/share/data/amber_u/
spce /grand/projects/TwinHostPath/software/nwchem/share/data/solvents/spce.rst
EOF

cp ${{SLURM_SUBMIT_DIR}}/*.frg ${{SCRATCH}}
cp ${{SLURM_SUBMIT_DIR}}/{complex} ${{SCRATCH}}
cp ${{SLURM_SUBMIT_DIR}}/prepare.nw ${{SCRATCH}}

mpiexec -hostfile $PBS_NODEFILE -n ${{NTOTRANKS}} -ppn ${{NRANKS_PER_NODE}} --depth=${{NDEPTH}} --cpu-bind core apptainer exec --bind $BINDS --workdir `pwd` $NWBIN nwchem prepare.nw > prepare.log

cleanup
"""

nwconstraint = "constrain  {: 10.6f}  {:5d}\n "
pyconstraint = "[{},{}],\n"
coordinates = "{}   {: 14.8f}   {: 14.8f}   {: 14.8f}\n"
pyprint = """print("{name}: {{:10.6f}}".format(q[{atom}]))\n"""
runsingularity = {}
runsingularity['Tahoma'] = "srun --mpi=pmi2 -N $SLURM_NNODES -n $SLURM_NPROCS apptainer exec --bind {scratch},$NWCHEM_BASIS_LIBRARY $NWBIN nwchem {name}.nw > {name}.log\n\n"
runsingularity['Frontier'] = "srun -N $SLURM_NNODES -n $SLURM_NPROCS apptainer exec --bind $BINDS --workdir `pwd` $NWBIN nwchem {name}.nw > {name}.log\n\n"
runsingularity['Polaris'] = "mpiexec -hostfile $PBS_NODEFILE -n ${{NTOTRANKS}} -ppn ${{NRANKS_PER_NODE}} --depth=${{NDEPTH}} --cpu-bind=core $NWBIN {name}.nw > {name}.log\n\n"
runsingularity_prefix = {}
runsingularity_prefix_tahoma = "srun --mpi=pmi2 -N $SLURM_NNODES -n $SLURM_NPROCS apptainer exec --bind {scratch},$NWCHEM_BASIS_LIBRARY $NWBIN nwchem"
runsingularity_prefix_frontier = "srun -N $SLURM_NNODES -n $SLURM_NPROCS apptainer exec --bind $BINDS --workdir `pwd` $NWBIN nwchem"
script_copy = {}
script_copy['slurm'] = "cp ${{SLURM_SUBMIT_DIR}}/{filename} . \n"
script_copy['pbs'] = "cp $PBS_O_WORKDIR/{filename} . \n"
