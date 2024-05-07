ions = """; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
"""

minim = """; minim.mdp - used as input into grompp to generate minim.tpr
title           = Minimization
{restraint}

; Parameters describing what to do, when to stop and what to save
integrator      = steep
emtol           = 10.0
emstep          = 0.01
nsteps          = 50000

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = cutoff
rvdw            = 1.2
rvdw-switch     = 1.0
vdw-modifier    = force-switch
pbc             = xyz
DispCorr        = no
"""

heating = """; heating.mdp - used as input to grompp to generate heating.tpr
title           = Heating
define          = -DPOSRES

integrator      = md
dt              = {timestep}
nsteps          = {nsteps}

; Bond parameters
continuation          = no
constraints           = h-bonds
constraint-algorithm  = lincs
lincs-iter            = 1
lincs-order           = 4

; Output control
nstxout-compressed    = 5000
nstenergy             = 500
nstlog                = 500
nstdisreout           = 0

; Neighbor searching
cutoff-scheme         = Verlet
nstlist               = 10
ns_type               = grid
rlist                 = 1.2
vdwtype               = cutoff
vdw-modifier          = force-switch
rcoulomb              = 1.2
rvdw                  = 1.2
rvdw-switch           = 1.0

; Electrostatics
coulombtype           = PME
pme_order             = 4
fourierspacing        = 0.16

; Temperature coupling
Tcoupl                = V-rescale
tc_grps               = Protein Non-Protein
tau_t                 = 0.1 0.1
ref_t                 = {temp} {temp}

; Pressure coupling
Pcoupl                = no

; Periodic Boundary conditions
pbc                   = xyz

; Velocity generation
gen-vel               = yes
gen-temp              = {temp}
gen-seed              = -1

; Dispersion Correction
DispCorr              = EnerPres

disre                 = simple
disre-weighting       = conservative
disre-fc              = 1000
disre-mixed           = yes

nstcomm               = 1
comm-mode             = Linear
comm-grps             = Protein Water_and_ions
"""

npt = """; npt.mdp - used as input to grompp to generate npt.tpr
title           = NPT equilibration
{restraint}

integrator      = md
dt              = {timestep}
nsteps          = {nsteps}

; Bond parameters
continuation          = yes
constraints           = h-bonds
constraint-algorithm  = lincs
lincs-iter            = 1
lincs-order           = 4

; Output control
nstxout-compressed    = 5000
nstenergy             = 500
nstlog                = 500
nstdisreout           = 0

; Neighbor searching
cutoff-scheme         = Verlet
nstlist               = 20
ns-type               = grid
rlist                 = 1.2
vdw-modifier          = force-switch
vdwtype               = cutoff
rcoulomb              = 1.2
rvdw                  = 1.2
rvdw-switch           = 1.0

; Electrostatics
coulombtype           = PME
pme-order             = 4
fourierspacing        = 0.16

; Temperature coupling
Tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
tau-t                 = 0.1 0.1
ref-t                 = {temp} {temp}

; Pressure coupling
Pcoupl                = C-rescale
pcoupltype            = isotropic
tau-p                 = 1.0
ref-p                 = 1.0
compressibility       = 4.5e-5
refcoord-scaling      = com

; Periodic Boundary conditions
pbc                   = xyz

; Velocity generation
gen-vel               = no

; Dispersion Correction
DispCorr              = EnerPres

; Distance Restraint
disre                 = simple
disre-weighting       = conservative
disre-fc              = 1000
"""

md = """
title           = MD simulation

integrator      = md
dt              = {timestep}
nsteps          = {nsteps}

; Bond parameters
continuation          = yes
constraints           = h-bonds
constraint-algorithm  = lincs
lincs-iter            = 1
lincs-order           = 4
comm-mode             = linear

; Output control
nstxout-compressed    = 5000
nstenergy             = 5000
nstlog                = 5000
nstxout               = 0
nstvout               = 0

; Neighbor searching
cutoff-scheme         = Verlet
nstlist               = 100
ns-type               = grid
rlist                 = 1.0
vdw-modifier          = force-switch
vdwtype               = cutoff
rcoulomb              = 1.2
rvdw                  = 1.2
rvdw-switch           = 1.0
verlet-buffer-tolerance = 2.0e-3

; Electrostatics
coulombtype           = PME
pme-order             = 4
fourierspacing        = 0.16

; Temperature coupling
Tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
tau-t                 = 0.1 0.1
ref-t                 = {temp} {temp}

; Pressure coupling
Pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau-p                 = 2.0
ref-p                 = 1.0
compressibility       = 4.5e-5

; Periodic Boundary conditions
pbc                   = xyz

; Velocity generation
gen-vel               = no

; Dispersion Correction
DispCorr              = EnerPres

; Distance Restraint
nstdisreout           = 0
disre                 = simple
disre-weighting       = conservative
disre-fc              = 1000
"""


SLURM_tail = """cleanup\n"""

slurm_header = {}

slurm_header['Tahoma'] = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --job-name={jname}
#SBATCH --error={jname}-%j.err
#SBATCH --output={jname}-%j.out
#SBATCH --partition={partition}
{gpu}

trap cleanup SIGINT SIGTERM SIGKILL SIGSEGV SIGCONT

cleanup()
{{
  cp /fast_scratch/*.gro ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.xtc ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.ndx ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.cpt ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.log ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.tpr ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.edr ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.ndx ${{SLURM_SUBMIT_DIR}}/{results} || :
  cp /fast_scratch/*.top ${{SLURM_SUBMIT_DIR}}/{results} || :
}}

source /etc/profile.d/modules.sh
module use /tahoma/emsle60550/modulefiles
module load gromacs/2023

mkdir -p ${{SLURM_SUBMIT_DIR}}/{results}

cd /fast_scratch
cp ${{SLURM_SUBMIT_DIR}}/{infile} .
cp ${{SLURM_SUBMIT_DIR}}/*.mdp .
cp -r ${{SLURM_SUBMIT_DIR}}/{ff} .

"""

cleanup = """
  cp {{scratch}}/*.gro ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.xtc ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.ndx ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.cpt ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.log ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.tpr ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.edr ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.ndx ${{SLURM_SUBMIT_DIR}}/{{results}} || :
  cp {{scratch}}/*.top ${{SLURM_SUBMIT_DIR}}/{{results}} || :
"""

slurm_header['AQE-LDRD'] = """#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --gpus-per-node={ngpus}
#SBATCH --cpus-per-task={nthreads}
#SBATCH --job-name={jname}
#SBATCH --get-user-env
#SBATCH --exclusive
#SBATCH --error={jname}-%j.err
#SBATCH --output={jname}-%j.out

# NVIDIA LIBRARIES
source /anfhome/.profile
export NVHPC_ROOT=/anfhome/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placehold/linux-almalinux8-x86_64_v3/gcc-8.5.0/nvhpc-23.7-7xowtxnqw4gn4fg2cqwvqynpy5qx32lj/Linux_x86_64/23.7
export LD_LIBRARY_PATH=${{NVHPC_ROOT}}/math_libs/12.2/lib64:${{LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH=${{NVHPC_ROOT}}/cuda/12.2/lib64:${{LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH=${{NVHPC_ROOT}}/comm_libs/12.2/nvshmem_cufftmp_compat/lib:${{LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH=${{NVHPC_ROOT}}/comm_libs/12.2/nccl/lib:${{LD_LIBRARY_PATH}}

NTASKS=$SLURM_NTASKS
CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
NGPUS=$((SLURM_GPUS_PER_NODE * SLURM_NNODES))

NTASKS=$SLURM_NTASKS
CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
NGPUS=$((SLURM_GPUS_PER_NODE * SLURM_NNODES))

export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK}}
export TMPDIR={scratch}
export APPTAINER_CACHEDIR=$TMPDIR
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_CUDA_GRAPH=1
export GMX_MAXBACKUP=-1
export UCX_POSIX_USE_PROC_LINK=n
export UCX_TLS=^cma
export UCX_LOG_LEVEL=ERROR
export HWLOC_HIDE_ERRORS=1
export APPTAINERENV_HWLOC_HIDE_ERRORS=1
export APPTAINERENV_UCX_LOG_LEVEL=ERROR
export APPTAINERENV_UCX_TLS=^cma
export APPTAINERENV_UCX_POSIX_USE_PROC_LINK=n
export APPTAINERENV_GMX_ENABLE_DIRECT_GPU_COMM=${{GMX_ENABLE_DIRECT_GPU_COMM}}
export APPTAINERENV_GMX_GPU_PME_DECOMPOSITION=${{GMX_GPU_PME_DECOMPOSITION}}
export APPTAINERENV_GMX_CUDA_GRAPH=${{GMX_CUDA_GRAPH}}
export APPTAINERENV_GMX_MAXBACKUP=${{GMX_MAXBACKUP}}
export APPTAINERENV_OMP_NUM_THREADS=${{OMP_NUM_THREADS}}
export APPTAINERENV_LD_LIBRARY_PATH="${{LD_LIBRARY_PATH}}:\$LD_LIBRARY_PATH"
myimage=/anfhome/daniel.mejia/sources/apptainer/gromacs_cufftmp_clean.simg

cat > rankfile1  <<EOF
rank 0=+n0 slot=0:0-4
rank 1=+n0 slot=0:5-9
rank 2=+n0 slot=0:10-14
rank 3=+n0 slot=0:15-19
EOF


cat > rankfile2  <<EOF
rank 0=+n0 slot=1:0-4
rank 1=+n0 slot=1:5-9
rank 2=+n0 slot=1:10-14
rank 3=+n0 slot=1:15-19
EOF

"""
