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


minimcg = """; minim.mdp - used as input into grompp to generate minim.tpr
title           = Minimization
{restraint}

; Parameters describing what to do, when to stop and what to save
integrator      = cg
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
{restraint}

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
tau-t                 = 1.0 1.0
ref-t                 = {temp} {temp}

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

nstcomm               = 100
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
tau-t                 = 1.0 1.0
ref-t                 = {temp} {temp}

; Pressure coupling
Pcoupl                = C-rescale
pcoupltype            = isotropic
tau-p                 = 5.0
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
tau-t                 = 1.0 1.0
ref-t                 = {temp} {temp}

; Pressure coupling
Pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau-p                 = 5.0
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

qlambdas = """
integrator = sd
dt         = 0.001
nsteps     = {nsteps};
nstlog     = 5000 ; update log file every 10.0 ps
nstxout    = 0
nstvout    = 0
nstenergy  = 1000 ; save energy every 10.0 ps
nstxout-compressed = 1000 ; save corrdinates every 10.0 ps

; Bond parameters
continuation            = yes       ; continuing from NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to All are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 8         ; also related to accuracy
comm-mode = linear

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
ns_type         = grid          ; Method to determine neighbor list (simple, grid)
coulombtype     = PME           ; Treatment of long range electrostatic interactions
fourierspacing  = 0.16          ; Grid spacing for FFT
cutoff-scheme   = Verlet
rlist           = 1.0           ; Cut-off for making neighbor list (short range forces)
verlet-buffer-tolerance   = 2e-03
rcoulomb        = 1.2
rvdw            = 1.2           ; long range Van der Waals cut-off
pbc             = xyz           ; Periodic Boundary Conditions (yes/no)

; Temp and Pressure Controls
tcoupl              = v-rescale ;berendsen ; nose-hoover              ; temperature coupling
ref-t               = 300 300
tc-grps             = Protein Water_and_ions
tau-t               = 1.0 1.0
pcoupl              = Parrinello-Rahman
compressibility     = 4.5e-5              ; Compressibility of water (Isothermal), bar^-1
pcoupl_type         = isotropic
tau-p               = 5.0
ref-p               = 1.0
gen_vel             = no                      ; generate initial velocities
gen_temp            = 300                      ; initial temperature

; Alchemistry
free_energy              = yes
delta_lambda             = 0
;
;couple-moltype          = ATC-A
couple-intramol         = no
couple-lambda0          = vdw-q
couple-lambda1          = vdw-q
init_lambda_state       = {lambda_state}
calc_lambda_neighbors   = 1
;                            0    1    2    3    4    5    6    7    8    9   10   11   12
vdw_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
coul_lambdas            = 0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00
; We are not transforming any bonded or restrained interactions
bonded_lambdas          = 0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00
restraint_lambdas       = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas            = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Not doing simulated temperting here
temperature_lambdas     = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = yes       ; yes for vdw??? linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10
"""


vdwlambdas = """
integrator = sd
dt         = 0.001
nsteps     = {nsteps};
nstlog     = 1000 ; update log file every 10.0 ps
nstxout    = 0
nstvout    = 0
nstenergy  = 1000 ; save energy every 10.0 ps
nstxout-compressed = 1000 ; save corrdinates every 10.0 ps

; Bond parameters
continuation            = yes       ; continuing from NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to All are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 8         ; also related to accuracy
comm-mode = linear

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
ns_type         = grid          ; Method to determine neighbor list (simple, grid)
coulombtype     = PME           ; Treatment of long range electrostatic interactions
fourierspacing  = 0.16          ; Grid spacing for FFT
cutoff-scheme   = Verlet
rlist           = 1.0           ; Cut-off for making neighbor list (short range forces)
verlet-buffer-tolerance   = 2e-03
rcoulomb        = 1.2
rvdw            = 1.2           ; long range Van der Waals cut-off
pbc             = xyz           ; Periodic Boundary Conditions (yes/no)

; Temp and Pressure Controls
tcoupl              = v-rescale ;berendsen ; nose-hoover              ; temperature coupling
ref-t               = 300 300
tc-grps             = Protein Water_and_ions
tau-t               = 1.0 1.0
pcoupl              = Parrinello-Rahman
compressibility     = 4.5e-5              ; Compressibility of water (Isothermal), bar^-1
pcoupl_type         = isotropic
tau-p               = 5.0
ref-p               = 1.0
gen_vel             = no                      ; generate initial velocities
gen_temp            = 300                      ; initial temperature

; Alchemistry
free_energy              = yes
delta_lambda             = 0
;
;couple-moltype          = ATC-A
couple-intramol         = vdw-q
couple-lambda0          = vdw-q
couple-lambda1          = none
init_lambda_state       = {lambda_state}
calc_lambda_neighbors   = 1
;                            0    1    2    3    4    5    6    7    8    9   10   11   12
vdw_lambdas             = 0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00
coul_lambdas            = 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
; We are not transforming any bonded or restrained interactions
bonded_lambdas          = 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
restraint_lambdas       = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
mass_lambdas            = 0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00
; Not doing simulated temperting here
temperature_lambdas     = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
; Options for the decoupling
sc-alpha                 = 0.5
sc-coul                  = yes       ; yes for vdw??? linear interpolation of Coulomb (none in this case)
sc-power                 = 1
sc-sigma                 = 0.3
nstdhdl                  = 10
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

# From the user
{user}

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

# From the user
{user}

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
#export GMX_CUDA_GRAPH=1
export GMX_MAXBACKUP=-1
export UCX_POSIX_USE_PROC_LINK=n
export UCX_TLS=^cma
export UCX_LOG_LEVEL=ERROR
export UCX_LOG_LEVEL_TRIGGER=ERROR
export UCX_RNDV_THRESH=8192
export HWLOC_HIDE_ERRORS=1
export APPTAINERENV_HWLOC_HIDE_ERRORS=1
export APPTAINERENV_UCX_LOG_LEVEL=${{UCX_LOG_LEVEL}}
export APPTAINERENV_UCX_LOG_LEVEL_TRIGGER=${{UCX_LOG_LEVEL_TRIGGER}}
export APPTAINERENV_UCX_TLS=${{UCX_TLS}}
export APPTAINERENV_UCX_POSIX_USE_PROC_LINK=${{UCX_POSIX_USE_PROC_LINK}}
export APPTAINERENV_UCX_RNDV_THRESH=${{UCX_RNDV_THRESH}}
export APPTAINERENV_GMX_ENABLE_DIRECT_GPU_COMM=${{GMX_ENABLE_DIRECT_GPU_COMM}}
export APPTAINERENV_GMX_GPU_PME_DECOMPOSITION=${{GMX_GPU_PME_DECOMPOSITION}}
#export APPTAINERENV_GMX_CUDA_GRAPH=${{GMX_CUDA_GRAPH}}
export APPTAINERENV_GMX_MAXBACKUP=${{GMX_MAXBACKUP}}
export APPTAINERENV_OMP_NUM_THREADS=${{OMP_NUM_THREADS}}
export APPTAINERENV_LD_LIBRARY_PATH="${{LD_LIBRARY_PATH}}:\$LD_LIBRARY_PATH"
myimage=/anfhome/shared/gromacs2023.4+plumed+cufftmp.simg

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

cat > rankfileA100 <<EOF
rank 0=+n0 slot=0:0-23
rank 1=+n0 slot=0:24-47
rank 2=+n0 slot=1:0-23
rank 3=+n0 slot=1:24-47
EOF


cat > rankfileH100 <<EOF
rank 0=+n0 slot=0:0-11
rank 1=+n0 slot=0:12-23
rank 2=+n0 slot=0:24-35
rank 3=+n0 slot=0:36-47
rank 4=+n0 slot=1:0-11
rank 5=+n0 slot=1:12-23
rank 6=+n0 slot=1:24-35
rank 7=+n0 slot=1:36-47
EOF

cat > rankfileH100_1 <<EOF
rank 0=+n0 slot=0:0-11
rank 1=+n0 slot=0:12-23
rank 2=+n0 slot=0:24-35
rank 3=+n0 slot=0:36-47
EOF

cat > rankfileH100_2 <<EOF
rank 0=+n0 slot=1:0-11
rank 1=+n0 slot=1:12-23
rank 2=+n0 slot=1:24-35
rank 3=+n0 slot=1:36-47
EOF
"""

SNC = """ N      -0.4157      14.01
 H       0.2719      1.008
CT       0.0213      12.01
H1       0.1124      1.008
CT      -0.1231      12.01
H1       0.1112      1.008
H1       0.1112      1.008
SH      -0.3119      32.06
HS       0.1933      1.008
DU       0.0000         16
 C       0.5973      12.01
 O      -0.5679         16"""

CSO = """ N      -0.4157       14.01
 H       0.2719       1.008
CT       0.0213       12.01
H1       0.1124       1.008
CT      -0.1231       12.01
H1       0.1112       1.008
H1       0.1112       1.008
SH      -0.3119       32.06
HS       0.1933       1.008
DU       0.0000       1.008
 C       0.5973       12.01
 O      -0.5679          16"""

CGL = """ N      -0.4157       14.01
 H       0.2719       1.008
CT       0.0213       12.01
H1       0.1124       1.008
CT      -0.1231       12.01
H1       0.1112       1.008
H1       0.1112       1.008
SH      -0.3119       32.06
DU       0.0000       14.01
DU       0.0000       1.008
DU       0.0000       1.008
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000          16
DU       0.0000          16
DU       0.0000       12.01
DU       0.0000          16
DU       0.0000       14.01
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       1.008
HS       0.1933       1.008
DU       0.0000       12.01
DU       0.0000          16
DU       0.0000       14.01
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000       1.008
DU       0.0000       1.008
DU       0.0000       12.01
DU       0.0000          16
DU       0.0000          16
 C       0.5973       12.01
 O      -0.5679          16"""

IYY = """ N      -0.4157       14.01
 H       0.27190     1.008
CT       0.02130     12.01
H1       0.11240     1.008
CT      -0.12310     12.01
H1       0.11120     1.008
H1       0.11120     1.008
SH      -0.31190     32.06
HS       0.19330     1.008
DU       0.00000     12.01
DU       0.00000     1.008
DU       0.00000     1.008
DU       0.00000     12.01
DU       0.00000     1.008
DU       0.00000     12.01
DU       0.00000        16
DU       0.00000        16
DU       0.00000     14.01
DU       0.00000     1.008
DU       0.00000     1.008
DU       0.00000     1.008
 C       0.59730     12.01
 O      -0.56790        16"""


update_topology = f"""#!/usr/bin/env python3
import os

with open("topol.top", "r") as topo:
  oldtopo = topo.readlines()

SNC = '''{SNC}'''
CSO = '''{CSO}'''
CGL = '''{CGL}'''
IYY = '''{IYY}'''

IYY_list = [
 [[2, 4, 7, 8],
   ["   0.0   1.39467   3     0.0   1.04600   3"]],
 [[4, 7, 8, 9],
   ["   0.0  14.64400   2     0.0  14.64400   2",
    "   0.0   2.51040   3     0.0   2.51040   3"]]
]

CSO_list = [
  [[2, 4, 7, 8], 
    ["180.0   8.02970   1   180.0   0.00000   1", 
     "  0.0   0.52192   2     0.0   0.00000   2", 
     "  0.0   4.95200   3     0.0   1.04600   3"]], 
  [[4, 7, 8, 9], 
    ["  0.0   2.34510   1     0.0   2.34510   1", 
     "  0.0   9.36040   2     0.0   9.36040   2"]]
]

SNC_list = [
  [[2, 4, 7, 8],
    ["180.0   0.23990   1   180.0   0.00000   1",
     "  0.0   1.47420   2     0.0   0.00000   2",  
     "  0.0   0.00000   3     0.0   1.04600   3"]],
  [[4, 7, 8, 9],
    ["  0.0   1.07690   2     0.0   1.07690   2",
     "180.0  27.29300   2   180.0  27.29300   2",
     "180.0   2.72140   3   180.0   2.72140   3",
     "180.0   0.50325   4   180.0   0.50325   4"]]
]

CGL_list = [
 [[2, 4, 7, 32],
   ["   0.0   1.39467   3     0.0   1.04600   3"]],
 [[4, 7, 32, 29],
   ["   0.0  14.64400   2     0.0  14.64400   2",
    "   0.0   2.51040   3     0.0   2.51040   3"]]
]

substitutions = []
with open("TItop.top", "w") as topo:
  dihedrals = False
  iline = 0
  while iline < len(oldtopo):
    fields = oldtopo[iline].split()

    if len(fields) < 1: 
      topo.write(oldtopo[iline])
      iline += 1
      continue

    if not dihedrals and "[ dihedrals ]" in oldtopo[iline]:
      dihedrals = True
      topo.write(oldtopo[iline])
      iline += 1
      continue

    if dihedrals:
      if fields[0] == ";" or len(substitutions) < 1:
        topo.write(oldtopo[iline])
        iline += 1
        continue
      
      found = False
      for isub in range(len(substitutions)):
        if [int(x) for x in fields[:4]] == substitutions[isub][0]:
          found = True
          _oldtopo = oldtopo[iline].strip("\\n")
          for value in substitutions[isub][1]:
            topo.write(f"{{_oldtopo}}  {{value}}\\n")
          substitutions.pop(isub)
          break
      if not found: topo.write(oldtopo[iline])
      iline += 1
      continue
    elif len(fields) < 4:
      topo.write(oldtopo[iline])
      iline += 1
      continue
    elif fields[3] not in ["SNC", "CSO", "CGL", "IYY"]: 
      topo.write(oldtopo[iline])
      iline += 1
      continue
    else:
      topo.write(oldtopo[iline])
      iline += 1

    ptm = eval(f"{{fields[3]}}").splitlines()
    _list = eval(f"{{fields[3]}}_list")
    _first = int(oldtopo[iline].split()[0])
    for item in _list:
      substitutions.append(tuple( ([_first+int(x) for x in item[0]], item[1]) ))

    for jline in range(len(ptm)):
      _oldtopo = oldtopo[iline].strip("\\n")
      _comment = _oldtopo.find(";")
      _oldtopo = _oldtopo[:_comment-3] if _comment > -1 else _oldtopo
      _ptmline = ptm[jline].strip("\\n")
      topo.write(f"{{_oldtopo}}     {{_ptmline}} \\n")
      iline += 1
"""
