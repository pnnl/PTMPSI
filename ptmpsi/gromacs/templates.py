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
rlist           = 1.0
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = cutoff
rvdw            = 1.0
pbc             = xyz
DispCorr        = no
"""

heating = """; heating.mdp - used as input to grompp to generate heating.tpr
title           = Heating
define          = -DPOSRES

integrator      = md
dt              = 0.002
nsteps          = 50000

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
rlist                 = 1.0
vdwtype               = cutoff
vdw-modifier          = potential-shift
rcoulomb              = 1.0
rvdw                  = 1.0

; Electrostatics
coulombtype           = PME
pme_order             = 4
fourierspacing        = 0.12

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
"""

npt = """; npt.mdp - used as input to grompp to generate npt.tpr
title           = NPT equilibration
{restraint}

integrator      = md
dt              = 0.002
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
rlist                 = 1.0
vdw-modifier          = potential-shift
vdwtype               = cutoff
rcoulomb              = 1.0
rvdw                  = 1.0

; Electrostatics
coulombtype           = PME
pme-order             = 4
fourierspacing        = 0.12

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
dt              = 0.002
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
vdw-modifier          = potential-shift
vdwtype               = cutoff
rcoulomb              = 1.0
rvdw                  = 1.0

; Electrostatics
coulombtype           = PME
pme-order             = 4
fourierspacing        = 0.12

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

; Distance Restraint ; We need to comment out this part as we do not use Dist. Restraint during the MD.
; nstdisreout           = 0
; disre                 = simple
; disre-weighting       = conservative
; disre-fc              = 1000
"""

SLURM_header = """#!/bin/bash
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

{omp}

echo -e "1\\n1" | gmx_mpi pdb2gmx -f {infile} -o step1.gro -merge all
gmx editconf -f step1.gro -o step2.gro -bt {bt} -d 1.0 -c
gmx solvate -cp step2.gro -cs spc216.gro -p topol.top -o step3.gro
gmx grompp -f ions.mdp -c step3.gro -p topol.top -o ions.tpr
echo SOL | gmx genion -s ions.tpr -o step4.gro -p topol.top -pname NA -nname CL {addion}
echo q | gmx make_ndx -f step4.gro -o index.ndx
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

SLURM_tail = """cleanup\n"""
