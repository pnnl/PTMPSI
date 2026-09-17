#!/bin/bash
#SBATCH --partition=batch
#SBATCH --account=bip258
#SBATCH --time=12:00:00
#SBATCH --nodes=1012
#SBATCH --job-name=MED4md_bundle
#SBATCH --error=md_bundle-%j.err
#SBATCH --output=md_bundle-%j.out

module use /path/to/apps/modulefiles
module load gromacs/2025.0
module load hwloc/2.11.1-gpu # Flux requires a GPU-enabled hwloc to see the GPUs
module load flux

echo $SLURM_JOB_ID > {prefix}{jobname}current.jobid
srun -N $SLURM_NNODES -n $SLURM_NNODES -c 56 --gpus-per-node=8 flux start ./md_bundle_flux.sh
