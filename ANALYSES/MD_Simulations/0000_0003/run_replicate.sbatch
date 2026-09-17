#!/bin/bash
#SBATCH --account=m4742
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4
#SBATCH --job-name=replicate_md
#SBATCH --error=replicate_%j.err
#SBATCH --output=replicate_%j.out

# One replicate: full equilibration + 100 ns production + postprocessing.
# Equilibration ~20 min, production ~24 h, postprocess ~10 min.

bash run_replicate_pipeline.sh
