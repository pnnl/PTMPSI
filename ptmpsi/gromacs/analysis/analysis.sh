#!/bin/bash
#SBATCH --job-name=analysis_job      # Job name
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=90           # Number of CPUs per task
#SBATCH --mem=100G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=analysis_%j.log     # Standard output and error log
#SBATCH --mail-type=END,FAIL         # Send email at job completion and failure
#SBATCH --mail-user=your-email@domain.com  # Where to send mail

source /path/to/your/virtualenv/activate  # Activate your virtual environment if necessary
module load gcc/13.2.0                   # Load necessary modules
module load intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx

# Run your Python script
srun python analysis.py
