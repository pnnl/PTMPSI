slurm_header = {}
slurm_body = {}


slurm_header["AQE-LDRD"] = """#!/bin/bash
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --job-name={jname}
#SBATCH --account={account}
#SBATCH --output={jname}-%j.out
#SBATCH --error={jname}-%j.out
# {scratch}
# {nthreads}
"""

slurm_header["Tahoma"] = """#!/bin/bash
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --job-name={jname}
#SBATCH --account={account}
#SBATCH --output={jname}-%j.out
#SBATCH --error={jname}-%j.out
# {scratch}
# {nthreads}
"""

slurm_header["Perlmutter"] = """#!/bin/bash
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --time={time}
#SBATCH --qos={partition}
#SBATCH --job-name={jname}
#SBATCH --account={account}
#SBATCH --output={jname}-%j.out
#SBATCH --error={jname}-%j.out
#SBATCH --licenses=scratch,cfs
"""

slurm_header["Frontier"] = """#!/bin/bash
#SBATCH --nodes={nnodes}
#SBATCH --cpus-per-task={ncpus}
#SBATCH --ntasks-per-node=1
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --job-name={jname}
#SBATCH --account={account}
#SBATCH --output={jname}-%j.out
#SBATCH --error={jname}-%j.out
#SBATCH --gpus-per-node={ngpus}
"""

slurm_body["Tahoma"] = """{header}
source /etc/profile.d/modules.sh
module load python/3.8.1

export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1

cd /big_scratch

# Create a Virtual Environment
if [ -d "venv" ]; then
  echo "Virtual environment already exists"
else
  python -m venv venv
fi
source venv/bin/activate

# Set proxy server
export https_proxy=http://proxy.emsl.pnl.gov:3128

# Remove previous alphafold leftovers
rm -rf alphafold

# Set useful envinroment variables
export DOWNLOAD_DIR={data_dir}
export ALPHAFOLD_DIR=/big_scratch/alphafold
export MODFILES_DIR=${{SLURM_SUBMIT_DIR}}
export TMPDIR=/big_scratch/alphafold_run
export ALPHAFOLD_VERSION={version}

mkdir -p alphafold

# Pull singularity container
{pull}

# Make tmp directory
mkdir -p $TMPDIR

python -m pip install --upgrade pip
python -m pip install absl-py spython

cp $MODFILES_DIR/*.fasta .
cp $MODFILES_DIR/run_singularity.py .

python ./run_singularity.py \\
            --data_dir=$DOWNLOAD_DIR \\
            --model_preset={model} \\
            --max_template_date={date} \\
            --fasta_paths={fasta_paths} \\
            --db_preset={dbs} \\
            --use_gpu={use_gpu} \\
            --enable_gpu_relax={enable_gpu_relax} \\
            {relax}

cp -rp $TMPDIR ${{SLURM_SUBMIT_DIR}}/AF_results.$SLURM_JOBID
"""

slurm_body["Perlmutter"] = """{header}
module load python/3.9

export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1

cd ${{SCRATCH}}

# Create a Virtual Environment
if [ -d "venv" ]; then
  echo "Virtual environment already exists"
else
  python -m venv venv
fi
source venv/bin/activate

# Remove previous alphafold leftovers
rm -rf alphafold

# Set useful envinroment variables
export DOWNLOAD_DIR={data_dir}
export ALPHAFOLD_DIR=${{SCRATCH}}/alphafold
export MODFILES_DIR=${{SLURM_SUBMIT_DIR}}
export TMPDIR=${{SCRATCH}}/alphafold_run
export ALPHAFOLD_VERSION={version}

mkdir -p alphafold

# Pull singularity container
{pull}

# Make tmp directory
mkdir -p $TMPDIR

python -m pip install --upgrade pip
python -m pip install absl-py spython

cp $MODFILES_DIR/*.fasta .
cp $MODFILES_DIR/run_singularity.py .

srun python ./run_singularity.py \\
            --data_dir=$DOWNLOAD_DIR \\
            --model_preset={model} \\
            --max_template_date={date} \\
            --fasta_paths={fasta_paths} \\
            --db_preset={dbs} \\
            --use_gpu={use_gpu} \\
            --enable_gpu_relax={enable_gpu_relax} \\
            {relax}

cp -rp $TMPDIR ${{SLURM_SUBMIT_DIR}}/AF_results.$SLURM_JOBID
"""
slurm_body["Frontier"] = """{header}

export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1

export SCRATCH="/lustre/orion/bip258/scratch/${{USER}}"

cd ${{SCRATCH}}

# Set proxy server
export https_proxy=http://proxy.emsl.pnl.gov:3128

# Remove previous alphafold leftovers
rm -rf alphafold

# Set useful envinroment variables
export DOWNLOAD_DIR={data_dir}
export ALPHAFOLD_DIR=${{SCRATCH}}/alphafold
export MODFILES_DIR=${{SLURM_SUBMIT_DIR}}
export TMPDIR=${{SCRATCH}}/alphafold_run
export ALPHAFOLD_VERSION={version}
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PATH=/ccs/proj/bip258/apps/alphafold/2.3.2/af-conda-env/hh-suite/install/bin:$PATH

mkdir -p alphafold

# Load experimental alphafold module
{pull}

# Make tmp directory
mkdir -p $TMPDIR

cp $MODFILES_DIR/*.fasta .
cp $MODFILES_DIR/run_singularity.py .

srun --gpus-per-node=8 -N1 -n1 -c56 python ./run_singularity.py \\
            --data_dir=$DOWNLOAD_DIR \\
            --model_preset={model} \\
            --max_template_date={date} \\
            --fasta_paths={fasta_paths} \\
            --db_preset={dbs} \\
            --use_gpu={use_gpu} \\
            --enable_gpu_relax={enable_gpu_relax} \\
            {relax}

cp -rp $TMPDIR ${{SLURM_SUBMIT_DIR}}/AF_results.$SLURM_JOBID
"""
