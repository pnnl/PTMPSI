import pkgutil
from datetime import datetime
from os.path import isfile, isdir

tahoma_datasets = "/tahoma/datasets/alphafold"
tahoma_datasets_231 = "/tahoma/datasets/alphafold-20230208"
edo_singularity = "oras://ghcr.io/edoapra/alphafold_singularity/alphafold:latest"
daniel_singularity = "oras://ghcr.io/dmejiar/alphafold_singularity/alphafold_v231:latest"

SLURM_template = """{header}
source /etc/profile.d/modules.sh
module load python

export PYTHONUNBUFFERED=1

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
{env}

mkdir -p alphafold

# Pull singularity container
{pull}

# Make tmp directory
mkdir -p $TMPDIR

python -m pip install --upgrade pip
python -m pip install absl-py spython

cp $MODFILES_DIR/*.fasta .
cp $MODFILES_DIR/run_singularity.py .

{command}

cp -rp $TMPDIR ${{SLURM_SUBMIT_DIR}}/AF_results.$SLURM_JOBID
"""

def gen_script(version="2.3.1"):
    filename = "run_singularity.py" if version == "2.2.4" else "run_singularity_231.py"
    with open("run_singularity.py","wb") as fh:
        data = pkgutil.get_data(__name__,filename)
        fh.write(data)
    return

def gen_header(nodes=1,account=None,partition="normal",time="180:00",jobname="AF",**kwargs):
    string = "#!/bin/bash\n"
    string += "#SBATCH -N {}\n".format(nodes)
    string += "#SBATCH -t {}\n".format(time)
    string += "#SBATCH -p {}\n".format(partition)
    string += "#SBATCH -J {}\n".format(jobname)
    _account = "XXXXXXXXXX" if account is None else account
    string += "#SBATCH -A {}\n".format(_account)
    string += "#SBATCH -o {}-%j.out\n".format(jobname)
    string += "#SBATCH -e {}-%j.out\n".format(jobname)
    string += "#SBATCH --gres=gpu:1"
    return string


def gen_envs(data_dir=tahoma_datasets_231, output_dir="/big_scratch/alphafold_run", version="2.3.1", **kwargs):
    # Check if data_dir exists
    if not isdir(data_dir):
        print("\t Warning: directory '{}' does not exists in this filesystem".format(data_dir))

    string  = f"export DOWNLOAD_DIR={data_dir}\n"
    string += "export ALPHAFOLD_DIR=/big_scratch/alphafold\n"
    string += "export MODFILES_DIR=${SLURM_SUBMIT_DIR}\n"
    string += f"export TMPDIR={output_dir}\n"
    string += f"export ALPHAFOLD_VERSION={version}\n"
    return string

def gen_pull(singularity=edo_singularity,**kwargs):
    # Check if Singularity container exists
    if isfile(singularity):
        print("\t Info: run script will use '{}' container".format(singularity))
        string = "ln -sf {} $ALPHAFOLD_DIR/alphafold.sif\n".format(singularity)
    else:
        _split = singularity.split(":")
        if _split[0] in ["oras", "library", "docker", "shub", "http", "https"]:
            print("\t Info: run script will pull container from '{}'".format(singularity))
            string = "/usr/bin/time -p singularity pull -F --name $ALPHAFOLD_DIR/alphafold.sif {}\n".format(singularity)
    return string
            

def gen_command(fasta,**kwargs):
    dbs = kwargs.get("db_preset","full_dbs")
    model = kwargs.get("model","monomer")
    _date = kwargs.get("date",datetime.today().strftime('%Y-%m-%d'))

    string  = """singularity run \\
   -B $(realpath $DOWNLOAD_DIR):/data \\
   -B $(realpath $DOWNLOAD_DIR/bfd):/data/bfd \\
   -B $(realpath $DOWNLOAD_DIR/mgnify):/data/mgnify \\
   -B $(realpath $DOWNLOAD_DIR/pdb70):/data/pdb70 \\
   -B $(realpath $DOWNLOAD_DIR/uniref30):/data/uniref30 \\
   -B $(realpath $DOWNLOAD_DIR/uniref90):/data/uniref90 \\
   -B .:/etc \\
   --pwd /app/alphafold \\
   --nv $ALPHAFOLD_DIR/alphafold.sif \\
   --data_dir=/data \\
   --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
   --pdb70_database_path=/data/pdb70/pdb70 \\
   --uniclust30_database_path=/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \\
   --uniref90_database_path=/data/uniref90/uniref90.fasta \\
   --mgnify_database_path=/data/mgnify/mgy_clusters_2018_12.fa \\
   --template_mmcif_dir=/data/pdb_mmcif/mmcif_files \\
   --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \\
   --max_template_date={} \\
   --model_preset={} \\
   --fasta_paths={} \\
   --db_preset=full_dbs \\
   --use_gpu_relax=True \\
   --output_dir=/big_scratch/alphafold_run
   """.format(_date,model,",".join(fasta))

    string = """python ./run_singularity.py \\
            --data_dir=$DOWNLOAD_DIR \\
            --model_preset={} \\
            --max_template_date={} \\
            --fasta_paths={},
            --db_preset={}""".format(model,_date,",".join(fasta),dbs)
    return string

def prediction(fasta,**kwargs):
    """
    Predicts the 3D structure from a given FASTA sequence.

    fasta could be a string with the FASTA sequence, or a path to a given file
    """

    version = kwargs.setdefault("version","2.3.1")
    singularity = kwargs.setdefault("singularity",daniel_singularity if version=="2.3.1" else edo_singularity)
    model = kwargs.setdefault("model","monomer")

    # Check if a FASTA sequence list was given
    if isinstance(fasta,list):
        _fasta = fasta
    else:
        _fasta = [fasta]

    # Put all consecutive strings in one file
    j = 0
    _nfiles = 0
    _files  = []
    _newfile = True
    for i,sequence in enumerate(_fasta):
        if not(isfile(sequence)):
            if _newfile:
                j = 0
                _newfile = False
                _nfiles += 1
                _filename = "temp{}.fasta".format(_nfiles)
                _files.append("./"+_filename)
                fh = open(_filename,"w")
            j += 1
            fh.write(">sequence_{}\n".format(j))
            _length = len(sequence)
            for k in range(0,_length,80):
                fh.write("{}\n".format(sequence[k:min(k+80,_length)]))
        else:
            if j>0:
                fh.close()
                j = 0
            _newfile = True
            _files.append("./"+sequence)
    if j>0: fh.close() 
    

    # Generate run_singularity script
    gen_script(version=version)

    # Generate SLURM script file
    header  = gen_header(**kwargs)
    env     = gen_envs(**kwargs)
    pull    = gen_pull(**kwargs)
    command = gen_command(_files,**kwargs)


    with open("alphafold.sbatch","w") as fh:
        fh.write(SLURM_template.format(header=header,pull=pull,command=command,env=env))

    print("\n\t Info: Make sure to copy alphafold.sbatch, run_singularity.py, your")
    print("\t       fasta files and all generated temp*.fasta files to a submission")
    print("\t       directory and execute `sbatch alphafold.sbatch`")
    return
