import pkgutil
from datetime import datetime
from os.path import isfile, isdir
from ptmpsi.slurm import Slurm
from ptmpsi.alphafold.templates import slurm_body

tahoma_datasets = "/tahoma/datasets/alphafold"
tahoma_datasets_231 = "/tahoma/datasets/alphafold-20230208"
edo_singularity = "oras://ghcr.io/edoapra/alphafold_singularity/alphafold:latest"
singularity_231 = "oras://ghcr.io/dmejiar/alphafold_singularity/alphafold_v231:latest"
singularity_232 = "oras://ghcr.io/dmejiar/alphafold_singularity/alphafold_v232:latest"

class AlphaFoldOptions:
    def __init__(self, fasta_paths, **kwargs):
        self.model     = kwargs.get("model","monomer")
        self.version   = kwargs.get("version","2.3.1")
        self.dbs       = kwargs.get("dbs","full_dbs")
        self.container = kwargs.get("container", None)
        self.use_gpu       = kwargs.get("use_gpu", True)
        self.gpu_relax     = kwargs.get("enable_gpu_relax", True)
        self.date          = kwargs.get("date",datetime.today().strftime('%Y-%m-%d'))
        self.run_relax     = kwargs.get("run_relax", True)
        self.models_to_relax = kwargs.get("models_to_relax", "all")
        self.num_multimer  = kwargs.get("num_multimer",5)
        self.xla_mem_fraction = kwargs.get("mem_fraction",4.0)
        self.unified_memory   = kwargs.get("unified_memory",1)
        self.fasta_paths      = fasta_paths
        self.data_dir         = kwargs.get("data_dir", None)

        if self.container is None:
            if self.version == "2.3.2":
                self.container = singularity_232
            elif self.version == "2.3.1":
                self.container = singularity_231
            else:
                self.container = edo_singularity

        if self.data_dir is None:
            if self.version in ["2.3.1","2.3.2"]:
                self.data_dir = tahoma_datasets_231
            else:
                self.data_dir = tahoma_datasets

        if self.run_relax is False and self.version == "2.3.2":
            self.models_to_relax = "none"
            
        return



def gen_script(options):
    if options.version == "2.3.1":
        filename = "run_singularity_231.py"
    elif options.version == "2.3.2":
        filename = "run_singularity_232.py"
    else:
        filename = "run_singularity.py"
    data = pkgutil.get_data(__name__,filename).decode('utf-8')
    with open("run_singularity.py","w") as fh:
        fh.write(data.format(options.unified_memory, options.xla_mem_fraction))
    return

def gen_pull(options):
    # Check if Singularity container exists
    if isfile(options.container):
        print("\t Info: run script will use '{}' container".format(container))
        string = "ln -sf {} $ALPHAFOLD_DIR/alphafold.sif\n".format(container)
    else:
        _split = options.container.split(":")
        if _split[0] in ["oras", "library", "docker", "shub", "http", "https"]:
            print("\t Info: run script will pull container from '{}'".format(options.container))
            string = "/usr/bin/time -p singularity pull -F --name $ALPHAFOLD_DIR/alphafold.sif {}\n".format(options.container)
    return string
            

def gen_fasta_paths(fasta):
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
    return _files

def prediction(fasta,**kwargs):
    """
    Predicts the 3D structure from a given FASTA sequence.

    fasta could be a string with the FASTA sequence, or a path to a given file
    """
    # Generate fasta paths
    fasta_paths = gen_fasta_paths(fasta)

    # Override default options with user's settings
    options = AlphaFoldOptions(fasta_paths, **kwargs)
    slurm   = Slurm("alphafold", **kwargs)
    if options.use_gpu == "True":
        slurm.header += "#SBATCH --gres=gpu:1\n"

    # Generate run_singularity script
    gen_script(options)
    pull    = gen_pull(options)

    if options.version == "2.3.2":
        relax = f"--models_to_relax={options.models_to_relax}"
    else:
        relax = f"--run_relax={options.run_relax}"

    with open("alphafold.sbatch","w") as fh:
        fh.write(slurm_body.format(header=slurm.header,
                 data_dir=options.data_dir, 
                 version=options.version,
                 pull=pull,
                 dbs=options.dbs,
                 use_gpu=options.use_gpu,
                 enable_gpu_relax=options.gpu_relax,
                 relax=relax,
                 fasta_paths=",".join(options.fasta_paths),
                 date=options.date,
                 model=options.model))

    print("\n\t Info: Make sure to copy alphafold.sbatch, run_singularity.py, your")
    print("\t       fasta files and all generated temp*.fasta files to a submission")
    print("\t       directory and execute `sbatch alphafold.sbatch`")
    return
