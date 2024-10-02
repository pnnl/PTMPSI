class Machine:
    def __init__(self, name=None, partitions=None, scratchdir=None, modules=None):
        self.name       = name
        self.partitions = partitions
        self.scratchdir = scratchdir
        self.modules    = modules
        self.account    = None

    def default_partition(self):
        if self.partitions is None:
            raise RuntimeError("Machine does not have a list of partitions")
        for partition in self.partitions:
            if self.partitions[partition].default == True:
                return self.partitions[partition]
        return self.partitions[0]

    def default_account(self):
        return "XXXX" if self.account is None else self.account
        

class Partition:
    def __init__(self, name=None, memory=4, ncpus=1, ngpus=0, maxtime=12, maxnode=1, options={}):
        self.ncpus   = ncpus
        self.ngpus   = ngpus
        self.memory  = memory
        self.name    = name
        self.maxtime = maxtime
        self.maxnode = maxnode
        self.default = False
        self.options = options


aqe_ldrd = Machine(
  name="AQE-LDRD",
  partitions={
    "basic"   : Partition(
                  name="bsc120c", 
                  memory=456, 
                  ncpus=120,        
                  ngpus=0, 
                  maxtime=0, 
                  maxnode=200
                  ),
    "standard": Partition(
                  name="std120c", 
                  memory=456, 
                  ncpus=120, 
                  ngpus=0 ,
                  maxtime=0, 
                  maxnode=10
                  ),
    "premium":  Partition(
                  name="prm40c8g", 
                  memory=672, 
                  ncpus=40, 
                  ngpus=8, 
                  maxtime=0, 
                  maxnode=8,
                  options = {
                    "gromacs": {
                      "mpirun": "mpirun -np 4 --rankfile rankfile1 --bind-to core",
                      "container": "apptainer exec --nv --bind /anfhome,/mnt,/etc,/sched,/run,/opt/hpcx $myimage",
                      "gmx": "gmx_mpi",
                      "gpu_id": "0123",
                      "ntasks": 8,
                      "nthreads": 5,
                      "nstlist": 200
                      }
                    }
                  ),
    "prm96c4g":  Partition(
                  name="prm96c4g", 
                  memory=672, 
                  ncpus=96, 
                  ngpus=4, 
                  maxtime=0, 
                  maxnode=3,
                  options = {
                    "gromacs": {
                      "mpirun": "mpirun -np 4 --bind-to core --rankfile rankfileA100",
                      "container": "apptainer exec --nv --bind /anfhome,/mnt,/etc,/sched,/run,/opt/hpcx $myimage",
                      "gmx": "gmx_mpi",
                      "gpu_id": "0123",
                      "ntasks": 4,
                      "nthreads": 24,
                      "nstlist": 200
                      }
                    }
                  ),
    },
    scratchdir="/mnt/scratch"
  )
aqe_ldrd.partitions["prm40c8g"] = aqe_ldrd.partitions["premium"]
aqe_ldrd.partitions["std120c"]  = aqe_ldrd.partitions["standard"]
aqe_ldrd.partitions["bsc120c"]  = aqe_ldrd.partitions["basic"]
aqe_ldrd.partitions["premium"].default = True



aqe_h100 = Machine(name="AQE-H100",
  partitions={
    "premium": Partition(
      name="prm96c8g", 
      memory=0, 
      ncpus=96,
      ngpus=8, 
      maxtime=0, 
      maxnode=28,
                  options = {
                    "gromacs": {
                      "mpirun": "mpirun -np 4 --bind-to core --rankfile rankfileH100",
                      "container": "apptainer exec --nv --bind /anfhome,/mnt,/etc,/sched,/run,/opt/hpcx $myimage",
                      "gmx": "gmx_mpi",
                      "gpu_id": "0123",
                      "ntasks": 4,
                      "nthreads": 12,
                      "nstlist": 300
                      }
                    }
      )
            
    },
    scratchdir="/mnt/scratch"
  )
aqe_h100.partitions["prm96c8g"] = aqe_h100.partitions["premium"]
aqe_h100.partitions["premium"].default = True

tahoma = Machine(name="Tahoma",
        partitions={
            "normal": Partition(name="normal", memory=0, ncpus=36,
                              ngpus=0, maxtime=72, maxnode=10),
            "analysis": Partition(name="analysis", memory=0, ncpus=36,
                              ngpus=1, maxtime=3, maxnode=1)
            },
        modules = {
            },
        scratchdir="/big_scratch")
tahoma.partitions["normal"].default = True


deception = Machine(name="Deception",
        partitions={
            "slurm": Partition(name="slurm", memory=0, ncpus=64,
                              ngpus=0, maxtime=168, maxnode=10),
            "short": Partition(name="short", memory=0, ncpus=64,
                              ngpus=0, maxtime=3, maxnode=7)
            },
        modules={"apptainer": "apptainer/1.2.5",
                 "gcc": "gcc/11.2.0",
                 "python": "python/3.11.5",
                 "openmpi": "openmpi/4.1.0"},
        scratchdir="/scratch")
deception.partitions["slurm"].default = True

perlmutter = Machine(name="Perlmutter",
        partitions={
            "regular": Partition(name="regular", memory=0, ncpus=64,
                              ngpus=0, maxtime=48, maxnode=10),
            "preempt": Partition(name="preempt", memory=0, ncpus=64,
                              ngpus=0, maxtime=48, maxnode=128),
            "shared": Partition(name="shared", memory=0, ncpus=64,
                              ngpus=0, maxtime=48, maxnode=0.5)
            },
        modules={"apptainer": "apptainer/1.2.5",
                 "gcc": "gcc/11.2.0",
                 "python": "cray-python/3.11.5",
                 "openmpi": "openmpi/5.0.3"},
        scratchdir="$SCRATCH")
perlmutter.partitions["regular"].default = True

frontier = Machine(name="Frontier",
        partitions={
            "batch": Partition(name="batch", memory=0, ncpus=56,
                              ngpus=8, maxtime=2, maxnode=8192,
                              options = {
                    "gromacs": {
                      "mpirun": "srun -N1 -n8 -c7 --gpus-per-task=1 --gpu-bind=closest",
                      "container": "",
                      "gmx": "gmx_mpi",
                      "ntasks": 8,
                      "nthreads": 7,
                      "nstlist": 300
                      }
                    }),
            "extended": Partition(name="extended", memory=0, ncpus=56,
                              ngpus=8, maxtime=24, maxnode=64,
                              options = {
                    "gromacs": {
                      "mpirun": "srun -N1 -n8 -c7 --gpus-per-task=1 --gpu-bind=closest",
                      "container": "",
                      "gmx": "gmx_mpi",
                      "ntasks": 8,
                      "nthreads": 7,
                      "nstlist": 300
                      }
                    }),
            },
        modules={"apptainer": "apptainer/1.2.5",
                 "gcc": "gcc/11.2.0",
                 "python": "cray-python/3.11.5",
                 "openmpi": "openmpi/5.0.3"},
        scratchdir="/lustre/orion/bip258/scratch/$USER")
frontier.partitions["batch"].default = True

polaris = Machine(name="Polaris",
        partitions={
            "debug": Partition(name="debug", memory=0, ncpus=32,
                              ngpus=4, maxtime=1, maxnode=2,
                              options = {
                    "gromacs": {
                      "mpirun": "mpiexec -hostfile $PBS_NODEFILE -n $NTOTRANKS -ppn $NRANKS_PER_NODE --depth=$NDEPTH --env OMP_NUM_THREADS=$NTHREADS --cpu-bind=core",
                      "container": "",
                      "gmx": "/home/wjiang/gromacs-2024.2/build/bin/gmx_mpi",
                      "gpu_tasks": "00112233",
                      "ntasks": 4,
                      "nthreads": 8,
                      "nstlist": 300
                      }
                    }),
            "prod": Partition(name="prod", memory=0, ncpus=32,
                              ngpus=4, maxtime=24, maxnode=496,
                              options = {
                    "gromacs": {
                      "mpirun": "mpiexec -hostfile $PBS_NODEFILE -n $NTOTRANKS -ppn $NRANKS_PER_NODE --depth=$NDEPTH --env OMP_NUM_THREADS=$NTHREADS --cpu-bind=core",
                      "container": "",
                      "gmx": "/home/wjiang/gromacs-2024.2/build/bin/gmx_mpi",
                      "gpu_tasks": "00112233",
                      "ntasks": 4,
                      "nthreads": 8,
                      "nstlist": 300
                      }
                    }),
            "preemptable": Partition(name="preemptable", memory=0, ncpus=32,
                              ngpus=4, maxtime=72, maxnode=10,
                              options = {
                    "gromacs": {
                      "mpirun": "mpiexec -hostfile $PBS_NODEFILE -n $NTOTRANKS -ppn $NRANKS_PER_NODE --depth=$NDEPTH --env OMP_NUM_THREADS=$NTHREADS --cpu-bind=core",
                      "container": "",
                      "gmx": "/home/wjiang/gromacs-2024.2/build/bin/gmx_mpi",
                      "gpu_tasks": "00112233",
                      "ntasks": 4,
                      "nthreads": 8,
                      "nstlist": 300
                      }
                    }),
            },
        modules={"apptainer": "apptainer/1.2.2",
                 "gcc": "gcc/11.2.0",
                 "python": "cray-python/3.11.5",
                 "openmpi": "openmpi/5.0.3"},
        scratchdir="/local/scratch")
polaris.partitions["prod"].default = True

machines = {'aqe_ldrd': aqe_ldrd,
            'aqe_h100': aqe_h100,
            'tahoma': tahoma,
            'deception': deception,
            'perlmutter': perlmutter,
            'frontier': frontier,
            'polaris': polaris
           }

default_machine = "tahoma"

class Slurm:
    def __init__(self, caller, **kwargs):
        if caller == "nwchem":
            from ptmpsi.nwchem.templates import slurm_header
        elif caller == "alphafold":
            from ptmpsi.alphafold.templates import slurm_header
        elif caller == "gromacs":
            from ptmpsi.gromacs.templates import slurm_header

        __machine = kwargs.pop("machine", default_machine)
        if isinstance(__machine, Machine):
            self.machine = __machine
        elif isinstance(__machine, str):
            self.machine = machines[__machine]
        else:
            self.machine = None
        if not isinstance(self.machine, Machine):
            raise KeyError(f"Machine is not an instance of the Machine class")

        self.partition = kwargs.pop("partition", self.machine.default_partition().name)
        _partition = None
        for k,v in self.machine.partitions.items():
            if self.partition == k:
                self.partition = v.name
                _partition = v
                break
            elif self.partition == v.name:
                _partition = v
                break
        if _partition is None:
            raise KeyError(f"Partition '{self.partition}' is not part of the Machine")

        self.ngpus = kwargs.pop("ngpus", _partition.ngpus)
        if self.ngpus > _partition.ngpus:
            raise KeyError(f"Partition '{self.partition}' does not have '{self.ngpus}' GPUs available")
        self.gpuline = f"#SBATCH --gres=gpu:{self.ngpus}" if self.ngpus > 0 else "" 

        self.ncpus = kwargs.pop("ntasks", _partition.ncpus)
        if self.ncpus > _partition.ncpus:
            raise KeyError(f"Partition '{self.partition}' does not have '{self.ncpus}' CPUs available")

        self.nthreads = kwargs.pop("nthreads", 1)
        if self.ncpus*self.nthreads > _partition.ncpus:
            raise KeyError(f"Partition '{self.partition}' can only use {_partition.ncpus} or less total threads (ntasks*nthreads)")

        _time = kwargs.pop("time", _partition.maxtime if _partition.maxtime > 0 else 144)
        if _partition.maxtime > 0 and _time > _partition.maxtime:
            raise KeyError(f"Partition '{self.partition}' has a maximum time policy of '{_partition.maxtime}' per job")
        self.time = f"{_time}:00:00"

        self.account = kwargs.pop("account", self.machine.default_account())

        self.nnodes  = kwargs.pop("nnodes", 1)
        if self.nnodes > _partition.maxnode:
            raise KeyError(f"Partition '{self.partition}' has only {_partition.maxnode} nodes available")

        self.scratch = kwargs.pop("scratch", self.machine.scratchdir)

        self.jobname = kwargs.pop("jobname", f"ptmpsi_{caller}")

        self.user    = kwargs.pop("slurm_header", "")

        self._header_template = slurm_header[self.machine.name]

        self.options_dictionary = {
                "partition": self.partition,
                "account"  : self.account,
                "time"     : self.time,
                "nnodes"   : self.nnodes,
                "ntasks"   : self.ncpus,
                "ncpus"    : self.ncpus,
                "nthreads" : self.nthreads,
                "jname"    : self.jobname,
                "scratch"  : self.scratch,
                "user"     : self.user
                }
        self.options_dictionary["ngpus"] = self.ngpus


        self.header = self._header_template.format(**self.options_dictionary)

        return


