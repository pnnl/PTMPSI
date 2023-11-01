class Slurm:
    def __init__(self, caller, **kwargs):
        if caller == "nwchem":
            from ptmpsi.nwchem.templates import slurm_header
        elif caller == "alphafold":
            from ptmpsi.alphafold.templates import slurm_header

        self.time      = kwargs.get("time","12:00:00")
        self.partition = kwargs.get("partition","normal")
        self.account   = kwargs.get("account","emsls60202")
        self.nnodes    = kwargs.get("nnodes",1)
        self.ntasks    = kwargs.get("ntasks",36)
        self.scratch   = kwargs.get("scratch","/big_scratch")
        self.jobname   = kwargs.get("jobname",f"ptmpsi_{caller}")            
        self.header    = slurm_header.format(
            account=self.account, time=self.time, nodes=self.nnodes,
            ntasks=self.ntasks, jname=self.jobname, partition=self.partition,
            np=self.nnodes*self.ntasks, scratch=self.scratch)
        return
