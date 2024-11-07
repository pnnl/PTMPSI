import errno
import os
from collections import defaultdict

def check_file_exists(filename, canbenone=False):
    #
    # Check if filename is None
    if filename is None: 
        if canbenone: 
            return
        else:
            raise KeyError("filename cannot be None")
    #
    # Check if file exists
    if not os.path.isfile(filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
    return

def totss(result):
    return result.count("H")+ result.count("E") + result.count("G") + result.count("I") + result.count("P") + result.count("B")

def totab(result):
    return result.count("H") + result.count("E")


def post_process(trajectories, structures, indices=None, prefix=None, suffix=None, interactive=False, maxprocs=90, group="Non-Water"):
    __trajectories = trajectories if isinstance(trajectories, list) else [trajectories]
    __structures   = structures   if isinstance(structures,   list) else [structures]

    # Special handling for indices list
    if indices is None:
        __indices = [None] * len(__trajectories)
    else:
        __indices = indices if isinstance(indices, list) else [indices]
    #
    # Check that lists have the same length
    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))
    #
    # Total number of trajectories to convert
    ntraj = len(__trajectories)
    #
    # Loop over all trajectories
    outputs = []
    for i in range(ntraj):
        #
        # Check that current set of files exists
        trajectory, structure, index = __trajectories[i], __structures[i], __indices[i] 
        check_file_exists(trajectory)
        check_file_exists(structure)
        check_file_exists(index, canbenone=True)
        #
        # All checks have passed
        # Now get the output filenames
        outputs.append("" if prefix is None else f"{prefix}_")
        outputs[-1] += f"{trajectory[:-4]}"
        outputs[-1] += "" if suffix is None else f"_{suffix}"
    #
    # Run interactively
    if interactive:
        import psutil
        import subprocess
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical/nprocs)
        __nprocs  = min(maxprocs, nprocs)
        #
        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")
        #
        commands = "source /anfhome/.profile "
        commands += "&& module  load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx " 
        commands += f"&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && " 
        #
        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids         = []
            for j in range(__nprocs):
                if i+j == ntraj: break
                trajectory, structure, index, output = __trajectories[i+j], __structures[i+j], __indices[i+j], outputs[i+j] 
                #
                # Extract the desired group from the first from of the trajectory 
                __commands = commands
                __commands += f"echo {group} | gmx trjconv -dump 0 -f {trajectory} -s {structure} -o {output}.gro "
                __commands += f"-n {index} && " if index is not None else "&& "
                #
                # Make the molecules in the desired group "Whole"
                __commands += f"echo {group} | gmx trjconv -f {trajectory} -s {structure} -o {output}_whole.xtc -pbc whole "
                __commands += f"-n {index} && " if index is not None else "&& "
                #
                # Remove jumps
                __commands += f"gmx trjconv -f {output}_whole.xtc -o {output}_whole_nojump.xtc -pbc nojump && "
                #
                # Remove rotations and translations
                __commands += f"echo 3 {group} | gmx trjconv -f {output}_whole_nojump.xtc -o {output}_whole_nojump_fit.xtc -fit rot+trans -s {output}.gro "
                __commands += f"-n {index} && " if index is not None else "&& "
                #
                # Generate new TPR file with desired group
                __commands += f"echo {group} | gmx convert-tpr -s {structure} -o {output}_{group}.tpr "
                __commands += f"-n {index}" if index is not None else ""
                #
                # Run subprocess and get current process id (PID)
                subprocesses.append(subprocess.Popen(_commands, shell=True))
                pids.append(subprocesses[-1].pid)
                #
                # Set CPU affinity to avoid oversubscription
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j*nskip])
                #
                # Sleep 2 seconds
                time.sleep(2)
            #
            # Wait for the current batch to finish
            exit_codes = [p.wait() for p in subprocesses]
    #
    # Generate a SLURM file
    else:
        pass
    return

def process_dssp(datfile, ref=0, doss="all"):
    """
    Takes a dssp.dat file, or a list of them, and process it to output the change
    in secondary structure with respect to a reference (taken to be the first dssp.dat
    file. The reference can be changed by passing the ref index.

    The doss parameter selects the configuration that are counted as secondary structure.
    With doss="ab", only Alpha Helices and Beta Sheets are counted ("H" and "E" labels)
    """
    __datfiles = datfile if isinstance(datfile, list) else [datfile]
    ntraj = len(__datfiles)

    # Select what secondary structure to count
    if doss == "all":
        countss = totss
        listss  = ["H", "E", "B", "P", "I", "G"]
    elif doss == "ab":
        countss = totab
        listss  = ["H", "E"]
    else:
        raise KeyError(f"doss value {doss} is not valid")
    
    results = []
    for i in range(0, ntraj):
        basename = __datfiles[i].strip("_dssp.dat")
        with open(__datfiles[i], "r") as fh:
            lines = [x for x in (line.strip() for line in fh) if x]
        nchars    = len(lines[0])
        nbreaks   = lines[0].count("=") 
        nresidues = nchars - nbreaks
        #
        dicts = []
        for k in range(nchars):
            dicts.append({"H": 0, "B": 0, "E": 0, 
                          "G": 0, "I": 0, "P": 0,
                          "S": 0, "T": 0, "~": 0,
                          "=": 0})
        #
        fh_summary = open(f"{basename}_summary.dat", "w")
        fh_summary.write("# Step          H        B        E        G        I        P        S       T        ~\n")
        for j in range(len(lines)):
            step = {"H": 0, "B": 0, "E": 0, 
                    "G": 0, "I": 0, "P": 0,
                    "S": 0, "T": 0, "~": 0,
                    "=": 0}
            for k in range(nchars):
                dicts[k][lines[j][k]] += 1
                step[lines[j][k]] += 1
            #    
            fh_summary.write(f"{j:10d} ")
            for ss in ["H", "B", "E", "G", "I", "P", "S", "T", "~"]:
                fraction = step[ss]/nresidues
                fh_summary.write(f"{fraction:8.4f} ")
            fh_summary.write("\n")
        fh_summary.close()
        #
        result = ""
        for k in range(nchars):
            maxval = 0
            for key, val in dicts[k].items():
                if val > maxval: 
                    char = key
                    maxval = val
            result += char
        #
        results.append(result)
        #
        with open(f"{basename}_result.out", "w") as fh:
            fh.write(result)
            fh.write("\n")
            fh.write(f"Alpha Helix: {result.count('H'):>5d}\n")
            fh.write(f" 3_10 Helix: {result.count('G'):>5d}\n")
            fh.write(f"   Pi Helix: {result.count('I'):>5d}\n")
            fh.write(f"Kappa Helix: {result.count('P'):>5d}\n")
            fh.write(f"Beta Sheets: {result.count('E'):>5d}\n")
            fh.write(f"Beta Bridge: {result.count('B'):>5d}\n")
            fh.write(f"    TotalAA: {nresidues:>5d}\n")
            fh.write(f"     Breaks: {nbreaks:>5d}\n\n\n")
            split = result.split("=")
            for ichain, chain in enumerate(split):
                fh.write(f"chain: {ichain:04d}\n")
                fh.write(chain)
                fh.write("\n")
                fh.write(f"Alpha Helix: {chain.count('H'):>5d}\n")
                fh.write(f" 3_10 Helix: {chain.count('G'):>5d}\n")
                fh.write(f"   Pi Helix: {chain.count('I'):>5d}\n")
                fh.write(f"Kappa Helix: {chain.count('P'):>5d}\n")
                fh.write(f"Beta Sheets: {chain.count('E'):>5d}\n")
                fh.write(f"Beta Bridge: {chain.count('B'):>5d}\n")
                fh.write(f"    TotalAA: {len(chain):>5d}\n\n\n") 
    #
    delspec = open("DELSPEC.out", "w")
    sdiff   = open("RESULT_sDiffRES.out", "w")
    #
    _numres = []
    for i in range(ntraj):
        nchars = len(results[i])
        _numres.append( [ 1 if results[i][j] in listss else 0 for j in range(nchars) ] )
    #   
    delspec.write("Case   totSS       SS+       SS-       Delta\n")
    sdiff.write(" Case   ")
    for i in range(ntraj):
        sdiff.write(f"{i:04d}   ")
    sdiff.write("\n")
    #
    for i in range(ntraj):
        ssp = 0
        ssm = 0
        ssa = countss(results[i])
        for j in range(nchars):
            _res = _numres[i][j] - _numres[ref][j]
            if _res > 0: ssp += 1
            if _res < 0: ssm += 1
        delspec.write(f"{i:04d}   {ssa:5d}     {ssp:5d}     {ssm:5d}      {ssa - countss(results[ref]):+5d}\n")  
    ibreak = 0
    for j in range(len(results[ref])):
        if results[ref][j] == "=": 
            ibreak += 1
            continue
        #   
        sdiff.write(f"{j-ibreak+1:5d} ")
        for i in range(ntraj):
            _res = _numres[i][j] - _numres[ref][j]
            if results[i][j] == "~" and results[ref][j] == "~":
                sdiff.write("    ~  ")
            else:
                sdiff.write(f"{_res:+5d}  ")
        sdiff.write("\n")
    delspec.close()
    sdiff.close()


def dssp(trajectory, structure, index=None, interactive=False, maxprocs=90, **kwargs):
    """
    Takes a trajectory file, or a list of trajectories, alongside their structure files (tpr)
    to run the dssp analysis in GROMACS.

    The script can be run in interactive mode, or it can be made to output a SLURM script to
    submit to a queue (TODO).

    By default, maxprocs processes are started at the same time.
    """
    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure,  list) else [structure]
    __indices      = index      if isinstance(index,      list) else [index]
    #
    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))
    #
    ntraj = len(__trajectories)
    #
    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])
    #
    hmode = kwargs.pop("hmode", "gromacs")
    hbond = kwargs.pop("hbond", "energy")
    ref   = kwargs.pop("ref", 0)
    doss  = kwargs.pop("doss", "all")
            
    # All checks have passed
    if interactive:
        import psutil
        import subprocess
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical/nprocs)
        __nprocs  = min(maxprocs, nprocs)
    
        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")
        #
        datfiles = []
        # Run DSSP on all trajectories
        commands = "source /anfhome/.profile "
        commands += "&& module  load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx " 
        commands += f"&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && gmx dssp -hmode {hmode} -hbond {hbond} -f " 
        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids         = []
            for j in range(__nprocs):
                if i+j == ntraj: break
                _commands = commands
                _commands += f"{__trajectories[i+j]} -s {__structures[i+j]} -o {i+j:04d}_dssp.dat"
                if __indices[i+j] is not None: 
                    _commands += "-n {__indices[i+j]}"
                datfiles.append(f"{i+j:04d}_dssp.dat")
                subprocesses.append(subprocess.Popen(_commands, shell=True))
                pids.append(subprocesses[-1].pid)
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j*nskip])
                time.sleep(2)
            exit_codes = [p.wait() for p in subprocesses]
        #
        # Process the Data
        #
        process_dssp(datfiles, ref=ref, doss=doss)
    else:
        """
        TODO
        """
        pass
 