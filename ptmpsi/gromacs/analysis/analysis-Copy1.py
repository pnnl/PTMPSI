import errno
import os
from collections import defaultdict
import time
import shutil  
import psutil
import subprocess
import numpy as np

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
    for i in range(0, ntraj):
        #
        # Check that current set of files exists
        trajectory, structure, index = __trajectories[i], __structures[i], __indices[i] 
        check_file_exists(trajectory)
        check_file_exists(structure)
        check_file_exists(index, canbenone=True)
        #
        # All checks have passed
        # Make sure the directory exists
        if prefix:
            os.makedirs(prefix, exist_ok=True)

        # Now get the output filenames
        base_output_name = os.path.splitext(os.path.basename(trajectory))[0]
        output = os.path.join(prefix, base_output_name) if prefix else base_output_name
        if suffix:
            output += f"_{suffix}"
        outputs.append(output)

    #
    # Run interactively
    if interactive:
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
                __commands += f"echo 3 0 | gmx trjconv -f {output}_whole_nojump.xtc -o {output}_whole_nojump_fit.xtc -fit rot+trans -s {output}.gro && "
                
                # Generate new TPR file with desired group
                __commands += f"echo '{group}' | gmx convert-tpr -s {structure} -o {output}_{group}.tpr "
                __commands += f"-n {index}" if index is not None else ""
                #
                # Run subprocess and get current process id (PID)
                subprocesses.append(subprocess.Popen(__commands, shell=True))
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

def process_dssp(datfile, trajectory, ref=0, doss="all"):
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
    for i in range(0, ntraj):
        nchars = len(results[i])
        _numres.append( [ 1 if results[i][j] in listss else 0 for j in range(nchars) ] )
    #   
    delspec.write("Case   totSS       SS+       SS-       Delta\n")
    sdiff.write(" Case   ")
    for i in range(0, ntraj):
        sdiff.write(f"{i:04d}   ")
    sdiff.write("\n")
    #
    for i in range(0, ntraj):
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
    
def construct_residue_path(file_path):
    """
    Remove 'postprocess' and the last two parts of the path, and append 'residuetypes.dat'.
    """
    # Split the path into parts
    path_parts = file_path.split('/')
    
    # Remove 'postprocess' if it exists
    if 'postprocess' in path_parts:
        path_parts.remove('postprocess')
    
    # Remove the last two parts (file name and the preceding directory)
    base_path = '/'.join(path_parts[:-2])
    
    # Append 'residuetypes.dat' to the remaining path
    residue_path = os.path.join(base_path, 'residuetypes.dat')
    
    return residue_path

def add_residue_file(file_path):
    """
    Add 'residuetypes.dat' to the working directory if it doesn't exist.
    """
    # Construct the path where 'residuetypes.dat' should be
    residue_path = construct_residue_path(file_path)
    print("residue path is:", residue_path)

    # Check if the file already exists
    working_dir = os.getcwd()
    parent_dir = os.path.dirname(working_dir) 
    destination_path = os.path.join(parent_dir, 'residuetypes.dat')
    if not os.path.exists(destination_path):
        # If not, copy it from the source path to working directory
        print(residue_path)
        shutil.copy(residue_path, working_dir)
        print("copied")
    else:
        print(f"'{residue_path}' already exists. Nothing to do.")

def extract_nested_directory(file_path):
    """
    Extract the second-to-last directory part from the given file path.
    Ensure it's a four-digit integer, including leading zeros.
    """
    # Split the path into parts using os.path.split
    dir_path, file_name = os.path.split(file_path)
    parent_dir, nested_dir = os.path.split(dir_path)
    
    try:
        # Convert the extracted part to an integer
        nested_dir_int = int(nested_dir)
        # Ensure the integer is four digits, including leading zeros
        if len(nested_dir) != 4 or not nested_dir.isdigit():
            raise ValueError
    except ValueError:
        raise ValueError(f"The extracted directory part '{nested_dir}' is not a four-digit integer.")
    
    return nested_dir_int

def dssp(trajectory, structure, output_dir, index=None, interactive=False, maxprocs=90, **kwargs):
    """
    Takes a trajectory file, or a list of trajectories, alongside their structure files (tpr)
    to run the dssp analysis in GROMACS.

    The script can be run in interactive mode, or it can be made to output a SLURM script to
    submit to a queue (TODO).

    By default, maxprocs processes are started at the same time.
    """
    #add_residue_file(trajectory)
    print("1) working directory is", os.getcwd()) 
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")
        
    original_directory= os.getcwd()
    os.chdir(output_dir)
    
    
    add_residue_file(trajectory)
    print("2) working directory is", os.getcwd()) 
    
    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    print(__trajectories)
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
                basename = os.path.basename(os.path.dirname(__trajectories[i + j]))
                _commands = commands
                _commands += f"{__trajectories[i+j]} -s {__structures[i+j]} -o {basename}_dssp.dat"
                if __indices[i+j] is not None: 
                    _commands += "-n {__indices[i+j]}"
                datfiles.append(f"{basename}_dssp.dat")
                subprocesses.append(subprocess.Popen(_commands, shell=True, cwd=output_dir))
                pids.append(subprocesses[-1].pid)
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j*nskip])
                time.sleep(2)
            exit_codes = [p.wait() for p in subprocesses]
        #
        # Process the Data
        #
        process_dssp(datfiles, trajectory, ref=ref, doss=doss)
        os.chdir(original_directory)
    else:
        """
        TODO
        """
        pass

def check_file_exists(filename, canbenone=False):
    if filename is None:
        if canbenone:
            return
        else:
            raise KeyError("filename cannot be None")
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    return

def process_rmsd(result_files, reference_file):
    ref_rmsd = []
    if not os.path.exists(reference_file):
        raise FileNotFoundError(f"Reference file {reference_file} not found.")
    
    with open(reference_file, "r") as ref:
        for line in ref:
            if not line.startswith("@") and not line.startswith("#"):
                ref_rmsd.append(float(line.split()[1]))

    delspec = open("DELSPEC_rmsd.out", "w")
    rms_diff = open("RESULT_rmsDiffRES.out", "w")

    delspec.write("Case   RMSD Total     Delta\n")
    rms_diff.write(" Case   ")

    for i, result_file in enumerate(result_files):
        if not os.path.exists(result_file):
            print(f"Warning: Result file {result_file} not found. Skipping.")
            continue

        rmsd_values = []
        with open(result_file, "r") as rf:
            for line in rf:
                if not line.startswith("@") and not line.startswith("#"):
                    rmsd_values.append(float(line.split()[1]))

        total_rmsd = sum(rmsd_values)
        ref_total_rmsd = sum(ref_rmsd)
        delta = total_rmsd - ref_total_rmsd

        delspec.write(f"{i:04d}   {total_rmsd:10.4f}     {delta:+10.4f}\n")

        rms_diff.write(f"{i:04d} ")
        for val in rmsd_values:
            rms_diff.write(f"{val:8.4f} ")
        rms_diff.write("\n")

    delspec.close()
    rms_diff.close()

def rmsd(trajectory, structure, output_dir, group=4, index=None, interactive=False ,maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    ref = kwargs.pop("ref", 0)

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                _commands = commands
                _commands += f"echo {group} 4 | gmx rms -s {__structures[i + j]} -f {__trajectories[i + j]} -o {i + j:04d}_rmsd.xvg"
                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"
                output_file = f"{i + j:04d}_rmsd.xvg"
                print(f"Executing command: {_commands}")
                print(f"Output file will be saved at: {os.path.join(output_dir, output_file)}")
                result_files.append(output_file)
                subprocesses.append(subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
                pids.append(subprocesses[-1].pid)
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j * nskip])
                time.sleep(2)
            exit_codes = [p.wait() for p in subprocesses]
            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    print(f"Command failed for file {result_files[exit_codes.index(ex)]}, checking logs...")
                    stdout, stderr = subprocesses[exit_codes.index(ex)].communicate()
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")

        process_rmsd(result_files, result_files[ref])
        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())
    else:
        """
        TODO
        """
        pass

def check_file_exists(filename, canbenone=False):
    if filename is None:
        if canbenone:
            return
        else:
            raise KeyError("filename cannot be None")
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    return

def process_rmsd(result_files, reference_file):
    ref_rmsd = []
    if not os.path.exists(reference_file):
        raise FileNotFoundError(f"Reference file {reference_file} not found.")
    
    with open(reference_file, "r") as ref:
        for line in ref:
            if not line.startswith("@") and not line.startswith("#"):
                ref_rmsd.append(float(line.split()[1]))

    delspec = open("DELSPEC_rmsd.out", "w")
    rms_diff = open("RESULT_rmsDiffRES.out", "w")

    delspec.write("Case   RMSD Total     Delta\n")
    rms_diff.write(" Case   ")

    for i, result_file in enumerate(result_files):
        if not os.path.exists(result_file):
            print(f"Warning: Result file {result_file} not found. Skipping.")
            continue

        rmsd_values = []
        with open(result_file, "r") as rf:
            for line in rf:
                if not line.startswith("@") and not line.startswith("#"):
                    rmsd_values.append(float(line.split()[1]))

        total_rmsd = sum(rmsd_values)
        ref_total_rmsd = sum(ref_rmsd)
        delta = total_rmsd - ref_total_rmsd

        delspec.write(f"{i:04d}   {total_rmsd:10.4f}     {delta:+10.4f}\n")

        rms_diff.write(f"{i:04d} ")
        for val in rmsd_values:
            rms_diff.write(f"{val:8.4f} ")
        rms_diff.write("\n")

    delspec.close()
    rms_diff.close()

def rmsd(trajectory, structure, output_dir, index=None, group=4, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    ref = kwargs.pop("ref", 0)

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                basename = os.path.basename(os.path.dirname(__trajectories[i + j]))
                print(basename)
                _commands = commands
                _commands += f"echo {group} 4 | gmx rms -s {__structures[i + j]} -f {__trajectories[i + j]} -o {basename}_rmsd.xvg"

                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"
                output_file = f"{basename}_rmsd.xvg"
                print(f"Executing command: {_commands}")
                print(f"Output file will be saved at: {os.path.join(output_dir, output_file)}")
                result_files.append(output_file)
                subprocesses.append(subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
                pids.append(subprocesses[-1].pid)
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j * nskip])
                time.sleep(2)
            exit_codes = [p.wait() for p in subprocesses]
            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    print(f"Command failed for file {result_files[exit_codes.index(ex)]}, checking logs...")
                    stdout, stderr = subprocesses[exit_codes.index(ex)].communicate()
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")

        process_rmsd(result_files, result_files[ref])
        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())
    else:
        """
        TODO
        """
        pass

def process_rmsf(result_files):
    delspec = open("DELSPEC_rmsf.out", "w")
    rms_diff = open("RESULT_rmsfDiffRES.out", "w")

    delspec.write("Case   RMSF Total\n")
    rms_diff.write(" Case   ")

    for i, (rmsf_file, bfactors_file) in enumerate(result_files):
        if not os.path.exists(rmsf_file):
            print(f"Warning: RMSF file {rmsf_file} not found. Skipping.")
            continue

        rmsf_values = []
        with open(rmsf_file, "r") as rf:
            for line in rf:
                if not line.startswith("@") and not line.startswith("#"):
                    rmsf_values.append(float(line.split()[1]))

        total_rmsf = sum(rmsf_values)

        delspec.write(f"{i:04d}   {total_rmsf:10.4f}\n")

        rms_diff.write(f"{i:04d} ")
        for val in rmsf_values:
            rms_diff.write(f"{val:8.4f} ")
        rms_diff.write("\n")

    delspec.close()
    rms_diff.close()

def rmsf(trajectory, structure, output_dir, index=None, group=4, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                basename = os.path.basename(os.path.dirname(__trajectories[i + j]))
                _commands = commands
                _commands += f"echo {group} 4 | gmx rmsf -f {__trajectories[i + j]} -s {__structures[i + j]} -o {basename}_rmsf.xvg -res -oq {basename}_bfactors.pdb -fit"
                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"
                rmsf_file = f"{basename}_rmsf.xvg"
                bfactors_file = f"{basename}_bfactors.pdb"
                print(f"Executing command: {_commands}")
                print(f"Output files will be saved at: {os.path.join(output_dir, rmsf_file)} and {os.path.join(output_dir, bfactors_file)}")
                result_files.append((rmsf_file, bfactors_file))
                subprocesses.append(subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
                pids.append(subprocesses[-1].pid)
                p = psutil.Process(pids[-1])
                p.cpu_affinity([j * nskip])
                time.sleep(2)
            exit_codes = [p.wait() for p in subprocesses]
            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    print(f"Command failed for files {result_files[exit_codes.index(ex)]}, checking logs...")
                    stdout, stderr = subprocesses[exit_codes.index(ex)].communicate()
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")

        process_rmsf(result_files)
        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())
    else:
        """
        TODO
        """
        pass

def process_gyration(result_files):
    delspec = open("DELSPEC_gyration.out", "w")
    gyr_diff = open("RESULT_gyrationDiffRES.out", "w")

    delspec.write("Case   Gyration Total\n")
    gyr_diff.write(" Case   ")

    for i, gyr_file in enumerate(result_files):
        if not os.path.exists(gyr_file):
            print(f"Warning: Gyration file {gyr_file} not found. Skipping.")
            continue

        gyr_values = []
        with open(gyr_file, "r") as gf:
            for line in gf:
                if not line.startswith("@") and not line.startswith("#"):
                    gyr_values.append(float(line.split()[1]))

        total_gyr = sum(gyr_values)

        delspec.write(f"{i:04d}   {total_gyr:10.4f}\n")

        gyr_diff.write(f"{i:04d} ")
        for val in gyr_values:
            gyr_diff.write(f"{val:8.4f} ")
        gyr_diff.write("\n")

    delspec.close()
    gyr_diff.close()

def gyration(trajectory, structure, output_dir, index=None, group=4, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                
                # Extract the directory name based on the trajectory
                basename_index = os.path.basename(os.path.dirname(__trajectories[i + j]))

                _commands = commands
                _commands += f"echo {group} | gmx gyrate -f {__trajectories[i + j]} -s {__structures[i + j]} -o {basename_index}_radgyr.xvg"
                
                
                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"
                
                gyr_file = f"{basename_index}_radgyr.xvg"

                print(f"Executing command: {_commands}")
                print(f"Output file will be saved at: {os.path.join(output_dir, gyr_file)}")

                result_files.append(gyr_file)

                process = subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocesses.append(process)
                pids.append(process.pid)

                # Set CPU affinity to avoid oversubscription
                p = psutil.Process(process.pid)
                p.cpu_affinity([j * nskip])

                # Sleep for 2 seconds to stagger commands
                time.sleep(2)
            
            # Wait for the current batch of subprocesses to finish
            exit_codes = [p.wait() for p in subprocesses]

            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    idx = exit_codes.index(ex)
                    process = subprocesses[idx]
                    stdout, stderr = process.communicate()
                    print(f"Command failed for file {result_files[idx]}, checking logs...")
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")
        
        process_gyration(result_files)
        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())
    else:
        """
        TODO
        """
        pass
        
def process_gyration(result_files):
    delspec = open("DELSPEC_gyration.out", "w")
    gyr_diff = open("RESULT_gyrationDiffRES.out", "w")

    delspec.write("Case   Gyration Total\n")
    gyr_diff.write(" Case   ")

    for i, gyr_file in enumerate(result_files):
        if not os.path.exists(gyr_file):
            print(f"Warning: Gyration file {gyr_file} not found. Skipping.")
            continue

        gyr_values = []
        with open(gyr_file, "r") as gf:
            for line in gf:
                if not line.startswith("@") and not line.startswith("#"):
                    gyr_values.append(float(line.split()[1]))

        total_gyr = sum(gyr_values)

        delspec.write(f"{i:04d}   {total_gyr:10.4f}\n")

        gyr_diff.write(f"{i:04d} ")
        for val in gyr_values:
            gyr_diff.write(f"{val:8.4f} ")
        gyr_diff.write("\n")

    delspec.close()
    gyr_diff.close()

def gyration(trajectory, structure, output_dir, index=None, group=4, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    assert(len(__trajectories) == len(__structures))
    assert(len(__trajectories) == len(__indices))

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                
                # Extract the directory name based on the trajectory
                basename_index = os.path.basename(os.path.dirname(__trajectories[i + j]))

                _commands = commands
                _commands += f"echo {group} | gmx gyrate -f {__trajectories[i + j]} -s {__structures[i + j]} -o {basename_index}_radgyr.xvg"
                
                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"
                
                gyr_file = f"{basename_index}_radgyr.xvg"

                print(f"Executing command: {_commands}")
                print(f"Output file will be saved at: {os.path.join(output_dir, gyr_file)}")

                result_files.append(gyr_file)

                process = subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocesses.append(process)
                pids.append(process.pid)

                # Set CPU affinity to avoid oversubscription
                p = psutil.Process(process.pid)
                p.cpu_affinity([j * nskip])

                # Sleep for 2 seconds to stagger commands
                time.sleep(2)
            
            # Wait for the current batch of subprocesses to finish
            exit_codes = [p.wait() for p in subprocesses]

            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    idx = exit_codes.index(ex)
                    process = subprocesses[idx]
                    stdout, stderr = process.communicate()
                    print(f"Command failed for file {result_files[idx]}, checking logs...")
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")
        
        process_gyration(result_files)
        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())

def contact_map_old(trajectory, structure, md_file, output_dir, index=None, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                
                basename_index = os.path.basename(os.path.dirname(__trajectories[i + j]))

                _commands = commands
                _commands += f"gmx pairdist -f {__trajectories[i + j]} -s {__structures[i + j]} -refgrouping res -selgrouping res -ref 0 -sel 0 -o {basename_index}_dist.xvg -cutoff 1.0 -dt 1000.0"
                
                if __indices[i + j] is not None:
                    _commands += f" -n {__indices[i + j]}"

                contact_file = f"{basename_index}_dist.xvg"

                print(f"Executing command: {_commands}")
                print(f"Output file will be saved at: {os.path.join(output_dir, contact_file)}")

                result_files.append(contact_file)

                process = subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocesses.append(process)
                pids.append(process.pid)

                p = psutil.Process(process.pid)
                p.cpu_affinity([j * nskip])

                time.sleep(2)
            
            exit_codes = [p.wait() for p in subprocesses]

            for ex in exit_codes:
                print(f"Process exited with code: {ex}")
                if ex != 0:
                    idx = exit_codes.index(ex)
                    process = subprocesses[idx]
                    stdout, stderr = process.communicate()
                    print(f"Command failed for file {result_files[idx]}, checking logs...")
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")

        os.chdir(original_directory)
        print("3) Changed back to original directory:", os.getcwd())
        
def get_number_of_residues(gro_filepath):
    """Extract the number of residues from the .gro file."""
    with open(gro_filepath, 'r', encoding='latin-1') as f:
        lines = f.readlines()
        if len(lines) < 2:
            raise ValueError(f"File {gro_filepath} doesn't contain enough lines to extract the number of residues.")
        # Extract the first integer from the second-to-last line
        second_to_last_line = lines[-2].strip().split()
        number_of_residues = int(''.join(filter(str.isdigit, second_to_last_line[0])))

        return number_of_residues

def get_frame_count(trajectory_file):
    """Get the total number of frames and print it
    some issues with this"""
    command = f"source /anfhome/.profile && source /anfhome/shared/qipd/gromacs/bin/GMXRC && gmx check -f {trajectory_file}"
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"Failed to check frame count for {trajectory_file}, error: {stderr.decode()}")
    for line in stdout.decode().split('\n'):
        if "Read n frames" in line:
            frame_count = int(line.split()[2])
            print(f"Total number of frames in {trajectory_file}: {frame_count}")
            return frame_count
    raise RuntimeError("Failed to parse frame count from gmx check output")

def contact_map(trajectory, structure, md_file, output_dir, index=None, interactive=False, maxprocs=90, **kwargs):
    print("1) Initial working directory is:", os.getcwd())
    if not os.path.exists(output_dir):
        print(f"Directory {output_dir} does not exist. Creating it now.")
        try:
            os.makedirs(output_dir)
            print(f"Directory {output_dir} created successfully.")
        except OSError as e:
            print(f"Failed to create {output_dir}. Error: {e}")
            raise
    else:
        print(f"Directory {output_dir} already exists.")

    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures   = structure  if isinstance(structure, list) else [structure]
    __indices      = index      if isinstance(index, list) else [index]

    ntraj = len(__trajectories)

    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    if interactive:
        nprocs    = psutil.cpu_count(logical=False)
        nlogical  = psutil.cpu_count()
        nskip     = int(nlogical / nprocs)
        __nprocs  = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        result_files = []
        commands = (
            "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
            "&& source /anfhome/shared/qipd/gromacs/bin/GMXRC && "
        )
        # Get number of residues from the postprocessed gro file
        nres = get_number_of_residues(md_file)
        print("number of residues is: ", nres)

        for i in range(0, ntraj, __nprocs):
            subprocesses = []
            pids = []
            matrix = np.zeros((nres, nres))
            num_frames = 0

            for j in range(__nprocs):
                if i + j == ntraj:
                    break
                
                basename_index = os.path.basename(os.path.dirname(__trajectories[i + j]))
                total_frames=3000
                while num_frames < total_frames:
                    print("num frames: ",num_frames)
                    begin, end = num_frames, num_frames+1000 #removed +1
                    stderr_file = os.path.join(output_dir, f"stderr_{basename_index}_{num_frames}.log")
                    stdout_file = os.path.join(output_dir, f"stdout_{basename_index}_{num_frames}.log")

                    _commands = commands
                    print(begin, end)
                    _commands += f"gmx pairdist -f {__trajectories[i + j]} -s {__structures[i + j]} -b {begin} -e {end} -refgrouping res -selgrouping res -ref 0 -sel 0 -o {basename_index}_dist_{num_frames}.xvg -cutoff 1.0 -dt 500"
                    
                    if __indices[i + j] is not None:
                        _commands += f" -n {__indices[i + j]}"

                    contact_file = f"{basename_index}_dist_{num_frames}.xvg"

                    print(f"Executing command: {_commands}")
                    print(f"Output file will be saved at: {os.path.join(output_dir, contact_file)}")

                    result_files.append(contact_file)
                    process = subprocess.Popen(_commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=output_dir)
                    subprocesses.append(process)
                    pids.append(process.pid)
                    
                    # Set CPU affinity
                    p = psutil.Process(pids[-1])
                    p.cpu_affinity([j * nskip])
                    time.sleep(2)
                    exit_codes = [p.wait() for p in subprocesses]
                    num_frames+=1000
                    
                    for ex in exit_codes:
                        print(f"Process exited with code: {ex}")
                        if ex != 0:
                            print(f"Command failed for file {result_files[exit_codes.index(ex)]}, checking logs...")
                            stdout, stderr = subprocesses[exit_codes.index(ex)].communicate()
                            print(f"stdout: {stdout.decode()}")
                            print(f"stderr: {stderr.decode()}")
        
                            # Process generated files
                            for process in subprocesses:
                                stdout, stderr = process.communicate()
                                if process.returncode == 0:
                                        basename_index = os.path.basename(os.path.dirname(__trajectories[i]))
                                        contact_file = f"{basename_index}_dist_{num_frames - 1}.xvg"
                        
    process_contact(result_files, nres)
    os.chdir(original_directory)
    print("3) Changed back to original directory:", os.getcwd())


def process_contact(contact_files):
    # Create a list to store results for each contact file
    results = []
    print(contact_files)
    print("processing contact files")
    
    # Process each contact file
    for contact_file in contact_files:
        with open(contact_file, 'r', encoding='latin-1') as file:
            # Skip the first 24 lines; adjust this if necessary based on comment lines
            lines = file.readlines()[24:] 

            # Assuming all lines have the same number of entries based on the first line after skipping
            first_line = lines[0]
            parts = first_line.split()

            # Calculate nres as the square root of the number of entries and rounding down
            nres = int(np.floor(np.sqrt(len(parts))))

            # Create a zero matrix of dimensions nres x nres
            contact_matrix = np.zeros((nres, nres))

            # Accumulate matrix information from each line
            for line in lines:
                parts = [float(entry) for entry in line.split()]
                
                if len(parts) < nres * nres:
                    # Adjust nres in cases of non-square number of parts
                    print(f"Skipping line due to insufficient data: {line}")
                    continue

                index = 0
                for i in range(nres):
                    for j in range(nres):
                        if index < len(parts) and parts[index] > 0.5:
                            contact_matrix[i, j] += 1
                        index += 1

            # Calculate the sum of entries in the contact matrix
            matrix_sum = np.sum(contact_matrix)
            
            # Calculate the average sum by nres*nres
            average_contact = matrix_sum / (nres * nres)
            
            # Append the result for the current file
            results.append(average_contact)

    # Return the list of results for each file
    return results



#for bash script
if __name__ == "__main__":
    # Example usage
    trajectory = "/path/to/your/trajectory/whole_nojump_fit.xtc"
    structure = "/path/to/your/structure/protein.tpr"
    output_dir = "/path/to/your/output/directory"
    contact_map([trajectory], [structure], output_dir)
