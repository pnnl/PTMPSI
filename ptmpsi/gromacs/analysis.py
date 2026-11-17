"""GROMACS trajectory post-processing for PTM-Psi runs.

Each analysis calls the matching gmx tool (trjconv, rms, rmsf, gyrate, dssp, pairdist) on one or
more trajectories and then summarizes the .xvg/.dat files in the output directory. The gmx
program is gmx_mpi on cluster modules (Deception) and gmx elsewhere.

PTM-Psi renames modified residues (SNC, CGL, CSO, ...). gmx only treats them as Protein when
PTM-Psi's residuetypes.dat is in the working directory, so every function stages it there first.
"""

import contextlib
import errno
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np

AQE_GMXRC = "/anfhome/shared/qipd/gromacs/bin/GMXRC"
AQE_PREAMBLE = (
    "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx "
    f"&& source {AQE_GMXRC}"
)

_PACKAGED_RESIDUETYPES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "forcefield", "residuetypes.dat")

def check_file_exists(filename, canbenone=False):
    if filename is None:
        if canbenone:
            return
        else:
            raise KeyError("filename cannot be None")
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    return

def totss(result):
    return result.count("H")+ result.count("E") + result.count("G") + result.count("I") + result.count("P") + result.count("B")

def totab(result):
    return result.count("H") + result.count("E")


def find_gmx(gmx=None):
    """Return the gmx executable: gmx if present (as before), else gmx_mpi (the Deception module)."""
    if gmx:
        found = shutil.which(gmx)
        if found is None:
            raise FileNotFoundError(f"GROMACS executable '{gmx}' was not found on PATH")
        return found
    for name in ("gmx", "gmx_mpi"):
        found = shutil.which(name)
        if found:
            return found
    if os.path.isfile(AQE_GMXRC):
        return "gmx"   # the AQE cluster only has gmx after sourcing GMXRC; _run does that
    raise FileNotFoundError("Neither 'gmx' nor 'gmx_mpi' is on PATH. Load the GROMACS module or install GROMACS")


def _available_cpus():
    try:
        return len(os.sched_getaffinity(0))   # the cores this Slurm job was given
    except AttributeError:
        return os.cpu_count() or 1


def _run(command, cwd, stdin=None):
    if command[0] == "gmx" and shutil.which("gmx") is None and os.path.isfile(AQE_GMXRC):
        command = ["bash", "-c", AQE_PREAMBLE + ' && exec "$@"', "bash", *command]
    result = subprocess.run(command, cwd=cwd, input=stdin or "", text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed with exit code {result.returncode}:\n{result.stderr[-2000:]}")
    return result


def _run_all(tasks, maxprocs):
    """Run callables concurrently. A failure is printed and the others carry on, as before."""
    if not tasks:
        return []
    failures = []
    with ThreadPoolExecutor(max_workers=max(1, min(maxprocs, len(tasks), _available_cpus()))) as pool:
        for future in [pool.submit(task) for task in tasks]:
            try:
                future.result()
            except Exception as error:
                failures.append(str(error))
                print(f"Command failed: {error}")
    return failures


@contextlib.contextmanager
def _working_dir(path):
    os.makedirs(path, exist_ok=True)
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _not_interactive(name):
    print(f"{name}(interactive=False) does nothing: the Slurm-script branch was never written. Pass interactive=True.")


def _inputs(trajectory, structure, index):
    trajectories = list(trajectory) if isinstance(trajectory, (list, tuple)) else [trajectory]
    structures = list(structure) if isinstance(structure, (list, tuple)) else [structure]
    indices = list(index) if isinstance(index, (list, tuple)) else [index] * len(trajectories)
    if not (len(trajectories) == len(structures) == len(indices)):
        raise ValueError("trajectory, structure and index must have the same number of entries")
    for t, s, n in zip(trajectories, structures, indices):
        check_file_exists(t)
        check_file_exists(s)
        check_file_exists(n, canbenone=True)
    absolute = lambda p: None if p is None else os.path.abspath(p)
    return [absolute(p) for p in trajectories], [absolute(p) for p in structures], [absolute(p) for p in indices]


def _labels(trajectories):
    """Name each output after the folder holding its trajectory (0001, 0002, ...)."""
    labels = [os.path.basename(os.path.dirname(t)) for t in trajectories]
    if len(set(labels)) != len(labels):
        labels = [f"{i:04d}" for i in range(len(trajectories))]
    return labels


def _index_args(index):
    return ["-n", index] if index is not None else []


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
        basename = __datfiles[i][:-len("_dssp.dat")] if __datfiles[i].endswith("_dssp.dat") else __datfiles[i]
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

def add_residue_file(file_path=None, destination=None):
    """Copy residuetypes.dat into `destination` (default: the working directory) if it is not there.

    The copy next to the run (see construct_residue_path) wins, because gen_ptm_files put it there;
    the copy shipped with PTM-Psi is the fallback.
    """
    target = os.path.join(destination or os.getcwd(), "residuetypes.dat")
    if os.path.exists(target):
        return target
    candidates = []
    if file_path is not None:
        candidates.append(construct_residue_path(os.path.abspath(file_path)))
    candidates.append(_PACKAGED_RESIDUETYPES)
    for candidate in candidates:
        if os.path.isfile(candidate):
            shutil.copy(candidate, target)
            return target
    raise FileNotFoundError(errno.ENOENT, "residuetypes.dat was not found", ", ".join(candidates))


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


def post_process(trajectories, structures, indices=None, prefix=None, suffix=None, interactive=False, maxprocs=90,
                 group="Non-Water", **kwargs):
    """Make `group` whole, remove periodic jumps and fit it, for each trajectory.

    Writes <output>.gro (frame 0), <output>_whole.xtc, <output>_whole_nojump.xtc,
    <output>_whole_nojump_fit.xtc and <output>_<group>.tpr, where <output> is the trajectory name
    inside `prefix` (plus _<suffix>). Returns the fitted trajectories. `fit_group` is the gmx group
    the fit uses (default 3, C-alpha). Raises RuntimeError if any command fails or an
    expected output is missing.
    """
    if not interactive:
        return _not_interactive("post_process")
    fit_group = kwargs.pop("fit_group", 3)
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectories, structures, indices)
    workdir = os.path.abspath(prefix) if prefix else os.getcwd()
    os.makedirs(workdir, exist_ok=True)
    outputs = []
    for trajectory in trajectories:
        output = os.path.join(workdir, os.path.splitext(os.path.basename(trajectory))[0])
        outputs.append(output + (f"_{suffix}" if suffix else ""))
    add_residue_file(trajectories[0], destination=workdir)

    def one(trajectory, structure, index, output):
        n = _index_args(index)
        commands = [
            ([gmx, "trjconv", "-dump", "0", "-f", trajectory, "-s", structure, "-o", f"{output}.gro", *n], f"{group}\n"),
            ([gmx, "trjconv", "-f", trajectory, "-s", structure, "-o", f"{output}_whole.xtc", "-pbc", "whole", *n], f"{group}\n"),
            ([gmx, "trjconv", "-f", f"{output}_whole.xtc", "-o", f"{output}_whole_nojump.xtc", "-pbc", "nojump"], "0\n"),
            ([gmx, "trjconv", "-f", f"{output}_whole_nojump.xtc", "-o", f"{output}_whole_nojump_fit.xtc",
              "-fit", "rot+trans", "-s", f"{output}.gro"], f"{fit_group}\n0\n"),
            ([gmx, "convert-tpr", "-s", structure, "-o", f"{output}_{group}.tpr", *n], f"{group}\n"),
        ]
        for command, stdin in commands:
            _run(command, workdir, stdin)
            expected = command[command.index("-o") + 1]
            if not os.path.isfile(expected):
                raise RuntimeError(f"GROMACS did not create expected output: {expected}")

    failures = _run_all(
        [partial(one, t, s, n, o) for t, s, n, o in zip(trajectories, structures, indices, outputs)], maxprocs
    )
    if failures:
        raise RuntimeError("Trajectory post-processing failed:\n" + "\n".join(failures))
    return [f"{output}_whole_nojump_fit.xtc" for output in outputs]


def dssp(trajectory, structure, output_dir, index=None, interactive=False, maxprocs=90, **kwargs):
    """Run gmx dssp on each trajectory, then summarize with process_dssp.

    Outputs <folder>_dssp.dat in `output_dir`, where <folder> is the folder holding the trajectory, and
    returns their paths.
    `hmode`, `hbond`, `ref` and `doss` are passed on as before.
    """
    if not interactive:
        return _not_interactive("dssp")
    hmode = kwargs.pop("hmode", "gromacs")
    hbond = kwargs.pop("hbond", "energy")
    ref = kwargs.pop("ref", 0)
    doss = kwargs.pop("doss", "all")
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectory, structure, index)
    labels = _labels(trajectories)
    output_dir = os.path.abspath(output_dir)
    with _working_dir(output_dir):
        add_residue_file(trajectories[0])
        datfiles = [f"{label}_dssp.dat" for label in labels]
        _run_all([partial(_run, [gmx, "dssp", "-hmode", hmode, "-hbond", hbond, "-f", t, "-s", s, "-o", out, *_index_args(n)], output_dir)
                  for t, s, n, out in zip(trajectories, structures, indices, datfiles)], maxprocs)
        process_dssp(datfiles, trajectories, ref=ref, doss=doss)
    return [os.path.join(output_dir, f) for f in datfiles]


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
    """gmx rms of `group` (default 4, Backbone) against the structure, fitted on `group`.

    Writes <folder>_rmsd.xvg per trajectory plus DELSPEC_rmsd.out and RESULT_rmsDiffRES.out, which
    compare each case with trajectory number `ref` (default 0).
    """
    if not interactive:
        return _not_interactive("rmsd")
    ref = kwargs.pop("ref", 0)
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectory, structure, index)
    labels = _labels(trajectories)
    output_dir = os.path.abspath(output_dir)
    with _working_dir(output_dir):
        add_residue_file(trajectories[0])
        result_files = [f"{label}_rmsd.xvg" for label in labels]
        _run_all([partial(_run, [gmx, "rms", "-s", s, "-f", t, "-o", out, *_index_args(n)], output_dir, f"{group} 4\n")
                  for t, s, n, out in zip(trajectories, structures, indices, result_files)], maxprocs)
        process_rmsd(result_files, result_files[ref])
    return [os.path.join(output_dir, f) for f in result_files]


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
    """gmx rmsf per residue of `group` (default 4, Backbone), fitted on `group`.

    Writes <folder>_rmsf.xvg and <folder>_bfactors.pdb per trajectory plus DELSPEC_rmsf.out and
    RESULT_rmsfDiffRES.out.
    """
    if not interactive:
        return _not_interactive("rmsf")
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectory, structure, index)
    labels = _labels(trajectories)
    output_dir = os.path.abspath(output_dir)
    with _working_dir(output_dir):
        add_residue_file(trajectories[0])
        result_files = [(f"{label}_rmsf.xvg", f"{label}_bfactors.pdb") for label in labels]
        _run_all([partial(_run, [gmx, "rmsf", "-f", t, "-s", s, "-o", xvg, "-res", "-oq", pdb, "-fit", *_index_args(n)],
                          output_dir, f"{group} 4\n")
                  for t, s, n, (xvg, pdb) in zip(trajectories, structures, indices, result_files)], maxprocs)
        process_rmsf(result_files)
    return [(os.path.join(output_dir, xvg), os.path.join(output_dir, pdb)) for xvg, pdb in result_files]


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
    """gmx gyrate (mass-weighted radius of gyration) of `group` (default 4, Backbone).

    Writes <folder>_radgyr.xvg per trajectory plus DELSPEC_gyration.out and RESULT_gyrationDiffRES.out.
    Pass group=1 for the whole protein.
    """
    if not interactive:
        return _not_interactive("gyration")
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectory, structure, index)
    labels = _labels(trajectories)
    output_dir = os.path.abspath(output_dir)
    with _working_dir(output_dir):
        add_residue_file(trajectories[0])
        result_files = [f"{label}_radgyr.xvg" for label in labels]
        _run_all([partial(_run, [gmx, "gyrate", "-f", t, "-s", s, "-o", out, *_index_args(n)], output_dir, f"{group}\n")
                  for t, s, n, out in zip(trajectories, structures, indices, result_files)], maxprocs)
        process_gyration(result_files)
    return [os.path.join(output_dir, f) for f in result_files]


def process_contact(contact_files, contact_cutoff=0.5):
    """Average the contact fractions of gmx pairdist files and plot them, as before.

    Each file holds a time column and then one minimum distance in nm per residue pair. A pair counts as
    in contact in a frame when it is closer than `contact_cutoff` nm; each file gives the fraction of its
    frames, and the files are averaged. Writes contact_fraction.dat and .png next to the first file, shows
    the plot, and returns the averaged matrix.
    """
    import matplotlib.pyplot as plt

    matrices = []
    print("Processing contact files")
    for contact_file in contact_files:
        rows = []
        with open(contact_file, "r", encoding="latin-1") as fh:
            for line in fh:
                if line.strip() and not line.startswith(("#", "@")):
                    rows.append(np.array(line.split()[1:], dtype=float))   # drop the time column
        if not rows:
            print(f"No data available in file: {contact_file}")
            continue
        data = np.vstack(rows)
        nres = int(round(np.sqrt(data.shape[1])))
        if nres * nres != data.shape[1]:
            raise ValueError(f"{contact_file}: {data.shape[1]} distances per frame is not a square number of residue pairs")
        matrices.append((data < contact_cutoff).mean(axis=0).reshape(nres, nres))

    if not matrices:
        print("No valid data matrices found.")
        return

    average_contact_matrix = sum(matrices) / len(matrices)
    stem = os.path.join(os.path.dirname(os.path.abspath(contact_files[0])), "contact_fraction")
    np.savetxt(f"{stem}.dat", average_contact_matrix, fmt="%.4f")

    plt.figure(figsize=(10, 8))
    plt.imshow(average_contact_matrix, cmap="Greens", interpolation="nearest")
    plt.colorbar(label="Fraction of Close Contact")
    plt.xlabel("Residue Index (Columns)")
    plt.ylabel("Residue Index (Rows)")
    plt.title("Fraction of Close Contact for Residues")
    plt.savefig(f"{stem}.png", dpi=150)
    plt.show()
    return average_contact_matrix


def contact_map(trajectory, structure, md_file, output_dir, index=None, interactive=False, maxprocs=90, **kwargs):
    """Residue-residue contact fractions from gmx pairdist, averaged over windows of the trajectory.

    As before, each trajectory is cut into 1000 ps windows (0-1000, 1000-2000, 2000-3000 ps), sampled
    every `dt` = 500 ps, and written to <folder>_dist_<begin>.xvg. `md_file`, `start_residue` and
    `end_residue` are accepted and not used. Optional keyword arguments: `dt`, `window` and `total_time`
    (ps), `selection` (gmx group or selection, default "0", the whole system), `contact_cutoff` (nm, default
    0.5) and `gmx`. Returns the averaged matrix.
    """
    if not interactive:
        return _not_interactive("contact_map")
    dt = kwargs.pop("dt", 500)
    window = kwargs.pop("window", 1000)
    total_time = kwargs.pop("total_time", 3000)
    selection = kwargs.pop("selection", "0")
    contact_cutoff = kwargs.pop("contact_cutoff", 0.5)
    gmx = find_gmx(kwargs.pop("gmx", None))
    trajectories, structures, indices = _inputs(trajectory, structure, index)
    labels = _labels(trajectories)
    output_dir = os.path.abspath(output_dir)
    with _working_dir(output_dir):
        add_residue_file(trajectories[0])
        tasks, result_files = [], []
        for t, s, n, label in zip(trajectories, structures, indices, labels):
            for begin in range(0, total_time, window):
                out = f"{label}_dist_{begin}.xvg"
                result_files.append(out)
                tasks.append(partial(_run, [gmx, "pairdist", "-f", t, "-s", s, "-b", str(begin), "-e", str(begin + window),
                                            "-refgrouping", "res", "-selgrouping", "res", "-ref", selection, "-sel", selection,
                                            "-o", out, "-cutoff", "1.0", "-dt", str(dt), *_index_args(n)], output_dir))
        _run_all(tasks, maxprocs)
        return process_contact([os.path.join(output_dir, f) for f in result_files if os.path.isfile(f)], contact_cutoff=contact_cutoff)
