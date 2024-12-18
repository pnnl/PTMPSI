import os
import psutil
import subprocess
import time
import numpy as np

def check_file_exists(filepath):
    """Check if a file exists, raise FileNotFoundError if not."""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File {filepath} does not exist.")

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
def contact_map(trajectory, structure, md_file, output_dir, index=None, group=0, interactive=False, maxprocs=90, **kwargs):
    """
    Generates contact maps.
    
    Args:
        trajectory (str or list): The trajectory file(s).
        structure (str or list): The structure file(s).
        md_file (str): The .gro file containing the number of residues.
        output_dir (str): The directory to write output files to.
        index (str or list, optional): Index file(s). If None, no index is used.
        interactive (bool, optional): Whether to use an interactive process allocation.
        maxprocs (int, optional): Max number of processors to use.
    """

    print("1) Initial working directory is:", os.getcwd())

    # Ensure the output directory exists
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

    # Change to the output directory
    print(f"Changing to output directory: {output_dir}")
    original_directory = os.getcwd()
    os.chdir(output_dir)
    print("2) Current working directory is:", os.getcwd())

    # Handle multi input cases
    __trajectories = trajectory if isinstance(trajectory, list) else [trajectory]
    __structures = structure if isinstance(structure, list) else [structure]
    __indices = index if isinstance(index, list) else [index]

    ntraj = len(__trajectories)
    
    # Validate all required files exist
    for i in range(ntraj):
        check_file_exists(__trajectories[i])
        check_file_exists(__structures[i])
        if __indices[i] is not None:
            check_file_exists(__indices[i])

    # Extract the number of residues from the md file
    nres = get_number_of_residues(md_file)

    # Interactive process handling
    if interactive:
        nprocs = psutil.cpu_count(logical=False)  # Available physical CPUs
        nlogical = psutil.cpu_count()
        nskip = int(nlogical / nprocs)
        __nprocs = min(maxprocs, nprocs)

        print(f"Found {nprocs} physical CPUs available")
        print(f"Will use {__nprocs} CPUs at the same time")

        commands = "source /anfhome/.profile && module load gcc/13.2.0 intel-oneapi-mkl/2023.2.0/gcc13.2.0-hpcx_mpi-ddupx"

        for i in range(ntraj):
            basename_index = os.path.basename(os.path.dirname(__trajectories[i]))

            matrix = np.zeros((nres, nres))
            num_frames = 0

            while True:
                begin, end = num_frames, num_frames + 1
                _commands = commands     
                _commands += f"gmx pairdist -f {__trajectories[i]} -s {__structures[i]} -b {begin} -e {end} -refgrouping res -selgrouping res -ref {group} -sel {group} -o {basename_index}_dist_{num_frames}.xvg -cutoff 1.0"
                
                if __indices[i] is not None:
                    _commands += f" -n {__indices[i]}"
                
                print(f"Executing command: {_commands}")
                process = subprocess.Popen(_commands, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()

                if process.returncode != 0:
                    print(f"Process exited with code: {process.returncode}")
                    print(f"Command failed for file {basename_index}_dist_{num_frames}.xvg, checking logs...")
                    print(f"stdout: {stdout.decode()}")
                    print(f"stderr: {stderr.decode()}")
                    break

                # Check if file was generated
                contact_file = f"{basename_index}_dist_{num_frames}.xvg"
                if not os.path.exists(os.path.join(output_dir, contact_file)):
                    break

                # Process file
                with open(contact_file, 'r', encoding='latin-1') as f:
                    for line in f.readlines():
                        distances = line.split()
                        for idx, dist in enumerate(distances):
                            row = idx // nres
                            col = idx % nres
                            if float(dist) < 0.5:
                                matrix[row, col] += 1
                
                num_frames += 1

            # Calculate the fraction of frames in close contact
            fraction_matrix = matrix / num_frames if num_frames > 0 else matrix
            np.savetxt(os.path.join(output_dir, f"{basename_index}_contact_map.dat"), fraction_matrix)

    # Change back to the original directory
    os.chdir(original_directory)
    print("3) Changed back to original directory:", os.getcwd())
