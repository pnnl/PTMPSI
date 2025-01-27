from argparse import ArgumentParser

import parsl

from parsl_config_polaris import gen_parsl_config
from parsl.app.app import bash_app
import os

@bash_app
def md_simulations(stdout=None, stderr=None, bash_str=''):
    return bash_str

# Protect the script from running twice.
#  See "Safe importing of main module" in Python multiprocessing docs
#  https://docs.python.org/3/library/multiprocessing.html#multiprocessing-programming
if __name__ == "__main__":
    # Get user instructions
    parser = ArgumentParser()
    parser.add_argument('--start-id', default=0, type=int, help='Starting ID of the simulation batch')
    parser.add_argument('--nodes-per-block', default=1, type=int, help='Number of nodes per block')
    parser.add_argument('--workers-per-node', default=1, type=int, help='Number of workers per node')
    parser.add_argument('--run-dir', default=None, type=str, help='Directory to run the simulation in')
    parser.add_argument('--worker-init', default="worker_init.sh", type=str, help='Worker initialization script')
    parser.add_argument('--joblist', default="parsl_jobs.txt", type=str, help='List of jobs to run')

    args = parser.parse_args()

    pbs_nodefile = os.environ['PBS_NODEFILE']
    with open(pbs_nodefile, 'r') as f:
        nodelist = f.readlines()
    num_nodes = len(nodelist)
    
    pbs_jobid = os.environ['PBS_JOBID']

    with open(args.worker_init, 'r') as f:
        parsl_worker_init = f.readlines()
    
    if args.run_dir is None:
        args.run_dir = os.path.join(os.getcwd(), "runinfo")

    # Load the configuration
    #  As a context manager so resources are shutdown on exit
    with parsl.load(gen_parsl_config(''.join(parsl_worker_init), num_nodes, args.workers_per_node, args.nodes_per_block, args.run_dir)) as cfg:
        
        with open(args.joblist, "r") as f:
            job_list = f.readlines()
        
        job_list = [job.strip() for job in job_list]

        if len(job_list) < args.start_id + num_nodes:
            raise ValueError(f"Not enough jobs in {args.joblist} to run on {num_nodes} nodes starting at {args.start_id}")  

        for i, job in enumerate(job_list):
            dir_name = os.path.dirname(job)
            with open(os.path.join(dir_name, "hostfile"), "w") as f:
                f.write(nodelist[i])

        # Spawn tasks
        futures = [
            md_simulations(stdout=f'{job_list[i]}_{pbs_jobid}_out.txt', stderr=f'{job_list[i]}_{pbs_jobid}_err.txt', bash_str=job_list[i])
            for i in range(args.start_id, args.start_id + num_nodes)
        ]

        # Retrieve task results with error handling
        for i, future in enumerate(futures, start=args.start_id):
            try:
                result = future.result()
                print(result)
            except Exception as e:
                print(f"An error occurred in job {i}: {e}")