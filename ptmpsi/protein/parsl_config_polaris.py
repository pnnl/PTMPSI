import os
from parsl.config import Config

# Use LocalProvider to launch workers within a submitted batch job:
from parsl.providers import LocalProvider
# The high throughput executor is for scaling to HPC systems:
from parsl.executors import HighThroughputExecutor
# The SimpleLauncher is the default launcher for the LocalProvider:
from parsl.launchers import SimpleLauncher
# address_by_interface is needed for the HighThroughputExecutor:
from parsl.addresses import address_by_interface
def gen_parsl_config(worker_init=None, num_nodes=1, workers_per_node=1, nodes_per_block=1, run_dir=None):
    if run_dir is None:
        run_dir = os.path.join(os.getcwd(), "runinfo")

    config = Config(
        executors=[
            HighThroughputExecutor(
                label="htex",
                heartbeat_period=15,
                heartbeat_threshold=120,
                worker_debug=False,
                max_workers_per_node=workers_per_node,
                address=address_by_interface("bond0"),
                cpu_affinity="block-reverse",
                prefetch_capacity=0,
                cores_per_worker=32,
                provider=LocalProvider(
                    # Number of nodes job
                    worker_init=worker_init,
                    nodes_per_block=nodes_per_block,
                    launcher=SimpleLauncher(),
                    init_blocks=num_nodes,
                    max_blocks=num_nodes,
                    parallelism=1.0,
                ),
            ),
        ],
        run_dir=run_dir,
    )
    return config
