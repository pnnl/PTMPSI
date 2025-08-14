#!/usr/bin/env bash
# run_in_gppenv.sh — wrapper to activate gppenv and then exec any Python script

source "/ccs/home/sama578/miniconda3/etc/profile.d/conda.sh"
conda activate gppenv
exec python "$@"

