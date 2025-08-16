#!/bin/bash
set -euo pipefail

# Ensure conda is available in this shell
source "$(conda info --base)/etc/profile.d/conda.sh"

# Pick the YAML file (adjust the path if you move it)
YML="./pca_frontier.yml"

if [[ ! -f "$YML" ]]; then
  echo "[ERROR] $YML not found. Put setup.sh in the same folder as pca_frontier.yml, or edit YML path."
  exit 1
fi

# Create env if missing (name comes from the YAML; it's 'gppenv' in your file)
if ! conda env list | awk '{print $1}' | grep -qx gppenv; then
  echo "[SETUP] Creating conda env 'gppenv' from $YML"
  conda env create -n gppenv -f "$YML"
else
  echo "[SETUP] Conda env 'gppenv' already exists — updating from $YML (safe)."
  conda env update -n gppenv -f "$YML"
fi

# Activate for interactive use / sanity checks
conda activate gppenv
exec python "$@"
which python
python -V

echo "[SETUP] gppenv is ready. Next:"
echo "  python generate_pipelines.py"
