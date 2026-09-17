#!/bin/bash
#SBATCH --job-name=med4_pca_95pct
#SBATCH --account=m4742
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=03:00:00
#SBATCH --licenses=scratch
#SBATCH --output=/scratch/path/to/PRK_SSG_RERUN/ensemble/logs/%x_%j.out
#SBATCH --error=/scratch/path/to/PRK_SSG_RERUN/ensemble/logs/%x_%j.err

set -euo pipefail

PYTHON=/path/to/.conda/envs/gppenv/bin/python
SCRIPTS_DIR=/scratch/path/to/scripts_perlmutter/PTMPSI/ptmpsi-PCA
PCA_SCRIPT=/scratch/path/to/PRK_SSG_RERUN/ensemble/run_pca_analysis.py
POST_BASE=/scratch/path/to/MED4_NAD/postprocess_final
REFERENCE_SOURCE=${POST_BASE}/REF
PCA_INPUT=${POST_BASE}/pca_input_ref_plus_numeric
REF_DIR=${PCA_INPUT}/REF
FEATURE_SOURCE_DIR=${POST_BASE}/suman_pca_fixed_90pct
OUTPUT_DIR=${POST_BASE}/suman_pca_fixed_fullpca_95pct
AFFECTED=/scratch/path/to/PRK_SSG_RERUN/MED4/PRK_SSG_RERUN_INTERACTIVE/affected_systems.txt

mkdir -p /scratch/path/to/PRK_SSG_RERUN/ensemble/logs
mkdir -p "${OUTPUT_DIR}" "${PCA_INPUT}"

export PYTHONPATH="${SCRIPTS_DIR}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/${USER}_matplotlib_${SLURM_JOB_ID}
mkdir -p "${MPLCONFIGDIR}"

for required in \
  "${PYTHON}" \
  "${PCA_SCRIPT}" \
  "${AFFECTED}" \
  "${REFERENCE_SOURCE}/solute_fixed.pdb" \
  "${REFERENCE_SOURCE}/solute_fit.xtc" \
  "${FEATURE_SOURCE_DIR}/feature_outputs"; do
  if [[ ! -e "${required}" ]]; then
    echo "Missing required input/cache: ${required}" >&2
    exit 1
  fi
done

"${PYTHON}" - <<'PY'
import MDAnalysis as mda
import sklearn
print(f"Using MDAnalysis {mda.__version__}")
print(f"Using scikit-learn {sklearn.__version__}")
PY

missing=0
checked=0
while read -r id; do
  [[ -z "${id}" ]] && continue
  checked=$((checked + 1))
  for required in "${POST_BASE}/${id}/solute_fixed.pdb" "${POST_BASE}/${id}/solute_fit.xtc"; do
    if [[ ! -f "${required}" ]]; then
      echo "Missing MED4 postprocess input for ${id}: ${required}" >&2
      missing=$((missing + 1))
    fi
  done
done < "${AFFECTED}"

echo "Checked affected rerun systems: ${checked}"
if [[ "${missing}" -ne 0 ]]; then
  echo "MED4 postprocess outputs incomplete: ${missing} missing files" >&2
  exit 1
fi

if [[ -f "${OUTPUT_DIR}/pca_distance_to_ref.csv" ]]; then
  echo "95 pct PCA output already exists: ${OUTPUT_DIR}/pca_distance_to_ref.csv"
  echo "Skipping PCA run."
  exit 0
fi

ln -sfn "${FEATURE_SOURCE_DIR}/feature_outputs" "${OUTPUT_DIR}/feature_outputs"

# Build a clean PCA input namespace with REF exactly once.
# MED4 REF and numeric 0003 are the same physical system; include REF and skip 0003.
find "${PCA_INPUT}" -mindepth 1 -maxdepth 1 -type l -delete
ln -sfn "${REFERENCE_SOURCE}" "${PCA_INPUT}/REF"
while IFS= read -r d; do
  sys=$(basename "${d}")
  if [[ "${sys}" == "0003" ]]; then
    continue
  fi
  if [[ "${sys}" =~ ^[0-9]{4}$ && -f "${d}/solute_fixed.pdb" && -f "${d}/solute_fit.xtc" ]]; then
    ln -sfn "${d}" "${PCA_INPUT}/${sys}"
  fi
done < <(find "${POST_BASE}" -mindepth 1 -maxdepth 1 -type d)

if [[ ! -f "${REF_DIR}/solute_fixed.pdb" || ! -f "${REF_DIR}/solute_fit.xtc" ]]; then
  echo "PCA reference is missing after symlink setup: ${REF_DIR}" >&2
  exit 1
fi

echo "PCA input systems including REF: $(find "${PCA_INPUT}" -mindepth 1 -maxdepth 1 -type l | wc -l)"
echo "Starting MED4 95 pct full-PCA run at $(date)"
echo "Base: ${PCA_INPUT}"
echo "Reference: ${REF_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Feature cache: ${FEATURE_SOURCE_DIR}/feature_outputs"

"${PYTHON}" "${PCA_SCRIPT}" \
  --base_dir "${PCA_INPUT}" \
  --ref_dir "${REF_DIR}" \
  --output_dir "${OUTPUT_DIR}" \
  --topology solute_fixed.pdb \
  --trajectories solute_fit.xtc \
  --step 1 \
  --pca_variance_threshold 0.95 \
  --pca_sensitivity_thresholds 0.80,0.90,0.95 \
  --n_clusters 3 \
  --per_system_jobs 32 \
  --dt_ns 0.1

for expected in \
  "${OUTPUT_DIR}/explained_variance.csv" \
  "${OUTPUT_DIR}/full_scree_variance.csv" \
  "${OUTPUT_DIR}/pca_scree_plot.pdf" \
  "${OUTPUT_DIR}/pca_scree_plot.png" \
  "${OUTPUT_DIR}/pca_scree_elbow_candidates.csv" \
  "${OUTPUT_DIR}/pca_distance_to_ref.csv" \
  "${OUTPUT_DIR}/pca_distance_to_ref_ranked.csv" \
  "${OUTPUT_DIR}/pca_distance_to_ref_95pct_ranked.csv" \
  "${OUTPUT_DIR}/pca_sensitivity_rank_stability.csv" \
  "${OUTPUT_DIR}/pca_sensitivity_recommendation.txt" \
  "${OUTPUT_DIR}/kinetic_observables_by_system.csv"; do
  if [[ ! -f "${expected}" ]]; then
    echo "Expected output was not created: ${expected}" >&2
    exit 1
  fi
done

echo "MED4 95 pct full-PCA run completed at $(date)"
echo "Output: ${OUTPUT_DIR}"
