#!/usr/bin/env bash
#SBATCH -A m4742
#SBATCH -J mmgbsa_final
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH -o /scratch/path/to/PRK_SSG_RERUN/ensemble/logs/mmgbsa_final_%x_%A_%a.out
#SBATCH -e /scratch/path/to/PRK_SSG_RERUN/ensemble/logs/mmgbsa_final_%x_%A_%a.err
#SBATCH --array=0-5

set -euo pipefail

ENSEMBLE="${ENSEMBLE:?Set ENSEMBLE=HM2 or ENSEMBLE=MED4 with sbatch --export}"

ROOT="/scratch/path/to/PRK_SSG_RERUN"
SELECT_ROOT="${ROOT}/ensemble/mmgbsa_selection"

case "${ENSEMBLE}" in
  HM2)
    BASE="/scratch/path/to/HM2_NAD/postprocess_final"
    OLD_BASE="/scratch/path/to/HM2_NAD/postprocess_new"
    RERUN_INPUT="${ROOT}/HM2/PRK_SSG_RERUN_INTERACTIVE"
    FORCEFIELD_SRC="${RERUN_INPUT}/amber99sb.ff"
    FRONTIER_TUPLES="/path/to/HM2/1tuples"
    LIST="${SELECT_ROOT}/hm2_systems_ready_for_new_mmgbsa.txt"
    TEMPLATE="/scratch/path/to/HM2_NAD/postprocess_new/REF/mmgbsa"
    ;;
  MED4)
    BASE="/scratch/path/to/MED4_NAD/postprocess_final"
    OLD_BASE="/scratch/path/to/MED4_NAD/postprocess_new"
    RERUN_INPUT="${ROOT}/MED4/PRK_SSG_RERUN_INTERACTIVE"
    FORCEFIELD_SRC="${RERUN_INPUT}/amber99sb.ff"
    FRONTIER_TUPLES="/path/to/MED4/1tuples"
    LIST="${SELECT_ROOT}/med4_systems_ready_for_new_mmgbsa.txt"
    TEMPLATE="/scratch/path/to/MED4_NAD/postprocess_new/REF/mmgbsa"
    ;;
  *)
    echo "ERROR: unsupported ENSEMBLE=${ENSEMBLE}" >&2
    exit 2
    ;;
esac

LIST="${LIST_OVERRIDE:-${LIST}}"

RANKS_PER_SYS=32
TASKS_PER_NODE=${SLURM_NTASKS:-128}
PER_NODE=$(( TASKS_PER_NODE / RANKS_PER_SYS ))
BASE_IDX=$(( SLURM_ARRAY_TASK_ID * PER_NODE ))

if [[ ! -f "${LIST}" ]]; then
  echo "ERROR: system list not found: ${LIST}" >&2
  exit 2
fi

if (( TASKS_PER_NODE % RANKS_PER_SYS != 0 )); then
  echo "ERROR: SLURM_NTASKS (${TASKS_PER_NODE}) not divisible by ${RANKS_PER_SYS}" >&2
  exit 2
fi

mapfile -t SYSTEMS < <(awk 'NF && $1 !~ /^#/ {print $1}' "${LIST}")
N_SYS=${#SYSTEMS[@]}

echo "JobID=${SLURM_JOB_ID} ArrayTask=${SLURM_ARRAY_TASK_ID} Ensemble=${ENSEMBLE} NodeList=${SLURM_NODELIST}"
echo "TASKS_PER_NODE=${TASKS_PER_NODE} RANKS_PER_SYS=${RANKS_PER_SYS} PER_NODE=${PER_NODE}"
echo "Systems total=${N_SYS} starting index=${BASE_IDX}"

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/.conda/envs/gmxMMPBSA

export PATH="/path/to/gromacs_local/gromacs-2025.0-install/bin:${PATH}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

find_topol () {
  local sys="$1"
  local candidate
  for candidate in \
    "${RERUN_INPUT}/${sys}/topol.top" \
    "${OLD_BASE}/${sys}/mmgbsa/topol.top" \
    "${OLD_BASE}/${sys}/topol.top" \
    "${FRONTIER_TUPLES}/${sys}/topol.top" \
    "${BASE}/${sys}/mmgbsa/topol.top"
  do
    if [[ -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

prepare_one () {
  local sys="$1"
  local src="${BASE}/${sys}"
  local workdir="${src}/mmgbsa"
  local topol

  topol="$(find_topol "${sys}")" || {
    echo "ERROR ${sys}: no topol.top source found" >&2
    return 2
  }

  for required in solute_fixed.pdb solute_fit.xtc solute.tpr; do
    if [[ ! -f "${src}/${required}" ]]; then
      echo "ERROR ${sys}: missing ${src}/${required}" >&2
      return 2
    fi
  done

  mkdir -p "${workdir}"
  cp "${TEMPLATE}/mmpbsa.in" "${workdir}/mmpbsa.in"
  cp "${TEMPLATE}/make_index_groups_withNAD.py" "${workdir}/make_index_groups_withNAD.py"
  if [[ "${topol}" != "${workdir}/topol.top" ]]; then
    cp "${topol}" "${workdir}/topol.top"
  fi
  if [[ ! -d "${FORCEFIELD_SRC}" ]]; then
    echo "ERROR ${sys}: missing custom force field directory ${FORCEFIELD_SRC}" >&2
    return 2
  fi
  if [[ -e "${workdir}/amber99sb.ff" || -L "${workdir}/amber99sb.ff" ]]; then
    rm -rf "${workdir}/amber99sb.ff"
  fi
  ln -s "${FORCEFIELD_SRC}" "${workdir}/amber99sb.ff"
  ln -sfn "../solute_fixed.pdb" "${workdir}/solute_fixed.pdb"
  ln -sfn "../solute_fit.xtc" "${workdir}/solute_fit.xtc"
  ln -sfn "../solute.tpr" "${workdir}/solute.tpr"

  (
    cd "${workdir}"
    python make_index_groups_withNAD.py solute_fixed.pdb mmpbsa.ndx
  )
}

run_one () {
  local sys="$1"
  local workdir="${BASE}/${sys}/mmgbsa"

  (
    prepare_one "${sys}"
    cd "${workdir}"

    if [[ -f FINAL_RESULTS_MMPBSA.csv && "${FORCE:-0}" != "1" ]]; then
      echo "SKIP ${sys}: FINAL_RESULTS_MMPBSA.csv exists; set FORCE=1 to overwrite"
      exit 0
    fi

    echo "RUN ${ENSEMBLE} ${sys} in ${workdir} using ${RANKS_PER_SYS} ranks"
    srun --exclusive --mpi=pmix --cpu-bind=cores -n "${RANKS_PER_SYS}" \
      gmx_MMPBSA -O \
        -i  mmpbsa.in \
        -cr solute_fixed.pdb \
        -cs solute.tpr \
        -ct solute_fit.xtc \
        -ci mmpbsa.ndx \
        -cg 5 4 \
        -cp topol.top \
        -o  FINAL_RESULTS_MMPBSA.dat \
        -eo FINAL_RESULTS_MMPBSA.csv \
        -nogui \
      > "mmpbsa_${sys}.log" 2> "mmpbsa_${sys}.err"
    echo "DONE ${ENSEMBLE} ${sys}"
  ) &
}

pids=()
rc=0
for ((i=0; i<PER_NODE; i++)); do
  idx=$(( BASE_IDX + i ))
  if (( idx < N_SYS )); then
    run_one "${SYSTEMS[$idx]}" || rc=$?
    pids+=($!)
  fi
done

for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    rc=1
  fi
done

exit "${rc}"
