#!/usr/bin/env bash
# Submit MM/GBSA arrays for the original top-5%-per-ensemble plan from the
# 90% structural-only PCA rankings.

set -euo pipefail

ROOT="/scratch/path/to/PRK_SSG_RERUN"
SBATCH_SCRIPT="${ROOT}/ensemble/run_mmgbsa.sbatch"
SELECT_ROOT="${ROOT}/ensemble/mmgbsa_selection"

sbatch \
  --array=0-11 \
  --export=ALL,ENSEMBLE=HM2,FORCE=1,LIST_OVERRIDE="${SELECT_ROOT}/hm2_top5pct_ready_for_new_mmgbsa.txt" \
  "${SBATCH_SCRIPT}"

sbatch \
  --array=0-11 \
  --export=ALL,ENSEMBLE=MED4,FORCE=1,LIST_OVERRIDE="${SELECT_ROOT}/med4_top5pct_ready_for_new_mmgbsa.txt" \
  "${SBATCH_SCRIPT}"
