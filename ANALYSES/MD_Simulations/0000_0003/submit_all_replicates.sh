#!/usr/bin/env bash
set -euo pipefail

# submit_all_replicates.sh — Prepare and optionally submit all 12 replicates.
# Usage: ./submit_all_replicates.sh [--submit]

SUBMIT=0
[[ "${1:-}" == "--submit" ]] && SUBMIT=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE="$(dirname "$SCRIPT_DIR")"
MDP_DIR="$BASE/mdp"
REF_VS_REF="/scratch/path/to/MED4_NAD_vs_HM2_NAD/NEW/REF_vs_REF"
NAD_AWARE_FF="/scratch/path/to/HM2_NAD/postprocess_new/0260/mmgbsa/amber99sb.ff"
GMX="/path/to/gromacs_local/gromacs-2025.0-install/bin/gmx"

# Topology sources per condition
declare -A TOP_SRC
TOP_SRC["med4_0000"]="/scratch/path/to/MED4_NAD_vs_HM2_NAD/NEW/REF_vs_REF/replicate_fix_L2/topology/med4_0000"
TOP_SRC["med4_0003"]="/scratch/path/to/MED4_NAD/postprocess_new/REF"
TOP_SRC["hm2_0000"]="/scratch/path/to/MED4_NAD_vs_HM2_NAD/NEW/REF_vs_REF/replicate_fix_L2/topology/hm2_0000"
TOP_SRC["hm2_0003"]="/scratch/path/to/HM2_NAD/postprocess_new/REF/mmpbsa_nad"

# Rep2 source (for extracting starting coords)
rep2_tpr() {
  local s="$1" sys="$2"
  local redox; [[ "$sys" == "0000" ]] && redox="red" || redox="ox"
  echo "${REF_VS_REF}/${s}/${s}_${redox}/rep2/md.tpr"
}

# Seeds: strain sys rep heating_seed fheating_seed
SEEDS=(
  "med4 0000 rep1 1992983281 1399811784"
  "med4 0000 rep2 1883869658 904763506"
  "med4 0000 rep3 632103819 861547114"
  "med4 0003 rep1 2017480395 1072361057"
  "med4 0003 rep2 10867561 410789386"
  "med4 0003 rep3 501868637 1851879791"
  "hm2 0000 rep1 1751167737 1588524991"
  "hm2 0000 rep2 1335661749 1034025992"
  "hm2 0000 rep3 1310592593 726501058"
  "hm2 0003 rep1 1080287867 1987670669"
  "hm2 0003 rep2 2065477894 143415490"
  "hm2 0003 rep3 2044689295 495353219"
)

echo "=== Replicate setup (12 reps) ==="
echo ""

for entry in "${SEEDS[@]}"; do
  read -r STRAIN SYS REP HSEED FSEED <<< "$entry"

  d="$BASE/$STRAIN/$SYS/$REP"
  top_src="${TOP_SRC[${STRAIN}_${SYS}]}"
  src_tpr=$(rep2_tpr "$STRAIN" "$SYS")

  echo "--- $STRAIN/$SYS/$REP (h_seed=$HSEED, f_seed=$FSEED) ---"

  # Skip if already set up
  if [[ -f "$d/md.tpr" || -f "$d/md.xtc" ]]; then
    echo "  Already set up or running, skipping."
    [[ $SUBMIT -eq 1 ]] && echo "  To resubmit: (cd $d && sbatch run_replicate.sbatch)"
    echo ""
    continue
  fi

  # Extract starting coords (PDB for higher precision)
  echo "  Extracting starting coords from rep2..."
  ${GMX} editconf -f "$src_tpr" -o "$d/start.pdb" > /dev/null 2>&1

  # Symlink topology + NAD-aware FF
  ln -sf "$top_src/topol.top" "$d/topol.top"
  ln -sf "$NAD_AWARE_FF" "$d/amber99sb.ff"

  # Build custom index file (Protein_NAD for tc-grps, matching the
  # T-coupling group structure used elsewhere in this project)
  ${GMX} make_ndx -f "$d/start.pdb" -o "$d/index.ndx" > /dev/null 2>&1 <<'NDX_EOF'
"Protein" | "NAD"
name 21 Protein_NAD
! "Protein_NAD"
name 22 non_Protein_NAD
q
NDX_EOF

  # Copy all mdp files, substituting seeds
  for f in "$MDP_DIR"/*.mdp; do
    cp "$f" "$d/"
  done
  sed -i "s/gen-seed.*=.*/gen-seed              = ${HSEED}/" "$d/heating.mdp"
  sed -i "s/gen-seed.*=.*/gen-seed              = ${FSEED}/" "$d/fheating.mdp"

  # Copy pipeline script + sbatch
  cp "$SCRIPT_DIR/run_replicate_pipeline.sh" "$d/"
  cp "$SCRIPT_DIR/run_replicate.sbatch" "$d/"
  chmod +x "$d/run_replicate_pipeline.sh"

  # Verify
  h_check=$(grep "gen-seed" "$d/heating.mdp" | awk '{print $NF}')
  f_check=$(grep "gen-seed" "$d/fheating.mdp" | awk '{print $NF}')
  echo "  Setup done. heating_seed=$h_check  fheating_seed=$f_check"

  if [[ $SUBMIT -eq 1 ]]; then
    out=$(cd "$d" && sbatch run_replicate.sbatch 2>&1)
    echo "  Submitted: $out"
  else
    echo "  To submit: (cd $d && sbatch run_replicate.sbatch)"
  fi
  echo ""
done

echo "=== Done ==="
[[ $SUBMIT -eq 0 ]] && echo "Run with --submit to queue all 12 jobs."
