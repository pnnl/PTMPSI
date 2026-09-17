#!/usr/bin/env bash
set -euo pipefail

# run_replicate_pipeline.sh — Full replicate pipeline for one rep.
# Runs: genrestr → minimization chain → heating(unique seed) → NPT →
#       minim → fheating(unique seed) → fnpt → production → postprocess
#
# Called from run_replicate.sbatch with cwd = the rep directory.
# Expects these files already in the rep dir (placed by submit_all_replicates.sh):
#   start.pdb, topol.top, amber99sb.ff (symlink), index.ndx,
#   all *.mdp files (with seeds substituted in heating.mdp + fheating.mdp)

# ----------------------------------------------------- GROMACS env --------
. /pscratch/sd/w/widm288/PRK_SSG_RERUN/setup_gmx.sh
export GMX_MAXBACKUP=-1
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export TMPDIR=${SCRATCH:-$PWD}

GMX_MPI="${GROMACS_HOME}/bin/gmx_mpi"
GMX="/path/to/gromacs_local/gromacs-2025.0-install/bin/gmx"
PY_BIN="/path/to/.conda/envs/gppenv/bin/python"
PCA_DIR="/scratch/path/to/PTMPSI/ptmpsi-PCA"
export PYTHONPATH="/scratch/path/to/PTMPSI/ptmpsi-PCA:${PYTHONPATH:-}"

REF_SYS_MED4="/scratch/path/to/MED4_NAD/postprocess_new/REF"
REF_SYS_HM2="/scratch/path/to/HM2_NAD/postprocess_new/REF"
WHICH_REF="${REF_SYS_MED4}"
[[ "$PWD" == *"/hm2/"* ]] && WHICH_REF="${REF_SYS_HM2}"

echo "[REPLICATE] Starting in $PWD at $(date)"

# ------------------------------------------ 1) position restraint itps ----
if [[ ! -f 1000p.itp ]]; then
  echo "[REPLICATE] Generating position restraint .itp files..."
  echo q | ${GMX} make_ndx -f start.pdb -o genrestr.ndx > /dev/null 2>&1
  for fc in 1000 500 100 50 10 5 1; do
    echo Protein | ${GMX} genrestr -f start.pdb -n genrestr.ndx \
        -o "${fc}p.itp" -fc "${fc}" "${fc}" "${fc}" > /dev/null 2>&1
  done
  echo "[REPLICATE] Position restraint .itp files created."
fi

# ------------------------------------------ 2) restrained minimization ----
prev_gro="start.pdb"
for stage in minim1000 minim500 minim100 minim50 minim10 minim5 minim1; do
  if [[ ! -f "${stage}.gro" ]]; then
    echo "[REPLICATE] ${stage}..."
    srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f "${stage}.mdp" \
        -c "${prev_gro}" -r "${prev_gro}" -p topol.top -n index.ndx -o "${stage}.tpr" > /dev/null 2>&1
    srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun -nb gpu -deffnm "${stage}" -v > /dev/null 2>&1
  fi
  prev_gro="${stage}.gro"
done

# ------------------------------------------ 3) heating (independence!) ----
if [[ ! -f heating.gro ]]; then
  echo "[REPLICATE] Heating (gen-vel=yes, unique seed)..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f heating.mdp \
      -c minim1.gro -r minim1.gro -p topol.top -n index.ndx -o heating.tpr > /dev/null 2>&1
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm heating 2>&1 | tail -3
fi

# ------------------------------------------ 4) NPT equilibration ---------
if [[ ! -f npt.gro ]]; then
  echo "[REPLICATE] NPT equilibration..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f npt.mdp \
      -c heating.gro -t heating.cpt -r heating.gro -p topol.top -n index.ndx -o npt.tpr > /dev/null 2>&1
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm npt 2>&1 | tail -3
fi

# ------------------------------------------ 5) unrestrained minimization --
if [[ ! -f minim.gro ]]; then
  echo "[REPLICATE] Unrestrained minimization..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f minim.mdp \
      -c npt.gro -p topol.top -n index.ndx -o minim.tpr > /dev/null 2>&1
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun -nb gpu -deffnm minim -v > /dev/null 2>&1
fi

# ------------------------------------------ 6) final heating -------------
if [[ ! -f fheating.gro ]]; then
  echo "[REPLICATE] Final heating (gen-vel=yes, unique seed)..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f fheating.mdp \
      -c minim.gro -p topol.top -n index.ndx -o fheating.tpr > /dev/null 2>&1
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm fheating 2>&1 | tail -3
fi

# ------------------------------------------ 7) final NPT ----------------
if [[ ! -f fnpt.gro ]]; then
  echo "[REPLICATE] Final NPT..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f fnpt.mdp \
      -c fheating.gro -t fheating.cpt -p topol.top -n index.ndx -o fnpt.tpr > /dev/null 2>&1
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm fnpt 2>&1 | tail -3
fi

# ------------------------------------------ 8) production tpr -----------
if [[ ! -f md.tpr ]]; then
  echo "[REPLICATE] Building production tpr..."
  srun -n 1 ${GMX_MPI} grompp -maxwarn 99 -f md.mdp \
      -c fnpt.gro -t fnpt.cpt -p topol.top -n index.ndx -o md.tpr > /dev/null 2>&1
fi

# ------------------------------------------ 9) production mdrun ----------
echo "[REPLICATE] Production (100 ns)..."
if [[ -f md.cpt ]]; then
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -cpi md.cpt -maxh 29.0 \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm md
else
  srun -n 4 -c 16 --cpu-bind=cores ${GMX_MPI} mdrun \
      -maxh 29.0 \
      -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm md
fi

# Check completion
nsteps_target=$(awk '$1 == "nsteps" {print $3; exit}' md.mdp)
last_step=$(grep -oE '^\s+[0-9]+\s+[0-9]' md.log | tail -1 | awk '{print $1}' || echo 0)
last_step=${last_step:-0}
echo "[REPLICATE] last step: ${last_step} / target: ${nsteps_target}"
if [[ "${last_step}" -lt "${nsteps_target}" ]]; then
  echo "[REPLICATE] Production not yet complete; resubmit to continue from md.cpt."
  exit 0
fi

# ------------------------------------------ 10) postprocessing -----------
echo "[REPLICATE] Postprocessing..."
if [[ ! -f solute_fit.xtc ]]; then
  ${GMX} make_ndx -f md.tpr -o solute_index.ndx <<'EOF'
"System" & ! "Water_and_ions"
q
EOF
  sed -i 's/System_&_!Water_and_ions/Solute/' solute_index.ndx
  echo Solute | ${GMX} convert-tpr -s md.tpr -n solute_index.ndx -o solute.tpr
  echo Solute | ${GMX} trjconv -f md.xtc -s solute.tpr -n solute_index.ndx -o solute_raw.xtc
  echo 0     | ${GMX} trjconv -f solute_raw.xtc -s solute.tpr -o step1.xtc -pbc nojump
  echo 0 0   | ${GMX} trjconv -f step1.xtc      -s solute.tpr -o step2.xtc -pbc cluster
  echo 0     | ${GMX} trjconv -f step2.xtc      -s solute.tpr -o step3.xtc -pbc atom -ur compact
  echo 0     | ${GMX} trjconv -f step3.xtc      -s solute.tpr -o reference.gro -dump 0
  echo 0     | ${GMX} trjconv -f solute_raw.xtc -s reference.gro -o md_fixed_nojump.xtc -pbc nojump
  echo 1 0   | ${GMX} trjconv -f md_fixed_nojump.xtc -s solute.tpr -fit rot+trans -o solute_fit.xtc
  echo 0     | ${GMX} trjconv -f solute_fit.xtc -s solute.tpr -dump 0 -o solute_0.gro
fi

if [[ ! -f solute_fixed.pdb ]]; then
  ${PY_BIN} - <<'PY'
import preprocess
preprocess.gro_to_pdb("solute_0.gro",   "solute_0_raw.pdb")
preprocess.add_chain_ids("solute_0_raw.pdb", "solute_0_chains.pdb")
preprocess.fix_pbc("solute_0_chains.pdb", "solute_fixed.pdb")
PY
fi

name="$(basename "$(dirname "$PWD")")_$(basename "$PWD")"
if [[ ! -f figures_rmsd_rg/*_rmsd.pdf ]]; then
  ${PY_BIN} ${PCA_DIR}/rmsd_rg.py \
      --top solute_fixed.pdb --xtc solute_fit.xtc \
      --ref_top "${WHICH_REF}/solute_fixed.pdb" --ref_xtc "${WHICH_REF}/solute_fit.xtc" \
      --skip 1 --window 10 --name "$name"
fi

if [[ ! -f figures_rmsf_secstruct/*_rmsf_secstruct.pdf ]]; then
  ${PY_BIN} ${PCA_DIR}/rmsf_secstruct.py \
      --top solute_fixed.pdb --xtc solute_fit.xtc \
      --ref_top "${WHICH_REF}/solute_fixed.pdb" --ref_xtc "${WHICH_REF}/solute_fit.xtc" \
      --all --name "$name"
fi

echo "[REPLICATE] COMPLETE at $(date)"
