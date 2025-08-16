#!/usr/bin/env python
# coding: utf-8
"""
generate_pipelines.py

Writes:
  <POST_BASE>/<ID>/md.sh          – per-system post-processing script (plain bash; no Flux/Slurm)
  <POST_BASE>/md_fluxjobs.txt     – list of per-system jobs (kept for continuity)
  <POST_BASE>/md_bundle_flux.sh   – LOCAL sequential runner:
                                     REF → each system; after each: run REF+SYS pairwise PCA;
                                     finally: JOIN pairwise CSVs and PLOT (no heavy global PCA)

# Legacy (preserved as comments; NOT emitted now):
#  <POST_BASE>/<ID>/md.sbatch      – Flux script per system
#  <POST_BASE>/md_bundle.sbatch    – Slurm wrapper that used to start Flux
#  <POST_BASE>/md_bundle_flux.sh   – Flux submitter (old behavior shown below)

This version calls the Python analysis scripts by absolute path from
SCRIPTS_DIR so you can keep them all in one place, and skips solute‐only
operations if they've already been run.

Execution order:
  REF (if present) → each numeric system (sequential)
  After each system OK: run PAIRWISE PCA (REF + that system) to see results early.
  Finally: join pairwise CSVs into POST_BASE/PC_score_ranking.csv and
           POST_BASE/combined_score_ranking.csv, then run plots.
"""

import os, stat
from pathlib import Path
import subprocess

# ── Location of Simulation Files ─────────────────────────────────────────────
#RAW_DIR     = Path("/lustre/orion/.../gromacs_sims/MED4/.../1tuples")
#RAW_DIR     = Path("/lustre/orion/.../gromacs_sims/HM2/.../1tuples")

# ── Prompt for where your raw simulation folders live ──
raw_path = input("Enter path to RAW_DIR (where your numeric subfolders live, e.g. /lustre/orion/.../1tuples): ").strip()
RAW_DIR  = Path(raw_path)

# ── Location of Reference Files ─────────────────────────────────────────────
#REF_DIR     = RAW_DIR.parent / "REF"
# Ask the user which folder to use as reference
ref_id = input("Enter reference system (REF or numeric ID, e.g. 0003): ").strip()
if ref_id.isdigit():
    # numeric subfolder under RAW_DIR
    REF_DIR = RAW_DIR / ref_id
else:
    # sibling directory next to RAW_DIR (e.g. “REF”)
    REF_DIR = RAW_DIR.parent / ref_id

# ── Location for Processing Files to be created ─────────────────────────────
#POST_BASE   = Path("/lustre/orion/.../gromacs_sims/postprocess/MED4/...")
#POST_BASE   = Path("/lustre/orion/.../gromacs_sims/postprocess/HM2/...")

# ── Prompt for where to write all post-processing ──
post_path = input("Enter path to POST_BASE (will be created if needed, e.g. /lustre/orion/.../postprocess/...): ").strip()
POST_BASE = Path(post_path)
POST_BASE.mkdir(parents=True, exist_ok=True)

# ── Setup Environment ───────────────────────────────────────────────────────
GROMACS     = "gmx_mpi"
#PYTHON      = "python"
PYTHON      = str(Path.home() / "miniconda3/envs/gppenv/bin/python")
SCRIPTS_DIR = Path("/anfhome/shared/qipd/frontier/pca_analysis_scripts")

# ── Per-system script template (plain bash; legacy Flux lines kept as comments) ─────────
MD_TEMPLATE = """#!/bin/bash
# Auto-generated post-processing for {sys}

#export SCRATCH="/lustre/orion/bip258/scratch/${{USER}}"
#module use /ccs/proj/bip258/apps/modulefiles
#module load gromacs/2025.0
source "/anfhome/shared/qipd/gromacs-2025.2/bin/GMXRC"
export PATH="/anfhome/shared/qipd/gromacs-2025.2/bin:${{PATH}}"

export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=7
#export TMPDIR=${{{{SCRATCH}}}}
export TMPDIR=${{SCRATCH:-/tmp}}
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_MAXBACKUP=-1
export UCX_POSIX_USE_PROC_LINK=n
export UCX_TLS=^cma
export UCX_LOG_LEVEL=ERROR
export UCX_RNDV_THRESH=8192
export HWLOC_HIDE_ERRORS=1
export SYCL_CACHE_PERSISTENT=1
export PYTHONPATH="{scripts}:${{PYTHONPATH}}"

# Original: flux resource list
# flux resource list
set -e

RAW="{raw}"
POST_SYS="{post_sys}"
REF_SYS="{ref_sys}"

mkdir -p "$POST_SYS"
cd "$POST_SYS"

# ── 0) CLEAN OLD OUTPUT ─────────────────────────────────────────────────────
#rm -rf figures figures_rmsd_rg figures_rmsf figures_rmsf_secstruct figures_contactmaps representative_structures

# ── 1) Solute extraction: only if not already done ─────────────────────────
if [ ! -f "solute_fit.xtc" ]; then
  echo "Solute extraction not found; running solute-only steps…"

  # 1) build Solute index
  # Original: flux run -n 1 {gmx} make_ndx -f "$RAW/md.tpr" -o index.ndx <<EOF
  {gmx} make_ndx -f "$RAW/md.tpr" -o index.ndx <<EOF
"System" & ! "Water_and_ions"
q
EOF

  sed -i 's/System_&_!Water_and_ions/Solute/' index.ndx

  # 2) Solute-only TPR & trajectory conversions ─────────────────────────────
  # Original: echo Solute | flux run -n 1 {gmx} convert-tpr -s "$RAW/md.tpr" -n index.ndx -o solute.tpr
  echo Solute | {gmx} convert-tpr -s "$RAW/md.tpr" -n index.ndx -o solute.tpr 
  
  # Original: echo Solute | flux run -n 1 {gmx} trjconv ...
  echo Solute | {gmx} trjconv -f "$RAW/md.xtc" -s solute.tpr -n index.ndx -o solute_raw.xtc
  # Original: echo 0 | flux run -n 1 {gmx} trjconv ...
  echo 0 | {gmx} trjconv -f solute_raw.xtc -s solute.tpr -pbc whole -o solute_whole.xtc
  echo 0 | {gmx} trjconv -f solute_whole.xtc -s solute.tpr -pbc nojump -o solute_nojump.xtc
  echo 0 0 | {gmx} trjconv -f solute_nojump.xtc -s solute.tpr -fit rot+trans -o solute_fit.xtc
  echo 0 | {gmx} trjconv -f solute_fit.xtc -s solute.tpr -dump 0 -o solute_0.gro
else
  echo "Solute-only steps already completed; skipping."
fi

# ── 2) Preprocess: GRO→PDB, add chains, fix PBC ────────────────────────────
if [ ! -f solute_fixed.pdb ]; then
  echo "Running preprocessing (GRO→PDB, add chains, fix PBC)…"
  # Original: flux run -n 1 {py} - <<'PY'
  {py} - <<'PY'
import preprocess
#sys.path.insert(0, "{scripts}")
preprocess.gro_to_pdb("solute_0.gro", "solute_0_raw.pdb")
preprocess.add_chain_ids("solute_0_raw.pdb", "solute_0_chains.pdb")
preprocess.fix_pbc("solute_0_chains.pdb",  "solute_fixed.pdb")
PY
else
  echo "Skipping preprocessing (solute_fixed.pdb exists)"
fi

# ── 3) RMSD/Rg (skip if done) ──────────────────────────────────────────────
# Original: flux run -n 1 {py} {scripts}/rmsd_rg.py --top ... --name {sys}
if [ ! -f figures_rmsd_rg/{sys}_rmsd.pdf ]; then
  {py} {scripts}/rmsd_rg.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --skip 1 --window 10 --name {sys} 
else
  echo "Skipping RMSD/Rg"
fi

# ── 4) RMSF/SecStruct (skip if done) ────────────────────────────────────────
# Original: flux run -n 1 {py} {scripts}/rmsf_secstruct.py --top ... --name {sys}
if [ ! -f figures_rmsf_secstruct/{sys}_rmsf_secstruct.pdf ]; then
  {py} {scripts}/rmsf_secstruct.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --all --name {sys}
else
  echo "Skipping RMSF/SecStruct"
fi

# ── 5) Contact maps (skip if done) ──────────────────────────────────────────
# Original: flux run -n 1 {py} {scripts}/contact_maps.py --top ... --name {sys}
if [ ! -f figures_contactmaps/{sys}_contact_maps.pdf ]; then
  {py} {scripts}/contact_maps.py --top solute_fixed.pdb --xtc solute_fit.xtc --cutoff 3.0 --name {sys}
else
  echo "Skipping Contact maps"
fi

echo "Finished post-processing {sys}"
"""

# ── System ordering: REF first (if present), then numerics ─────────────────
def system_ids():
    ids = ["REF"] if REF_DIR.is_dir() else []
    ref_child = REF_DIR.name if REF_DIR.parent == RAW_DIR else None
    others = [
        d.name for d in RAW_DIR.iterdir()
        if d.is_dir() and d.name.isdigit() and d.name != ref_child
    ]
    ids += sorted(others)
    return ids

def write_per_system(sysid):
    raw_path = REF_DIR if sysid == "REF" else RAW_DIR / sysid
    post_sys = POST_BASE / sysid
    ref_sys  = POST_BASE / "REF"
    post_sys.mkdir(parents=True, exist_ok=True)

    script = MD_TEMPLATE.format(
        raw=raw_path, post_sys=post_sys, ref_sys=ref_sys,
        sys=sysid, gmx=GROMACS, py=PYTHON, scripts=SCRIPTS_DIR
    )

    # Plain bash runnable (no sbatch)
    sh = post_sys / "md.sh"
    sh.write_text(script)
    sh.chmod(stat.S_IRWXU)
    return sh.relative_to(POST_BASE)

def main():
    # 1) Write per-system scripts (plain bash; legacy Flux lines kept as comments inside)
    rel_paths = [write_per_system(s) for s in system_ids()]

    # 2) md_fluxjobs.txt (kept for continuity)
    (POST_BASE / "md_fluxjobs.txt").write_text("\n".join(map(str, rel_paths)) + "\n")

    # 3) md_bundle_flux.sh – LOCAL sequential runner:
    #    REF → each system; after each: pairwise PCA in POST_BASE/pairs/<SYS>;
    #    finally: JOIN all pairwise CSVs into POST_BASE/*_ranking.csv and PLOT.
    flux_sh = POST_BASE / "md_bundle_flux.sh"
    flux_sh.write_text(f"""#!/bin/bash

# ==== CONFIG (edit if needed) ==============================================
RUN_BASE="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
RESUME=1                     # 1 = skip jobs that already have .done marker
STOP_ON_FAIL=0               # 1 = stop immediately on any failure
LOG_DIR="$RUN_BASE/__runner_logs"      # runner logs (under POST_BASE)
PAIRS_BASE="$RUN_BASE/pairs"           # pairwise outputs live here (REF+SYS)
# ==========================================================================

mkdir -p "$LOG_DIR" "$PAIRS_BASE"

# Legacy Flux submitter (not used now; kept for reference):
# ---------------------------------------------------------------------------
# flux resource list
# while IFS= read -r line; do
#   dir=$(dirname "$line")
#   base=$(basename "$line")
#   flux batch -N 1 -t 11.8h --gpus-per-slot=8 --cores-per-slot=56 -x --cwd="$PWD/$dir" --output="$base.flux.out" --error="$base.flux.err" "$PWD/$line"
# done < md_fluxjobs.txt
# flux queue drain
# ---------------------------------------------------------------------------

# Load job list
mapfile -t JOBS < md_fluxjobs.txt
TOTAL=${{#JOBS[@]}}
if (( TOTAL == 0 )); then
  echo "No jobs found in md_fluxjobs.txt"; exit 1
fi

# Helpers
ts() {{ date +%Y-%m-%dT%H:%M:%S; }}
sec() {{ date +%s; }}
fmt_dur() {{ local T=$1; printf "%d:%02d:%02d" $((T/3600)) $(((T%3600)/60)) $((T%60)); }}

START=$(sec)
LAST_DUR=0
OK=0
FAIL=0
FAILED_LIST=()

i=0
for rel in "${{JOBS[@]}}"; do
  ((i++))
  dir=$(dirname "$rel")
  base=$(basename "$rel")
  abs_dir="$PWD/$dir"
  job_name="$(basename "$dir")"
  mark="$abs_dir/.done"

  # Resume?
  if (( RESUME==1 )) && [[ -f "$mark" ]]; then
    elapsed=$(( $(sec) - START ))
    echo "[SKIP] [$i/$TOTAL | $(fmt_dur $elapsed) elapsed]  $rel  (.done exists)"
    ((OK++))
  else
    elapsed=$(( $(sec) - START ))
    echo "[RUN ] [$i/$TOTAL | $(fmt_dur $elapsed) elapsed | $(fmt_dur $LAST_DUR) last]  $rel"

    job_log="$abs_dir/${{base%.sh}}.run.log"
    mkdir -p "$abs_dir"
    cd "$abs_dir"

    BEG=$(sec)
    if bash "$base" >"$job_log" 2>&1; then
      LAST_DUR=$(( $(sec) - BEG ))
      echo "[ OK ] $job_name ($(fmt_dur $LAST_DUR))"
      echo "$(ts) OK   $rel $(fmt_dur $LAST_DUR)" | tee -a "$LOG_DIR/runner.log" >/dev/null || true
      touch "$mark"
      ((OK++))
    else
      LAST_DUR=$(( $(sec) - BEG ))
      echo "[FAIL] $job_name ($(fmt_dur $LAST_DUR)) — see $job_log"
      echo "$(ts) FAIL $rel $(fmt_dur $LAST_DUR)" | tee -a "$LOG_DIR/runner.log" >/dev/null || true
      FAILED_LIST+=("$rel")
      ((FAIL++))
      if (( STOP_ON_FAIL==1 )); then
        echo "Stopping due to STOP_ON_FAIL=1"; exit 2
      fi
    fi

    # ---------- Pairwise PCA: REF + this system ----------
    if [[ "$job_name" != "REF" ]]; then
      pair_dir="$PAIRS_BASE/$job_name"
      mkdir -p "$pair_dir"
      # Symlink REF and system so the analysis can *read* from them
      ln -sfn "$RUN_BASE/REF"        "$pair_dir/REF"
      ln -sfn "$RUN_BASE/$job_name"  "$pair_dir/$job_name"

      # Create the real WRITE destination inside the system dir
      sys_pca_dir="$RUN_BASE/$job_name/figures_pca"
      mkdir -p "$sys_pca_dir"

      echo "[PAIR] REF + $job_name  →  $pair_dir (writing to $sys_pca_dir)"
      {PYTHON} {SCRIPTS_DIR}/run_pca_analysis.py \
        --base_dir     "$pair_dir" \
        --ref_dir      "$pair_dir/REF" \
        --output_dir   "$sys_pca_dir" \
        --topology     solute_fixed.pdb \
        --trajectories solute_fit.xtc \
        --step 10 --n_comp 5 --n_clusters 3 --dt_ns 0.1
    fi

    cd - >/dev/null
  fi
done

# Summary
ELAPSED=$(( $(sec) - START ))
echo
echo "==================== SUMMARY ===================="
echo "  Total jobs : $TOTAL"
echo "  Succeeded  : $OK"
echo "  Failed     : $FAIL"
echo "  Elapsed    : $(fmt_dur $ELAPSED)"
if (( FAIL>0 )); then
  echo "  Failed list:"
  for x in "${{FAILED_LIST[@]}}"; do echo "   - $x"; done
fi
echo "================================================="
echo

# ---------- JOIN per-system CSVs from figures_pca ----------
echo "[JOIN] Building {POST_BASE}/PC_score_ranking.csv and combined_score_ranking.csv from */figures_pca …"

JOIN_PC="{POST_BASE}/PC_score_ranking.csv"
JOIN_CB="{POST_BASE}/combined_score_ranking.csv"

# find the first available per-system CSVs to copy headers from
first_pc=""
first_cb=""
shopt -s nullglob
for p in "{POST_BASE}"/*/figures_pca/PC_score_ranking.csv; do first_pc="$p"; break; done
for p in "{POST_BASE}"/*/figures_pca/PC_loadings.csv; do first_cb="$p"; break; done
shopt -u nullglob

if [[ -z "$first_pc" || -z "$first_cb" ]]; then
  echo "No per-system figures_pca CSVs found; nothing to join."; exit 1
fi

# write headers
head -n 1 "$first_pc" > "$JOIN_PC"
echo "system,$(head -n 1 "$first_cb")" > "$JOIN_CB"

# append rows where first column equals the system id (for PC_score_ranking.csv)
shopt -s nullglob
for csv in "{POST_BASE}"/*/figures_pca/PC_score_ranking.csv; do
  sys=$(basename "$(dirname "$(dirname "$csv")")")   # .../<SYS>/figures_pca/...
  awk -F',' -v s="$sys" 'NR>1 && $1==s {{print $0}}' "$csv" >> "$JOIN_PC"
done

# build combined_score_ranking.csv from PC_loadings.csv, prefixing the system column
for csv in "{POST_BASE}"/*/figures_pca/PC_loadings.csv; do
  sys=$(basename "$(dirname "$(dirname "$csv")")")
  awk -F',' -v s="$sys" 'NR>1 {{print s "," $0}}' "$csv" >> "$JOIN_CB"
done
shopt -u nullglob

# ---------- PLOTS over the joined CSVs ----------
echo "[PLOTS] Generating combined plots from joined results…"
{PYTHON} {SCRIPTS_DIR}/plot_pca_selected_systems.py --base_dir {POST_BASE} || true

# Optional additional (non-PCA) plots: expects a selected_systems.csv
if [[ -f "{POST_BASE}/selected_systems.csv" ]]; then
  systems=$(tail -1 {POST_BASE}/selected_systems.csv | tr ',' ' ')
  {PYTHON} {SCRIPTS_DIR}/plot_nonpca_selected_systems.py \\
    --base_dir {POST_BASE} \\
    --output_dir {POST_BASE}/combined_analysis \\
    --systems $systems \\
    --skip 1 --dt 0.10 --mode all || true
fi

echo "[DONE] All steps complete."
""")
    flux_sh.chmod(stat.S_IRWXU)

    # 4) (Legacy) md_bundle.sbatch — preserved here as a reference (NOT written now)
    # -------------------------------------------------------------------------
    # n = len(rel_paths)
    # bundle = POST_BASE / "md_bundle.sbatch"
    # bundle.write_text(f\"\"\"#!/bin/bash
    # #SBATCH --partition=batch
    # #SBATCH --account=bip258
    # #SBATCH --time=12:00:00
    # #SBATCH --nodes={{n}}
    # #SBATCH --job-name=md_bundle
    # #SBATCH --output=md_bundle-%j.out
    # #SBATCH --error=md_bundle-%j.err
    #
    # module use /ccs/proj/bip258/apps/modulefiles
    # module load gromacs/2025.0
    # module load hwloc/2.11.1-gpu # Flux requires a GPU-enabled hwloc
    # module load flux
    #
    # srun -N $SLURM_NNODES -n $SLURM_NNODES -c 56 --gpus-per-node=8 flux start ./md_bundle_flux.sh
    # \"\"\")
    # bundle.chmod(stat.S_IRWXU)
    # -------------------------------------------------------------------------

    print(f"\nGenerated per-system scripts under: {POST_BASE}")
    print(f"Job list: {POST_BASE/'md_fluxjobs.txt'}")
    print(f"Sequential runner: {POST_BASE/'md_bundle_flux.sh'}")
    print("\nTo run:")
    print(f"  cd {POST_BASE}")
    print( "  bash md_bundle_flux.sh")

if __name__ == "__main__":
    main()
