#!/usr/bin/env python
# coding: utf-8
"""
generate_pipelines.py

Writes:
  <POST_BASE>/<ID>/md.sbatch      – Flux script per system
  <POST_BASE>/md_fluxjobs.txt
  <POST_BASE>/md_bundle_flux.sh
  <POST_BASE>/md_bundle.sbatch

This version calls the Python analysis scripts by absolute path from
SCRIPTS_DIR so you can keep them all in one place, and skips solute‐only
operations if they've already been run.
"""

import os, stat
from pathlib import Path
import subprocess

# ── Location of Simulation Files ──────────────────────────────────────────────────────
#RAW_DIR     = Path("/lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/MED4/ae816d19-d80c-4446-9f60-a07cfd79ca3c/1tuples")
#RAW_DIR     = Path("/lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/HM2/ad8ac6d3-c270-4bf7-b8ff-f05c3a2c30d8/1tuples")

# ── Prompt for where your raw simulation folders live ──
raw_path = input("Enter path to RAW_DIR (where your numeric subfolders live, e.g. /lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/HM2/ad8ac6d3-c270-4bf7-b8ff-f05c3a2c30d8/1tuples)): ").strip()
RAW_DIR  = Path(raw_path)

# ── Location of Reference Files ──────────────────────────────────────────────────────
#REF_DIR     = RAW_DIR.parent / "REF"
# Ask the user which folder to use as reference
ref_id = input("Enter reference system (REF or numeric ID, e.g. 0003): ").strip()
if ref_id.isdigit():
    # numeric subfolder under RAW_DIR
    REF_DIR = RAW_DIR / ref_id
else:
    # sibling directory next to RAW_DIR (e.g. “REF”)
    REF_DIR = RAW_DIR.parent / ref_id

# ── Location for Processing Files to be created ──────────────────────────────────────────────────────
#POST_BASE   = Path("/lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/postprocess/MED4/ae816d19-d80c-4446-9f60-a07cfd79ca3c")
#POST_BASE   = Path("/lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/postprocess/HM2/ad8ac6d3-c270-4bf7-b8ff-f05c3a2c30d8")

# ── Prompt for where to write all post-processing ──
post_path = input("Enter path to POST_BASE (will be created if needed, e.g. /lustre/orion/bip258/proj-shared/ptmpsi-production-runs/gromacs_sims/postprocess/HM2/ad8ac6d3-c270-4bf7-b8ff-f05c3a2c30d8): ").strip()
POST_BASE = Path(post_path)
# make sure the output directory exists
POST_BASE.mkdir(parents=True, exist_ok=True)


# ── Setup Environment ──────────────────────────────────────────────────────
GROMACS     = "gmx_mpi"
#PYTHON      = "python"
PYTHON      = str(Path.home() / "miniconda3/envs/gppenv/bin/python")
SCRIPTS_DIR = Path("/lustre/orion/bip258/proj-shared/scripts")

MD_TEMPLATE = """#!/bin/bash
# Auto-generated post-processing for {sys}

export SCRATCH="/lustre/orion/bip258/scratch/${{USER}}"
module use /ccs/proj/bip258/apps/modulefiles
module load gromacs/2025.0


export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=7
export TMPDIR=${{SCRATCH}}
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_MAXBACKUP=-1
export UCX_POSIX_USE_PROC_LINK=n
export UCX_TLS=^cma
export UCX_LOG_LEVEL=ERROR
export UCX_RNDV_THRESH=8192
export HWLOC_HIDE_ERRORS=1
export SYCL_CACHE_PERSISTENT=1
export PYTHONPATH="{scripts}":${{PYTHONPATH}}

flux resource list
module list
set -e

RAW="{raw}"
POST_SYS="{post_sys}"
REF_SYS="{ref_sys}"

mkdir -p "$POST_SYS"
cd "$POST_SYS"

# ── 0) CLEAN OLD OUTPUT ──────────────────────────────────────────────────────
#rm -rf figures figures_rmsd_rg figures_rmsf figures_rmsf_secstruct figures_contactmaps representative_structures

# ── 1) Solute extraction: only if not already done ─────────────────────────
if [ ! -f "solute_fit.xtc" ]; then
  echo "Solute extraction not found; running solute-only steps…"

  # 1) build Solute index
  flux run -n 1 {gmx} make_ndx -f "$RAW/md.tpr" -o index.ndx <<EOF
"System" & ! "Water_and_ions"
q
EOF

  sed -i 's/System_&_!Water_and_ions/Solute/' index.ndx

  # 2) Solute-only TPR & trajectory conversions ──────────────────────────────────────────────────────
  echo Solute | flux run -n 1 {gmx} convert-tpr -s "$RAW/md.tpr" -n index.ndx -o solute.tpr 
  
  echo Solute | flux run -n 1 {gmx} trjconv -f "$RAW/md.xtc" -s solute.tpr -n index.ndx -o solute_raw.xtc
  echo 0 | flux run -n 1 {gmx} trjconv -f solute_raw.xtc -s solute.tpr -pbc whole -o solute_whole.xtc
  echo 0 | flux run -n 1 {gmx} trjconv -f solute_whole.xtc -s solute.tpr -pbc nojump -o solute_nojump.xtc
  echo 0 0 | flux run -n 1 {gmx} trjconv -f solute_nojump.xtc -s solute.tpr -fit rot+trans -o solute_fit.xtc
  echo 0 | flux run -n 1 {gmx} trjconv -f solute_fit.xtc -s solute.tpr -dump 0 -o solute_0.gro
else
  echo "Solute-only steps already completed; skipping."
fi

# ── 2) Preprocess: GRO→PDB, add chains, fix PBC ──────────────────────────────────────────────────────
if [ ! -f solute_fixed.pdb ]; then
  echo "Running preprocessing (GRO→PDB, add chains, fix PBC)…"
  flux run -n 1 {py} - <<'PY'
import preprocess
#sys.path.insert(0, "{scripts}")
preprocess.gro_to_pdb("solute_0.gro", "solute_0_raw.pdb")
preprocess.add_chain_ids("solute_0_raw.pdb", "solute_0_chains.pdb")
preprocess.fix_pbc("solute_0_chains.pdb",  "solute_fixed.pdb")
PY
else
  echo "Skipping preprocessing (solute_fixed.pdb exists)"
fi

# ── 3) RMSD/Rg (skip if done) ──────────────────────────────────────────────────────
#flux run -n 1 {py} {scripts}/rmsd_rg.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --skip 1 --window 10 --name {sys}
if [ ! -f figures_rmsd_rg/{sys}_rmsd.pdf ]; then
  flux run -n 1 {py} {scripts}/rmsd_rg.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --skip 1 --window 10 --name {sys} 
else
  echo "Skipping RMSD/Rg"
fi

# ── 4) RMSF/SecStruct (skip if done) ──────────────────────────────────────────────────────
#flux run -n 1 {py} {scripts}/rmsf_secstruct.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --all --name {sys}
if [ ! -f figures_rmsf_secstruct/{sys}_rmsf_secstruct.pdf ]; then
  flux run -n 1 {py} {scripts}/rmsf_secstruct.py --top solute_fixed.pdb --xtc solute_fit.xtc --ref_top "$REF_SYS/solute_fixed.pdb" --ref_xtc "$REF_SYS/solute_fit.xtc" --all --name {sys}
else
  echo "Skipping RMSF/SecStruct"
fi

# ── 5) Contact maps (skip if done) ──────────────────────────────────────────────────────
#flux run -n 1 {py} {scripts}/contact_maps.py --top solute_fixed.pdb --xtc solute_fit.xtc --cutoff 3.0 --name {sys}
if [ ! -f figures_contactmaps/{sys}_contact_maps.pdf ]; then
  flux run -n 1 {py} {scripts}/contact_maps.py --top solute_fixed.pdb --xtc solute_fit.xtc --cutoff 3.0 --name {sys}
else
  echo "Skipping Contact maps"
fi

echo "Finished post-processing {sys}"
"""

def system_ids():
    # 1) Always start with “REF” if that folder exists
    ids = ["REF"] if REF_DIR.is_dir() else []

    # 2) If REF_DIR lives under RAW_DIR, note its name so we can skip it
    ref_child = REF_DIR.name if REF_DIR.parent == RAW_DIR else None

    # 3) List all numeric subfolders except the one we’re using as REF
    others = [
        d.name
        for d in RAW_DIR.iterdir()
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
        raw=raw_path,
        post_sys=post_sys,
        ref_sys=ref_sys,
        sys=sysid,
        gmx=GROMACS,
        py=PYTHON,
        scripts=SCRIPTS_DIR
    )

    sb = post_sys / "md.sbatch"
    sb.write_text(script)
    sb.chmod(stat.S_IRWXU)
    return sb.relative_to(POST_BASE)

def main():
    # 1) Write per-system Flux SBATCH scripts
    rel_paths = [write_per_system(s) for s in system_ids()]

    # 2) md_fluxjobs.txt
    (POST_BASE / "md_fluxjobs.txt").write_text("\n".join(map(str, rel_paths)) + "\n")

    # 3) md_bundle_flux.sh – submit all, wait, then merge & plot
    flux_sh = POST_BASE / "md_bundle_flux.sh"
    flux_sh.write_text(f"""#!/bin/bash
flux resource list
while IFS= read -r line; do
  dir=$(dirname "$line")
  base=$(basename "$line")
  flux batch -N 1 -t 11.8h --gpus-per-slot=8 --cores-per-slot=56 -x --cwd="$PWD/$dir" --output="$base.flux.out" --error="$base.flux.err" "$PWD/$line"
done < md_fluxjobs.txt

# wait for all jobs to finish
flux queue drain

if [ ! -f combined_score_ranking.csv ]; then
  {PYTHON} {SCRIPTS_DIR}/run_pca_analysis.py --base_dir {POST_BASE} --ref_dir {POST_BASE}/REF --output_dir {POST_BASE} --topology solute_fixed.pdb --trajectories solute_fit.xtc --step 1 --n_comp 5 --n_clusters 3 --profile --dt_ns 0.1 
  # ── Single Script ──
# {PYTHON} {SCRIPTS_DIR}/unified_global_pca_analysis.py --base_dir {POST_BASE} --ref_dir {POST_BASE}/REF --output_dir {POST_BASE} --topology solute_fixed.pdb --trajectories solute_fit.xtc --step 1 --n_comp 5 --n_clusters 3 --profile --dt_ns 0.1 
else
  echo "Global PCA analysis already done; skipping."
fi


{PYTHON} {SCRIPTS_DIR}/plot_pca_selected_systems.py --base_dir {POST_BASE}

# grab the final systems line
systems=$(tail -1 {POST_BASE}/selected_systems.csv | tr ',' ' ')

# feed into your multi-system plotting
{PYTHON} {SCRIPTS_DIR}/plot_nonpca_selected_systems.py --base_dir {POST_BASE} --output_dir {POST_BASE}/combined_analysis --systems $systems --skip 1 --dt 0.10 --mode all

########################################################################################################################

""")
    flux_sh.chmod(stat.S_IRWXU)

    # 4) md_bundle.sbatch (for legacy Slurm)
    n = len(rel_paths)
    bundle = POST_BASE / "md_bundle.sbatch"
    bundle.write_text(f"""#!/bin/bash
#SBATCH --partition=batch
#SBATCH --account=bip258
#SBATCH --time=12:00:00
#SBATCH --nodes={n}
#SBATCH --job-name=md_bundle
#SBATCH --output=md_bundle-%j.out
#SBATCH --error=md_bundle-%j.err

module use /ccs/proj/bip258/apps/modulefiles
module load gromacs/2025.0
module load hwloc/2.11.1-gpu # Flux requires a GPU-enabled hwloc to see the GPUs
module load flux

srun -N $SLURM_NNODES -n $SLURM_NNODES -c 56 --gpus-per-node=8 flux start ./md_bundle_flux.sh
""")
    bundle.chmod(stat.S_IRWXU)

    # 5) ensure analysis directories exist
    #(POST_BASE / "score_analysis").mkdir(exist_ok=True)
    #(POST_BASE / "combined_analysis").mkdir(exist_ok=True)

    print(f"Change directory to: {POST_BASE} and sbatch md_bundle.sbatch")

if __name__ == "__main__":
    main()
