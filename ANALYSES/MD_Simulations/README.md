# MD Simulations

Scripts and input files used to run all GROMACS molecular dynamics simulations for the
GAP2·CP12·PRK dark complex: the focused replicate set (reduced and oxidized reference systems)
and the full 1,012-system PTM combinatorial set. All production runs are 100 ns.

## Subdirectories

### `0000_0003/`

Scripts and MDP files for the **reduced (0000) and oxidized (0003) reference replicate
simulations** across both host (MED4) and viral (HM2) dark complexes — four conditions
(MED4 ox, MED4 red, HM2 ox, HM2 red) × 3 replicates = 12 runs total. Each replicate undergoes
its own full equilibration chain (restrained minimization cascade → two-pass heating/NPT →
unrestrained minimization → production) with unique random seeds to ensure independent sampling.
Runs are targeted at NERSC Perlmutter.

Key scripts:

- **`submit_all_replicates.sh`** — Prepares (and optionally submits with `--submit`) all 12
  replicate directories: sets up topology symlinks, builds index groups, substitutes per-replicate
  `gen-seed` values into `heating.mdp` / `fheating.mdp`, and copies in the pipeline scripts.
- **`run_replicate.sbatch`** — SLURM submission script for a single replicate (1 node, 4 GPUs,
  30 h); calls `run_replicate_pipeline.sh`.
- **`run_replicate_pipeline.sh`** — Full per-replicate pipeline from restrained minimization
  through production and postprocessing (PBC fix, fit, RMSD/Rg/RMSF/secondary-structure figures).
  Stage-completion is checked by output-file existence, making it safely resumable.

MDP files cover the complete equilibration chain: `minim1000.mdp` … `minim1.mdp` (restrained
minimization cascade), `minim.mdp` (unrestrained), `heating.mdp`, `npt.mdp`, `fheating.mdp`,
`fnpt.mdp`, and `md.mdp` (100 ns production, 50,000,000 steps at 2 fs).

### `PTMs/`

Scripts and input structures for generating and running the **full 1,012-system PTM
combinatorial set** (46 GAP2 states × 11 PRK states × 2 CP12 states) for both MED4 and HM2.
All systems are submitted via `ptmpsi`'s `Protein.gen_ptm_files` with `bundling=True` and
dispatched at scale on OLCF Frontier using a nested Flux scheduler inside a single large SLURM
allocation.

Key scripts:

- **`generate_ptm_combinations_med4.py`** — Enumerates every PTM combination for MED4 and
  calls `gen_ptm_files`; uses `med4_dark_complex_modeller_CYX_250101.pdb` + `NAD_DMR.pdb`.
- **`generate_ptm_combinations_hm2.py`** — HM2 counterpart; uses
  `hm2_dark_complex_modeller_CYX_250101.pdb` + `NAD_HM2.pdb`.
- **`bundle.sbatch`** / **`bundle_flux.sh`** — Initial SLURM+Flux dispatch: requests all
  1,012 nodes in one allocation, submits each system's setup/equilibration job as a 1-node
  Flux batch job, then chains further production rounds via `md_bundle.sbatch` until every
  system reaches the target trajectory length.
- **`md_bundle.sbatch`** / **`md_bundle_flux.sh`** — Production-continuation dispatch; reads
  `md_fluxjobs.txt` and per-system `estruns.txt` to determine and submit the required number
  of additional MD rounds.

## Usage

```bash
# --- 0000/0003 replicates (Perlmutter) ---
cd 0000_0003
./submit_all_replicates.sh          # prepare all 12 replicate directories only
./submit_all_replicates.sh --submit # prepare and submit all 12 SLURM jobs

# --- PTM combinatorial set (Frontier) ---
cd PTMs
python generate_ptm_combinations_med4.py   # generates MED4 run dirs + fluxjobs.txt
python generate_ptm_combinations_hm2.py    # generates HM2 run dirs + fluxjobs.txt
sbatch bundle.sbatch                        # submit all 1,012 systems
```

## Dependencies

- GROMACS (MPI build `gmx_mpi` for production; `gmx` for setup/postprocessing)
- Python 3 with `ptmpsi` installed (PTM combinatorial set only)
- SLURM (both clusters); Flux scheduler (Frontier only)
- Python 3 with the `preprocess`, `rmsd_rg`, and `rmsf_secstruct` modules from this project's
  PCA package (used in the 0000/0003 postprocessing stage)

## Notes

- `0000_0003/` scripts are hardcoded to NERSC Perlmutter scratch paths; update the
  topology-source, force-field, and reference-structure path constants in
  `submit_all_replicates.sh` and `run_replicate_pipeline.sh` before running elsewhere.
- `PTMs/` scripts are written for OLCF Frontier; the `machine`, `account`, `partition`, and
  `gromacs_directory` constants are Frontier-specific and require updating for Perlmutter.
