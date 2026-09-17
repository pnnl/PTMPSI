# 0000/0003 Replicate Simulation Scripts

Scripts used to run the MED4 and HM2 dark-complex replicate simulations for the focused CP12
analysis: system 0003 (oxidized reference, CP12 and PRK regulatory cysteines modeled as
disulfides) and system 0000 (reduced reference, same cysteines as free thiols, no additional
PTMs), giving four conditions (MED4 ox, MED4 red, HM2 ox, HM2 red) with three independent
replicates each (12 runs total).

## Replicate design

Each replicate starts from the same solvated structure per condition and diverges through its own
full equilibration chain before production:

```
start.pdb
  -> genrestr (position-restraint .itp files, force constants 1000..1)
  -> minim1000 -> minim500 -> minim100 -> minim50 -> minim10 -> minim5 -> minim1  (restrained minimization)
  -> heating   (gen-vel = yes, unique gen-seed per replicate)   <- independence begins here
  -> npt       (continuation = yes)
  -> minim     (unrestrained minimization)
  -> fheating  (gen-vel = yes, second unique gen-seed per replicate)
  -> fnpt      (continuation = yes)
  -> md        (production, continuation = yes, gen-vel = no, 100 ns / 50,000,000 steps @ 2 fs)
  -> postprocessing (solute extraction, PBC fix, fit, RMSD/Rg/RMSF)
```

Each replicate is assigned two unique random seeds (heating and final heating), so thermostat,
barostat, and box dimensions are all re-equilibrated independently per replicate rather than
inherited from a shared equilibration.

Production parameters (integrator=md, dt=0.002, T=300 K/V-rescale, P=1 bar/Parrinello-Rahman,
PME rcoulomb=1.2, force-switch VDW rvdw=1.2, LINCS order=4/iter=1, DispCorr=EnerPres) are held
identical across all 12 replicates.

## Scripts

- **`submit_all_replicates.sh`** — Prepares (and optionally submits, with `--submit`) all 12
  replicate directories: extracts starting coordinates, symlinks topology and force field, builds
  the custom `Protein_NAD`/`non_Protein_NAD` index groups, copies in the mdp set with the
  per-replicate `gen-seed` values substituted into `heating.mdp` and `fheating.mdp`, and copies in
  `run_replicate_pipeline.sh` + `run_replicate.sbatch`.
- **`run_replicate.sbatch`** — SLURM submission script for one replicate (1 node, 4 GPUs, 30 h);
  invokes `run_replicate_pipeline.sh`.
- **`run_replicate_pipeline.sh`** — The full per-replicate pipeline: restrained minimization chain,
  heating, NPT, unrestrained minimization, final heating, final NPT, production, and
  postprocessing (trajectory PBC/fit correction, chain-ID assignment, RMSD/Rg and RMSF/secondary-
  structure figure generation). Stage-completion is checked by output-file existence, so it's
  resumable; production also resumes from a checkpoint (`md.cpt`) if present.

## MDP files

`minim1000.mdp` ... `minim1.mdp`, `minim.mdp` (restrained/unrestrained minimization chain),
`heating.mdp`, `npt.mdp`, `fheating.mdp`, `fnpt.mdp` (two-pass heating/NPT equilibration), and
`md.mdp` (100 ns production). `gen-seed` in `heating.mdp` and `fheating.mdp` is a placeholder
substituted per replicate by `submit_all_replicates.sh`.

## Usage

```bash
./submit_all_replicates.sh          # prepare all 12 rep directories only
./submit_all_replicates.sh --submit # prepare and submit all 12 jobs
```

## Dependencies

- GROMACS (MPI build for production, `gmx_mpi`; a separate `gmx` binary for setup/postprocessing)
- Python 3 with the `preprocess`, `rmsd_rg`, and `rmsf_secstruct` modules from this project's PCA
  package (used only in the postprocessing stage)

## Notes

- Paths are hardcoded to this project's NERSC Perlmutter scratch layout. Reusing these scripts
  elsewhere requires updating the topology-source, force-field, and reference-structure path
  constants in `submit_all_replicates.sh` and `run_replicate_pipeline.sh`.
