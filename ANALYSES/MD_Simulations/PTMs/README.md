# PTM Combinatorial System Generator

Generates the full combinatorial PTM system set for the GAP2·CP12·PRK dark complex: every
combination of PTM state across the GAP2 (4 sites × 8 chains), PRK (2 sites × 4 chains), and CP12
(disulfide state × 4 chains) positions, yielding the full 1,012-system set (46 GAP2 states × 11
PRK states × 2 CP12 states) submitted via `ptmpsi`'s `Protein.gen_ptm_files`.

## Scripts

- **`generate_ptm_combinations_med4.py`** — MED4 (host) version: input PDB
  `med4_dark_complex_modeller_CYX_250101.pdb`, `NAD_DMR.pdb` cofactor, CP12 cysteines 60/69,
  output path `gromacs_sims/MED4_NAD`.
- **`generate_ptm_combinations_hm2.py`** — HM2 (viral) version: input PDB
  `hm2_dark_complex_modeller_CYX_250101.pdb`, `NAD_HM2.pdb` cofactor, CP12 cysteines 56/65,
  output path `gromacs_sims/HM2_NAD`.

Both build `full_combo_list` (one tuple of `(chain:residue, chemistry)` pairs per system) by
enumerating:
- **GAP2** — reduction (`None`) or nitrosylation/sulfenylation/glutathionylation
  (SNO/SOH/SSG) at each of Cys20/78/157/161, across all non-empty subsets of the 4 sites, applied
  identically across all 8 GAP2 chains.
- **PRK** — reduction, disulfide (`SS`), or SNO/SOH/SSG at Cys201/225, across all non-empty
  subsets of the 2 sites, applied identically across all 4 PRK chains.
- **CP12** — reduced or disulfide-bonded state, applied identically across all 4 CP12 chains.
  (CP12 carries no chemistry PTM in this model — only the redox/disulfide axis.)

The resulting combinatorial list is submitted as one batch via `protein.gen_ptm_files(...)` with
`bundling=True`.

## Usage

```bash
python generate_ptm_combinations_med4.py
python generate_ptm_combinations_hm2.py
```

No CLI arguments; all paths, machine/account/partition settings, and simulation parameters
(`lenmd`, `lennvt`, `lennpt`, `time`) are module-level constants (see Notes).

## Job submission / bundling scripts (Frontier)

`ptmpsi`'s `gen_ptm_files(...)` (called above) writes one run directory per generated system plus
a `fluxjobs.txt` manifest. The following scripts dispatch and run all 1,012 per-system GROMACS
jobs on OLCF Frontier, one compute node per system, using a nested Flux scheduler launched inside
a single large SLURM allocation:

- **`bundle.sbatch`** — Requests all 1,012 nodes in one SLURM job, loads the GROMACS/hwloc/Flux
  modules, starts a Flux instance spanning the allocation, and runs `bundle_flux.sh` inside it.
- **`bundle_flux.sh`** — Reads `fluxjobs.txt` (one per-system job script path per line) and submits
  each as a 1-node Flux batch job (23.8 h walltime, 8 GPUs, 56 cores), then waits for the whole
  queue to drain. Afterward it reads `md_fluxjobs.txt` and each system's `estruns.txt` (estimated
  additional runs needed to reach the target trajectory length) to determine how many more rounds
  are required, and chains further `md_bundle.sbatch` submissions via `sbatch --dependency=afterok`
  until every system's production run is complete.
- **`md_bundle.sbatch`** — Same shape as `bundle.sbatch` (1,012 nodes), starts Flux and runs
  `md_bundle_flux.sh`; this is the job resubmitted by `bundle_flux.sh` for each additional round.
- **`md_bundle_flux.sh`** — Reads `md_fluxjobs.txt` and submits each system's MD-continuation job
  as a 1-node Flux batch job, then drains the queue.

### Usage

```bash
sbatch bundle.sbatch
```

Requires `fluxjobs.txt` and (once the initial round completes) `md_fluxjobs.txt` with per-system
`estruns.txt` files to already exist in the working directory, as produced by the generator step
above.

## Dependencies

- Python 3 with `ptmpsi` (this project's PTM-application package) installed and importable
- The corresponding input structure and cofactor PDB in the working directory (see Scripts above)
- GROMACS 2025.0, a GPU-enabled `hwloc`, and Flux (all loaded as OLCF Frontier modules)

## Notes

- Unlike most scripts in this repo, these were written for **OLCF Frontier**, not NERSC
  Perlmutter — the `gromacs_directory`, `machine`, `account`, and `partition` constants are
  Frontier-specific and would need updating to run on Perlmutter.
