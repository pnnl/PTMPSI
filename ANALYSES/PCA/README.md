# PCA Scripts

This directory contains the driver script and SLURM launchers used to run the per-ensemble
structural PCA on `<ENSEMBLE>_NAD/postprocess_final/` trajectories, including the primary
structural-only run and the full-feature sensitivity/validation reruns used to check it.

`run_pca_analysis.py` extracts per-system structural/dynamic descriptors, runs PCA (with a
multi-threshold sensitivity sweep), KMeans clustering, MFPT-based kinetic-network analysis, and
ranks non-REF systems by variance-normalized PCA distance to REF. It depends on the
`pca_config`, `pca_selection_features`, `pca_clustering`, and `pca_plotting` modules from the
general PTMPSI PCA package (imported via a hardcoded `sys.path` entry — see Notes).

## Pipeline order

1. **`run_hm2_pca_features.sbatch`** / **`run_med4_pca_features.sbatch`** — First pass per
   ensemble: extracts and caches per-system PCA features (`feature_outputs/`) into
   `suman_pca_fixed_90pct`. Run once after trajectory corrections/postprocessing land in
   `postprocess_final`; later runs reuse this feature cache instead of re-extracting.
2. **`run_hm2_pca.sbatch`** / **`run_med4_pca.sbatch`** — **Primary run.** Reuses the cached
   features from step 1, excludes `PTM_`-prefixed features (`--exclude_feature_prefixes PTM_`),
   and writes the canonical structural-only PCA outputs to `suman_pca_structural_only_90pct`
   (`pca_distance_to_ref_90pct_ranked.csv`, `PC_loadings.csv`, `explained_variance.csv`,
   `full_scree_variance.csv`, `pca_sensitivity_*`, per-system cluster/KTN/top-5-PC figures). This
   is the PCA referenced throughout the project and manuscript.
3. **`run_hm2_pca_fullfeature_95pct.sbatch`** / **`run_med4_pca_fullfeature_95pct.sbatch`** and
   **`run_hm2_pca_fullfeature_sensitivity.sbatch`** / **`run_med4_pca_fullfeature_sensitivity.sbatch`**
   — Sensitivity/validation reruns of the same script with PTM features *included* (no
   `--exclude_feature_prefixes`), at a fixed 95% variance threshold and across an 80/90/95%
   sensitivity sweep respectively. These are comparison runs, not the primary structural-only
   path — used to check that the structural-only exclusion doesn't distort the ranking.

Each sbatch script builds a `pca_input_ref_plus_numeric/` symlink tree (REF + all numeric system
dirs with both `solute_fixed.pdb` and `solute_fit.xtc`, skipping `0003` since it duplicates REF)
before invoking `run_pca_analysis.py`, and verifies expected output files exist afterward.

## `run_pca_analysis.py` key CLI arguments

| Argument | Purpose |
|---|---|
| `--base_dir` / `--ref_dir` / `--output_dir` | Input system-symlink tree, REF directory, output directory |
| `--topology` / `--trajectories` | Per-system structure/trajectory filenames (`solute_fixed.pdb`, `solute_fit.xtc`) |
| `--step` | Frame stride |
| `--pca_variance_threshold` | Primary cumulative-variance cutoff for PC retention |
| `--pca_sensitivity_thresholds` | Comma-separated thresholds (e.g. `0.80,0.90,0.95`) for the sensitivity sweep and rank-stability check |
| `--exclude_feature_prefixes` | Feature-name prefixes to drop before PCA (e.g. `PTM_` for the structural-only run) |
| `--n_clusters` | KMeans cluster count |
| `--per_system_jobs` | Joblib parallelism per system |
| `--dt_ns` | Trajectory time step in ns, for kinetic-network/MFPT calculations |

Run `run_pca_analysis.py --help` for the complete argument list (radii/cutoffs for PTM, NAD, and
interface contacts, DSSP window, plotting flags, profiling).

## Dependencies

- Python environment `gppenv` (conda), with MDAnalysis, scikit-learn, joblib, pandas, matplotlib
- The `pca_config` / `pca_selection_features` / `pca_clustering` / `pca_plotting` modules from the
  PTMPSI PCA package

## Notes

- Paths are hardcoded to this project's NERSC Perlmutter scratch layout, including the
  `sys.path` entry that locates the PTMPSI PCA package modules. Reusing these scripts elsewhere
  requires updating those paths.
- Only the primary + sensitivity/validation launchers are included here; the general-purpose
  PTMPSI PCA package (feature extraction, clustering, and plotting library code) lives separately
  and is not duplicated in this directory.
