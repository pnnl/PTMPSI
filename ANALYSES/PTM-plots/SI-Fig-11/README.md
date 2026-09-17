# SI Fig 11 — Combined 8-Panel MED4/HM2 Figure

This directory contains the plotting script for the combined 8-panel supplementary figure
comparing MED4 (host) and HM2 (viral CP12) ensembles across four structural/dynamic metrics.

## Script: `plot_combined_8panel.py`

Publication-quality figure, 2 columns (MED4, HM2) × 4 rows:

| Row | Metric | Panel (MED4 / HM2) |
|---|---|---|
| 1 | RMSD vs. time | A / B |
| 2 | Per-residue RMSF heatmap | C / D |
| 3 | Radius of gyration (Rg) vs. time | E / F |
| 4 | Secondary-structure composition | G / H |

The script reads only pre-computed cache `.npz` files — it does not touch trajectories directly.

### Input cache structure

`--med4_cache` / `--hm2_cache` must each point to a directory containing:

```
rmsd_*.npz         # systems, x_list (time), y_list (RMSD) per system
rg_*.npz           # same structure, Rg
rmsf_*.npz         # arr (per-residue RMSF matrix), labs0 (residue/column labels)
secstruct_*.npz    # valid_systems, frac_arr (per-system secondary-structure fractions)
```

The script globs for the first `<metric>_*.npz` match in each cache directory.

### Usage

```bash
python plot_combined_8panel.py \
  --med4_cache MED4_NAD/postprocess_final/cache_matrices \
  --hm2_cache  HM2_NAD/postprocess_final/cache_matrices \
  --out_dir combined_8panel
```

**Key arguments:**

| Argument | Default | Description |
|---|---|---|
| `--med4_cache` / `--hm2_cache` | required | Cache directories with the `.npz` files above |
| `--out_dir` | `combined_8panel` | Output directory |
| `--nonref_alpha` | `0.25` | Line opacity for non-REF systems |
| `--nonref_lw` | `0.8` | Line width for non-REF systems |
| `--ref_lw` | `1.0` | Line width for the REF system |

### Output

Writes `combined_8panel.pdf` and `combined_8panel.tiff` (600 dpi) to `--out_dir`.

## Dependencies

- Python 3 with `numpy`, `matplotlib`

## Notes

- Paths in the usage example are hardcoded to this project's NERSC Perlmutter scratch layout;
  point `--med4_cache`/`--hm2_cache` at the corresponding cache directories in your own layout.
