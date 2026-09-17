# RMSF-bin

## `interface_rmsf_by_region.py`

Computes and plots CP12 Cα RMSF broken down by interface region (PRK-binding, GAP2-binding,
Non-interface) across all four dark-complex conditions (MED4 ox/red, HM2 ox/red) and three
replicates each.

### Outputs

| File | Description |
|---|---|
| `section31_strengthening/interface_rmsf.csv` | Per-residue RMSF table (strain, redox, rep, chain, resid, region, rmsf) |
| `FINAL_FIGURES/interface_rmsf_by_region.pdf/png/tiff` | Bar chart of mean RMSF ± SEM per region, oxidized and reduced panels side by side |

Summary statistics and Welch t-tests (MED4 vs HM2, per redox state and region) are printed to
stdout.

### Usage

```bash
python interface_rmsf_by_region.py
```

No arguments. Edit the constants at the top of the script to change paths or settings:

| Constant | Default | Description |
|---|---|---|
| `BASE_DIR` | `/path/to/simualtions` | Root directory containing per-strain simulation output |
| `OUT_DIR` | `BASE_DIR/CP12_focused/section31_strengthening` | CSV output directory |
| `FIG_DIR` | `BASE_DIR/CP12_focused/FINAL_FIGURES` | Figure output directory |
| `REPS` | `["rep1","rep2","rep3"]` | Replicate subdirectory names |
| `STRIDE` | `10` | Frame stride for RMSF calculation |

### Expected input layout

```
<BASE_DIR>/<strain_dir>/<sys_id>/<rep>/solute_fixed.pdb
<BASE_DIR>/<strain_dir>/<sys_id>/<rep>/solute_fit.xtc
```

where `strain_dir` is `MED4_NAD` or `HM2_NAD`, `sys_id` is `0003` (ox) or `0000` (red), and
`rep` is one of the values in `REPS`. Missing topology/trajectory files are skipped with a
warning.

### Dependencies

- Python 3
- `MDAnalysis`
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
