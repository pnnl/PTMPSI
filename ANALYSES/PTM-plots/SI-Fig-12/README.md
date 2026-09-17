# SI Fig 12 — PTM Chain Count vs. PCA Distance (Combined 2+3 Panel)

Combined figure: a full-width boxplot of total PTM count vs. PCA distance-to-REF (panel a) stacked
above a 3-panel per-subunit breakdown (CP12 | PRK | GAP2, panel b), for MED4 vs. P-HM2.

## Scripts

### `plot_ptm_dist_box_2plus3.py` (entry point)

Builds the combined 2-row figure and writes
`fig_ptm_chaincount_vs_distance_box_combined_2plus3.{pdf,png,tiff}`. It does **not** load data or
draw boxes itself — it imports `plot_ptm_chaincount_vs_distance.py` as a module and reuses its
`load()`, `count_vals()`, and `_draw_combined_box()` functions, then applies styling
on top (font sizes, box/whisker linewidths, opaque box edges, zero-baseline y-axis, compact
Spearman-ρ annotations).

**This script cannot run standalone** — it does `sys.path.insert(0, ".../PRK_SSG_RERUN/ensemble")`
and `import plot_ptm_chaincount_vs_distance as P`, so both files must be run from the same
directory (or the `sys.path.insert` line updated to point at wherever
`plot_ptm_chaincount_vs_distance.py` lives).

### `plot_ptm_chaincount_vs_distance.py` (dependency)

Generates the full family of PTM-count-vs-PCA-distance figures this repo's `audit_figure_refresh`
output is built from: per-grain (total / GAP2 / PRK / CP12) scatter, box, violin, and MED4/P-HM2
paired "combined" box/violin, plus a 3-panel combined layout and a Spearman-ρ forest plot summary.
Only its `load()`, `count_vals()`, `_draw_combined_box()`, and shared constants (`ENS_COLORS`,
`ENS_DISPLAY`, `PANELS`, `DIST_COL`, `OUT`) are used by `plot_ptm_dist_box_2plus3.py`; running this
script's own `main()` regenerates the standalone figure family (not the SI-Fig-12 combined figure).

**Counts (chain-level):**
- `count_total` — all 120 chain-residue PTM slots (8 GAP2 chains × 4 sites + 4 PRK chains × 2
  sites), range 0–40.
- `count_gap2` — GAP2-only PTM slots, multiples of 8, range 0–32.
- `count_prk` — PRK-only PTM slots, multiples of 4, range 0–8.
- `count_cp12` — CP12 disulfide-bonded chain count (0 or 4; the CP12 redox axis is all-or-nothing
  across its 4 chains, and is not part of the 120-feature PTM set).

## Input data

Both scripts read, per ensemble (`MED4`, `HM2`), from `<ENSEMBLE>_NAD/postprocess_final/`:

- `pca_distance_to_ref_90pct_ranked.csv` — PCA rank/distance table
- `ptm_occupancy_binary_chainlevel.csv` — per-chain PTM occupancy
- `cp12_redox_summary.csv` — per-system CP12 disulfide count

## Usage

```bash
python plot_ptm_dist_box_2plus3.py
```

No CLI arguments; input/output paths are set as module-level constants (see Notes).

## Dependencies

- Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`

## Notes

- Paths are hardcoded to this project's NERSC Perlmutter scratch layout (both the input data
  paths in `plot_ptm_chaincount_vs_distance.py` and the `sys.path.insert`/output directory in
  `plot_ptm_dist_box_2plus3.py`). Reusing these scripts elsewhere requires updating those paths.
