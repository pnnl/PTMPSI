# SI Fig 14 — Overlap-12 Kinetics Dumbbell (MED4 vs P-HM2)

Dumbbell plot comparing kinetic metrics between MED4 (host) and P-HM2 (viral) for the 12
"overlap-12" systems — the systems that fall in the top-5%-per-ensemble PCA selection under
*both* ensembles. Since these 12 carry identical GAP2+PRK PTM combinations in both CP12
backgrounds, they form matched pairs differing only in CP12 origin (host vs. viral), isolating the
effect of the CP12 substitution on kinetics.

## Script: `plot_overlap12_kinetics_dumbbell.py`

Self-contained — no other local script dependencies. Writes
`fig_overlap12_kinetics_dumbbell.{pdf,png,tiff}` to `audit_figure_refresh/kinetics_vs_pca/`.

Two panels, one row per system, dots joined by a connector line:

- **Left — MFPT:** log2(mean first-passage time / REF), MED4 vs. P-HM2
- **Right — τ (tau):** log2(intra-cluster residence time / REF), MED4 vs. P-HM2

Systems are ordered top-to-bottom by P-HM2 τ shift (most-accelerated first). A left-pointing
connector in the τ panel indicates P-HM2 acceleration relative to MED4 for that system.

## Input data

`combined_kinetics_ref_relative.csv`
— precomputed REF-relative log2 ratios (`log2_mean_mfpt_ns_ratio_vs_ref`,
`log2_mean_tau_ns_ratio_vs_ref`) per system per ensemble, sourced from the canonical
`kinetic_summary_structural_only_90pct` kinetics summary.

The 12 overlap systems are hardcoded in the script: `0273 0332 0334 0620 0925 0928 0929 0935 0956
0997 1005 1006`.

## Usage

```bash
python plot_overlap12_kinetics_dumbbell.py
```

No CLI arguments; input/output paths are set as module-level constants (see Notes).

## Dependencies

- Python 3 with `numpy`, `pandas`, `matplotlib`

## Notes

- Paths are hardcoded to this project's NERSC Perlmutter scratch layout. Reusing this script
  elsewhere requires updating the `SCRATCH`/`KIN`/`OUT` path constants.
