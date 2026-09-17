# SI Fig 13 — Combined PTM Enrichment Figure

Combined figure: PTM-type fold-enrichment barplot (panel a) stacked above the per-site PTM
enrichment dotplot (panel b), faceted MED4 | P-HM2.

## Script: `plot_ptm_enrichment_combined.py`

Self-contained — no other local script dependencies. Writes
`fig_ptm_enrichment_combined.{pdf,png,tiff}`.

**Panel (a) — fold enrichment by PTM type:**
log2(mean chain-level PTM count in the top-5% systems / mean count in the background), per PTM
type (SSG/SOH/SNO, displayed names for the Amber CGL/CSO/SNC codes), faceted MED4 | P-HM2, with
BH-q-value significance stars per bar.

**Panel (b) — per-site enrichment dotplot:**
one row per PTM site (`<type>-<subunit> Cys<residue>`), x = PTM ratio in the top-5% set
(count / 51), dot size = count in top-5%, dot color = −log10(BH q-value), faceted MED4 | P-HM2
with a shared colorbar and a count-size legend.

## Input data

Per ensemble (`MED4`, `HM2`), reads from `<ENSEMBLE>_NAD/postprocess_final/ptm_enrichment/`:

- `ptm_occupancy_binary.csv` — per-system chain-level PTM occupancy + `in_top5pct` flag (panel a)
- `enrichment_by_type.csv` — BH-adjusted q-values per PTM type (panel a)
- `enrichment_per_feature.csv` — per-site top-5% counts and BH-adjusted q-values (panel b)

## Usage

```bash
python plot_ptm_enrichment_combined.py
```

No CLI arguments; input/output paths are set as module-level constants (see Notes).

## Dependencies

- Python 3 with `numpy`, `pandas`, `matplotlib`

## Notes

- Paths are hardcoded to this project's NERSC Perlmutter scratch layout. Reusing this script
  elsewhere requires updating the `SCRATCH`/`OUT` path constants.
