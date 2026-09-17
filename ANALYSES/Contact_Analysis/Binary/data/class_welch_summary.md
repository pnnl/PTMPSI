# Class-stratified Welch test on log τ_on (CP12 ↔ partner, MED4 vs HM2)

For each (partner × redox × chemistry class), Welch t-test on the log10(mean_on_ps) distribution comparing MED4 vs HM2 (pooled across CP12 chains × replicates × pair instances). Bootstrap 95% CI on the median host difference is the conservative anchor of significance; BH-FDR corrects across the full 20-test family.

Interpretation guide:
- `log_diff_hm2_minus_med4 > 0` means HM2 has longer dwell than MED4 in that cell.
- `log_diff_ps_factor = 10^Δlog` gives the multiplicative dwell shift (e.g., 1.30 ⇒ HM2 contacts last 30% longer; 0.77 ⇒ MED4 contacts last 30% longer).
- A test is flagged `sig_q05 = True` iff BH-q < 0.05 AND the bootstrap CI on Δlog-dwell excludes zero.

## Significant cells (BH-FDR q < 0.05 AND CI excludes 0)

- none —
## All cells

| partner | redox | class | n_med4 | n_hm2 | median MED4 (ps) | median HM2 (ps) | Δlog (HM2−MED4) | 95% CI | dwell factor | d | Welch p | BH q | sig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| prk | ox | hydrophobic | 220 | 312 | 1480 | 1768 | +0.077 | [-0.057, +0.232] | 1.19 | +0.11 | 0.211 | 0.469 | n |
| prk | ox | polar_hbond | 1585 | 1657 | 1462 | 1500 | +0.011 | [-0.063, +0.071] | 1.03 | +0.00 | 0.945 | 0.974 | n |
| prk | ox | salt_bridge | 417 | 303 | 2200 | 2621 | +0.076 | [-0.090, +0.215] | 1.19 | +0.13 | 0.0849 | 0.243 | n |
| prk | ox | cation_pi | 83 | 54 | 1157 | 1545 | +0.126 | [-0.206, +0.351] | 1.34 | +0.14 | 0.438 | 0.875 | n |
| prk | ox | mixed | 1174 | 970 | 1000 | 1100 | +0.041 | [-0.013, +0.100] | 1.10 | +0.08 | 0.0527 | 0.181 | n |
| prk | red | hydrophobic | 205 | 328 | 1286 | 1798 | +0.146 | [-0.016, +0.266] | 1.40 | +0.23 | 0.00692 | 0.0529 | n |
| prk | red | polar_hbond | 1506 | 1701 | 1300 | 1500 | +0.062 | [+0.010, +0.115] | 1.15 | +0.09 | 0.00793 | 0.0529 | n |
| prk | red | salt_bridge | 426 | 304 | 2195 | 2347 | +0.029 | [-0.103, +0.183] | 1.07 | +0.01 | 0.913 | 0.974 | n |
| prk | red | cation_pi | 84 | 47 | 913 | 1580 | +0.238 | [-0.150, +0.501] | 1.73 | +0.03 | 0.845 | 0.974 | n |
| prk | red | mixed | 1188 | 1018 | 1100 | 1188 | +0.033 | [-0.039, +0.083] | 1.08 | +0.02 | 0.685 | 0.974 | n |
| gap2 | ox | hydrophobic | 123 | 109 | 2020 | 1513 | -0.125 | [-0.325, +0.104] | 0.75 | -0.09 | 0.5 | 0.908 | n |
| gap2 | ox | polar_hbond | 797 | 784 | 1900 | 1856 | -0.010 | [-0.111, +0.082] | 0.98 | -0.07 | 0.151 | 0.377 | n |
| gap2 | ox | salt_bridge | 242 | 201 | 2177 | 2264 | +0.017 | [-0.209, +0.273] | 1.04 | -0.02 | 0.871 | 0.974 | n |
| gap2 | ox | cation_pi | 76 | 53 | 3068 | 5286 | +0.236 | [-0.291, +0.678] | 1.72 | +0.07 | 0.7 | 0.974 | n |
| gap2 | ox | mixed | 640 | 584 | 1537 | 1475 | -0.018 | [-0.146, +0.102] | 0.96 | -0.11 | 0.0475 | 0.181 | n |
| gap2 | red | hydrophobic | 117 | 104 | 1160 | 1014 | -0.058 | [-0.238, +0.116] | 0.87 | -0.02 | 0.906 | 0.974 | n |
| gap2 | red | polar_hbond | 766 | 788 | 2062 | 1800 | -0.059 | [-0.190, +0.052] | 0.87 | -0.14 | 0.00611 | 0.0529 | n |
| gap2 | red | salt_bridge | 251 | 194 | 1822 | 2275 | +0.096 | [-0.100, +0.265] | 1.25 | -0.00 | 0.974 | 0.974 | n |
| gap2 | red | cation_pi | 67 | 48 | 3943 | 9079 | +0.362 | [-0.057, +1.073] | 2.30 | +0.37 | 0.0542 | 0.181 | n |
| gap2 | red | mixed | 605 | 597 | 1383 | 1420 | +0.011 | [-0.093, +0.136] | 1.03 | -0.02 | 0.667 | 0.974 | n |
