# Host & redox per-HMM tests on log10(τ_on)

## PRK host effect (MED4 vs HM2)
- HMM columns tested per redox state: **103**
- significant after BH-FDR q<0.05 AND CI excludes 0: **4**
- top significant positions:

| HMM col | redox | n_med4 | n_hm2 | log-dwell HM2-MED4 | 95% CI | BH q |
|---|---|---|---|---|---|---|
| 45 | red | 100 | 82 | -0.43 | [-0.55, -0.16] | 6.11e-06 |
| 34 | red | 164 | 107 | +0.51 | [+0.30, +0.87] | 3.35e-05 |
| 23 | red | 69 | 64 | +0.29 | [+0.05, +0.55] | 0.00725 |
| 34 | ox | 130 | 116 | +0.36 | [+0.04, +0.65] | 0.0282 |

## PRK redox effect (ox vs red)
- HMM columns tested per host: **109**
- significant after BH-FDR q<0.05 AND CI excludes 0: **0**

## PRK convergence flags
- HMM columns evaluated: 63
- flagged 'undersampled' (within-CV > 1.5 and between-range < within-spread): **43**
- flagged columns: 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16, 17, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 54

## GAP2 host effect (MED4 vs HM2)
- HMM columns tested per redox state: **50**
- significant after BH-FDR q<0.05 AND CI excludes 0: **3**
- top significant positions:

| HMM col | redox | n_med4 | n_hm2 | log-dwell HM2-MED4 | 95% CI | BH q |
|---|---|---|---|---|---|---|
| 47 | ox | 50 | 15 | -0.51 | [-0.80, -0.15] | 0.00898 |
| 51 | red | 63 | 83 | +0.35 | [+0.09, +0.53] | 0.0359 |
| 66 | ox | 27 | 31 | -1.06 | [-1.71, -0.31] | 0.0359 |

## GAP2 redox effect (ox vs red)
- HMM columns tested per host: **51**
- significant after BH-FDR q<0.05 AND CI excludes 0: **1**
- top significant positions:

| HMM col | host | n_ox | n_red | log-dwell red-ox | 95% CI | BH q |
|---|---|---|---|---|---|---|
| 66 | med4 | 27 | 40 | -1.22 | [-1.84, -0.29] | 0.0299 |

## GAP2 convergence flags
- HMM columns evaluated: 36
- flagged 'undersampled' (within-CV > 1.5 and between-range < within-spread): **19**
- flagged columns: 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 67, 69, 70

