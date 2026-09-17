# Spatial-neighborhood partition of CP12 HMM-column MI coupling networks

CP12 HMM columns are partitioned into three structural/functional zones:

- **Block N (cols 1–22)** divergent N-terminus; substituted-position cluster
- **Block M (cols 23–50)** structured PRK-binding middle
- **Block C (cols 51–70)** C-terminal regulatory region (regulatory cys at 57, 66)

For each (condition × partner) MI matrix, the upper triangle is decomposed into three intra-block sums (within N, M, C) and three cross-block sums (N↔M, N↔C, M↔C). Edges with MI ≥ 0.05 bits are the same as those reported in the main coupling-graph analysis (Main Figure 5B).

## Total MI per condition × partner (sum of upper triangle, bits)

| partner | condition | intra N | intra M | intra C | cross N-M | cross N-C | cross M-C | total |
|---|---|---|---|---|---|---|---|---|
| prk | med4_ox | 1.70 | 0.77 | 2.13 | 2.30 | 4.37 | 2.35 | 13.61 |
| gap2 | med4_ox | 0.63 | 1.29 | 0.17 | 2.75 | 0.42 | 1.06 | 6.31 |
| prk | med4_red | 3.20 | 0.60 | 3.74 | 2.77 | 3.51 | 2.77 | 16.58 |
| gap2 | med4_red | 0.00 | 0.54 | 0.43 | 0.09 | 0.03 | 1.20 | 2.29 |
| prk | hm2_ox | 3.79 | 1.00 | 14.34 | 2.04 | 7.84 | 2.84 | 31.85 |
| gap2 | hm2_ox | 0.00 | 0.89 | 0.21 | 0.00 | -0.00 | 0.75 | 1.84 |
| prk | hm2_red | 2.72 | 1.52 | 7.03 | 4.84 | 7.38 | 4.92 | 28.41 |
| gap2 | hm2_red | 0.00 | 1.19 | 0.39 | 0.29 | 0.10 | 1.31 | 3.28 |

## Host difference per block (MED4 − HM2)

Positive value ⇒ MED4 has more coupling in that block; negative ⇒ HM2 has more.

| partner | redox | Δ intra N | Δ intra M | Δ intra C | Δ cross N-M | Δ cross N-C | Δ cross M-C | Δ total |
|---|---|---|---|---|---|---|---|---|
| prk | ox | -2.09 | -0.24 | -12.21 | +0.26 | -3.46 | -0.48 | -18.23 |
| prk | red | +0.48 | -0.92 | -3.28 | -2.07 | -3.87 | -2.15 | -11.82 |
| gap2 | ox | +0.63 | +0.41 | -0.05 | +2.75 | +0.42 | +0.31 | +4.47 |
| gap2 | red | -0.00 | -0.65 | +0.05 | -0.20 | -0.07 | -0.12 | -0.98 |

