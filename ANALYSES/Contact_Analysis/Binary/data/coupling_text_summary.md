# Contact-coupling network similarity (MI-thresholded HMM-column graphs)

- MI threshold for edge inclusion: 0.05 bits

## PRK coupling graph similarity (Jaccard on edges)

| condition A | condition B | edges A | edges B | shared | Jaccard | Frobenius | Spectral |
|---|---|---|---|---|---|---|---|
| med4_ox | med4_ox | 83 | 83 | 83 | 1.000 | 0.000 | 0.000 |
| med4_ox | med4_red | 83 | 96 | 17 | 0.105 | 1.977 | 0.385 |
| med4_ox | hm2_ox | 83 | 153 | 25 | 0.118 | 3.648 | 2.407 |
| med4_ox | hm2_red | 83 | 167 | 24 | 0.106 | 2.739 | 1.487 |
| med4_red | med4_red | 96 | 96 | 96 | 1.000 | 0.000 | 0.000 |
| med4_red | hm2_ox | 96 | 153 | 30 | 0.137 | 3.732 | 2.212 |
| med4_red | hm2_red | 96 | 167 | 27 | 0.114 | 3.010 | 1.220 |
| hm2_ox | hm2_ox | 153 | 153 | 153 | 1.000 | 0.000 | 0.000 |
| hm2_ox | hm2_red | 153 | 167 | 51 | 0.190 | 3.975 | 1.175 |
| hm2_red | hm2_red | 167 | 167 | 167 | 1.000 | 0.000 | 0.000 |

## GAP2 coupling graph similarity (Jaccard on edges)

| condition A | condition B | edges A | edges B | shared | Jaccard | Frobenius | Spectral |
|---|---|---|---|---|---|---|---|
| med4_ox | med4_ox | 33 | 33 | 33 | 1.000 | 0.000 | 0.000 |
| med4_ox | med4_red | 33 | 14 | 5 | 0.119 | 1.577 | 0.946 |
| med4_ox | hm2_ox | 33 | 13 | 1 | 0.022 | 1.610 | 0.983 |
| med4_ox | hm2_red | 33 | 22 | 3 | 0.058 | 1.667 | 0.725 |
| med4_red | med4_red | 14 | 14 | 14 | 1.000 | 0.000 | 0.000 |
| med4_red | hm2_ox | 14 | 13 | 2 | 0.080 | 0.736 | 0.094 |
| med4_red | hm2_red | 14 | 22 | 5 | 0.161 | 0.888 | 0.272 |
| hm2_ox | hm2_ox | 13 | 13 | 13 | 1.000 | 0.000 | 0.000 |
| hm2_ox | hm2_red | 13 | 22 | 2 | 0.061 | 0.926 | 0.323 |
| hm2_red | hm2_red | 22 | 22 | 22 | 1.000 | 0.000 | 0.000 |

