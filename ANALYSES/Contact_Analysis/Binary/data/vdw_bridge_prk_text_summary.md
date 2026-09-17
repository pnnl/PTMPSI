# Δdwell vs Δ(per-residue VDW) — CP12 vs PRK, MED4 − HM2

VDW source: direct trajectory L-J computation (cutoff 10 Å, stride 4, frames 200-end), using amber99sb σ/ε from `mmgbsa_cp12_prk/COM.prmtop`. This is the equivalent of gmx_MMPBSA's idecomp=2 DELTAS VDW computed without re-running MM/PBSA.

**Validation:** the same direct method applied to CP12 ↔ GAP2 has Pearson r = +0.908 (n = 270 HMM-col × condition cells) against the existing gmx_MMPBSA idecomp=2 numbers from `mmgbsa_cp12_gap2_nad/`.

## OX
| Δ component | n_all | r_all | p_all | n_subs | r_subs | p_subs |
|---|---|---|---|---|---|---|
| vdw_mean | 57 | -0.294 | 0.0266 | 31 | -0.744 | 1.59e-06 |

## RED
| Δ component | n_all | r_all | p_all | n_subs | r_subs | p_subs |
|---|---|---|---|---|---|---|
| vdw_mean | 56 | -0.451 | 0.000487 | 30 | -0.511 | 0.00392 |

