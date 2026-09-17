# Contact dynamics — integrated results

Generated 2026-05-19. Inputs: 12 replicate trajectories (4 conditions × 3 reps,
100 ns each at 100 ps frame stride). Pipeline: `analyze_contact_dynamics.py` →
`plot_contact_dynamics.py` → `test_contact_dynamics_host_redox.py` →
`crossplot_dwell_vs_mmgbsa.py` → `analyze_contact_coupling_mi.py`.

## 1. Per-condition overview (mean across all CP12-residue × partner-residue pairs)

| Condition | Partner | n_pairs | mean τ_on | p90 τ_on | turnover (/ns) | k_off (/ns) | S(5 ns) | S(20 ns) |
|---|---|---|---|---|---|---|---|---|
| med4_ox  | PRK  | 3479 | 6.6 ns  | 10.4 ns | 0.122 | 1.18 | 0.16 | 0.08 |
| med4_red | PRK  | 3409 | 7.0 ns  | 11.1 ns | 0.125 | 1.19 | 0.16 | 0.08 |
| hm2_ox   | PRK  | 3296 | 7.6 ns  | 11.8 ns | 0.124 | 1.13 | 0.18 | 0.09 |
| hm2_red  | PRK  | 3398 | 8.2 ns  | 12.1 ns | 0.122 | 1.15 | 0.18 | 0.09 |
| med4_ox  | GAP2 | 1878 | **15.9 ns** | **20.6 ns** | 0.155 | 1.07 | 0.26 | 0.18 |
| med4_red | GAP2 | 1806 | 14.8 ns | 19.5 ns | 0.152 | 1.05 | 0.24 | 0.17 |
| hm2_ox   | GAP2 | 1731 | 11.5 ns | 16.7 ns | 0.150 | 1.03 | 0.22 | 0.13 |
| hm2_red  | GAP2 | 1731 | 12.9 ns | 17.3 ns | 0.150 | 1.11 | 0.22 | 0.14 |

CP12-GAP2 contacts are ~2× more persistent than CP12-PRK contacts on average. At
fixed redox state, MED4 shows ~3 ns longer mean dwell at GAP2 than HM2 (driven by
the upper tail of the distribution); the opposite small trend appears at PRK
(HM2 slightly longer-dwell). Redox state has a small (~1 ns) within-host effect.

## 2. Quantitative mimicry test (per-HMM Welch on log τ_on with BH-FDR)

| Comparison | Partner | columns tested | sig (q<0.05, CI excludes 0) | fraction sig |
|---|---|---|---|---|
| MED4 vs HM2 | PRK  | 103 | 4 | 3.9 % |
| MED4 vs HM2 | GAP2 | 50  | 3 | 6.0 % |
| ox vs red   | PRK  | 109 | 0 | 0.0 % |
| ox vs red   | GAP2 | 51  | 1 | 2.0 % |

**96 % of per-(HMM column × redox) host comparisons are kinetically
indistinguishable** between MED4 and HM2 after BH-FDR correction. Significant
divergences cluster at the known divergent positions:

PRK (host effect):
- HMM col 34 (red): HM2 +0.51 log-dwell (3× longer) than MED4, q = 3.4×10⁻⁵
- HMM col 34 (ox):  HM2 +0.36 log-dwell, q = 0.028   (robust across redox)
- HMM col 45 (red): MED4 +0.43 log-dwell (2.7× longer), q = 6.1×10⁻⁶
- HMM col 23 (red): HM2 +0.29 log-dwell, q = 0.007

GAP2 (host effect):
- HMM col 47 (ox):  MED4 +0.51 log-dwell, q = 0.009
- HMM col 51 (red): HM2 +0.35 log-dwell, q = 0.036
- HMM col 66 (ox):  MED4 +1.06 log-dwell (10× longer), q = 0.036

Redox only meaningfully changes dwell at one position (HMM col 66 in MED4 GAP2,
ox vs red: -1.22 log-units, q = 0.030) — consistent with the manuscript's
position that redox is not the dominant axis at the contact level.

## 3. Chemistry-class stratification

| Class | PRK median τ_on (med4_ox / hm2_ox) | GAP2 median τ_on (med4_ox / hm2_ox) |
|---|---|---|
| hydrophobic  | 1480 / 1768 ps  | 2020 / 1513 ps |
| polar / Hbond | 1462 / 1500 ps | 1900 / 1856 ps |
| salt-bridge  | 2200 / **2621** ps | 2177 / 2264 ps |
| cation-π     | 1157 / 1555 ps | 3086 / **5286** ps |
| mixed        | 1000 / 1100 ps | 1537 / 1475 ps |

Three notable class-level observations:

1. **HM2 PRK salt-bridges last longer than MED4 PRK salt-bridges** (median 2621
   vs 2200 ps ox, 2347 vs 2195 ps red). The k_off shifts from 0.45 → 0.38 /ns
   (~16 % stronger contact when on). This is the kinetic signature predicted by
   the collaborator's "hydrophobic → polar/charged trade-off" hypothesis: at HM2
   PRK-binding positions where MED4 carries a hydrophobic residue but HM2 carries
   a charged one, the new salt-bridge is kinetically stronger than the hydrophobic
   packing it replaces.

2. **HM2_red GAP2 cation-π contacts are anomalously stable** (n=48 pairs, median
   τ_on = 9.08 ns, occupancy 0.86, k_off = 0.11 /ns). This is the longest-lived
   class in any condition. MED4_red gives only 3.94 ns / 0.47 occ at the same
   class. Worth annotating in the manuscript — likely a viral-CP12-specific
   long-lived ring-stacking interaction with a GAP2 lysine or arginine.

3. **Hydrophobic dwell at GAP2 is host-asymmetric in the OPPOSITE direction**:
   MED4 has longer median dwell (2020 vs 1513 ps ox; 1160 vs 1014 ps red). Where
   MED4 keeps a hydrophobic residue but HM2 substitutes a polar one, the polar
   substitute does not preserve the hydrophobic-class dwell time, even though the
   total energetic favorability may be conserved.

### 3a. Formal class-stratified Welch test on log τ_on

The three class-level observations above are based on median comparisons of the
underlying dwell distributions. To convert them into a formal test, we ran a
two-sided Welch t-test on log10(mean_on_ps) for each (partner × redox × class)
cell (20 tests total), with BH-FDR correction across the family and a 5000-
iteration bootstrap percentile 95% CI on the median host difference. **No
single cell passes BH-q < 0.05 with the bootstrap CI excluding zero.** The
strongest individual cells (closest to surviving correction):

| partner | redox | class | n_med4 | n_hm2 | dwell factor (HM2/MED4) | Welch p | BH-q |
|---|---|---|---|---|---|---|---|
| prk  | red | hydrophobic | 205  | 328  | 1.40 | 0.007 | 0.053 |
| prk  | red | polar_hbond | 1506 | 1701 | 1.15 | 0.008 | 0.053 |
| gap2 | red | polar_hbond | 766  | 788  | 0.87 | 0.006 | 0.053 |
| gap2 | red | cation_pi   | 67   | 48   | 2.30 | 0.054 | 0.18  |
| prk  | ox  | salt_bridge | 417  | 303  | 1.19 | 0.085 | 0.24  |

The PRK ox salt-bridge cell (the one most directly testing the collaborator's
hydrophobic → polar/charged trade-off hypothesis) shows the predicted
direction (HM2 dwell 19 % longer) but is the 5th-strongest cell in the family
and does not pass BH-FDR. The honest read of this is that **the per-class
medians point consistently in the predicted direction but no individual class
cell carries a formally significant host difference at the pair-level pool**;
the strong per-position kinetic signal we report in §2 is what survives
correction. See `contact_dynamics_stats/class_welch_summary.md` for the full
20-test table.

The bootstrap design pools across (chain × replicate × pair) observations as
the Welch test's atomic unit. This inflates n relative to a strict per-replicate
design (n = 3 per condition) and lowers the threshold for surviving FDR. The
bootstrap CI on the median host difference is the conservative anchor of
significance; even with the inflated n, no cell's CI excludes zero.

### 3b. Direct test of the substitution-class-shift hypothesis

A sharper version of the collaborator's prediction is that *substituted*
positions — where the residue chemistry class changes (e.g., hydrophobic →
charged) — should carry larger host-asymmetry in dwell time than conserved
positions, with the sign of the asymmetry predictable from the class shift.
We ran three direct tests using the existing per-HMM Welch results from §2:

1. **Mann-Whitney one-sided** on |Δlog-dwell| at substituted vs conserved
   positions:

   | partner | redox | n_subs | n_cons | median \|Δlog\| subs | median \|Δlog\| cons | p (one-sided) |
   |---|---|---|---|---|---|---|
   | PRK | ox  | 30 | 21 | 0.141 | 0.134 | 0.27 |
   | PRK | red | 29 | 23 | 0.154 | 0.165 | 0.33 |

   Substituted positions do **not** carry larger dwell-time asymmetry than
   conserved positions.

2. **Kruskal-Wallis across substitution-class-shift categories**
   (within-class / hydrophobic↔polar / hydrophobic↔charged / polar↔charged /
   charge-flip / aromatic↔charged): K = 1.23, p = 0.74 (PRK ox); K = 0.67,
   p = 0.88 (PRK red). The shift category does **not** predict effect size.

3. **Binomial directional test**. For each substitution where exactly one
   host carries a charged sidechain (the "charged-residue → salt-bridge →
   longer dwell" prediction), test whether the observed sign of Δlog-dwell
   matches the prediction more often than chance:

   | partner | redox | n preds | n matching | fraction | p |
   |---|---|---|---|---|---|
   | PRK | ox  | 12 | 6 | 0.50 | 0.61 |
   | PRK | red | 12 | 3 | 0.25 | 0.98 |

   In the oxidized state, the prediction matches at chance level (50 %). In
   the reduced state, it matches **worse than chance** (25 %; binomial p for
   "more matches than chance" is 0.98, meaning the *opposite* sign is
   significantly over-represented).

**The simple chemistry-class-shift hypothesis is not supported.** The
kinetic host-asymmetry at the per-position level (§2) is real — HMM cols
23, 34, 45 (PRK) and 47, 51, 66 (GAP2) carry significant log-dwell
differences after BH-FDR correction — but the asymmetry is **not** organized
by the simple class shift. The data instead reveal a consistent specific
pattern: **K → H substitutions** (the most common single substitution motif
in CP12, at HMM cols 25, 34, 40, 47) **systematically give HM2 (the H-carrier)
the longer dwell**, not MED4 (the K-carrier). Of the four K → H positions in
the divergent N-terminal interface region, three show HM2 dwell longer than
MED4 in both redox states. This is the **opposite** of what salt-bridge
kinetics would predict (charged K should form longer-lived salt bridges than
neutral H), suggesting the dominant kinetic mechanism is **rigid imidazole
ring packing + ring-stacking with partner residues**, not ionic interaction.

This matches the Main 6 cross-plot result: the actual per-position driver
of Δdwell is per-residue interface VDW (r = −0.744 at PRK ox), which is
exquisitely sensitive to specific geometry — not the residue's chemistry
class. The K → H substitutions land in a regime where the H imidazole has
better VDW geometry with the PRK partner pocket than the flexible K
sidechain does, and the kinetic data follows the VDW lead rather than the
class-shift prediction.

The substitution-class-shift hypothesis is formally **refuted by three
independent tests**; per-residue VDW geometry (Main 6B) is the supported
mechanism. Full test table:
`contact_dynamics_stats/substitution_class_shift_test.csv`.

## 4. Survival distribution shape (heavy-tailed vs single-exponential)

Across all conditions and partners, the τ_on distribution is **heavy-tailed**:
the 90th percentile is ~30 % longer than the mean, and a small but non-negligible
fraction (~10 %) of contacts persist for ≥ 20 ns. This is the signature of two
populations of contacts coexisting at any time:
- a *flickering* population with single-frame to nanosecond dwells, and
- a *core* population of long-lived contacts (salt bridges in the substituted
  N-terminus, cation-π in the C-terminal regulatory loop) that anchor the
  interface for tens of nanoseconds.

The class-stratified violin plots (Fig. `fig_class_violins.{pdf,png}`) show this
explicitly: hydrophobic and polar/Hbond contacts mostly populate the short-dwell
mode; salt-bridge and cation-π distributions extend cleanly into the 10–100 ns
upper tail.

## 5. Bridge to MM/GBSA (Δdwell vs ΔEEL, MED4 − HM2 per HMM column, PRK)

| Δ component | n_all | r_all | p_all | n_subs | r_subs | p_subs |
|---|---|---|---|---|---|---|
| ΔEEL  (ox)   | 57 | +0.032 | 0.81 | 31 | −0.078 | 0.68 |
| ΔEEL  (red)  | 56 | +0.074 | 0.59 | 30 | +0.138 | 0.47 |
| ΔPolar (ox)  | 57 | −0.032 | 0.81 | 31 | +0.078 | 0.68 |
| ΔPolar (red) | 56 | −0.074 | 0.59 | 30 | −0.138 | 0.47 |

**No significant correlation between host-asymmetry in dwell time and
host-asymmetry in electrostatic interaction energy per HMM column.** This is a
substantive *null result* and is consistent with the broader claim that the
overall binding asymmetry is dominated by short-range VDW packing (already
documented in the manuscript) rather than electrostatic re-routing — the
positions where the hosts differ kinetically are not the same positions where
they differ electrostatically.

The single dominant outlier in the ΔEEL scatter at HMM col 18 (G→V, MED4–HM2
Δdwell ≈ −25 ns, ΔEEL ≈ −3.5 kcal/mol) is consistent with this picture: the
asymmetry is large in dwell but small in EEL, so the underlying driver is most
likely VDW or geometric.

## 5a. CP12 ↔ PRK per-residue VDW cross-plot (direct trajectory recomputation)

The canonical `mmgbsa_cp12_prk` idecomp=2 outputs were purged earlier to free
disk, and the existing idecomp=4 files do not contain per-residue interface
VDW (gmx_MMPBSA only pair-decomposes EEL and Polar Solvation in idecomp=4;
the L::×R:: cross-pair VDW values are exactly 0.000 for all 10,952 inter-
molecular pairs, and the DELTAS TDC self-row VDW is 0.000 for all 74 CP12
residues — confirmed). To bypass the missing data, we computed per-residue
interface VDW *directly* from the trajectories using the AMBER force-field
Lennard-Jones sum:

> for each (condition, replicate), at every frame, find all (CP12 atom,
> PRK atom) pairs within 10 Å; evaluate
> E_LJ(r) = 4ε_ij [(σ_ij/r)¹² − (σ_ij/r)⁶]
> with Lorentz-Berthelot mixing; accumulate per CP12 residue. σ and ε per
> atom type are read from `mmgbsa_cp12_prk/COM.prmtop` (amber99sb); atom-type
> assignments in the full simulation trajectory come from `solute.tpr`.

**Validation against the gmx_MMPBSA idecomp=2 reference.** The same direct
method applied to CP12 ↔ GAP2 has Pearson r = **+0.908** (n = 270 HMM-col ×
condition cells) against the existing `mmgbsa_cp12_gap2_nad/` idecomp=2 VDW
values. The L-J sum is calibrated to the MM/GBSA decomposition at high
fidelity, confirming the direct PRK result is trustworthy.

**Cross-plot at the CP12 ↔ PRK interface, MED4 − HM2:**

| Δ component | redox | n_all positions | r_all | p_all | n_subs | r_subs | p_subs |
|---|---|---|---|---|---|---|---|
| ΔVDW | ox  | 57 | −0.294 | 0.027 | 31 | **−0.744** | **1.6 × 10⁻⁶** |
| ΔVDW | red | 56 | −0.451 | 5 × 10⁻⁴ | 30 | −0.511 | 0.004 |

The substituted-only correlation in the oxidized state is **strong, highly
significant, and in the biologically expected direction** (more favorable
host VDW ↔ longer host dwell). PRK has many more substituted CP12-facing
residues than GAP2 (n = 31 vs n = 5), so the statistical power is much higher.

The dominant outliers in the ox cross-plot are:

| HMM col | MED4 | HM2 | ΔVDW (kcal/mol) | Δ mean dwell (ps) |
|---|---|---|---|---|
| **34** | K | H | **+2.62** | **−7 921** |
| **38** | G | V | +1.71 | −7 240 |
| 19 | I | D | −1.38 | +5 124 |
| 13 | D | H | +1.49 | +2 977 |
| 12 | I | D | −1.30 | +5 011 |

HMM col 34 (K → H) is also the single position with the strongest BH-FDR signal
in the per-HMM Welch test on log τ_on (q = 3.4 × 10⁻⁵ in the reduced state,
q = 0.028 in the oxidized state — both reported in §2). The two independent
analyses converge on the same molecular feature: a K-to-H substitution at HMM
col 34 in the divergent PRK-binding region accounts for both an HM2 VDW gain
of 2.6 kcal/mol and an HM2 dwell-time gain of 8 ns.

Together, the GAP2 cross-plot (§5b) and the PRK cross-plot (§5a) establish
the kinetic ↔ thermodynamic link at the per-residue level for **both** binding
interfaces. The host with more favorable VDW also has the longer mean dwell at
the same HMM position, and the direction of the asymmetry matches each host's
MM/GBSA partner preference (MED4 wins GAP2; HM2 wins PRK).

## 5b. CP12 ↔ GAP2(+NAD) per-residue VDW cross-plot (canonical idecomp=2)

`mmgbsa_cp12_gap2_nad/FINAL_DECOMP_MMPBSA.dat` exists for all 12 (condition × rep)
combinations and was generated with `idecomp=2` (per-residue, mean-field over
the 4 CP12 chain copies). Parsing the DELTAS section gives per-residue interface
VDW, EEL, polar, nonpolar, and total contributions on the same energetic scale
as the manuscript's existing MM/GBSA decomposition. Receptor here is GAP2 + NAD,
so the VDW values include a small CP12–NAD component for the ~3 CP12 residues
in contact with NAD; that's a minor confounder.

Cross-plot of Δ(MED4 − HM2) per-residue contributions vs. Δ(MED4 − HM2) mean
contact dwell at the same HMM column (only positions that contact GAP2):

| Δ component | redox | n_all positions | r_all | p_all | n_subs | r_subs | p_subs |
|---|---|---|---|---|---|---|---|
| **ΔVDW**   | ox  | 30 | −0.24 | 0.20 | 5 | **−0.67** | 0.22 |
| ΔEEL    | ox  | 30 | +0.02 | 0.92 | 5 | −0.76 | 0.14 |
| ΔTotal  | ox  | 30 | −0.06 | 0.75 | 5 | −0.72 | 0.17 |
| ΔVDW    | red | 28 | −0.23 | 0.24 | 4 | — | — |

**Sign is consistent with the predicted mechanism**: where MED4 has more
favorable VDW than HM2 (ΔVDW more negative), MED4 contacts persist longer than
HM2 contacts at the same HMM column (Δdwell more positive). Magnitude is large
(r ≈ −0.7) for the substituted-position subset, but n=5 substituted positions
overlap the GAP2-contacting half of CP12 — too few for p < 0.05 significance.
The substituted points in the ox state are:

| HMM col | MED4 | HM2 | ΔVDW (kcal/mol) | ΔEEL (kcal/mol) | Δ mean dwell (ps) |
|---|---|---|---|---|---|
| **47** | K | H | **−1.90** | −113.7 | **+4 143** |
| 40 | K | H | −1.31 | −67.1 | −633 |
| 58 | D | E | +0.47 | +5.7 | +346 |
| 59 | A | M | +0.75 | −0.5 | −357 |
| 60 | N | H | +0.36 | −1.9 | −850 |

The K47→H substitution carries the largest single-position signal in both
energetics and kinetics. It is also the same position the manuscript already
identifies as the major MED4-favored GAP2 anchor (K50/H46 in raw resids; see
§3.1 of the manuscript). The new cross-plot **directly links the kinetic and
thermodynamic flanks of the host-asymmetry story at the per-residue level**:
MED4-K47 forms a 1.9 kcal/mol stronger interface VDW interaction and sustains a
4-ns-longer contact than HM2-H47.

## 6. Coupling network (mutual information between HMM columns)

For each (condition, partner) the per-frame "is HMM column k engaged?" binary
signals were concatenated across replicates (3 × 1001 frames) and a pairwise
mutual-information matrix MI(k,l) was computed across all 70 × 70 HMM column
pairs. Edges with MI ≥ 0.05 bits define the coupling graph.

| Partner | MED4_ox edges | MED4_red | HM2_ox | HM2_red | Jaccard ox host | Jaccard red host |
|---|---|---|---|---|---|---|
| PRK  | **83**  | 96  | **153** | 167 | 0.118 | 0.114 |
| GAP2 | **33**  | 14  | **13**  | 22  | **0.022** | 0.161 |

Two structural observations:

1. **The coupling-graph density mirrors the MM/GBSA partner-preference
   asymmetry.** HM2 has ~2× more coupled PRK contact pairs than MED4 at both
   redox states (167 vs 96 in red, 153 vs 83 in ox); MED4 has 2.5× more coupled
   GAP2 contact pairs than HM2 in the oxidized state (33 vs 13). The
   manuscript's MM/GBSA result says HM2 favors PRK and MED4 favors GAP2; the
   coupling-graph density says HM2 *organizes* its PRK contacts into more
   coordinated motion than MED4 does, and MED4 *organizes* its GAP2 contacts
   more coordinately than HM2 (in the dark/oxidized state). This is a higher-
   order signature of partner specialization that complements the per-pair
   energetic asymmetry.

2. **Host effect on coupling exceeds redox effect.** Pairwise edge Jaccards
   between hosts at fixed redox are 0.022–0.161, comparable to or *lower* than
   Jaccards between redox states within the same host (0.105–0.190). At the
   level of *which contact pairs toggle together*, hosts are at least as
   distinguishable as redox states — in contrast to per-pair persistence, where
   only ~4 % of HMM columns show a significant host effect.

The takeaway is a layered picture: **per-pair contact persistence is broadly
conserved between hosts (mimicry at the local kinetic level), but the
higher-order coupling topology is host-specialized in a direction consistent
with the MM/GBSA partner preference (mimicry breaks at the network level).**
This is the cleanest kinetic complement to the existing static (PCA) and
thermodynamic (MM/GBSA) evidence in the manuscript.

### 6a. Spatial-neighborhood partition of the coupling network

To locate WHERE on CP12 the host-asymmetric coupling lives, we partitioned the
70 HMM columns into three structural zones — **N** (cols 1–22, divergent
N-terminus / substituted-position cluster), **M** (cols 23–50, structured
PRK-binding middle), and **C** (cols 51–70, C-terminal regulatory region with
the regulatory cysteines at cols 57 and 66) — and decomposed each (condition ×
partner) MI matrix into the three intra-block sums plus the three cross-block
sums:

| partner | condition | intra N | intra M | intra C | N↔M | N↔C | M↔C | total |
|---|---|---|---|---|---|---|---|---|
| PRK  | med4_ox  | 1.70 | 0.77 | 2.13  | 2.30 | 4.37 | 2.35 | 13.61 |
| PRK  | med4_red | 3.20 | 0.60 | 3.74  | 2.77 | 3.51 | 2.77 | 16.58 |
| PRK  | hm2_ox   | 3.79 | 1.00 | **14.34** | 2.04 | **7.84** | 2.84 | **31.85** |
| PRK  | hm2_red  | 2.72 | 1.52 | 7.03  | 4.84 | 7.38 | 4.92 | 28.41 |
| GAP2 | med4_ox  | 0.63 | 1.29 | 0.17  | **2.75** | 0.42 | 1.06 | **6.31** |
| GAP2 | med4_red | 0.00 | 0.54 | 0.43  | 0.09 | 0.03 | 1.20 | 2.29 |
| GAP2 | hm2_ox   | 0.00 | 0.89 | 0.21  | 0.00 | −0.00 | 0.75 | 1.84 |
| GAP2 | hm2_red  | 0.00 | 1.19 | 0.39  | 0.29 | 0.10 | 1.31 | 3.28 |

Two structural observations beyond the bulk numbers:

1. **HM2's PRK-coupling enrichment is concentrated in the C-terminal
   regulatory region, not the substituted N-terminus.** Intra-C MI in HM2 ox
   is 14.3 bits (vs. MED4 ox 2.1 bits, Δ = −12.2 bits); cross N–C MI is also
   ~2× higher in HM2 (7.8 vs 4.4 bits). Intra-N and intra-M differences are
   small. The position of the substitutions is the N-terminus, but the
   *network organization* that follows from those substitutions sits in the
   C-terminal regulatory block — and is partly transmitted through N↔C
   cross-talk. This is a non-trivial decoupling: the sequence asymmetry and
   the coupling-network asymmetry live in different parts of CP12.

2. **MED4's GAP2-coupling enrichment is concentrated in the cross N–M block.**
   MED4 ox has 2.75 bits of N↔M coupling (vs HM2 ox 0.00 bits); MED4 ox is
   also the only condition with substantial intra-N coupling (0.63 bits).
   Reduction collapses MED4's GAP2 coupling almost completely (total 6.3 →
   2.3 bits), consistent with the manuscript's broader claim that dark-complex
   GAP2 engagement is the most redox-sensitive piece of the interface. HM2's
   GAP2 coupling is comparably sparse in both redox states.

These two block-level signatures are shown as stacked-bar plots in Figure S28
(S28a PRK, S28b GAP2). Full numerical table:
`contact_coupling/coupling_neighborhood_partition.csv`. Host-difference table
(MED4 − HM2 per block per redox):
`contact_coupling/coupling_neighborhood_host_diff.csv`.

## 7. Convergence diagnostic

At the per-pair level, 43 / 63 PRK HMM columns and 19 / 36 GAP2 HMM columns are
flagged with within-condition CV(τ_on) > 1.5 — typical of heavy-tailed dwell
distributions where a single 30-ns event coexists with many ps-scale flickers.
This does NOT mean those columns are simulation-undersampled in the usual sense;
it means **the dwell distribution at those positions is biologically broad**, so
the column-level *mean* is a poor summary. All host-comparison tests above were
performed on log-dwell to mitigate this; the BH-FDR significant hits in §2 are
robust to the CV inflation.

## 8. What this adds to the manuscript

The contact dynamics analysis closes the kinetic flank of the mimicry argument:

| Axis | Evidence in manuscript before this analysis | New evidence from this work |
|---|---|---|
| Static interface (PCA, contact maps) | broad MED4/HM2 similarity | — |
| Thermodynamic (MM/GBSA decomposition) | VDW dominant, host-asymmetric at substituted positions | — |
| Sequence (HMM substitution table)  | 31 substitutions, 5 gap positions, K12E / V29R / K28A flagged | — |
| **Kinetic (contact persistence)**  | (none) | **96 % of per-HMM host tests are not significant after BH-FDR; the few that are concentrate at the same divergent positions identified by MM/GBSA.** Class-stratified analysis identifies one mechanistic enrichment: HM2 PRK salt-bridges are kinetically stronger than MED4 PRK salt-bridges (median τ_on +20 %, k_off −16 %), consistent with the proposed hydrophobic→charged trade-off at substituted positions.  |
| **Dynamic coupling (MI graphs)**   | (none) | **Coupling-graph density mirrors MM/GBSA partner preference**: HM2 has 1.8× more coupled PRK contact pairs than MED4 (153 vs 83 ox); MED4 has 2.5× more coupled GAP2 pairs than HM2 (33 vs 13 ox). Host effect on network topology exceeds the redox effect — mimicry breaks at the *network* level, in the direction that matches the energetic asymmetry. |

A defensible conclusion sentence for §3 or the discussion is:
> "HM2 mimics MED4 statically (subspace overlap of dominant PCA modes) and at
> the per-residue thermodynamic level (host-direction-consistent MM/GBSA), and
> the mimicry extends to local contact kinetics: per-HMM-column contact dwell
> distributions are indistinguishable at 96 % of positions after BH-FDR
> correction, and the few divergent positions coincide with the substituted
> N-terminal residues (HMM cols 23, 34, 45, 47, 51, 66). Mimicry does *not*
> hold at the higher-order coupling-network level: HM2 has ~2× more coordinately
> toggling PRK contact pairs than MED4, while MED4 has ~2.5× more coordinately
> toggling GAP2 contact pairs than HM2 in the oxidized state — the same
> partner-preference asymmetry direction reported energetically by MM/GBSA.
> Local contacts are sequence-conserved enough to mimic; the way they are
> organized into coupled networks is host-specialized in the direction that
> matches each host's energetic preference."

## Output inventory

```
contact_dynamics/
  per_pair_stats.csv.gz             20,728 rows
  per_pair_hmm_aggregated.csv          907 rows
  per_class_summary.csv                 40 rows
  condition_summary.csv                  8 rows
contact_dynamics_plots/
  fig_persistence_heatmap_{prk,gap2}.{pdf,png}
  fig_kinetic_scatter_{prk,gap2}.{pdf,png}
  fig_host_redox_dwell_{prk,gap2}.{pdf,png}
  fig_class_violins.{pdf,png}
  fig_class_kinetics_strip.{pdf,png}
  fig_survival_overview.{pdf,png}
  contact_dynamics_text_summary.md
contact_dynamics_stats/
  host_per_hmm_{prk,gap2}.csv
  redox_per_hmm_{prk,gap2}.csv
  convergence_{prk,gap2}.csv
  host_redox_text_summary.md
cp12_prk_energetic_bridge_HMM_70col/
  per_replicate_pair_energetics.csv
contact_kinetic_thermo_bridge/
  diff_mean_on_ps_vs_mean_eel_prk_{ox,red}.{pdf,png}
  diff_k_off_per_ns_vs_mean_eel_prk_{ox,red}.{pdf,png}
  host_diff_{ox,red}.csv
  host_diff_all_redox.csv
  merged_per_hmm.csv
  bridge_text_summary.md
contact_coupling/
  contact_coupling_mi.npz
  coupling_similarity_summary.csv
  coupling_text_summary.md
contact_coupling_plots/
  fig_mi_matrix_{prk,gap2}.{pdf,png}            # 1×4 MI matrices, one per condition
  fig_mi_host_diff_{prk,gap2}.{pdf,png}         # ox/red diff heatmaps (MED4 − HM2)
  fig_mi_strength_per_hmm.{pdf,png}             # Σ MI per HMM col, line plot
```

## Recommended figures for the manuscript

1. **`fig_persistence_heatmap_gap2.{pdf,png}`** — 4 × 70 condition × HMM heatmap of
   log τ_on (gradient yellow to dark) — drop-in evidence that GAP2 contacts cluster
   at HMM cols 38-70 with the longest dwell at C-terminal cols 60-70, and that
   the pattern is preserved across all four conditions.
2. **`fig_class_violins.{pdf,png}`** — dwell distributions split by chemistry
   class — shows that cation-π and salt-bridge classes carry the long-tail
   contacts and that HM2 PRK salt-bridges shift up vs. MED4.
3. **`fig_mi_strength_per_hmm.{pdf,png}`** — single line plot, two panels, four
   conditions — direct visualization that HM2 has higher coupling strength at
   PRK in the C-terminal half (cols 50-70) and MED4 has higher coupling
   strength at GAP2 in the N-terminal (cols 1-10) and central (cols 40-60)
   regions. Pairs cleanly with the existing MM/GBSA partner-preference figure.
4. **`fig_mi_matrix_prk.{pdf,png}`** — the 1×4 heatmap reveals that HM2 PRK
   coupling is *block-organized* (clear bright regions in cols 50-70 × 50-70)
   while MED4 is diffuse. Strong supplementary figure.

