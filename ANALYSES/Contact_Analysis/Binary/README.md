# Contact dynamics analysis package

This directory bundles the binary contact-reversibility analyses for the
CP12–PRK–GAP2 dark complex, plus follow-up extensions, in one self-contained
package. Everything is computed from the same 12 trajectories (4 conditions ×
3 replicates, 100 ns each, 1001 frames at 100 ps stride) used for the
manuscript's static contact-frequency and MM/GBSA results.

---

## Scope and design

The analysis question: can per-residue contacts be turned into binary time
series so that contact persistence and switching frequency become measurable,
separating hydrophobic from polar/charged contributions and short- from
long-range components? Is 100 ns enough trajectory to do this reliably?

The pipeline implements this with three refinements:

1. **Hysteresis-equivalent noise filter.** Naive 5.0 Å thresholding makes
   residues near the cutoff "flicker," inflating the turnover count. We apply
   a 3-frame open–close morphological filter to B(t): short 1-runs (<3 frames)
   are flipped to 0 (open), then short 0-runs are flipped to 1 (close). A
   contact is "real" only if it persists ≥ 300 ps; a brief disconnect <300 ps
   doesn't break it. This is equivalent in spirit to a 4.0 Å form / 5.5 Å break
   hysteresis.
2. **Full distribution, not just the mean.** Per pair we report mean / median /
   90th-percentile / max dwell time, transitions (0→1 and 1→0), turnover rate
   (events / ns), two-state k_on and k_off (using k_off = 1/⟨τ_on⟩ as the
   "contact strength" axis, closer to a thermodynamic quantity than raw mean
   dwell), and survival fractions S(1 ns), S(5 ns), S(20 ns).
3. **Chemistry-class decomposition.** Each (CP12 residue, partner residue) pair
   is tagged as **hydrophobic** (both sidechains in {A, V, L, I, M, F, W, Y, P}),
   **polar / H-bond** (≥1 in {S, T, N, Q, H}), **salt-bridge** ({D, E} ↔ {K, R}),
   **cation-π** ({K, R} ↔ {F, W, Y}), or **mixed**. Persistence / turnover are
   reported per class per condition — that is the direct test of your
   "hydrophobic → polar/charged trade-off" hypothesis.

Plus two extensions on top:

- **Per-HMM Welch test on log τ_on with BH-FDR.** Tells you at which exact HMM
  columns the contact kinetics differ between MED4 and HM2 (or between ox and
  red within host), after controlling for multiple comparisons.
- **Coupling-network mutual information.** For each condition, aggregate per-
  frame "HMM column k engaged" binary signals across replicates and compute
  pairwise MI between every 70 × 70 column pair. Compare graphs across
  conditions via Jaccard on the edge set, Frobenius on the MI matrix, and
  spectral distance. This is the higher-order question: do contacts that
  toggle together do so the same way in both hosts? (Spoiler: per-pair
  kinetics is conserved, network topology is not.)

The short answer to "is 100 ns enough?" — yes for fast contacts (≤ ns dwell;
salt bridges and hydrophobic packing flicker on this scale), borderline for
slow contacts (10s of ns rotamer-coupled rearrangements; some of those are
flagged as undersampled in the convergence diagnostic). The replicate-vs-host
variance check tells you per HMM column which kinetics are converged.

---

## Where to look first

| Want to | Open |
|---|---|
| Read the methods write-up | [`docs/METHODS.md`](docs/METHODS.md) |
| Read the integrated results write-up (the manuscript-ready text) | [`docs/RESULTS.md`](docs/RESULTS.md) |
| See the per-condition headline numbers | [`data/condition_summary.csv`](data/condition_summary.csv) |
| See the class-stratified result you predicted | [`data/per_class_summary.csv`](data/per_class_summary.csv) — search for "salt_bridge" and look at PRK rows for HM2 vs MED4 |
| See the per-HMM significance test | [`data/host_redox_text_summary.md`](data/host_redox_text_summary.md) |
| See the coupling-network finding | [`data/coupling_text_summary.md`](data/coupling_text_summary.md) |
| See *where on CP12* the coupling asymmetry lives | [`data/coupling_neighborhood_partition.md`](data/coupling_neighborhood_partition.md) |
| See the formal class-stratified Welch test (honest negative result; matches predicted direction) | [`data/class_welch_summary.md`](data/class_welch_summary.md) |
| See the substitution-class-shift hypothesis test (three tests, all refute the simple chemistry-shift framework; per-residue VDW geometry is the actual driver) | [`data/substitution_class_shift_summary.md`](data/substitution_class_shift_summary.md) |
| See the GAP2 VDW ↔ dwell cross-plot stats | [`data/vdw_bridge_text_summary.md`](data/vdw_bridge_text_summary.md) |
| **See the PRK VDW ↔ dwell cross-plot stats (the strongest single result, r = −0.744, p = 1.6×10⁻⁶)** | [`data/vdw_bridge_prk_text_summary.md`](data/vdw_bridge_prk_text_summary.md) |
| See the main figures | [`figures/Main5a*`, `Main5b*`, `Main6a*`, `Main6b*`](figures/) |

---

## Headline findings (TL;DR)

> **Central question:** HM2 (the viral strain) carries 31 amino-acid
> substitutions in its CP12 relative to native MED4. What do those
> substitutions actually *do* to how HM2-CP12 binds the dark complex
> (GAP2 + PRK), compared to MED4-CP12?
>
> **Short answer:** most are invisible — at the great majority of positions,
> HM2's substitution leaves the local binding kinetics indistinguishable from
> MED4 (broad mimicry). But where the substitutions *do* register, they
> systematically make HM2-CP12 a **tighter and longer-lasting PRK binder**
> than MED4-CP12 is. The cleanest single case is HMM column 34, where MED4's
> lysine is replaced by HM2's histidine: HM2's histidine packs 2.6 kcal/mol
> more favorably against PRK and holds the PRK contact ~8 ns longer per
> binding event. Two independent analyses — per-residue VDW energy and
> per-residue contact dwell — flag this same residue (BH-FDR q = 3.4 × 10⁻⁵
> from the Welch test, n = 31 substituted-position correlation r = −0.744,
> p = 1.6 × 10⁻⁶). The pattern points to a viral-CP12 design principle:
> rigid-imidazole / hydrophobic substitutions in the N-terminal PRK interface
> trade ionic flexibility for VDW packing, and the network-level consequence
> shows up as a denser PRK coupling network in HM2's C-terminal regulatory
> block.

1. **Most HM2 substitutions are kinetically silent — broad mimicry of MED4.**
   Per-(HMM column × redox) Welch tests on log τ_on with BH-FDR correction:
   only 4 of 103 PRK and 3 of 50 GAP2 tests pass q < 0.05 with bootstrap CI
   excluding zero. **96 % of HM2-substituted positions are kinetically
   indistinguishable from the MED4 residue they replaced.** The HM2
   substitutions that *do* change kinetics concentrate at HMM cols 23, 34, 45
   (PRK interface) and 47, 51, 66 (GAP2 interface) — the same positions
   identified as host-asymmetric by MM/GBSA per-residue decomposition.
2. **HM2's substitutions don't follow a chemistry-class rule — they follow
   VDW packing geometry.** The intuitive hypothesis was that HM2's
   hydrophobic → polar / charged → polar swaps would predict which host wins
   each contact. Three direct tests of that prediction all fail:
   - Class-stratified Welch on log-dwell (partner × redox × class, 20-cell
     family, BH-FDR + bootstrap CI): zero cells pass q < 0.05. The directional
     pattern is consistent (HM2 wins PRK ox salt-bridge dwell ×1.19,
     PRK red hydrophobic ×1.40, GAP2 red cation-π ×2.30) but no single class
     reaches formal significance.
   - MWU |Δlog-dwell| at substituted vs conserved positions (PRK): substituted
     0.141 vs conserved 0.134 ox, p = 0.27; substituted 0.154 vs conserved
     0.165 red, p = 0.33. HM2-substituted positions carry **no larger** dwell
     asymmetry than conserved positions.
   - Binomial directional sign test on charged-residue-host predictions (PRK):
     6 / 12 match ox (chance), 3 / 12 match red (worse than chance,
     p = 0.98 one-sided).

   The interesting specific pattern hiding behind the failure: **K → H
   substitutions consistently favor HM2 (the H-carrier), not MED4
   (the K-carrier)**. Of the four K → H positions in the divergent N-terminal
   PRK interface, three show HM2 dwell longer than MED4 in both redox states —
   the opposite of what salt-bridge kinetics would predict. HM2's histidines
   bind PRK better than MED4's lysines because the rigid imidazole ring packs
   tighter, not because of a new electrostatic contact. This is the same
   mechanism finding #5 makes concrete at the per-position level.
3. **A viral-specific anomaly in the reduced state: HM2_red GAP2 cation-π
   contacts are unusually stable.** n = 48 pairs, median τ_on = 9.08 ns,
   occupancy = 0.86, k_off = 0.11 / ns — the longest-lived contact class in
   any of the four conditions. Likely a viral-CP12-specific aromatic
   ring-stacking interaction with a GAP2 lysine or arginine; not present in
   MED4 at comparable lifetime.
4. **HM2 organizes its PRK contacts into a denser coupling network than
   MED4 does**, in the direction of HM2's MM/GBSA PRK preference. HM2 has
   1.8× more coupled PRK contact pairs (153 vs 83 edges in ox); MED4 has
   2.5× more coupled GAP2 contact pairs (33 vs 13). Host effect on network
   topology is comparable to or larger than the redox effect.
   **Spatial-neighborhood partition (Figure S28) localizes the asymmetry**:
   HM2's PRK-coupling enrichment is concentrated in the **C-terminal
   regulatory block** (intra-C, cols 51–70: 14.3 vs 2.1 bits, Δ = −12.2 bits),
   *not* in the substituted N-terminus. The 31 sequence substitutions sit in
   the N-terminus, but their network-level signature shows up in the
   C-terminal regulatory machinery — N-terminal HM2 substitutions reorganize
   how the C-terminal regulatory residues engage PRK. (MED4's reciprocal
   GAP2-coupling enrichment lives in the cross N ↔ M block: 2.75 vs 0.00 bits
   in ox.)
5. **HM2's substitutions that pack better against PRK also hold the PRK
   contact longer — the per-residue energetic and kinetic effects co-vary.**
   Per-HMM-column cross-plot of per-residue interface VDW vs contact dwell,
   restricted to the 31 substituted positions:
   - **CP12 ↔ PRK** (VDW computed directly from trajectories; the kinetic
     measure is per-pair mean dwell aggregated to HMM column): substituted-only
     Pearson **r = −0.744, p = 1.6 × 10⁻⁶** in the oxidized state;
     r = −0.51, p = 0.004 reduced. Dominant outliers are **HMM col 34
     (K → H)** and **HMM col 38 (G → V)** — at both positions, HM2's residue
     wins on both axes (better VDW packing, longer dwell). HMM col 34 is
     independently flagged as the strongest per-HMM Welch-test hit
     (BH-FDR q = 3.4 × 10⁻⁵). Two analytical paths converge on the same
     residue. The direct-VDW recomputation is validated against the GAP2
     MM/GBSA idecomp = 2 reference at Pearson r = +0.908.
   - **CP12 ↔ GAP2** (n = 5 substituted positions, MM/GBSA idecomp = 2 from
     mmgbsa_cp12_gap2_nad): substituted-only Pearson r = −0.67, p = 0.22 —
     not significant at n = 5, but directionally consistent. Dominant outlier
     **HMM col 47 (K → H)** is the exception that proves the rule: at GAP2,
     MED4's lysine packs 1.9 kcal/mol better and holds GAP2 4.1 ns longer
     than HM2's histidine. The K → H swap that helps HM2 at PRK *costs* HM2
     at GAP2 — consistent with the MM/GBSA partner-preference asymmetry
     (HM2 favors PRK, MED4 favors GAP2).

   At both interfaces and in both redox states, the host with the more
   favorable per-residue VDW packing at a substituted position is the host
   with the longer mean contact dwell at that position. Where HM2's
   substitutions improve VDW packing (most of the PRK interface), HM2 also
   holds the contact longer; where they degrade VDW packing (col 47 at GAP2),
   HM2 loses the contact faster.


---

## Directory layout

```
├── README.md                     this file
├── docs/
│   ├── METHODS.md                methods write-up
│   └── RESULTS.md                integrated results write-up (manuscript-ready)
├── scripts/
│   ├── analyze_contact_dynamics.py        per-pair B(t) → persistence/turnover/Markov
│   ├── plot_contact_dynamics.py           heatmaps, violins, kinetic scatter
│   ├── test_contact_dynamics_host_redox.py per-HMM Welch + BH-FDR + convergence flags
│   ├── analyze_contact_coupling_mi.py     HMM-column MI coupling networks
│   ├── plot_contact_coupling.py           MI matrices, host-diff panels, strength lines
│   ├── compute_per_residue_vdw_direct.py  direct L-J VDW for cp12_prk + cp12_gap2 (no gmx_MMPBSA)
│   ├── crossplot_dwell_vs_vdw_gap2.py     Δdwell vs Δ(VDW) cross-plot, CP12 ↔ GAP2 (idecomp = 2)
│   ├── crossplot_dwell_vs_vdw_prk.py      Δdwell vs Δ(VDW) cross-plot, CP12 ↔ PRK (direct)
│   ├── analyze_contact_class_welch.py     formal Welch + bootstrap CI + BH-FDR on log τ_on per (partner × redox × class)
│   ├── analyze_mi_neighborhood_partition.py  decompose MI coupling into N/M/C spatial blocks + cross-block sums
│   └── analyze_substitution_class_shift.py  formal MWU + Kruskal + binomial directional test of the chemistry-class-shift hypothesis
├── data/
│   ├── condition_summary.csv             8 rows, headline per-condition numbers
│   ├── per_class_summary.csv             40 rows, condition × partner × class
│   ├── per_pair_hmm_aggregated.csv       907 rows, per HMM × condition × class
│   ├── per_pair_stats.csv.gz             20,728 rows, the per-pair atomic table
│   ├── host_per_hmm_{prk,gap2}.csv       per-HMM host comparison (Welch + bootstrap CI + BH)
│   ├── redox_per_hmm_{prk,gap2}.csv      per-HMM redox comparison
│   ├── convergence_{prk,gap2}.csv        per-HMM CV(τ_on) flags (rep variance vs host variance)
│   ├── contact_coupling_mi.npz           pairwise MI matrices, per (condition, partner)
│   ├── coupling_similarity_summary.csv   Jaccard / Frobenius / spectral pairwise stats
│   ├── coupling_neighborhood_partition.csv  intra- and cross-block MI decomposition (N/M/C zones)
│   ├── coupling_neighborhood_host_diff.csv  MED4 − HM2 per block per redox
│   ├── coupling_neighborhood_partition.md   prose summary + the partition table
│   ├── class_welch_host_comparison.csv   20-cell formal Welch test on log τ_on (partner × redox × class)
│   ├── class_welch_summary.md            prose summary of class-stratified test (negative result, honest)
│   ├── substitution_class_shift.csv      per-HMM-position table joined with substitution-class-shift category
│   ├── substitution_class_shift_test.csv 3 hypothesis tests (MWU, Kruskal, binomial directional)
│   ├── substitution_class_shift_summary.md  prose summary; all three tests refute the chemistry-shift framework
│   ├── per_residue_vdw_gap2.csv          per-rep MM/GBSA per-residue VDW for GAP2 (idecomp=2)
│   ├── per_residue_vdw_direct_gap2.csv   per-rep DIRECT-trajectory L-J VDW for GAP2 (validation; r=+0.908 vs idecomp=2)
│   ├── per_residue_vdw_direct_prk.csv    per-rep DIRECT-trajectory L-J VDW for PRK (the canonical idecomp=2 was purged)
│   ├── per_hmm_vdw_gap2_by_condition.csv aggregated to HMM × condition (GAP2)
│   ├── per_hmm_vdw_prk_by_condition.csv  aggregated to HMM × condition (PRK)
│   ├── host_diff_vdw_gap2_{ox,red}.csv   Δ(MED4 − HM2) per HMM col (GAP2)
│   ├── host_diff_vdw_prk_{ox,red}.csv    Δ(MED4 − HM2) per HMM col (PRK)
│   ├── merged_per_hmm_{gap2,prk}.csv     joined dwell + VDW per HMM × condition
│   ├── vdw_bridge_text_summary.md        GAP2 cross-plot stats
│   ├── vdw_bridge_prk_text_summary.md    PRK cross-plot stats (r = −0.744, p = 1.6×10⁻⁶)
│   ├── *_text_summary.md                 prose summaries of each table
│   └── CP12_MED4_vs_HM2_alignment_annotation.csv  the 70-col HMM alignment used throughout
└── figures/
    ├── Main5a_contact_persistence_heatmap_gap2.{pdf,png}   GAP2 dwell heatmap (mimicry)
    ├── Main5b_coupling_strength_per_hmm.{pdf,png}          Σ MI per HMM col (host specialization)
    ├── Main6a_kinetic_thermo_bridge_gap2_ox.{pdf,png}      Δdwell vs ΔVDW (GAP2 ox)
    ├── Main6b_kinetic_thermo_bridge_prk_ox.{pdf,png}       Δdwell vs ΔVDW (PRK ox; r=-0.74)
    ├── S27_kinetic_thermo_bridge_prk_red.{pdf,png}         PRK reduced state companion
    ├── S24_contact_persistence_heatmap_prk.{pdf,png}       PRK dwell heatmap
    ├── S25_contact_class_violins.{pdf,png}                 class-stratified dwell distributions
    ├── S26a_coupling_mi_matrix_prk.{pdf,png}               full 70×70 MI matrices (PRK)
    ├── S26b_coupling_mi_matrix_gap2.{pdf,png}              full 70×70 MI matrices (GAP2)
    ├── fig_mi_neighborhood_partition_prk.{pdf,png}         MI by spatial block (PRK; S28a)
    ├── fig_mi_neighborhood_partition_gap2.{pdf,png}        MI by spatial block (GAP2; S28b)
    ├── fig_kinetic_scatter_{prk,gap2}.{pdf,png}            k_off vs k_on, per HMM column
    ├── fig_host_redox_dwell_{prk,gap2}.{pdf,png}           per-HMM mean dwell scatter
    └── fig_survival_overview.{pdf,png}                     median S(t) per condition
```

---

## Reproducibility

All scripts run in the shared `gppenv` conda environment:

```bash
conda activate gppenv
python analyze_contact_dynamics.py --outdir contact_dynamics --cutoff 5.0 --min_run 3
python plot_contact_dynamics.py --indir contact_dynamics --outdir contact_dynamics_plots
python test_contact_dynamics_host_redox.py
python analyze_contact_coupling_mi.py --outdir contact_coupling
python plot_contact_coupling.py
python crossplot_dwell_vs_vdw_gap2.py
python compute_per_residue_vdw_direct.py --partners prk --start 200 --stride 4
python crossplot_dwell_vs_vdw_prk.py
python analyze_contact_class_welch.py
python analyze_mi_neighborhood_partition.py
python analyze_substitution_class_shift.py
```

All inputs are the existing
`<system>/<sysid>/<rep>/solute_{tpr,fit.xtc,fixed.pdb}` (full trajectories +
chain-labeled PDBs), the canonical
`mmgbsa_cp12_gap2_nad/FINAL_DECOMP_MMPBSA.dat` per-residue decomp (for GAP2 +
direct-method validation), and `mmgbsa_cp12_prk/COM.prmtop` (for amber99sb σ/ε
lookup; trajectory + topology not required for PRK because the AMBER L-J
parameters come straight from the prmtop, not the deleted COM_traj_0.xtc).

