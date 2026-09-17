# Binary contact dynamics — methods and results scaffold

## Methods (to drop into MANUSCRIPT.md §2)

For each replicate (n=3 per condition, 4 conditions: med4_ox, med4_red, hm2_ox, hm2_red),
we extracted a per-(CP12_residue, partner_residue) binary contact indicator
B(t) ∈ {0,1} from the fitted production trajectory (`solute_fit.xtc`, 1001 frames at
100 ps/frame). A contact is defined by the existence of any heavy-atom pair within
5.0 Å between the CP12 residue and the partner residue (PRK or GAP2) in frame t,
computed with MDAnalysis (v2.x) FastNS under periodic boundary conditions.

To suppress thermal flicker near the cutoff, B(t) was post-processed with a
3-frame (300 ps) morphological open–close filter: short 1-runs of length < 3
frames are reassigned to 0 (opening), then short 0-gaps of length < 3 frames are
reassigned to 1 (closing). The 300 ps minimum is conservative relative to typical
sidechain rotamer relaxation (~100 ps) and exceeds rotational tumbling between
adjacent frames. The filter is equivalent in spirit to a 4.0 Å form / 5.5 Å break
hysteresis but operates on the cleaner pre-existing 5.0 Å indicator.

For every CP12 chain copy (C/H/M/P) and partner-residue partner that ever made
contact in any frame, the cleaned B(t) was run-length-encoded to extract:
- **occupancy** = ⟨B(t)⟩ — fraction of frames with contact (matches the
  pre-existing per-residue contact-frequency metric in `interface_results/redox/`).
- **τ_on (mean, median, p90, max)** — distribution of consecutive-1-run lengths,
  in frames × 100 ps/frame.
- **τ_off (mean, median, p90)** — distribution of consecutive-0-run lengths.
- **turnover rate** = (0→1 + 1→0 transitions) / total time, in events per ns.
- **k_off** = 1 / mean_τ_on (1/ns) — two-state Markov off-rate ("contact strength":
  lower k_off ⇒ stronger contact).
- **k_on** = 1 / mean_τ_off (1/ns).
- **Survival fractions S(1 ns), S(5 ns), S(20 ns)** — fraction of observed 1-runs
  that persist for at least 10, 50, or 200 frames.

Each (CP12_residue, partner_residue) pair was assigned a contact chemistry class:
- **hydrophobic** — both sidechains in {A, V, L, I, M, F, W, Y, P};
- **salt-bridge** — pair contains both a {D, E} and a {K, R};
- **cation-π** — pair contains both an aromatic {F, W, Y} and a {K, R};
- **polar/H-bond** — at least one residue in {S, T, N, Q, H} and no salt-bridge or
  cation-π already assigned;
- **mixed** — otherwise.

Per-(CP12 HMM column, partner, contact class) statistics were aggregated by
averaging across the 4 CP12 chain copies × 3 replicates × all matching partner
residues for that condition; standard errors are reported with n = number of
(chain, replicate, partner-residue) observations contributing to each cell.
Bootstrap 95% CIs on median dwell time were obtained from 2000 resampling
iterations across the underlying pair-level rows.

Replicate convergence was checked by comparing the within-condition coefficient
of variation of mean τ_on against the between-condition spread; HMM columns where
the within-condition CV exceeded the between-condition spread are flagged as
kinetically undersampled in the supporting figures.

## Results scaffold (drop into §3.x)

(Numbers populated from `contact_dynamics_plots/contact_dynamics_text_summary.md`
after the analysis run completes.)

> **Headline finding placeholder:** Across the four conditions, mean per-residue
> contact dwell at the CP12 interface ranges from ~X ps (transient electrostatic
> flickers) to ~Y ns (stable salt bridges and core hydrophobic packing). Median
> τ_on per HMM column is within 95% bootstrap CI between MED4 and HM2 at fixed
> redox state for [N] / 70 columns (host mimicry holds kinetically). Where they
> differ, the divergent positions concentrate in [region], coincident with
> physicochemical-class changes in the HMM substitution table.

> **Class-stratified test of "hydrophobic ↔ polar/charged trade-off":** At HMM
> columns where MED4 carries a hydrophobic residue but HM2 carries a polar or
> charged residue (N=…), the median τ_on shifts from … ps (MED4) to … ps (HM2),
> and the contact distribution moves from the **hydrophobic** class shape
> (broad, heavy-tailed) toward the **salt-bridge** class shape (narrower, lower
> k_off when on). This is the kinetic signature predicted by the collaborator's
> short-range / long-range component-balance argument.

> **Cross-axis check:** For the same HMM-substituted positions, Δ(mean τ_on)
> correlates with Δ(per-residue VDW)+Δ(per-residue EEL) from MM/GBSA
> decomposition (Pearson r = …, n = 31). This couples the kinetic and
> thermodynamic evidence at the per-residue level.

> **Convergence:** 100 ns × 3 replicates is sufficient (within-condition CV of
> median τ_on < between-condition spread) for [M] / 70 HMM columns; [70-M] are
> flagged as undersampled and excluded from the host-comparison test.

## Output inventory

```
contact_dynamics/
  per_pair_stats.csv.gz          # one row per (condition,rep,partner,CP12_chain,
                                  #              CP12_resid,partner_chain,partner_resid)
  per_pair_hmm_aggregated.csv    # rolled up to (condition, hmm_col, partner, class)
  per_class_summary.csv          # rolled up to (condition, partner, class)
  condition_summary.csv          # rolled up to (condition, partner)

contact_dynamics_plots/
  fig_persistence_heatmap_<partner>.{pdf,png}
  fig_kinetic_scatter_<partner>.{pdf,png}
  fig_host_redox_dwell_<partner>.{pdf,png}
  fig_class_violins.{pdf,png}
  fig_class_kinetics_strip.{pdf,png}
  fig_survival_overview.{pdf,png}
  contact_dynamics_text_summary.md
```
