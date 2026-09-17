# Substitution-class-shift hypothesis test

**Hypothesis**: the host-asymmetry in contact kinetics is concentrated at HMM positions where the substitution changes the residue chemistry class (e.g., hydrophobic↔charged); positions where the substitution preserves class (e.g., K↔R) and conserved positions should show smaller dwell-time differences between hosts.

Operationalized as three tests on the per-HMM log-dwell host comparison (from the existing `host_per_hmm_{prk,gap2}.csv`):

1. **MWU substituted vs conserved**: one-sided Mann-Whitney on |Δlog-dwell(HM2 − MED4)|. Expects substituted > conserved.
2. **Kruskal-Wallis across shift categories**: within substituted positions only, do the different substitution-shift categories carry different effect sizes?
3. **Binomial directional test**: for substitutions where one host is predicted to be the dwell winner (charged residue → predicted host = the charged-residue host; salt bridges are kinetically stronger), is the observed sign of Δlog-dwell more often the predicted sign than chance?

## Test 1: MWU substituted vs conserved (|Δlog τ_on|)

| partner | redox | n_subs | n_cons | med abs Δlog subs | med abs Δlog cons | U | p (one-sided) |
|---|---|---|---|---|---|---|---|
| prk | ox | 30 | 21 | 0.141 | 0.134 | 347 | 0.273 |
| prk | red | 29 | 23 | 0.154 | 0.165 | 358 | 0.329 |

## Test 2: Kruskal-Wallis across substitution-class shift categories

| partner | redox | categories (n) | n total | K | p (two-sided) |
|---|---|---|---|---|---|
| prk | ox | hydrophobic↔charged(4),hydrophobic↔polar(3),polar↔charged(8),within-class(11) | 26 | 1.23 | 0.745 |
| prk | red | hydrophobic↔charged(4),hydrophobic↔polar(3),polar↔charged(8),within-class(10) | 25 | 0.67 | 0.88 |

## Test 3: Binomial directional test (predicted-winner sign match)

Substitutions where exactly one host carries a charged residue: predict that host to be the dwell winner (charged-residue host → longer dwell, via salt-bridge kinetics). Test whether observed sign of Δlog matches prediction more often than chance.

| partner | redox | n preds | n matching | fraction matching | binomial p |
|---|---|---|---|---|---|
| prk | ox | 12 | 6 | 0.50 | 0.613 |
| prk | red | 12 | 3 | 0.25 | 0.981 |

## Per-substitution detail (predicted-winner cases)

### PRK OX
col6:Q>E→hm2(✗)
- col8:D>H→med4(✓)
- col10:S>E→hm2(✓)
- col12:I>D→hm2(✗)
- col16:K>L→med4(✓)
- col25:K>A→med4(✓)
- col26:V>R→hm2(✓)
- col30:T>K→hm2(✓)
- col34:K>H→med4(✗)
- col40:K>H→med4(✗)
- col41:D>N→med4(✗)
- col47:K>H→med4(✗)

### PRK RED
col6:Q>E→hm2(✗)
- col8:D>H→med4(✗)
- col10:S>E→hm2(✗)
- col12:I>D→hm2(✗)
- col16:K>L→med4(✓)
- col25:K>A→med4(✗)
- col26:V>R→hm2(✓)
- col30:T>K→hm2(✓)
- col34:K>H→med4(✗)
- col40:K>H→med4(✗)
- col41:D>N→med4(✗)
- col47:K>H→med4(✗)

## Substitution-class shift breakdown

Shift-category counts across the 70 HMM columns:

- **conserved**: 34
- **within-class**: 12
- **polar↔charged**: 8
- **gap**: 5
- **hydrophobic↔charged**: 4
- **hydrophobic↔polar**: 3
- **charge-flip**: 2
- **other**: 2

## Substituted positions per shift category (annotated table)

| HMM col | MED4 | HM2 | shift category | predicted winner |
|---|---|---|---|---|
| 4 | H | T | within-class | — |
| 6 | Q | E | polar↔charged | hm2 |
| 7 | K | E | charge-flip | — |
| 8 | D | H | polar↔charged | med4 |
| 9 | Q | I | hydrophobic↔polar | — |
| 10 | S | E | polar↔charged | hm2 |
| 11 | E | K | charge-flip | — |
| 12 | I | D | hydrophobic↔charged | hm2 |
| 13 | A | H | hydrophobic↔polar | — |
| 14 | S | A | hydrophobic↔polar | — |
| 15 | A | I | within-class | — |
| 16 | K | L | hydrophobic↔charged | med4 |
| 17 | E | D | within-class | — |
| 18 | E | D | within-class | — |
| 21 | G | T | other | — |
| 22 | N | S | within-class | — |
| 23 | L | P | within-class | — |
| 24 | P | A | within-class | — |
| 25 | K | A | hydrophobic↔charged | med4 |
| 26 | V | R | hydrophobic↔charged | hm2 |
| 29 | L | Y | within-class | — |
| 30 | T | K | polar↔charged | hm2 |
| 34 | K | H | polar↔charged | med4 |
| 35 | D | E | within-class | — |
| 38 | G | V | other | — |
| 40 | K | H | polar↔charged | med4 |
| 41 | D | N | polar↔charged | med4 |
| 47 | K | H | polar↔charged | med4 |
| 58 | D | E | within-class | — |
| 59 | A | M | within-class | — |
| 60 | N | H | within-class | — |
