#!/usr/bin/env python3
"""
Substitution-class-shift hypothesis test.

Collaborator's hypothesis: at HMM positions where the substitution changes the
residue chemistry class (e.g., MED4-hydrophobic → HM2-charged), the contact
kinetics at that position should shift in the direction predicted by the new
class (e.g., longer dwell, lower k_off — characteristic of salt bridges).

Operationalization:
  1. For each HMM column, classify the (MED4_aa, HM2_aa) pair into a
     substitution-class-shift category:
         conserved              same amino acid
         within-class           different AA but same chemistry class
         hydrophobic↔polar     hydrophobic ↔ polar/H-bond
         hydrophobic↔charged   hydrophobic ↔ {D,E,K,R}
         polar↔charged         polar/Hbond ↔ charged
         charge-flip            +charge ↔ −charge (K/R ↔ D/E)
         charged↔aromatic-π    {K,R} ↔ {F,W,Y}
         other                  mismatched or gap
  2. Read the existing per-HMM log-dwell host comparison from
     `contact_dynamics_stats/host_per_hmm_{prk,gap2}.csv`. Each row is a
     (HMM col × redox) cell with log_dwell_diff_hm2_med4, bootstrap CI, etc.
  3. Mann-Whitney one-sided test (alternative='greater') of |Δlog-dwell| at
     substituted positions vs conserved positions. Expect substituted > conserved.
  4. Kruskal-Wallis across shift categories within substituted positions.
  5. Per-category direction-of-effect test: for substitutions where the
     "expected winner" is one host (e.g., hydrophobic→charged: HM2 acquires
     salt-bridge → expected HM2 longer dwell → expected Δlog>0), test whether
     Δlog-dwell is more often the predicted sign than the opposite.

Outputs:
    contact_dynamics_stats/substitution_class_shift.csv
    contact_dynamics_stats/substitution_class_shift_summary.md
    contact_dynamics_stats/substitution_class_shift_test.csv  (test statistics)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HYDROPHOBIC = set("AVLIMFWYP")
POLAR_HBOND = set("STNQHC")  # H counted as polar here; class checks separately for cation-π
POS_CHARGED = set("KR")
NEG_CHARGED = set("DE")
AROMATIC    = set("FWY")


def aa_class(aa: str) -> str:
    if aa in ("-", "", "X"):
        return "gap"
    if aa in POS_CHARGED:
        return "pos"
    if aa in NEG_CHARGED:
        return "neg"
    if aa in HYDROPHOBIC:
        return "hydrophobic"
    if aa in POLAR_HBOND:
        return "polar"
    return "other"


def shift_category(m: str, h: str) -> str:
    if m == h:
        return "conserved"
    cm = aa_class(m); ch = aa_class(h)
    if cm == "gap" or ch == "gap":
        return "gap"
    if cm == ch:
        return "within-class"
    # mismatched classes
    s = {cm, ch}
    if s == {"hydrophobic", "polar"}:
        return "hydrophobic↔polar"
    if s == {"hydrophobic", "pos"} or s == {"hydrophobic", "neg"}:
        return "hydrophobic↔charged"
    if s == {"polar", "pos"} or s == {"polar", "neg"}:
        return "polar↔charged"
    if s == {"pos", "neg"}:
        return "charge-flip"
    # cation-π: aromatic (in hydrophobic) ↔ pos; this overlaps hydrophobic↔charged
    # but flagged separately if we detect aromatic specifically
    if (m in AROMATIC and h in POS_CHARGED) or (h in AROMATIC and m in POS_CHARGED):
        return "aromatic↔charged"
    return "other"


# For directional test: which host is "predicted to have longer dwell"?
# Hypothesis: a residue in the salt-bridge / cation-π class makes longer contacts
# than a hydrophobic-class residue at the same position (because the new charged
# residue forms persistent salt bridges; the original hydrophobic packing was
# weaker).
def predicted_winner(m: str, h: str) -> str | None:
    """Return 'med4', 'hm2', or None for no prediction."""
    if m == h:
        return None
    cm = aa_class(m); ch = aa_class(h)
    if cm == "gap" or ch == "gap":
        return None
    # If only one is charged, predict THAT host to have longer dwell
    m_charged = cm in ("pos", "neg")
    h_charged = ch in ("pos", "neg")
    if m_charged and not h_charged:
        return "med4"
    if h_charged and not m_charged:
        return "hm2"
    # If only one is hydrophobic and the other is polar (no charge), no strong prediction
    return None


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    if mask.sum() == 0:
        return q
    pv = p[mask]; n = pv.size
    order = np.argsort(pv); ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    qraw = pv * n / ranks
    qs = qraw[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    qf = np.empty_like(qs); qf[order] = np.clip(qs, 0, 1)
    q[mask] = qf
    return q


def load_alignment(aln_csv: Path) -> pd.DataFrame:
    aln = pd.read_csv(aln_csv).rename(columns={"aln_col_1based": "hmm_col"})
    aln["med4_aa"] = aln["med4_aa"].fillna("-")
    aln["hm2_aa"] = aln["hm2_aa"].fillna("-")
    aln["shift_category"] = aln.apply(lambda r: shift_category(r["med4_aa"], r["hm2_aa"]), axis=1)
    aln["predicted_winner"] = aln.apply(lambda r: predicted_winner(r["med4_aa"], r["hm2_aa"]), axis=1)
    aln["is_substituted"] = (aln["med4_aa"] != aln["hm2_aa"]) & (aln["med4_aa"] != "-") & (aln["hm2_aa"] != "-")
    return aln


def per_position_data(stats_dir: Path, aln: pd.DataFrame) -> pd.DataFrame:
    """Join per-HMM Welch results with substitution metadata."""
    rows = []
    for partner in ("prk", "gap2"):
        p = stats_dir / f"host_per_hmm_{partner}.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p)
        df["partner"] = partner
        rows.append(df)
    if not rows:
        raise SystemExit("no host_per_hmm CSVs found")
    df = pd.concat(rows, ignore_index=True)
    df = df.merge(aln[["hmm_col", "med4_aa", "hm2_aa", "shift_category",
                        "predicted_winner", "is_substituted"]],
                  on="hmm_col", how="left")
    df["abs_log_diff"] = df["log_dwell_diff_hm2_med4"].abs()
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats_dir", default="contact_dynamics_stats")
    ap.add_argument("--alignment", default="CP12_MED4_vs_HM2_alignment_annotation.csv")
    ap.add_argument("--outdir", default="contact_dynamics_stats")
    args = ap.parse_args()

    stats_dir = Path(args.stats_dir).resolve()
    outdir = Path(args.outdir).resolve()
    aln = load_alignment(Path(args.alignment).resolve())
    df = per_position_data(stats_dir, aln)
    df.to_csv(outdir / "substitution_class_shift.csv", index=False)

    # Test 1: Mann-Whitney one-sided, |Δlog| substituted vs conserved (per partner × redox)
    tests = []
    for partner in ("prk", "gap2"):
        for redox in ("ox", "red"):
            sub = df[(df["partner"] == partner) & (df["redox"] == redox)].copy()
            if sub.empty:
                continue
            subs = sub.loc[sub["is_substituted"], "abs_log_diff"].dropna().to_numpy()
            cons = sub.loc[~sub["is_substituted"] & sub["med4_aa"].ne("-") & sub["hm2_aa"].ne("-"),
                           "abs_log_diff"].dropna().to_numpy()
            if subs.size < 5 or cons.size < 5:
                continue
            mwu = stats.mannwhitneyu(subs, cons, alternative="greater")
            tests.append({
                "test": "MWU_substituted_vs_conserved_abs_log_dwell",
                "partner": partner, "redox": redox,
                "n_substituted": int(subs.size), "n_conserved": int(cons.size),
                "median_abs_log_dwell_substituted": float(np.median(subs)),
                "median_abs_log_dwell_conserved":   float(np.median(cons)),
                "U_stat": float(mwu.statistic), "p_one_sided": float(mwu.pvalue),
            })

    # Test 2: Kruskal-Wallis across shift categories (within substituted positions)
    for partner in ("prk", "gap2"):
        for redox in ("ox", "red"):
            sub = df[(df["partner"] == partner) & (df["redox"] == redox) & df["is_substituted"]].copy()
            if sub.empty:
                continue
            groups: list[np.ndarray] = []
            labels: list[str] = []
            for cat, g in sub.groupby("shift_category"):
                vals = g["abs_log_diff"].dropna().to_numpy()
                if vals.size >= 3:
                    groups.append(vals); labels.append(cat)
            if len(groups) < 2:
                continue
            k, p = stats.kruskal(*groups)
            tests.append({
                "test": "Kruskal_across_shift_categories_within_substituted",
                "partner": partner, "redox": redox,
                "categories_tested": ",".join(f"{l}({len(g)})" for l, g in zip(labels, groups)),
                "n_substituted_with_data": int(sum(len(g) for g in groups)),
                "K_stat": float(k), "p_two_sided": float(p),
            })

    # Test 3: Directional sign test on predicted-winner substitutions
    # For substitutions where one host is predicted to have longer dwell, test whether
    # Δlog-dwell sign matches prediction more often than chance (binomial p)
    for partner in ("prk", "gap2"):
        for redox in ("ox", "red"):
            sub = df[(df["partner"] == partner) & (df["redox"] == redox)
                      & df["is_substituted"] & df["predicted_winner"].notna()
                      & df["log_dwell_diff_hm2_med4"].notna()].copy()
            if sub.empty:
                continue
            # Δlog-dwell = log(hm2) - log(med4). Positive = hm2 longer dwell.
            # If predicted_winner='hm2', the prediction is Δlog > 0; if 'med4', Δlog < 0.
            sub["predicted_sign"] = np.where(sub["predicted_winner"] == "hm2", 1, -1)
            sub["observed_sign"] = np.sign(sub["log_dwell_diff_hm2_med4"])
            sub["matches"] = sub["predicted_sign"] == sub["observed_sign"]
            n = int(len(sub))
            k = int(sub["matches"].sum())
            if n >= 3:
                # binomial p, one-sided greater-than-chance
                binom = stats.binomtest(k, n, p=0.5, alternative="greater")
                tests.append({
                    "test": "Binomial_directional_predicted_winner_match",
                    "partner": partner, "redox": redox,
                    "n_substitutions_with_prediction": n,
                    "n_matching_prediction": k,
                    "fraction_matching": k / n if n else float("nan"),
                    "p_one_sided": float(binom.pvalue),
                    "predicted_winners": ",".join(
                        f"col{r.hmm_col}:{r.med4_aa}>{r.hm2_aa}→{r.predicted_winner}({'✓' if r.matches else '✗'})"
                        for r in sub.itertuples()),
                })

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(outdir / "substitution_class_shift_test.csv", index=False)

    # Markdown summary
    lines = ["# Substitution-class-shift hypothesis test",
             "",
             "**Hypothesis**: the host-asymmetry in contact kinetics is concentrated at "
             "HMM positions where the substitution changes the residue chemistry class "
             "(e.g., hydrophobic↔charged); positions where the substitution preserves "
             "class (e.g., K↔R) and conserved positions should show smaller dwell-time "
             "differences between hosts.",
             "",
             "Operationalized as three tests on the per-HMM log-dwell host comparison "
             "(from the existing `host_per_hmm_{prk,gap2}.csv`):",
             "",
             "1. **MWU substituted vs conserved**: one-sided Mann-Whitney on "
             "|Δlog-dwell(HM2 − MED4)|. Expects substituted > conserved.",
             "2. **Kruskal-Wallis across shift categories**: within substituted "
             "positions only, do the different substitution-shift categories carry "
             "different effect sizes?",
             "3. **Binomial directional test**: for substitutions where one host is "
             "predicted to be the dwell winner (charged residue → predicted host = "
             "the charged-residue host; salt bridges are kinetically stronger), is the "
             "observed sign of Δlog-dwell more often the predicted sign than chance?",
             "",
             "## Test 1: MWU substituted vs conserved (|Δlog τ_on|)",
             "",
             "| partner | redox | n_subs | n_cons | med abs Δlog subs | med abs Δlog cons | U | p (one-sided) |",
             "|---|---|---|---|---|---|---|---|"]
    t1 = [t for t in tests if t["test"] == "MWU_substituted_vs_conserved_abs_log_dwell"]
    for t in t1:
        lines.append(
            f"| {t['partner']} | {t['redox']} | {t['n_substituted']} | {t['n_conserved']} | "
            f"{t['median_abs_log_dwell_substituted']:.3f} | "
            f"{t['median_abs_log_dwell_conserved']:.3f} | "
            f"{t['U_stat']:.0f} | {t['p_one_sided']:.3g} |"
        )
    lines.append("")

    lines += ["## Test 2: Kruskal-Wallis across substitution-class shift categories",
              "",
              "| partner | redox | categories (n) | n total | K | p (two-sided) |",
              "|---|---|---|---|---|---|"]
    t2 = [t for t in tests if t["test"].startswith("Kruskal_across_shift")]
    for t in t2:
        lines.append(
            f"| {t['partner']} | {t['redox']} | {t['categories_tested']} | "
            f"{t['n_substituted_with_data']} | {t['K_stat']:.2f} | {t['p_two_sided']:.3g} |"
        )
    lines.append("")

    lines += ["## Test 3: Binomial directional test (predicted-winner sign match)",
              "",
              "Substitutions where exactly one host carries a charged residue: predict "
              "that host to be the dwell winner (charged-residue host → longer dwell, "
              "via salt-bridge kinetics). Test whether observed sign of Δlog matches "
              "prediction more often than chance.",
              "",
              "| partner | redox | n preds | n matching | fraction matching | binomial p |",
              "|---|---|---|---|---|---|"]
    t3 = [t for t in tests if t["test"].startswith("Binomial_directional")]
    for t in t3:
        lines.append(
            f"| {t['partner']} | {t['redox']} | {t['n_substitutions_with_prediction']} | "
            f"{t['n_matching_prediction']} | {t['fraction_matching']:.2f} | "
            f"{t['p_one_sided']:.3g} |"
        )
    lines.append("")

    lines += ["## Per-substitution detail (predicted-winner cases)",
              ""]
    for t in t3:
        lines.append(f"### {t['partner'].upper()} {t['redox'].upper()}")
        lines.append(t["predicted_winners"].replace(",", "\n- "))
        lines.append("")

    lines += ["## Substitution-class shift breakdown",
              "",
              "Shift-category counts across the 70 HMM columns:",
              ""]
    cat_counts = aln.groupby("shift_category").size().sort_values(ascending=False)
    for cat, n in cat_counts.items():
        lines.append(f"- **{cat}**: {n}")
    lines.append("")
    lines += ["## Substituted positions per shift category (annotated table)", ""]
    sub_rows = aln[aln["is_substituted"] & (aln["shift_category"] != "gap")].copy()
    lines.append("| HMM col | MED4 | HM2 | shift category | predicted winner |")
    lines.append("|---|---|---|---|---|")
    for _, r in sub_rows.iterrows():
        pred = r["predicted_winner"] if pd.notna(r["predicted_winner"]) else "—"
        lines.append(f"| {int(r['hmm_col'])} | {r['med4_aa']} | {r['hm2_aa']} | "
                     f"{r['shift_category']} | {pred} |")

    (outdir / "substitution_class_shift_summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {outdir / 'substitution_class_shift_summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
