#!/usr/bin/env python3
"""
Class-stratified Welch tests of CP12 ↔ partner contact dwell time, MED4 vs HM2.

For each (partner ∈ {prk, gap2}) × (redox ∈ {ox, red}) × (chemistry class ∈
{hydrophobic, polar_hbond, salt_bridge, cation_pi, mixed}), run a two-sided
Welch t-test on log10(mean_on_ps) comparing all MED4 pairs to all HM2 pairs.
Report:
  - n_med4, n_hm2
  - median dwell per host (ps)
  - log-dwell mean and SD per host
  - Welch t, p, df
  - bootstrap-percentile 95% CI on median(HM2) − median(MED4)
  - effect size d = (mean_log_HM2 − mean_log_MED4) / pooled_sd
BH-FDR correction across the full 20-test family.

Outputs:
    contact_dynamics_stats/class_welch_host_comparison.csv
    contact_dynamics_stats/class_welch_summary.md

The pool-across-pairs design is the same as the existing per-HMM Welch test
(host_per_hmm_*.csv) — each (chain × rep × cp12_resid × partner_resid) row is
an independent observation. This inflates n compared with a strict
per-replicate design (n=3 per condition), so the bootstrap CI and BH-FDR are
the conservative anchors of significance.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PARTNERS = ["prk", "gap2"]
REDOXES = ["ox", "red"]
CLASSES = ["hydrophobic", "polar_hbond", "salt_bridge", "cation_pi", "mixed"]


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, NaN-aware."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    if mask.sum() == 0:
        return q
    pv = p[mask]
    n = pv.size
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    qraw = pv * n / ranks
    qs = qraw[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    qf = np.empty_like(qs)
    qf[order] = np.clip(qs, 0, 1)
    q[mask] = qf
    return q


def welch_block(df: pd.DataFrame, partner: str, redox: str, contact_class: str,
                 n_boot: int = 5000, rng: np.random.Generator | None = None) -> dict:
    """Welch test on log10(mean_on_ps) for one (partner, redox, class) cell."""
    if rng is None:
        rng = np.random.default_rng(0)
    sub = df[
        (df["partner"] == partner)
        & (df["condition"].str.endswith("_" + redox))
        & (df["contact_class"] == contact_class)
        & df["mean_on_ps"].notna()
        & (df["mean_on_ps"] > 0)
    ].copy()
    sub["log_dwell"] = np.log10(sub["mean_on_ps"])
    sub["host"] = sub["condition"].str[:4].str.rstrip("_")

    m = sub.loc[sub["host"] == "med4", "log_dwell"].to_numpy()
    h = sub.loc[sub["host"] == "hm2",  "log_dwell"].to_numpy()
    n_m, n_h = m.size, h.size
    if n_m < 5 or n_h < 5:
        return dict(
            partner=partner, redox=redox, contact_class=contact_class,
            n_med4=int(n_m), n_hm2=int(n_h),
            med4_median_ps=float("nan"), hm2_median_ps=float("nan"),
            med4_log_mean=float("nan"), hm2_log_mean=float("nan"),
            welch_t=float("nan"), welch_p=float("nan"), welch_df=float("nan"),
            log_diff_hm2_minus_med4=float("nan"),
            ci_lo_log=float("nan"), ci_hi_log=float("nan"),
            cohens_d=float("nan"),
            log_diff_ps_factor=float("nan"),
        )

    t_stat = stats.ttest_ind(h, m, equal_var=False, nan_policy="omit")
    # bootstrap CI on median(log_h) - median(log_m)
    diffs = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        mb = rng.choice(m, size=n_m)
        hb = rng.choice(h, size=n_h)
        diffs[i] = np.median(hb) - np.median(mb)
    lo, hi = np.percentile(diffs, [2.5, 97.5])

    # Cohen's d (pooled SD across hosts)
    sd_pool = np.sqrt(((n_m - 1) * np.var(m, ddof=1) + (n_h - 1) * np.var(h, ddof=1))
                       / max(n_m + n_h - 2, 1))
    d = (np.mean(h) - np.mean(m)) / sd_pool if sd_pool > 0 else float("nan")

    return dict(
        partner=partner, redox=redox, contact_class=contact_class,
        n_med4=int(n_m), n_hm2=int(n_h),
        med4_median_ps=float(10 ** np.median(m)),
        hm2_median_ps =float(10 ** np.median(h)),
        med4_log_mean=float(np.mean(m)),
        hm2_log_mean =float(np.mean(h)),
        welch_t=float(t_stat.statistic),
        welch_p=float(t_stat.pvalue),
        welch_df=float(t_stat.df) if hasattr(t_stat, "df") else float("nan"),
        log_diff_hm2_minus_med4=float(np.median(h) - np.median(m)),
        ci_lo_log=float(lo), ci_hi_log=float(hi),
        cohens_d=float(d),
        # multiplicative factor on raw dwell time (10^Δlog)
        log_diff_ps_factor=float(10 ** (np.median(h) - np.median(m))),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per_pair",
                    default="contact_dynamics/per_pair_stats.csv.gz")
    ap.add_argument("--outdir", default="contact_dynamics_stats")
    args = ap.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.per_pair)
    print(f"loaded {len(df):,} per-pair rows", flush=True)

    rows = []
    rng = np.random.default_rng(42)
    for partner in PARTNERS:
        for redox in REDOXES:
            for cls in CLASSES:
                row = welch_block(df, partner, redox, cls, rng=rng)
                rows.append(row)

    out = pd.DataFrame(rows)
    out["bh_q"] = bh_fdr(out["welch_p"].to_numpy())
    # significant if BH-q < 0.05 AND bootstrap CI excludes zero
    out["sig_q05"] = (out["bh_q"] < 0.05) & ((out["ci_lo_log"] > 0) | (out["ci_hi_log"] < 0))
    out.to_csv(outdir / "class_welch_host_comparison.csv", index=False)

    # markdown summary
    lines = ["# Class-stratified Welch test on log τ_on (CP12 ↔ partner, MED4 vs HM2)",
             "",
             "For each (partner × redox × chemistry class), Welch t-test on the log10(mean_on_ps) "
             "distribution comparing MED4 vs HM2 (pooled across CP12 chains × replicates × pair "
             "instances). Bootstrap 95% CI on the median host difference is the conservative "
             "anchor of significance; BH-FDR corrects across the full 20-test family.",
             "",
             "Interpretation guide:",
             "- `log_diff_hm2_minus_med4 > 0` means HM2 has longer dwell than MED4 in that cell.",
             "- `log_diff_ps_factor = 10^Δlog` gives the multiplicative dwell shift "
             "(e.g., 1.30 ⇒ HM2 contacts last 30% longer; 0.77 ⇒ MED4 contacts last 30% longer).",
             "- A test is flagged `sig_q05 = True` iff BH-q < 0.05 AND the bootstrap CI on "
             "Δlog-dwell excludes zero.",
             "",
             "## Significant cells (BH-FDR q < 0.05 AND CI excludes 0)",
             ""]
    sig = out[out["sig_q05"]].sort_values(["partner", "bh_q"]).copy()
    if sig.empty:
        lines.append("- none —")
    else:
        lines.append("| partner | redox | class | n_med4 | n_hm2 | "
                     "median MED4 (ps) | median HM2 (ps) | Δlog (HM2−MED4) | "
                     "95% CI | dwell factor (HM2/MED4) | d | BH q |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for _, r in sig.iterrows():
            lines.append(
                f"| {r['partner']} | {r['redox']} | {r['contact_class']} | "
                f"{r['n_med4']} | {r['n_hm2']} | "
                f"{r['med4_median_ps']:.0f} | {r['hm2_median_ps']:.0f} | "
                f"{r['log_diff_hm2_minus_med4']:+.3f} | "
                f"[{r['ci_lo_log']:+.3f}, {r['ci_hi_log']:+.3f}] | "
                f"{r['log_diff_ps_factor']:.2f} | {r['cohens_d']:+.2f} | "
                f"{r['bh_q']:.3g} |"
            )
        lines.append("")

    lines += ["## All cells", ""]
    lines.append("| partner | redox | class | n_med4 | n_hm2 | "
                 "median MED4 (ps) | median HM2 (ps) | Δlog (HM2−MED4) | "
                 "95% CI | dwell factor | d | Welch p | BH q | sig |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in out.iterrows():
        sig_mark = "**Y**" if r["sig_q05"] else "n"
        if not np.isfinite(r["welch_p"]):
            lines.append(
                f"| {r['partner']} | {r['redox']} | {r['contact_class']} | "
                f"{r['n_med4']} | {r['n_hm2']} | — | — | — | — | — | — | — | — | n |"
            )
            continue
        lines.append(
            f"| {r['partner']} | {r['redox']} | {r['contact_class']} | "
            f"{r['n_med4']} | {r['n_hm2']} | "
            f"{r['med4_median_ps']:.0f} | {r['hm2_median_ps']:.0f} | "
            f"{r['log_diff_hm2_minus_med4']:+.3f} | "
            f"[{r['ci_lo_log']:+.3f}, {r['ci_hi_log']:+.3f}] | "
            f"{r['log_diff_ps_factor']:.2f} | {r['cohens_d']:+.2f} | "
            f"{r['welch_p']:.3g} | {r['bh_q']:.3g} | {sig_mark} |"
        )
    (outdir / "class_welch_summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {outdir / 'class_welch_summary.md'}")


if __name__ == "__main__":
    main()
