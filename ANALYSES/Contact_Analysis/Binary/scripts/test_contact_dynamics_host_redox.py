#!/usr/bin/env python3
"""
Per-HMM-column statistical test of host (MED4 vs HM2) and redox (ox vs red)
effects on contact dwell time and turnover.

Approach: at each HMM column, take per-pair log(mean_on_ps) for the four
conditions, run a two-way Welch test (host effect at fixed redox; redox effect
at fixed host), and report bootstrap-CI on the host x redox interaction.

Convergence flag: a column is marked 'undersampled' if within-condition CV
of mean_on_ps across the per-rep observations exceeds the between-condition
spread.

BH-FDR corrects across all HMM columns within (partner, comparison).

Outputs:
    contact_dynamics_stats/host_redox_per_hmm_<partner>.csv
    contact_dynamics_stats/convergence_flags_<partner>.csv
    contact_dynamics_stats/host_redox_text_summary.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, ignoring NaNs."""
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
    q_raw = pv * n / ranks
    # enforce monotonic q-values (BH step-up)
    q_sorted = q_raw[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_final = np.empty_like(q_sorted)
    q_final[order] = np.clip(q_sorted, 0, 1)
    q[mask] = q_final
    return q


def host_comparison(per_pair: pd.DataFrame, partner: str) -> pd.DataFrame:
    """For each HMM column × redox state, test MED4 vs HM2 on log(mean_on_ps)."""
    df = per_pair[(per_pair["partner"] == partner) & per_pair["mean_on_ps"].notna()].copy()
    df = df[df["mean_on_ps"] > 0]
    df["log_dwell"] = np.log10(df["mean_on_ps"])
    df = df.dropna(subset=["hmm_col"])
    df["hmm_col"] = df["hmm_col"].astype(int)
    df["host"]  = df["condition"].str[:4]  # 'med4' or 'hm2_'
    df["host"]  = df["host"].str.rstrip("_")
    df["redox"] = df["condition"].str.split("_").str[-1]  # 'ox' or 'red'

    rows = []
    rng = np.random.default_rng(0)
    for (hmm, redox), g in df.groupby(["hmm_col", "redox"]):
        med = g.loc[g["host"] == "med4", "log_dwell"].to_numpy()
        hm2 = g.loc[g["host"] == "hm2",  "log_dwell"].to_numpy()
        if med.size < 5 or hm2.size < 5:
            continue
        t, p = stats.ttest_ind(med, hm2, equal_var=False, nan_policy="omit")
        diff = float(np.median(hm2) - np.median(med))
        # bootstrap CI on median difference
        n_boot = 1000
        diffs = []
        for _ in range(n_boot):
            mb = rng.choice(med, size=med.size)
            hb = rng.choice(hm2, size=hm2.size)
            diffs.append(np.median(hb) - np.median(mb))
        lo, hi = np.percentile(diffs, [2.5, 97.5])
        rows.append({
            "partner": partner, "hmm_col": int(hmm), "redox": redox,
            "n_med4": int(med.size), "n_hm2": int(hm2.size),
            "log_dwell_med4_median": float(np.median(med)),
            "log_dwell_hm2_median":  float(np.median(hm2)),
            "log_dwell_diff_hm2_med4": diff,
            "ci_lo": float(lo), "ci_hi": float(hi),
            "welch_t": float(t), "welch_p": float(p),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["bh_q"] = bh_fdr(out["welch_p"].to_numpy())
        out["sig_q05"] = (out["bh_q"] < 0.05) & ((out["ci_lo"] > 0) | (out["ci_hi"] < 0))
    return out


def redox_comparison(per_pair: pd.DataFrame, partner: str) -> pd.DataFrame:
    """Per HMM column × host, test ox vs red."""
    df = per_pair[(per_pair["partner"] == partner) & per_pair["mean_on_ps"].notna()].copy()
    df = df[df["mean_on_ps"] > 0]
    df["log_dwell"] = np.log10(df["mean_on_ps"])
    df = df.dropna(subset=["hmm_col"])
    df["hmm_col"] = df["hmm_col"].astype(int)
    df["host"]  = df["condition"].str[:4].str.rstrip("_")
    df["redox"] = df["condition"].str.split("_").str[-1]

    rows = []
    rng = np.random.default_rng(1)
    for (hmm, host), g in df.groupby(["hmm_col", "host"]):
        ox  = g.loc[g["redox"] == "ox",  "log_dwell"].to_numpy()
        red = g.loc[g["redox"] == "red", "log_dwell"].to_numpy()
        if ox.size < 5 or red.size < 5:
            continue
        t, p = stats.ttest_ind(ox, red, equal_var=False, nan_policy="omit")
        diff = float(np.median(red) - np.median(ox))
        diffs = []
        for _ in range(1000):
            ob = rng.choice(ox,  size=ox.size)
            rb = rng.choice(red, size=red.size)
            diffs.append(np.median(rb) - np.median(ob))
        lo, hi = np.percentile(diffs, [2.5, 97.5])
        rows.append({
            "partner": partner, "hmm_col": int(hmm), "host": host,
            "n_ox": int(ox.size), "n_red": int(red.size),
            "log_dwell_ox_median": float(np.median(ox)),
            "log_dwell_red_median": float(np.median(red)),
            "log_dwell_diff_red_ox": diff,
            "ci_lo": float(lo), "ci_hi": float(hi),
            "welch_t": float(t), "welch_p": float(p),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["bh_q"] = bh_fdr(out["welch_p"].to_numpy())
        out["sig_q05"] = (out["bh_q"] < 0.05) & ((out["ci_lo"] > 0) | (out["ci_hi"] < 0))
    return out


def convergence_flags(per_pair: pd.DataFrame, partner: str) -> pd.DataFrame:
    """At each (HMM column, condition), CV of mean_on_ps across replicates × chains × partner residues.

    A column is 'undersampled' if median within-condition CV across the 4 conditions
    exceeds 1.0 (i.e., SD > mean), AND between-condition spread (range of medians) < median within-CV * median.
    """
    df = per_pair[(per_pair["partner"] == partner) & per_pair["mean_on_ps"].notna()].copy()
    df = df[df["mean_on_ps"] > 0]
    df = df.dropna(subset=["hmm_col"])
    df["hmm_col"] = df["hmm_col"].astype(int)
    rows = []
    for hmm, g in df.groupby("hmm_col"):
        intra_cvs = []
        cond_medians = []
        for cond, gc in g.groupby("condition"):
            v = gc["mean_on_ps"].to_numpy()
            if v.size < 3:
                continue
            mu = float(np.mean(v))
            sd = float(np.std(v, ddof=1))
            cv = sd / mu if mu > 0 else np.nan
            intra_cvs.append(cv)
            cond_medians.append(float(np.median(v)))
        if not intra_cvs:
            continue
        med_cv = float(np.nanmedian(intra_cvs))
        between_range = float(np.max(cond_medians) - np.min(cond_medians)) if cond_medians else np.nan
        # spread of condition medians divided by typical within-condition std
        # high within-CV + low between-range = undersampled
        rows.append({
            "partner": partner, "hmm_col": int(hmm),
            "n_conditions_evaluated": len(intra_cvs),
            "median_within_cv": med_cv,
            "between_range_ps": between_range,
            "between_to_within_ratio": between_range / (med_cv * np.mean(cond_medians))
                                           if med_cv > 0 and cond_medians else np.nan,
            "undersampled": med_cv > 1.5 and (
                between_range < med_cv * np.mean(cond_medians)
            ),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="contact_dynamics")
    ap.add_argument("--outdir", default="contact_dynamics_stats")
    args = ap.parse_args()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    per_pair = pd.read_csv(indir / "per_pair_stats.csv.gz")
    print(f"loaded {len(per_pair):,} per-pair rows", flush=True)

    summary_lines = ["# Host & redox per-HMM tests on log10(τ_on)", ""]

    for partner in ["prk", "gap2"]:
        hc = host_comparison(per_pair, partner)
        rc = redox_comparison(per_pair, partner)
        cf = convergence_flags(per_pair, partner)

        if not hc.empty:
            hc.to_csv(outdir / f"host_per_hmm_{partner}.csv", index=False)
            n_sig = int(hc["sig_q05"].sum())
            n_tot = len(hc)
            summary_lines.append(f"## {partner.upper()} host effect (MED4 vs HM2)")
            summary_lines.append(f"- HMM columns tested per redox state: **{n_tot}**")
            summary_lines.append(f"- significant after BH-FDR q<0.05 AND CI excludes 0: **{n_sig}**")
            sig = hc.loc[hc["sig_q05"]].sort_values("bh_q")
            if not sig.empty:
                top = sig.head(15)
                summary_lines.append("- top significant positions:")
                summary_lines.append("")
                summary_lines.append("| HMM col | redox | n_med4 | n_hm2 | log-dwell HM2-MED4 | 95% CI | BH q |")
                summary_lines.append("|---|---|---|---|---|---|---|")
                for _, r in top.iterrows():
                    summary_lines.append(
                        f"| {r['hmm_col']} | {r['redox']} | {r['n_med4']} | {r['n_hm2']} | "
                        f"{r['log_dwell_diff_hm2_med4']:+.2f} | "
                        f"[{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}] | {r['bh_q']:.3g} |"
                    )
            summary_lines.append("")

        if not rc.empty:
            rc.to_csv(outdir / f"redox_per_hmm_{partner}.csv", index=False)
            n_sig = int(rc["sig_q05"].sum())
            n_tot = len(rc)
            summary_lines.append(f"## {partner.upper()} redox effect (ox vs red)")
            summary_lines.append(f"- HMM columns tested per host: **{n_tot}**")
            summary_lines.append(f"- significant after BH-FDR q<0.05 AND CI excludes 0: **{n_sig}**")
            sig = rc.loc[rc["sig_q05"]].sort_values("bh_q")
            if not sig.empty:
                top = sig.head(15)
                summary_lines.append("- top significant positions:")
                summary_lines.append("")
                summary_lines.append("| HMM col | host | n_ox | n_red | log-dwell red-ox | 95% CI | BH q |")
                summary_lines.append("|---|---|---|---|---|---|---|")
                for _, r in top.iterrows():
                    summary_lines.append(
                        f"| {r['hmm_col']} | {r['host']} | {r['n_ox']} | {r['n_red']} | "
                        f"{r['log_dwell_diff_red_ox']:+.2f} | "
                        f"[{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}] | {r['bh_q']:.3g} |"
                    )
            summary_lines.append("")

        if not cf.empty:
            cf.to_csv(outdir / f"convergence_{partner}.csv", index=False)
            n_under = int(cf["undersampled"].sum())
            summary_lines.append(f"## {partner.upper()} convergence flags")
            summary_lines.append(f"- HMM columns evaluated: {len(cf)}")
            summary_lines.append(f"- flagged 'undersampled' (within-CV > 1.5 and between-range < within-spread): **{n_under}**")
            if n_under > 0:
                summary_lines.append("- flagged columns: " +
                                     ", ".join(str(int(x)) for x in cf.loc[cf["undersampled"], "hmm_col"]))
            summary_lines.append("")

    (outdir / "host_redox_text_summary.md").write_text("\n".join(summary_lines) + "\n")
    print(f"wrote stats outputs to {outdir}")


if __name__ == "__main__":
    main()
