#!/usr/bin/env python3
"""
Per-residue VDW cross-plot for CP12 against GAP2 (with NAD context).

Uses the canonical idecomp=2 per-residue decomposition files from
`mmgbsa_cp12_gap2_nad/FINAL_DECOMP_MMPBSA.dat` — full 12-replicate coverage.

For each (condition, rep) parses the DELTAS section, extracts per-CP12-residue
mean VDW (chain-averaged by gmx_MMPBSA across CP12 chains C/H/M/P), maps to
HMM column, then computes the host-difference (MED4 − HM2) per HMM col and
correlates it with the host-difference in contact dwell time from the binary
contact dynamics analysis.

Output:
    contact_kinetic_thermo_bridge/per_residue_vdw_gap2.csv
    contact_kinetic_thermo_bridge/diff_mean_on_ps_vs_vdw_gap2_{ox,red}.{pdf,png}
    contact_kinetic_thermo_bridge/diff_k_off_per_ns_vs_vdw_gap2_{ox,red}.{pdf,png}
    contact_kinetic_thermo_bridge/vdw_bridge_text_summary.md
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


CONDITIONS = {
    "MED4_NAD/0003": "med4_ox",
    "MED4_NAD/0000": "med4_red",
    "HM2_NAD/0003":  "hm2_ox",
    "HM2_NAD/0000":  "hm2_red",
}


DAT_PATH_TEMPLATE = "{root}/{system}/{sysid}/{rep}/mmgbsa_cp12_gap2_nad/FINAL_DECOMP_MMPBSA.dat"


def parse_idecomp2_deltas(path: Path) -> pd.DataFrame:
    """Return per-CP12-residue (L::) DELTAS rows with mean/sd/sem for each component.

    Output columns:
        side ('L' or 'R'), resname, resid,
        internal_mean, vdw_mean, eel_mean, polar_mean, nonpolar_mean, total_mean,
        vdw_sd, eel_sd, total_sd  (selected sds for downstream uncertainty)
    """
    rows = []
    in_deltas = False
    in_total = False
    row_re = re.compile(r"^([LR])::([A-Z]{3}):(\d+),(.*)$")
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s == "DELTAS:":
                in_deltas = True
                in_total = False
                continue
            if in_deltas and s == "Total Energy Decomposition:":
                in_total = True
                continue
            if in_deltas and s in ("Sidechain Energy Decomposition:",
                                    "Backbone Energy Decomposition:"):
                in_total = False
                continue
            if not (in_deltas and in_total):
                continue
            m = row_re.match(line)
            if not m:
                continue
            side, resname, resid_str, rest = m.groups()
            parts = rest.split(",")
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            if len(vals) < 18:
                continue
            # 6 columns × (avg, sd, sem). Order: Internal, VDW, EEL, Polar, NonPolar, TOTAL.
            internal_mean = vals[0]
            vdw_mean, vdw_sd = vals[3], vals[4]
            eel_mean, eel_sd = vals[6], vals[7]
            polar_mean = vals[9]
            nonpolar_mean = vals[12]
            total_mean, total_sd = vals[15], vals[16]
            rows.append({
                "side": side, "resname": resname, "resid": int(resid_str),
                "internal_mean": internal_mean,
                "vdw_mean": vdw_mean, "vdw_sd": vdw_sd,
                "eel_mean": eel_mean, "eel_sd": eel_sd,
                "polar_mean": polar_mean,
                "nonpolar_mean": nonpolar_mean,
                "total_mean": total_mean, "total_sd": total_sd,
            })
    return pd.DataFrame(rows)


def gather(root: Path, reps: list[str]) -> pd.DataFrame:
    """Load all 12 (condition, rep) decomps, return per-rep CP12-side per-residue rows."""
    out = []
    for sys_path, cond in CONDITIONS.items():
        system, sysid = sys_path.split("/")
        for rep in reps:
            path = Path(DAT_PATH_TEMPLATE.format(root=root, system=system,
                                                  sysid=sysid, rep=rep))
            if not path.exists():
                print(f"[miss] {cond}/{rep}: {path}")
                continue
            df = parse_idecomp2_deltas(path)
            cp12 = df[df["side"] == "L"].copy()
            cp12["condition"] = cond
            cp12["rep"] = rep
            out.append(cp12)
            print(f"[ok ] {cond}/{rep}: {len(cp12)} CP12 residues parsed")
    if not out:
        raise SystemExit("No decomp files found")
    return pd.concat(out, ignore_index=True)


def load_alignment(aln_csv: Path) -> tuple[dict[int, int], dict[int, int], pd.DataFrame]:
    df = pd.read_csv(aln_csv)
    med4_map = {int(p): int(c) for p, c in zip(df["med4_pos"], df["aln_col_1based"])
                if pd.notna(p)}
    hm2_map  = {int(p): int(c) for p, c in zip(df["hm2_pos"], df["aln_col_1based"])
                if pd.notna(p)}
    return med4_map, hm2_map, df


def merge_with_alignment(per_res: pd.DataFrame,
                          med4_map: dict[int, int],
                          hm2_map: dict[int, int]) -> pd.DataFrame:
    out = per_res.copy()
    out["hmm_col"] = np.nan
    for cond in out["condition"].unique():
        mask = out["condition"] == cond
        host_map = med4_map if cond.startswith("med4") else hm2_map
        out.loc[mask, "hmm_col"] = out.loc[mask, "resid"].map(host_map)
    return out


def aggregate_by_hmm(per_res: pd.DataFrame) -> pd.DataFrame:
    df = per_res.dropna(subset=["hmm_col"]).copy()
    df["hmm_col"] = df["hmm_col"].astype(int)
    # mean per (condition, hmm_col) across replicates (and across CP12 chains - already
    # collapsed by gmx_MMPBSA)
    keys = ["condition", "hmm_col"]
    out = (df.groupby(keys, as_index=False)
             .agg(vdw_mean=("vdw_mean", "mean"),
                  vdw_sem=("vdw_mean", lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0),
                  eel_mean=("eel_mean", "mean"),
                  total_mean=("total_mean", "mean"),
                  n_reps=("rep", "nunique")))
    return out


def host_diff(by_hmm: pd.DataFrame, value_col: str, redox: str) -> pd.DataFrame:
    sel = by_hmm[by_hmm["condition"].str.endswith("_" + redox)].copy()
    sel["host"] = sel["condition"].str[:4].str.rstrip("_")
    piv = sel.pivot(index="hmm_col", columns="host", values=value_col)
    if "med4" not in piv.columns or "hm2" not in piv.columns:
        return pd.DataFrame()
    piv["diff_med4_minus_hm2"] = piv["med4"] - piv["hm2"]
    piv = piv.reset_index().rename(columns={
        "med4": f"{value_col}_med4",
        "hm2":  f"{value_col}_hm2",
        "diff_med4_minus_hm2": f"diff_{value_col}_med4_minus_hm2",
    })
    return piv


def load_dwell_gap2(per_pair_hmm_path: Path) -> pd.DataFrame:
    df = pd.read_csv(per_pair_hmm_path)
    df = df[df["partner"] == "gap2"].copy()
    if df.empty:
        raise SystemExit("No GAP2 dwell data in per_pair_hmm_aggregated.csv")
    out = (df.groupby(["condition", "hmm_col"], as_index=False)
             .agg(mean_on_ps=("mean_on_ps_mean", "mean"),
                  median_on_ps=("median_on_ps_mean", "mean"),
                  k_off_per_ns=("k_off_per_ns_mean", "mean"),
                  occupancy=("occupancy_mean", "mean")))
    out["hmm_col"] = out["hmm_col"].astype(int)
    return out


def scatter_plot(merged: pd.DataFrame, energy_col: str, dwell_col: str,
                 redox: str, aln_df: pd.DataFrame, outdir: Path) -> None:
    df = merged.dropna(subset=[f"diff_{energy_col}_med4_minus_hm2",
                                f"diff_{dwell_col}_med4_minus_hm2"]).copy()
    if df.empty:
        return
    if "med4_aa" not in df.columns:
        sub = aln_df[["aln_col_1based", "med4_aa", "hm2_aa"]].rename(columns={"aln_col_1based": "hmm_col"})
        df = df.merge(sub, on="hmm_col", how="left")
    if "is_substituted" not in df.columns:
        df["is_substituted"] = (df["med4_aa"] != df["hm2_aa"]) & df["med4_aa"].ne("-") & df["hm2_aa"].ne("-")

    x = df[f"diff_{energy_col}_med4_minus_hm2"].to_numpy()
    y = df[f"diff_{dwell_col}_med4_minus_hm2"].to_numpy()
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.axhline(0, color="0.6", lw=0.5)
    ax.axvline(0, color="0.6", lw=0.5)
    conserved = df[~df["is_substituted"] & mask]
    substituted = df[df["is_substituted"] & mask]
    ax.scatter(conserved[f"diff_{energy_col}_med4_minus_hm2"],
               conserved[f"diff_{dwell_col}_med4_minus_hm2"],
               s=30, c="#bbbbbb", edgecolor="white", lw=0.5, alpha=0.75,
               label=f"conserved (n={len(conserved)})")
    ax.scatter(substituted[f"diff_{energy_col}_med4_minus_hm2"],
               substituted[f"diff_{dwell_col}_med4_minus_hm2"],
               s=44, c="#d62728", edgecolor="white", lw=0.5, alpha=0.9,
               label=f"substituted (n={len(substituted)})")
    for _, r in substituted.iterrows():
        ax.annotate(f"{int(r['hmm_col'])} ({r['med4_aa']}>{r['hm2_aa']})",
                    (r[f"diff_{energy_col}_med4_minus_hm2"],
                     r[f"diff_{dwell_col}_med4_minus_hm2"]),
                    fontsize=6, alpha=0.7, xytext=(2, 2), textcoords="offset points")

    rho_all, p_all = stats.pearsonr(x[mask], y[mask])
    sx = substituted[f"diff_{energy_col}_med4_minus_hm2"].to_numpy()
    sy = substituted[f"diff_{dwell_col}_med4_minus_hm2"].to_numpy()
    sm = np.isfinite(sx) & np.isfinite(sy)
    if sm.sum() >= 5:
        rho_sub, p_sub = stats.pearsonr(sx[sm], sy[sm])
        title = (f"CP12 vs GAP2(+NAD), {redox.upper()}: "
                 f"Δ{energy_col} vs Δ{dwell_col} (MED4 − HM2)\n"
                 f"all r={rho_all:+.2f} (p={p_all:.2g}), "
                 f"subs r={rho_sub:+.2f} (p={p_sub:.2g})")
    else:
        title = f"r_all={rho_all:+.2f} (p={p_all:.2g})"
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(f"Δ {energy_col} (kcal/mol)  MED4 − HM2")
    ax.set_ylabel(f"Δ {dwell_col}  MED4 − HM2")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fname = f"diff_{dwell_col}_vs_{energy_col}_gap2_{redox}"
    fig.savefig(outdir / f"{fname}.pdf")
    fig.savefig(outdir / f"{fname}.png", dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/home/scott.widmann/DOE/Perlmutter")
    ap.add_argument("--alignment_csv",
                    default="CP12_MED4_vs_HM2_alignment_annotation.csv")
    ap.add_argument("--per_hmm",
                    default="contact_dynamics/per_pair_hmm_aggregated.csv")
    ap.add_argument("--outdir", default="contact_kinetic_thermo_bridge")
    ap.add_argument("--reps", nargs="+", default=["rep1", "rep2", "rep3"])
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    per_res = gather(Path(args.root).resolve(), args.reps)
    per_res.to_csv(outdir / "per_residue_vdw_gap2.csv", index=False)

    med4_map, hm2_map, aln_df = load_alignment(Path(args.alignment_csv).resolve())
    per_res_hmm = merge_with_alignment(per_res, med4_map, hm2_map)
    by_hmm = aggregate_by_hmm(per_res_hmm)
    by_hmm.to_csv(outdir / "per_hmm_vdw_gap2_by_condition.csv", index=False)

    dwell = load_dwell_gap2(Path(args.per_hmm).resolve())

    merged_by_hmm = by_hmm.merge(dwell, on=["condition", "hmm_col"], how="inner")
    merged_by_hmm.to_csv(outdir / "merged_per_hmm_gap2.csv", index=False)
    print(f"merged rows: {len(merged_by_hmm)}")

    # build host-diff blocks
    text_lines = ["# Δdwell vs Δ(per-residue VDW/EEL/total) — CP12 vs GAP2(+NAD), MED4 − HM2", ""]
    for redox in ["ox", "red"]:
        rows = [host_diff(merged_by_hmm, "vdw_mean", redox),
                host_diff(merged_by_hmm, "eel_mean", redox),
                host_diff(merged_by_hmm, "total_mean", redox),
                host_diff(merged_by_hmm, "mean_on_ps", redox),
                host_diff(merged_by_hmm, "k_off_per_ns", redox),
                host_diff(merged_by_hmm, "occupancy", redox)]
        if not rows[0].empty:
            block = rows[0]
            for r in rows[1:]:
                if not r.empty:
                    block = block.merge(r, on="hmm_col", how="outer")
            block["redox"] = redox
            block.to_csv(outdir / f"host_diff_vdw_gap2_{redox}.csv", index=False)

            text_lines.append(f"## {redox.upper()}")
            text_lines.append("| Δ component | n_all | r_all | p_all | n_subs | r_subs | p_subs |")
            text_lines.append("|---|---|---|---|---|---|---|")
            sub = aln_df[["aln_col_1based", "med4_aa", "hm2_aa"]].rename(columns={"aln_col_1based": "hmm_col"})
            block_ann = block.merge(sub, on="hmm_col", how="left")
            block_ann["is_substituted"] = (block_ann["med4_aa"] != block_ann["hm2_aa"]) & block_ann["med4_aa"].ne("-") & block_ann["hm2_aa"].ne("-")
            for ec in ["vdw_mean", "eel_mean", "total_mean"]:
                xc = f"diff_{ec}_med4_minus_hm2"
                yc = "diff_mean_on_ps_med4_minus_hm2"
                if xc not in block_ann.columns or yc not in block_ann.columns:
                    continue
                x = block_ann[xc].to_numpy()
                y = block_ann[yc].to_numpy()
                m = np.isfinite(x) & np.isfinite(y)
                if m.sum() < 5:
                    text_lines.append(f"| {ec} | {int(m.sum())} | — | — | — | — | — |")
                    continue
                r_all, p_all = stats.pearsonr(x[m], y[m])
                ssub = block_ann[block_ann["is_substituted"] & m]
                sx = ssub[xc].to_numpy(); sy = ssub[yc].to_numpy()
                sm = np.isfinite(sx) & np.isfinite(sy)
                if sm.sum() >= 5:
                    r_sub, p_sub = stats.pearsonr(sx[sm], sy[sm])
                    text_lines.append(f"| {ec} | {int(m.sum())} | {r_all:+.3f} | {p_all:.3g} | "
                                      f"{int(sm.sum())} | {r_sub:+.3f} | {p_sub:.3g} |")
                else:
                    text_lines.append(f"| {ec} | {int(m.sum())} | {r_all:+.3f} | {p_all:.3g} | "
                                      f"{int(sm.sum())} | — | — |")
            text_lines.append("")

            # plot
            for ec in ["vdw_mean", "eel_mean", "total_mean"]:
                scatter_plot(block_ann, ec, "mean_on_ps", redox, aln_df, outdir)
                scatter_plot(block_ann, ec, "k_off_per_ns", redox, aln_df, outdir)

    (outdir / "vdw_bridge_text_summary.md").write_text("\n".join(text_lines) + "\n")
    print(f"wrote {outdir / 'vdw_bridge_text_summary.md'}")


if __name__ == "__main__":
    main()
