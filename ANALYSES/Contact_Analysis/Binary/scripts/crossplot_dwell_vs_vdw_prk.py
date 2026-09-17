#!/usr/bin/env python3
"""
Per-residue VDW cross-plot for CP12 ↔ PRK using direct trajectory-computed VDW
(via `compute_per_residue_vdw_direct.py`). Mirrors `crossplot_dwell_vs_vdw_gap2.py`.

Inputs:
    contact_kinetic_thermo_bridge/per_residue_vdw_direct_prk.csv
    contact_dynamics/per_pair_hmm_aggregated.csv
    CP12_MED4_vs_HM2_alignment_annotation.csv

Outputs:
    contact_kinetic_thermo_bridge/per_hmm_vdw_prk_by_condition.csv
    contact_kinetic_thermo_bridge/host_diff_vdw_prk_{ox,red}.csv
    contact_kinetic_thermo_bridge/diff_mean_on_ps_vs_vdw_mean_prk_{ox,red}.{pdf,png}
    contact_kinetic_thermo_bridge/diff_k_off_per_ns_vs_vdw_mean_prk_{ox,red}.{pdf,png}
    contact_kinetic_thermo_bridge/vdw_bridge_prk_text_summary.md

The direct VDW data also includes GAP2 for validation against the existing
gmx_MMPBSA idecomp=2 numbers — that cross-check is performed and reported.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def load_alignment(aln_csv: Path) -> tuple[dict[int, int], dict[int, int], pd.DataFrame]:
    df = pd.read_csv(aln_csv)
    med4_map = {int(p): int(c) for p, c in zip(df["med4_pos"], df["aln_col_1based"])
                if pd.notna(p)}
    hm2_map  = {int(p): int(c) for p, c in zip(df["hm2_pos"], df["aln_col_1based"])
                if pd.notna(p)}
    return med4_map, hm2_map, df


def aggregate_direct_vdw(direct_csv: Path, med4_map: dict[int, int],
                          hm2_map: dict[int, int]) -> pd.DataFrame:
    """Average direct-VDW across 4 CP12 chains × 3 replicates per condition, by HMM col."""
    df = pd.read_csv(direct_csv)
    df["hmm_col"] = np.nan
    med_mask = df["condition"].str.startswith("med4")
    hm_mask  = df["condition"].str.startswith("hm2")
    df.loc[med_mask, "hmm_col"] = df.loc[med_mask, "resid"].map(med4_map)
    df.loc[hm_mask,  "hmm_col"] = df.loc[hm_mask,  "resid"].map(hm2_map)
    df = df.dropna(subset=["hmm_col"]).copy()
    df["hmm_col"] = df["hmm_col"].astype(int)
    # mean across chain × rep within each (condition, hmm_col)
    out = (df.groupby(["condition", "hmm_col"], as_index=False)
             .agg(vdw_mean=("vdw_mean_kcal_mol", "mean"),
                  vdw_sem =("vdw_mean_kcal_mol",
                             lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0),
                  n_obs=("vdw_mean_kcal_mol", "count")))
    return out


def load_dwell(per_pair_hmm_path: Path, partner: str) -> pd.DataFrame:
    df = pd.read_csv(per_pair_hmm_path)
    df = df[df["partner"] == partner].copy()
    if df.empty:
        raise SystemExit(f"no dwell data for partner={partner}")
    out = (df.groupby(["condition", "hmm_col"], as_index=False)
             .agg(mean_on_ps=("mean_on_ps_mean", "mean"),
                  median_on_ps=("median_on_ps_mean", "mean"),
                  k_off_per_ns=("k_off_per_ns_mean", "mean"),
                  occupancy=("occupancy_mean", "mean")))
    out["hmm_col"] = out["hmm_col"].astype(int)
    return out


def host_diff(merged: pd.DataFrame, value_col: str, redox: str) -> pd.DataFrame:
    sel = merged[merged["condition"].str.endswith("_" + redox)].copy()
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


def validate_gap2(direct_csv: Path, mm_idecomp_csv: Path,
                   med4_map: dict[int, int], hm2_map: dict[int, int]) -> tuple[float, int]:
    """Pearson r between direct-VDW and gmx_MMPBSA idecomp=2 VDW for GAP2, per (condition, HMM col).

    The idecomp=2 source is `contact_kinetic_thermo_bridge/per_residue_vdw_gap2.csv`
    (output of `crossplot_dwell_vs_vdw_gap2.py`).
    """
    direct = aggregate_direct_vdw(direct_csv, med4_map, hm2_map)
    mm = pd.read_csv(mm_idecomp_csv)
    # mm reference (per_residue_vdw_gap2.csv from crossplot_dwell_vs_vdw_gap2.py)
    # has columns: side, resname, resid, internal_mean, vdw_mean, ..., condition, rep
    # Filter to CP12 ligand side only.
    if "side" in mm.columns:
        mm = mm[mm["side"] == "L"].copy()
    rep_col = "rep" if "rep" in mm.columns else "replicate"
    per_rep = (mm.groupby(["condition", rep_col, "resid"], as_index=False)
                 .agg(vdw_mean=("vdw_mean", "mean")))

    per_rep_agg = per_rep.copy()
    per_rep_agg["hmm_col"] = np.nan
    med_mask = per_rep_agg["condition"].str.startswith("med4")
    hm_mask  = per_rep_agg["condition"].str.startswith("hm2")
    per_rep_agg.loc[med_mask, "hmm_col"] = per_rep_agg.loc[med_mask, "resid"].map(med4_map)
    per_rep_agg.loc[hm_mask,  "hmm_col"] = per_rep_agg.loc[hm_mask,  "resid"].map(hm2_map)
    per_rep_agg = per_rep_agg.dropna(subset=["hmm_col"]).copy()
    per_rep_agg["hmm_col"] = per_rep_agg["hmm_col"].astype(int)
    mm_by_hmm = (per_rep_agg.groupby(["condition", "hmm_col"], as_index=False)
                              .agg(vdw_mm=("vdw_mean", "mean")))
    merged = direct.merge(mm_by_hmm, on=["condition", "hmm_col"], how="inner")
    x = merged["vdw_mean"].to_numpy()
    y = merged["vdw_mm"].to_numpy()
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return float("nan"), 0
    r, _ = stats.pearsonr(x[mask], y[mask])
    return r, int(mask.sum())


def scatter_plot(merged: pd.DataFrame, energy_col: str, dwell_col: str,
                  redox: str, aln_df: pd.DataFrame, outdir: Path,
                  partner: str = "prk") -> None:
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

    fig, ax = plt.subplots(figsize=(5.8, 5.3))
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
        title = (f"CP12 vs {partner.upper()}, {redox.upper()}: "
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
    fname = f"diff_{dwell_col}_vs_{energy_col}_{partner}_{redox}"
    fig.savefig(outdir / f"{fname}.pdf")
    fig.savefig(outdir / f"{fname}.png", dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--direct_csv",
                    default="contact_kinetic_thermo_bridge/per_residue_vdw_direct_prk.csv")
    ap.add_argument("--per_hmm",
                    default="contact_dynamics/per_pair_hmm_aggregated.csv")
    ap.add_argument("--alignment_csv",
                    default="CP12_MED4_vs_HM2_alignment_annotation.csv")
    ap.add_argument("--validate_gap2_csv",
                    default="contact_kinetic_thermo_bridge/per_replicate_pair_energetics.csv",
                    help="Optional: gmx_MMPBSA per-replicate pair energetics CSV from "
                         "analyze_cp12_prk_idecomp4_bridge.py for a cross-validation check.")
    ap.add_argument("--direct_gap2_csv",
                    default="contact_kinetic_thermo_bridge/per_residue_vdw_direct_gap2.csv",
                    help="Direct-VDW GAP2 CSV for validation against mmgbsa idecomp=2.")
    ap.add_argument("--mm_gap2_csv",
                    default="contact_kinetic_thermo_bridge/per_residue_vdw_gap2.csv",
                    help="gmx_MMPBSA idecomp=2 per-residue VDW for GAP2 (validation reference).")
    ap.add_argument("--outdir", default="contact_kinetic_thermo_bridge")
    ap.add_argument("--partner", default="prk")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    med4_map, hm2_map, aln_df = load_alignment(Path(args.alignment_csv).resolve())

    # Validate direct-VDW against gmx_MMPBSA idecomp=2 numbers on the GAP2 partner
    direct_gap2 = Path(args.direct_gap2_csv).resolve()
    mm_gap2 = Path(args.mm_gap2_csv).resolve()
    val_r, val_n = float("nan"), 0
    if direct_gap2.exists() and mm_gap2.exists():
        val_r, val_n = validate_gap2(direct_gap2, mm_gap2, med4_map, hm2_map)
        print(f"[validate] direct VDW vs gmx_MMPBSA idecomp=2 (GAP2): "
              f"Pearson r = {val_r:+.3f}, n = {val_n}")

    direct_prk_path = Path(args.direct_csv).resolve()
    if not direct_prk_path.exists():
        raise SystemExit(f"missing input: {direct_prk_path}")
    direct_prk = aggregate_direct_vdw(direct_prk_path, med4_map, hm2_map)
    direct_prk.to_csv(outdir / f"per_hmm_vdw_{args.partner}_by_condition.csv", index=False)

    dwell = load_dwell(Path(args.per_hmm).resolve(), args.partner)
    merged = direct_prk.merge(dwell, on=["condition", "hmm_col"], how="inner")
    merged.to_csv(outdir / f"merged_per_hmm_{args.partner}.csv", index=False)
    print(f"merged rows: {len(merged)}")

    text_lines = [
        f"# Δdwell vs Δ(per-residue VDW) — CP12 vs {args.partner.upper()}, MED4 − HM2",
        "",
        f"VDW source: direct trajectory L-J computation (cutoff 10 Å, stride 4, frames 200-end), "
        f"using amber99sb σ/ε from `mmgbsa_cp12_prk/COM.prmtop`. This is the equivalent of "
        f"gmx_MMPBSA's idecomp=2 DELTAS VDW computed without re-running MM/PBSA.",
        "",
    ]
    if val_n > 0:
        text_lines.append(
            f"**Validation:** the same direct method applied to CP12 ↔ GAP2 has "
            f"Pearson r = {val_r:+.3f} (n = {val_n} HMM-col × condition cells) against the "
            f"existing gmx_MMPBSA idecomp=2 numbers from `mmgbsa_cp12_gap2_nad/`."
        )
        text_lines.append("")

    for redox in ["ox", "red"]:
        rows = [host_diff(merged, "vdw_mean", redox),
                host_diff(merged, "mean_on_ps", redox),
                host_diff(merged, "k_off_per_ns", redox),
                host_diff(merged, "occupancy", redox)]
        if rows[0].empty:
            continue
        block = rows[0]
        for r in rows[1:]:
            if not r.empty:
                block = block.merge(r, on="hmm_col", how="outer")
        block["redox"] = redox
        block.to_csv(outdir / f"host_diff_vdw_{args.partner}_{redox}.csv", index=False)

        text_lines.append(f"## {redox.upper()}")
        sub = aln_df[["aln_col_1based", "med4_aa", "hm2_aa"]].rename(columns={"aln_col_1based": "hmm_col"})
        block_ann = block.merge(sub, on="hmm_col", how="left")
        block_ann["is_substituted"] = (block_ann["med4_aa"] != block_ann["hm2_aa"]) & block_ann["med4_aa"].ne("-") & block_ann["hm2_aa"].ne("-")

        text_lines.append("| Δ component | n_all | r_all | p_all | n_subs | r_subs | p_subs |")
        text_lines.append("|---|---|---|---|---|---|---|")
        for ec in ["vdw_mean"]:
            xc = f"diff_{ec}_med4_minus_hm2"
            yc = "diff_mean_on_ps_med4_minus_hm2"
            x = block_ann[xc].to_numpy()
            y = block_ann[yc].to_numpy()
            m = np.isfinite(x) & np.isfinite(y)
            r_all, p_all = stats.pearsonr(x[m], y[m]) if m.sum() >= 5 else (np.nan, np.nan)
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
        scatter_plot(block_ann, "vdw_mean", "mean_on_ps", redox, aln_df, outdir,
                      partner=args.partner)
        scatter_plot(block_ann, "vdw_mean", "k_off_per_ns", redox, aln_df, outdir,
                      partner=args.partner)

    (outdir / f"vdw_bridge_{args.partner}_text_summary.md").write_text("\n".join(text_lines) + "\n")
    print(f"wrote {outdir / ('vdw_bridge_' + args.partner + '_text_summary.md')}")


if __name__ == "__main__":
    main()
