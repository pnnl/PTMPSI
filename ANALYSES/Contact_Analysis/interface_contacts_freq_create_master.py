#!/usr/bin/env python3
"""Build the master CSV used by the CP12 contact-frequency heatmap figures.

Pipeline (per condition × partner × residue):

  per-rep per-chain freq   from interface_contacts_freq.py
        ↓ average across CP12 chains (C,H,M,P)
  per-rep chain_avg
        ↓ mean across reps (rep1, rep2, rep3)
  mean_freq_chain_avg, sem_freq_chain_avg

Then for each (host, redox, HMM aln column):
  - join the prk, gap2, receptor chain-rep averages
  - count partner residues per HMM column with freq >= PAIR_THRESH
    (from analyze_cp12_contact_pairs.py output)
  - classify by INTERFACE_THRESH on the chain-rep-averaged prk and gap2 freqs

Inputs:
  --interface_root <dir>   produced by interface_contacts_freq.py
                            expects <dir>/<condition>/<rep>/interface_resfreq_cp12_side_{partner}.csv
                            with --no-reps: <dir>/<condition>/interface_resfreq_cp12_side_{partner}.csv
  --pairs_dir_med4 <dir>   produced by analyze_cp12_contact_pairs.py for MED4
                            expects <dir>/differential_pairs_{prk,gap2}.csv
  --pairs_dir_hm2 <dir>    produced by analyze_cp12_contact_pairs.py for HM2
                            expects <dir>/differential_pairs_{prk,gap2}.csv
  --alignment_csv <file>   CP12 MED4↔HM2 HMM alignment annotation

Output:
  --out <file>             single master CSV used by the heatmap scripts

Run from MALIO/ (or anywhere). Example:

  python interface_contacts_freq_create_master.py \\
    --interface_root ../interface_results/redox_cutoff_5p0A \\
    --pairs_dir ../cp12_contact_pairs_5p0A \\
    --alignment_csv ../CP12_MED4_vs_HM2_alignment_annotation.csv \\
    --out ../FINAL_FIGURES/cp12_contact_freq_heatmap_data_master.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PARTNERS   = ["prk", "gap2", "receptor"]
REPS       = ["rep1", "rep2", "rep3"]
N_HMM      = 70
PAIR_THRESH      = 0.10   # partner-count threshold
INTERFACE_THRESH = 0.50   # interface class threshold


def _host_from_path(cond_path: Path, interface_root: Path) -> str:
    """Infer host by checking all path components between interface_root and cond_path."""
    for part in cond_path.relative_to(interface_root).parts:
        pl = part.lower()
        if "med4" in pl:
            return "med4"
        if "hm2" in pl:
            return "hm2"
    return "hm2"  # fallback


def chain_and_rep_average(interface_root: Path, no_reps: bool = False) -> pd.DataFrame:
    """Average per-rep per-chain freq across chains, then mean+SEM across reps."""
    # Discover condition directories: support both flat (<root>/<cond>/) and
    # one-level-nested (<root>/<host>/<cond>/) layouts by searching for the
    # directories that actually contain the expected CSV files.
    partner0 = PARTNERS[0]
    csv_glob = (f"interface_resfreq_cp12_side_{partner0}.csv"
                if no_reps else f"*/interface_resfreq_cp12_side_{partner0}.csv")
    cond_paths = sorted({p.parent for p in interface_root.glob(f"**/{csv_glob}")})
    print(f"  discovered condition paths: {[str(p) for p in cond_paths]}")
    rows = []
    for cond_path in cond_paths:
        cond = cond_path.name
        host = _host_from_path(cond_path, interface_root)
        for partner in PARTNERS:
            if no_reps:
                f = cond_path / f"interface_resfreq_cp12_side_{partner}.csv"
                per_rep = [pd.read_csv(f).assign(rep="rep1")]
            else:
                per_rep = []
                for rep in REPS:
                    f = cond_path / rep / f"interface_resfreq_cp12_side_{partner}.csv"
                    per_rep.append(pd.read_csv(f).assign(rep=rep))
            all_rep = pd.concat(per_rep, ignore_index=True)

            # within each rep: mean freq across the 4 CP12 chains
            chain_avg = (all_rep.groupby(["rep", "resid", "resname"], as_index=False)
                                 .agg(chain_avg=("freq", "mean"),
                                      chain_std=("freq", "std"),
                                      n_chains=("freq", "count")))
            chain_avg["chain_sem_within_rep"] = chain_avg["chain_std"] / np.sqrt(chain_avg["n_chains"])

            # across reps: mean and SEM of the chain-averaged values
            sem3 = lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0
            agg = (chain_avg.groupby(["resid", "resname"], as_index=False)
                            .agg(mean_freq_chain_avg=("chain_avg", "mean"),
                                 sem_freq_chain_avg=("chain_avg", sem3),
                                 mean_sem_within_chain=("chain_sem_within_rep", "mean"),
                                 n_chains_present=("n_chains", "max")))
            agg.insert(0, "condition", cond)
            agg.insert(1, "host", host)
            agg.insert(2, "partner", partner)
            rows.append(agg)
    return pd.concat(rows, ignore_index=True)


def build_master(chain_rep_df: pd.DataFrame, pairs_dirs: dict, aln_csv: Path) -> pd.DataFrame:
    aln = pd.read_csv(aln_csv)
    aln["aln_col_1based"] = aln["aln_col_1based"].astype(int)

    # CP12 residue → HMM column maps, per host
    med4_map = {int(r.med4_pos): int(r.aln_col_1based) for _, r in aln.iterrows() if pd.notna(r.med4_pos)}
    hm2_map  = {int(r.hm2_pos):  int(r.aln_col_1based) for _, r in aln.iterrows() if pd.notna(r.hm2_pos)}
    chain_rep_df = chain_rep_df.copy()
    chain_rep_df["aln_col_1based"] = chain_rep_df.apply(
        lambda r: (med4_map if r["host"] == "med4" else hm2_map).get(int(r["resid"])), axis=1)

    # preload differential-pairs tables keyed by host
    diff = {host: {"prk":  pd.read_csv(d / "differential_pairs_prk.csv"),
                   "gap2": pd.read_csv(d / "differential_pairs_gap2.csv")}
            for host, d in pairs_dirs.items()}

    out_rows = []
    cond_host_pairs = (chain_rep_df[["condition", "host"]]
                       .drop_duplicates()
                       .sort_values(["host", "condition"])
                       .itertuples(index=False))
    for cond, host in cond_host_pairs:
        redox = cond[len(host) + 1:] if cond.startswith(host + "_") else cond
        host_diff = diff[host]
        for aln_col in range(1, N_HMM + 1):
            arow      = aln[aln["aln_col_1based"] == aln_col].iloc[0]
            med4_aa   = arow.get("med4_aa", "-") or "-"
            hm2_aa    = arow.get("hm2_aa",  "-") or "-"
            host_pos  = arow.get("med4_pos") if host == "med4" else arow.get("hm2_pos")
            host_aa   = med4_aa if host == "med4" else hm2_aa
            is_gap    = pd.isna(host_pos)

            prk_f = prk_s = gap_f = gap_s = rec_f = rec_s = np.nan
            if not is_gap:
                sub = chain_rep_df[(chain_rep_df.condition == cond) &
                                    (chain_rep_df.host == host) &
                                    (chain_rep_df.resid == int(host_pos))]
                pick = lambda p: sub[sub.partner == p]
                for p, ftarget, starget in [("prk", "prk_f", "prk_s"),
                                              ("gap2","gap_f","gap_s"),
                                              ("receptor","rec_f","rec_s")]:
                    ps = pick(p)
                    if not ps.empty:
                        vals = (float(ps.mean_freq_chain_avg.iloc[0]),
                                float(ps.sem_freq_chain_avg.iloc[0]))
                        if p == "prk":      prk_f, prk_s = vals
                        elif p == "gap2":   gap_f, gap_s = vals
                        else:               rec_f, rec_s = vals

            if is_gap:
                prk_cnt = gap_cnt = 0
            else:
                fc = f"{cond}_freq"
                prk_cnt = int(((host_diff["prk"]["cp12_aln_pos"] == aln_col)  & (host_diff["prk"][fc]  >= PAIR_THRESH)).sum())
                gap_cnt = int(((host_diff["gap2"]["cp12_aln_pos"] == aln_col) & (host_diff["gap2"][fc] >= PAIR_THRESH)).sum())

            if is_gap:
                cls = "Gap"
            else:
                prk_above = (not np.isnan(prk_f)) and prk_f >= INTERFACE_THRESH
                gap_above = (not np.isnan(gap_f)) and gap_f >= INTERFACE_THRESH
                cls = ("Both" if prk_above and gap_above
                       else "PRK" if prk_above
                       else "GAP2" if gap_above
                       else "Neither")

            out_rows.append({
                "host": host, "redox": redox, "condition": cond,
                "aln_col_1based": aln_col,
                "med4_aa": med4_aa, "hm2_aa": hm2_aa,
                "med4_pos": arow.get("med4_pos"), "hm2_pos": arow.get("hm2_pos"),
                "host_aa": host_aa, "host_resid": host_pos,
                "is_gap_for_host": is_gap,
                "prk_freq_chain_rep_avg":      prk_f, "prk_freq_chain_sem":      prk_s,
                "gap2_freq_chain_rep_avg":     gap_f, "gap2_freq_chain_sem":     gap_s,
                "receptor_freq_chain_rep_avg": rec_f, "receptor_freq_chain_sem": rec_s,
                "prk_partner_count":  prk_cnt,
                "gap2_partner_count": gap_cnt,
                "interface_class":    cls,
            })
    return pd.DataFrame(out_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interface_root", required=True, type=Path,
                    help="dir containing <condition>/<rep>/interface_resfreq_cp12_side_*.csv")
    ap.add_argument("--no-reps", action="store_true",
                    help="CSVs live directly at <interface_root>/<condition>/ with no rep subdirectory")
    ap.add_argument("--pairs_dir_med4", required=True, type=Path,
                    help="dir containing differential_pairs_{prk,gap2}.csv for MED4")
    ap.add_argument("--pairs_dir_hm2", required=True, type=Path,
                    help="dir containing differential_pairs_{prk,gap2}.csv for HM2")
    ap.add_argument("--alignment_csv", required=True, type=Path,
                    help="CP12_MED4_vs_HM2_alignment_annotation.csv")
    ap.add_argument("--out", required=True, type=Path,
                    help="output master CSV path")
    args = ap.parse_args()

    print(f"[1/2] averaging per-chain × per-rep freqs from {args.interface_root}")
    chain_rep_df = chain_and_rep_average(args.interface_root, no_reps=args.no_reps)

    pairs_dirs = {"med4": args.pairs_dir_med4, "hm2": args.pairs_dir_hm2}
    print(f"[2/2] building master CSV with partner counts (med4: {args.pairs_dir_med4}, hm2: {args.pairs_dir_hm2})")
    master = build_master(chain_rep_df, pairs_dirs, args.alignment_csv)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(master)} rows)")


if __name__ == "__main__":
    main()
