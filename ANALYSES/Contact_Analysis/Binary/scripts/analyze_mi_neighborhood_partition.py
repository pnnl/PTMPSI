#!/usr/bin/env python3
"""
Spatial-neighborhood partition of the CP12 HMM-column MI coupling network.

Partitions the 70 HMM columns into three blocks chosen to match the structural
and functional zones of CP12:
  - Block N (1–22)   the divergent N-terminus; the substituted-position cluster.
  - Block M (23–50)  the structured PRK-binding middle (the brightest CP12-PRK
                     dwell band around HMM cols 27–40).
  - Block C (51–70)  the C-terminal regulatory region (regulatory cysteines at
                     HMM cols 57, 66; the longest-lived GAP2 contacts).

For each (condition, partner) MI matrix from analyze_contact_coupling_mi.py:
  intra_N = sum_{i<j in N}     MI[i, j]
  intra_M = sum_{i<j in M}     MI[i, j]
  intra_C = sum_{i<j in C}     MI[i, j]
  cross_NM = sum_{i in N, j in M} MI[i, j]
  cross_NC = sum_{i in N, j in C} MI[i, j]
  cross_MC = sum_{i in M, j in C} MI[i, j]
total = intra_N + intra_M + intra_C + cross_NM + cross_NC + cross_MC
fraction of total per block is reported, and host-difference (MED4 − HM2) is
computed per redox state.

Outputs:
    contact_coupling/coupling_neighborhood_partition.csv
    contact_coupling/coupling_neighborhood_partition.md
    contact_coupling_plots/fig_mi_neighborhood_partition_{prk,gap2}.{pdf,png}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_ORDER = ["med4_ox", "med4_red", "hm2_ox", "hm2_red"]
CONDITION_COLORS = {
    "med4_ox":  "#1f77b4",
    "med4_red": "#aec7e8",
    "hm2_ox":   "#d62728",
    "hm2_red":  "#ff9896",
}

BLOCKS = {
    "N (cols 1-22)":   list(range(1, 23)),      # divergent N-terminus
    "M (cols 23-50)":  list(range(23, 51)),     # structured PRK-binding middle
    "C (cols 51-70)":  list(range(51, 71)),     # C-terminal regulatory
}
BLOCK_ORDER = list(BLOCKS.keys())


def partition_mi(mi: np.ndarray, blocks: dict[str, list[int]]) -> dict[str, float]:
    """Return intra_X and cross_XY sums for each block pair."""
    # mi is 70x70 indexed 0-based; blocks use 1-based HMM cols
    out: dict[str, float] = {}
    block_idx = {name: np.array([c - 1 for c in cols if 1 <= c <= mi.shape[0]])
                  for name, cols in blocks.items()}
    names = list(block_idx)
    for i, ni in enumerate(names):
        idx_i = block_idx[ni]
        for j in range(i, len(names)):
            nj = names[j]
            idx_j = block_idx[nj]
            sub = mi[np.ix_(idx_i, idx_j)]
            if i == j:
                # upper triangle (no double counting)
                key = f"intra__{ni}"
                out[key] = float(np.triu(sub, k=1).sum())
            else:
                key = f"cross__{ni}__{nj}"
                out[key] = float(sub.sum())
    out["total"] = float(np.triu(mi, k=1).sum())
    return out


def make_table(npz_path: Path) -> pd.DataFrame:
    arr = np.load(npz_path)
    rows = []
    for cond in CONDITION_ORDER:
        for partner in ("prk", "gap2"):
            key = f"mi_{cond}_{partner}"
            if key not in arr.files:
                continue
            mi = arr[key]
            stats = partition_mi(mi, BLOCKS)
            row = dict(condition=cond, partner=partner, **stats)
            rows.append(row)
    return pd.DataFrame(rows)


def host_diff_block(df: pd.DataFrame, redox: str, partner: str) -> dict:
    sel = df[(df["condition"].str.endswith("_" + redox)) & (df["partner"] == partner)].copy()
    sel["host"] = sel["condition"].str[:4].str.rstrip("_")
    out = {}
    for col in [c for c in sel.columns if c.startswith(("intra__", "cross__", "total"))]:
        med = sel.loc[sel["host"] == "med4", col]
        hm  = sel.loc[sel["host"] == "hm2",  col]
        if med.empty or hm.empty:
            out[col + "__diff_med4_minus_hm2"] = float("nan")
        else:
            out[col + "__diff_med4_minus_hm2"] = float(med.iloc[0] - hm.iloc[0])
    return out


def make_stacked_bar(df: pd.DataFrame, partner: str, outdir: Path) -> None:
    sub = df[df["partner"] == partner].set_index("condition").reindex(CONDITION_ORDER)
    # 6 block contributions, in display order: intra_N, intra_M, intra_C, cross_NM, cross_NC, cross_MC
    blocks_in_order = [
        "intra__N (cols 1-22)",
        "intra__M (cols 23-50)",
        "intra__C (cols 51-70)",
        "cross__N (cols 1-22)__M (cols 23-50)",
        "cross__N (cols 1-22)__C (cols 51-70)",
        "cross__M (cols 23-50)__C (cols 51-70)",
    ]
    display_labels = [
        "intra N (1-22)",
        "intra M (23-50)",
        "intra C (51-70)",
        "cross N-M",
        "cross N-C",
        "cross M-C",
    ]
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    x = np.arange(len(sub))
    bottoms = np.zeros(len(sub))
    for b, lbl, c in zip(blocks_in_order, display_labels, colors):
        vals = sub[b].to_numpy()
        ax.bar(x, vals, bottom=bottoms, label=lbl, color=c, edgecolor="white", lw=0.4)
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels(sub.index.tolist())
    ax.set_ylabel("Σ MI per block (bits)")
    ax.set_title(f"CP12 ↔ {partner.upper()} coupling decomposed by spatial block")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.12))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"fig_mi_neighborhood_partition_{partner}.pdf")
    fig.savefig(outdir / f"fig_mi_neighborhood_partition_{partner}.png", dpi=300)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="contact_coupling/contact_coupling_mi.npz")
    ap.add_argument("--outdir", default="contact_coupling")
    ap.add_argument("--plot_outdir", default="contact_coupling_plots")
    args = ap.parse_args()

    npz_path = Path(args.npz).resolve()
    outdir = Path(args.outdir).resolve()
    plot_outdir = Path(args.plot_outdir).resolve()

    df = make_table(npz_path)
    df.to_csv(outdir / "coupling_neighborhood_partition.csv", index=False)

    # host-diff tables per redox per partner
    host_rows = []
    for partner in ("prk", "gap2"):
        for redox in ("ox", "red"):
            hd = host_diff_block(df, redox, partner)
            row = dict(partner=partner, redox=redox, **hd)
            host_rows.append(row)
    host_df = pd.DataFrame(host_rows)
    host_df.to_csv(outdir / "coupling_neighborhood_host_diff.csv", index=False)

    # markdown summary
    lines = ["# Spatial-neighborhood partition of CP12 HMM-column MI coupling networks",
             "",
             "CP12 HMM columns are partitioned into three structural/functional zones:",
             "",
             "- **Block N (cols 1–22)** divergent N-terminus; substituted-position cluster",
             "- **Block M (cols 23–50)** structured PRK-binding middle",
             "- **Block C (cols 51–70)** C-terminal regulatory region (regulatory cys at 57, 66)",
             "",
             "For each (condition × partner) MI matrix, the upper triangle is decomposed into "
             "three intra-block sums (within N, M, C) and three cross-block sums (N↔M, N↔C, M↔C). "
             "Edges with MI ≥ 0.05 bits are the same as those reported in the main coupling-graph "
             "analysis (Main Figure 5B).",
             "",
             "## Total MI per condition × partner (sum of upper triangle, bits)",
             "",
             "| partner | condition | intra N | intra M | intra C | cross N-M | cross N-C | cross M-C | total |",
             "|---|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(
            f"| {r['partner']} | {r['condition']} | "
            f"{r['intra__N (cols 1-22)']:.2f} | "
            f"{r['intra__M (cols 23-50)']:.2f} | "
            f"{r['intra__C (cols 51-70)']:.2f} | "
            f"{r['cross__N (cols 1-22)__M (cols 23-50)']:.2f} | "
            f"{r['cross__N (cols 1-22)__C (cols 51-70)']:.2f} | "
            f"{r['cross__M (cols 23-50)__C (cols 51-70)']:.2f} | "
            f"{r['total']:.2f} |"
        )
    lines.append("")
    lines += ["## Host difference per block (MED4 − HM2)",
              "",
              "Positive value ⇒ MED4 has more coupling in that block; negative ⇒ HM2 has more.",
              "",
              "| partner | redox | Δ intra N | Δ intra M | Δ intra C | Δ cross N-M | Δ cross N-C | Δ cross M-C | Δ total |",
              "|---|---|---|---|---|---|---|---|---|"]
    for _, r in host_df.iterrows():
        keys = [
            "intra__N (cols 1-22)__diff_med4_minus_hm2",
            "intra__M (cols 23-50)__diff_med4_minus_hm2",
            "intra__C (cols 51-70)__diff_med4_minus_hm2",
            "cross__N (cols 1-22)__M (cols 23-50)__diff_med4_minus_hm2",
            "cross__N (cols 1-22)__C (cols 51-70)__diff_med4_minus_hm2",
            "cross__M (cols 23-50)__C (cols 51-70)__diff_med4_minus_hm2",
            "total__diff_med4_minus_hm2",
        ]
        vals = [r[k] for k in keys]
        lines.append(
            f"| {r['partner']} | {r['redox']} | "
            + " | ".join(f"{v:+.2f}" for v in vals) + " |"
        )
    lines.append("")
    (outdir / "coupling_neighborhood_partition.md").write_text("\n".join(lines) + "\n")

    # plots
    for partner in ("prk", "gap2"):
        make_stacked_bar(df, partner, plot_outdir)

    print(f"wrote {outdir / 'coupling_neighborhood_partition.csv'}")
    print(f"wrote {outdir / 'coupling_neighborhood_partition.md'}")
    print(f"wrote {plot_outdir / 'fig_mi_neighborhood_partition_*.{pdf,png}'}")


if __name__ == "__main__":
    main()
