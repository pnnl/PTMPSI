"""
Four-panel contact frequency correlation plot.
Rows: CP12-PRK (top) and CP12-GAP (bottom)
Cols: Oxidized (left) and Reduced (right)
X-axis: MED4 contact frequency
Y-axis: HM2 contact frequency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats

# ── Load & prep ───────────────────────────────────────────────────────────────
df = pd.read_csv("FINAL_FIGURES/cp12_contact_freq_heatmap_data_master.csv")

# Drop rows that are gaps in either host
df = df[df["is_gap_for_host"] == False].copy()

def get_merged(redox_label):
    """Merge MED4 and HM2 rows for a given redox condition on alignment column."""
    med4 = df[df["condition"] == f"med4_{redox_label}"].copy()
    hm2  = df[df["condition"] == f"hm2_{redox_label}"].copy()
    merged = pd.merge(
        med4[["aln_col_1based", "prk_freq_chain_rep_avg", "gap2_freq_chain_rep_avg",
              "interface_class", "med4_aa", "hm2_aa"]],
        hm2[["aln_col_1based",  "prk_freq_chain_rep_avg", "gap2_freq_chain_rep_avg",
              "interface_class"]],
        on="aln_col_1based",
        suffixes=("_med4", "_hm2"),
    )
    return merged

ox  = get_merged("ox")
red = get_merged("red")

# ── Colour map for interface class ────────────────────────────────────────────
class_palette = {
    "PRK":     "#2166ac",   # blue
    "GAP2":    "#d6604d",   # red
    "Neither": "#999999",   # grey
    "Gap":     "#cccccc",   # light grey (shouldn't appear after filter)
}
# Use MED4 classification for colouring
def class_colors(series):
    return [class_palette.get(c, "#aaaaaa") for c in series]

# ── Helper: scatter one panel ─────────────────────────────────────────────────
def add_panel(ax, x, y, class_col, x_label=None, y_label=None, title=None, color="#4393c3"):
    # Drop NaN pairs
    mask = x.notna() & y.notna()
    xv, yv, cv = x[mask].values, y[mask].values, class_col[mask]

    ax.scatter(xv, yv, color=color, s=30, alpha=0.75, linewidths=0.3,
               edgecolors="white", zorder=3)

    # 1:1 diagonal reference line
    lim = max(xv.max(), yv.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="black", lw=0.8, ls="--",
            alpha=0.5, zorder=2)

    # Linear regression
    if len(xv) >= 3:
        slope, intercept, r, p, se = stats.linregress(xv, yv)
        xfit = np.linspace(0, lim, 200)
        yfit = slope * xfit + intercept
        ax.plot(xfit, yfit, color="dimgray", lw=1.4, ls="-", zorder=4)
        p_str = f"p < 0.001" if p < 0.001 else f"p = {p:.3f}"
        ax.text(0.97, 0.05, f"r = {r:.2f}\n{p_str}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="dimgray",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7))

    ax.set_xlim(-0.02, lim)
    ax.set_ylim(-0.02, lim)
    ax.set_aspect("equal")

    ax.set_xlabel(x_label or "MED4 contact frequency", fontsize=9)
    ax.set_ylabel(y_label or "HM2 contact frequency", fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold", pad=5)

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
fig.subplots_adjust(hspace=0.38, wspace=0.38)

# Top-left  : PRK, oxidized
add_panel(axes[0, 0],
          ox["prk_freq_chain_rep_avg_med4"],
          ox["prk_freq_chain_rep_avg_hm2"],
          ox["interface_class_med4"],
          x_label="MED4 CP12–PRK contact freq.",
          y_label="HM2 CP12–PRK contact freq.",
          title="CP12–PRK  |  Oxidized",
          color="u")

# Top-right : PRK, reduced
add_panel(axes[0, 1],
          red["prk_freq_chain_rep_avg_med4"],
          red["prk_freq_chain_rep_avg_hm2"],
          red["interface_class_med4"],
          x_label="MED4 CP12–PRK contact freq.",
          y_label="HM2 CP12–PRK contact freq.",
          title="CP12–PRK  |  Reduced",
          color="#E691E8")

# Bottom-left : GAP, oxidized
add_panel(axes[1, 0],
          ox["gap2_freq_chain_rep_avg_med4"],
          ox["gap2_freq_chain_rep_avg_hm2"],
          ox["interface_class_med4"],
          x_label="MED4 CP12–GAP contact freq.",
          y_label="HM2 CP12–GAP contact freq.",
          title="CP12–GAP  |  Oxidized",
          color="darkgray")

# Bottom-right: GAP, reduced
add_panel(axes[1, 1],
          red["gap2_freq_chain_rep_avg_med4"],
          red["gap2_freq_chain_rep_avg_hm2"],
          red["interface_class_med4"],
          x_label="MED4 CP12–GAP contact freq.",
          y_label="HM2 CP12–GAP contact freq.",
          title="CP12–GAP  |  Reduced",
          color="darkgray")


fig.suptitle("CP12 contact frequency correlation (5Å cutoff): MED4 vs. HM2",
             fontsize=12, fontweight="bold", y=1.01)

# ── Save ──────────────────────────────────────────────────────────────────────
out = "FINAL_FIGURES/cp12_contact_freq_correlation.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"Saved → {out}")
plt.close()
