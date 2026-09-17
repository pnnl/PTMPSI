"""
Plot physico-chemical changes between two aligned sequences from a FASTA file.
Computes per-position changes in hydrophobicity, volume, and charge (HM - MED).
Gaps (. or -) are treated as absent residue (property = 0).
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ── Property tables ─────────────────────────────────────────────────────────

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

# Amino acid van der Waals / side-chain volumes (Å³), Pontius et al.
VOLUME = {
    'A':  88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
    'Q': 143.8, 'E': 138.4, 'G':  60.1, 'H': 153.2, 'I': 166.7,
    'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
    'S':  89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0,
}

# Formal charge at physiological pH (~7.4)
CHARGE = {
    'A':  0, 'R':  1, 'N':  0, 'D': -1, 'C':  0,
    'Q':  0, 'E': -1, 'G':  0, 'H':  0.1, 'I':  0,
    'L':  0, 'K':  1, 'M':  0, 'F':  0, 'P':  0,
    'S':  0, 'T':  0, 'W':  0, 'Y':  0, 'V':  0,
}

GAP_CHARS = set('.-')


def is_gap(c):
    return c in GAP_CHARS


def prop(table, aa):
    """Return property value; 0 for gaps or unknown characters."""
    return table.get(aa.upper(), 0) if not is_gap(aa) else 0


# ── FASTA parser ─────────────────────────────────────────────────────────────

def parse_fasta(path):
    sequences = {}
    current_name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                current_name = line[1:].strip()
                sequences[current_name] = ''
            elif current_name is not None:
                sequences[current_name] += line
    return sequences


# ── Load alignment ────────────────────────────────────────────────────────────

fasta_path = Path(__file__).parent.parent / 'Data' / 'HMM_seq_align.fasta'
seqs = parse_fasta(fasta_path)
names = list(seqs.keys())

if len(seqs) < 2:
    raise ValueError(f"Need at least 2 sequences; found {len(seqs)}.")

name1, name2 = names[0], names[1]
seq1, seq2 = seqs[name1], seqs[name2]

if len(seq1) != len(seq2):
    raise ValueError(
        f"Aligned sequences have different lengths: {len(seq1)} vs {len(seq2)}. "
        "Ensure the input is a proper multiple-sequence alignment."
    )

# ── Compute per-position differences (seq2 − seq1) ──────────────────────────

positions = []       # alignment column index (1-based display label)
labels = []          # "aa1 / aa2" strings
hyd_changes = []
vol_changes = []
charge_changes = []
is_different = []    # True if the two residues differ

for i, (a1, a2) in enumerate(zip(seq1, seq2)):
    # skip columns where either sequence has a gap (insertion/deletion)
    if is_gap(a1) or is_gap(a2):
        continue

    positions.append(i + 1)
    labels.append(f"{a1.upper()} / {a2.upper()}")

    hyd_changes.append(prop(HYDROPHOBICITY, a2) - prop(HYDROPHOBICITY, a1))
    vol_changes.append(prop(VOLUME, a2) - prop(VOLUME, a1))
    charge_changes.append(prop(CHARGE, a2) - prop(CHARGE, a1))
    is_different.append(a1.upper() != a2.upper())

n = len(positions)
y = np.arange(n)

# ── Averages over changed positions only ────────────────────────────────────

def mean_changed(values, diffs):
    changed_vals = [v for v, d in zip(values, diffs) if d]
    return np.mean(changed_vals) if changed_vals else 0.0

avg_vol = mean_changed(vol_changes, is_different)
avg_hyd = mean_changed(hyd_changes, is_different)
avg_chg = mean_changed(charge_changes, is_different)

print(f"Average Δ (changed residues only):")
print(f"  Volume:          {avg_vol:+.2f} Å³")
print(f"  Hydrophobicity:  {avg_hyd:+.3f} (Kyte-Doolittle)")
print(f"  Charge:          {avg_chg:+.3f} (formal)")

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(11, max(8, n * 0.22)))
ax_vol, ax_hyd, ax_chg = axes

bar_pos = '#2196F3'   # blue – positive change
bar_neg = '#F44336'   # red  – negative change
bar_zero = '#B0BEC5'  # grey – no change


def color_bars(ax, values):
    for bar, v in zip(ax.patches, values):
        if v > 0:
            bar.set_color(bar_pos)
        elif v < 0:
            bar.set_color(bar_neg)
        else:
            bar.set_color(bar_zero)


for ax, data, xlabel, avg in [
    (ax_vol, vol_changes,    'Volume Change (Å³)',                  avg_vol),
    (ax_hyd, hyd_changes,    'Hydrophobicity Change\n(Kyte-Doolittle)', avg_hyd),
    (ax_chg, charge_changes, 'Charge Change\n(Formal, pH 7.4)',     avg_chg),
]:
    ax.barh(y, data, color=bar_pos)
    color_bars(ax, data)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.axvline(avg, color='black', linewidth=1.2, linestyle='--',
               label=f'Mean Δ = {avg:+.2f}')
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    ax.set_yticks([])
    ax.set_xlabel(xlabel, fontsize=9)
    ax.legend(fontsize=7.5, loc='lower right')
    ax.invert_yaxis()

# Sequence labels on left side of volume plot
label_x = ax_vol.get_xlim()[0]
for j, (label, diff) in enumerate(zip(labels, is_different)):
    color = 'teal' if diff else 'orange'
    ax_vol.text(
        label_x - 5, j, label,
        va='center', ha='right', fontsize=6.5, color=color,
        fontfamily='monospace'
    )

# Legend for label colours
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='teal',   lw=2, label='Changed'),
    Line2D([0], [0], color='orange', lw=2, label='Identical'),
    Line2D([0], [0], color=bar_pos,  lw=4, label=f'+Δ ({name2} > {name1})'),
    Line2D([0], [0], color=bar_neg,  lw=4, label=f'−Δ ({name2} < {name1})'),
]
fig.legend(
    handles=legend_elements, loc='lower center',
    ncol=4, fontsize=8, bbox_to_anchor=(0.5, 0.01)
)

plt.suptitle(
    f'Physico-Chemical Changes\n MED4 vs P-HM2 (HMM alignment)',
    fontsize=13, y=0.99
)
plt.tight_layout(rect=[0, 0.06, 1, 0.97])

out = Path(__file__).parent / 'physchem_alignment_plot.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved → {out}")
plt.show()
