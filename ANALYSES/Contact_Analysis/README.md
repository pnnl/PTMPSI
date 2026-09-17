# Contact Analysis Scripts

This directory contains Python scripts for analyzing protein-protein interface contacts and contact rewiring patterns from molecular dynamics (MD) trajectories, specifically focused on CP12 protein interactions with PRK and GAP2 partners across different redox conditions and host strains (MED4 and HM2).

## Overview

The Contact_Analysis workflow consists of two main scripts that work together to:

1. **Compute contact frequencies** between residues across MD trajectory replicates
2. **Analyze pairwise rewiring patterns** comparing different conditions (host × redox state)
3. **Generate master datasets** for visualization of contact frequency heatmaps

## Script Descriptions

### `analyze_cp12_contact_pairs.py`

Computes residue-pair contact frequencies and differential contact rewiring maps comparing conditions.

**What it does:**
- Extracts all atomic coordinates from MD trajectories
- Identifies spatial contacts between CP12 residues and partner proteins (PRK, GAP2)
- Aggregates pairwise contact frequencies across replicates
- Collapses multiple CP12 chains to canonical alignment positions
- Computes differential contact frequencies between conditions (host × redox)

**Outputs:**
- Per-replicate pair frequencies: `cp12_contact_pairs/<condition>/rep{N}/contact_pairs_{partner}.csv`
- Aggregated per-condition: `cp12_contact_pairs/<condition>/<condition>_contact_pairs_agg_{partner}.csv`
- Collapsed to alignment coordinates: `cp12_contact_pairs/<condition>/<condition>_contact_pairs_collapsed_{partner}.csv`
- Differential tables: `cp12_contact_pairs/differential_pairs_{partner}.csv`
- Top rewired pairs: `cp12_contact_pairs/top_rewired_pairs_{partner}.csv`

**Key parameters:**
- `--cutoff`: Heavy-atom distance threshold (default: 6.0 Å)
- `--stride`: Frame sampling interval (default: 10 frames)
- `--systems`: Condition labels to process (default: med4_ox, med4_red, hm2_ox, hm2_red)

---

### `interface_contacts_freq_create_master.py`

Aggregates contact frequency data and generates a master CSV for visualization.

**What it does:**
- Loads per-replicate contact frequencies from interface analysis
- Averages frequencies across CP12 chains (C, H, M, P) within each replicate
- Computes mean and standard error across replicates
- Maps residues to HMM alignment columns
- Joins with differential pair data to classify interface contacts
- Generates a unified master CSV for heatmap figures

**Outputs:**
- Master heatmap data: Single CSV file specified by `--out`

**Key parameters:**
- `--interface_root`: Directory containing per-condition/per-replicate contact frequency data
- `--pairs_dir_med4`: Directory with differential pairs data for MED4
- `--pairs_dir_hm2`: Directory with differential pairs data for HM2
- `--alignment_csv`: CP12 sequence alignment annotation file
- `--out`: Output path for master CSV

**Thresholds:**
- `PAIR_THRESH`: Minimum frequency for partner residues (default: 0.10)
- `INTERFACE_THRESH`: Classification threshold for interface contacts (default: 0.50)

---

## Input File Structure

### Trajectory Input Structure

Both scripts require MD trajectory files organized as follows:

```
<trajectory_root>/
├── med4/
│   ├── rep1/
│   │   ├── solute.tpr              # GROMACS topology
│   │   ├── solute_fit.xtc          # MD trajectory (fitted)
│   │   └── solute_fixed.pdb        # Reference structure with chain IDs
│   ├── rep2/
│   │   └── (same files)
│   └── rep3/
│       └── (same files)
└── hm2/
    ├── rep1/
    │   ├── solute.tpr
    │   ├── solute_fit.xtc
    │   └── solute_fixed.pdb
    ├── rep2/
    └── rep3/
```

**File descriptions:**
- **solute.tpr**: GROMACS binary topology file containing molecular definition
- **solute_fit.xtc**: GROMACS trajectory file with atomic coordinates aligned to reference
- **solute_fixed.pdb**: PDB reference structure (must have chain IDs in ATOM records)

### Contact Frequency Input Format

When running `interface_contacts_freq_create_master.py`, expected input CSVs from `--interface_root` have structure:

```
<interface_root>/
├── med4_ox/
│   ├── rep1/
│   │   ├── interface_resfreq_cp12_side_prk.csv
│   │   ├── interface_resfreq_cp12_side_gap2.csv
│   │   └── interface_resfreq_cp12_side_receptor.csv
│   ├── rep2/
│   └── rep3/
├── med4_red/
│   └── (same structure)
├── hm2_ox/
└── hm2_red/
```

**CSV columns** in `interface_resfreq_cp12_side_*.csv`:
```
resid          : Residue ID (integer)
resname        : Residue name (3-letter code, e.g., ALA)
chain          : Chain identifier (C, H, M, or P for CP12)
freq           : Contact frequency (0.0–1.0)
```

### Differential Pairs Input Format

Expected input from `--pairs_dir_*`:

```
<pairs_dir>/
├── differential_pairs_prk.csv
├── differential_pairs_gap2.csv
└── (optional) top_rewired_pairs_*.csv
```

**CSV columns** in `differential_pairs_*.csv`:
```
cp12_aln_pos          : CP12 position in HMM alignment (1-based)
partner_resid         : Partner residue ID
partner_resname       : Partner residue name
cp12_resname          : CP12 residue name
med4_ox_freq          : Contact frequency in MED4 oxidized condition
med4_red_freq         : Contact frequency in MED4 reduced condition
hm2_ox_freq           : Contact frequency in HM2 oxidized condition
hm2_red_freq          : Contact frequency in HM2 reduced condition
delta_med4_ox_vs_med4_red  : Frequency difference (redox delta)
delta_med4_ox_vs_hm2_ox    : Frequency difference (strain delta)
... (additional delta columns)
```

### Alignment Annotation Input

The `--alignment_csv` file maps residues between hosts and HMM alignment positions:

```
aln_col_1based   : Alignment column (1-based, 1–70)
med4_pos         : MED4 residue ID (integer or NA)
med4_aa          : MED4 residue name
hm2_pos          : HM2 residue ID (integer or NA)
hm2_aa           : HM2 residue name
conservation     : (optional) conservation score
```

**Example:**
```
aln_col_1based,med4_pos,med4_aa,hm2_pos,hm2_aa
1,1,M,1,M
2,2,D,2,D
3,,,,
4,3,T,3,T
```

---

## Usage

### Step 1: Run Contact Pair Analysis

```bash
cd /path/to/analysis/root
python analyze_cp12_contact_pairs.py \
  --alignment_csv CP12_MED4_vs_HM2_alignment_annotation.csv \
  --outdir cp12_contact_pairs \
  --cutoff 6.0 \
  --stride 10 \
  --systems med4_ox med4_red hm2_ox hm2_red
```

**Expected runtime:** ~1–2 hours per system (depending on trajectory length and frame stride)

### Step 2: Create Master CSV for Visualization

```bash
python interface_contacts_freq_create_master.py \
  --interface_root /path/to/interface_results \
  --pairs_dir_med4 ./cp12_contact_pairs/med4 \
  --pairs_dir_hm2 ./cp12_contact_pairs/hm2 \
  --alignment_csv CP12_MED4_vs_HM2_alignment_annotation.csv \
  --out cp12_contact_freq_heatmap_data_master.csv
```

---

## Data Processing Pipeline

The workflow processes data through the following stages:

### `analyze_cp12_contact_pairs.py` Pipeline:

```
Trajectories (solute.tpr, solute_fit.xtc, solute_fixed.pdb)
    ↓
Per-replicate contact detection (FastNS spatial neighbor search)
    ↓
Per-rep contact_pairs_{partner}.csv
    ↓
Aggregate across replicates (mean ± SEM)
    ↓
contact_pairs_agg_{partner}.csv
    ↓
Collapse across CP12 chains (C/H/M/P) using alignment
    ↓
contact_pairs_collapsed_{partner}.csv (alignment coordinates)
    ↓
Compute pairwise deltas across conditions
    ↓
differential_pairs_{partner}.csv
```

### `interface_contacts_freq_create_master.py` Pipeline:

```
interface_resfreq_cp12_side_*.csv (per-rep per-chain)
    ↓
Average across CP12 chains (C/H/M/P)
    ↓
Compute mean ± SEM across replicates
    ↓
differential_pairs_*.csv
    ↓
Map residues to HMM alignment columns
    ↓
Join and classify contacts
    ↓
cp12_contact_freq_heatmap_data_master.csv
```

---

## Technical Details

### Chain Mapping

- **CP12 chains:** C, H, M, P (4 identical copies in the complex)
- **PRK chains:** B, G, J, L
- **GAP2 chains:** A, D, E, F, I, K, N, O
- **NAD chains:** Q (optional, used in some analyses)

Chain assignments are inferred from the PDB reference structure's ATOM records and validated against trajectory topology.

### Contact Definition

A **contact** between two residues is counted when any heavy atom from residue 1 falls within the cutoff distance from any heavy atom of residue 2. Each residue pair is counted at most once per frame (binary contact per frame).

**Contact frequency** = (number of frames with contact) / (total frames analyzed)

### Alignment Coordinates

CP12 is a 70-residue protein. Residues from multiple chains (C/H/M/P) are mapped to shared HMM alignment positions (1–70) to enable cross-replicate and cross-condition comparisons.

---

## Dependencies

- **Python 3.7+**
- **MDAnalysis** (trajectory analysis)
- **pandas** (data manipulation)
- **numpy** (numerical computing)

Install via:
```bash
pip install MDAnalysis pandas numpy
```

---

## Troubleshooting

### Missing trajectory files
**Error:** `FileNotFoundError: Missing: /path/to/solute.tpr`

**Solution:** Ensure trajectory structure matches the expected layout. Check that all three files exist in each replicate directory.

### Chain ID mismatch
**Error:** `ValueError: Atom count mismatch` or `Atom key multiset mismatch`

**Solution:** Verify that `solute_fixed.pdb` was generated from the same `solute.tpr` and trajectory. Regenerate if necessary.

### Condition not found
**Error:** `discovered condition paths: []`

**Solution:** Check that interface CSV files exist in the `--interface_root` directory. Verify filenames match expected pattern: `interface_resfreq_cp12_side_{prk,gap2,receptor}.csv`

### Missing alignment CSV
**Error:** `FileNotFoundError: CP12_MED4_vs_HM2_alignment_annotation.csv`

**Solution:** Provide the alignment annotation file via `--alignment_csv` or ensure it exists in the current working directory.

---

## Output File Descriptions

### From `analyze_cp12_contact_pairs.py`

| Filename | Rows | Columns | Purpose |
|----------|------|---------|---------|
| `contact_pairs_{partner}.csv` | Per-rep × per-pair | cp12_chain, cp12_resid, partner_resid, freq, ... | Raw pairwise frequencies per replicate |
| `contact_pairs_agg_{partner}.csv` | Per-pair | cp12_chain, cp12_resid, partner_resid, mean_freq, sem_freq, n_reps | Aggregated across replicates |
| `contact_pairs_collapsed_{partner}.csv` | Per-(aln_pos × partner_resid) | cp12_aln_pos, partner_resid, max_freq | Collapsed to alignment coordinates |
| `differential_pairs_{partner}.csv` | Per-(aln_pos × partner_resid) | aln_pos, partner_resid, med4_ox_freq, ..., delta_*, ... | All pairwise condition deltas |

### From `interface_contacts_freq_create_master.py`

| Filename | Rows | Columns | Purpose |
|----------|------|---------|---------|
| `cp12_contact_freq_heatmap_data_master.csv` | Per-(condition × host × aln_col) | condition, host, aln_col, resid, resname, prk_f, prk_s, gap_f, gap_s, rec_f, rec_s, ... | Master table for visualization |

---

## References

- MDAnalysis documentation: https://docs.mdanalysis.org/
- GROMACS file formats: https://manual.gromacs.org/
