# SASA-BSA Analysis

This directory contains a Python script for computing solvent-accessible surface area (SASA) and buried surface area (BSA) for the CP12-partner protein complex using GROMACS.

## Overview

`complex_sasa_gmx.py` analyzes the buried surface area between CP12 and its binding partners (PRK and GAP2/GAPDH) across different redox conditions and host strains (MED4 and HM2). The script:

1. **Computes SASA** for multiple protein configurations:
   - Complete protein complex
   - Individual components (CP12, PRK, GAP2) in isolation
   - Components within the complex
   - Partner pair complexes (PRK+GAP2)
   - Binary complexes (CP12+PRK, CP12+GAP2)

2. **Calculates buried surface area (BSA)** using the thermodynamic formula:
   ```
   BSA(A,B) = SASA(A_isolated) + SASA(B_isolated) - SASA(A+B_complex)
   ```

3. **Supports both static and dynamic analysis**:
   - Static mode: analyzes a single frame (PDB structure or trajectory snapshot)
   - Dynamic mode: analyzes full trajectories to compute time-averaged SASA

4. **Optionally maps CP12 residues to HMM alignment columns** for cross-strain comparisons

---

## Script: `complex_sasa_gmx.py`

### What it does

- Reads PDB structures to extract chain information
- Creates GROMACS index files for different protein components
- Runs `gmx sasa` for all specified configurations
- Parses GROMACS XVG output files into CSV tables
- Computes buried surface area (BSA) from isolated and complex SASA values
- Optionally aggregates per-residue SASA from atom-level data

### Dependencies

- **GROMACS** (installed and accessible as `gmx` or via `--gmx-bin`)
- **Python 3.7+**
- **pandas** (data manipulation)
- **numpy** (numerical computing)

Install Python dependencies:
```bash
pip install pandas numpy
```

---

## Input File Structure

The script expects MD trajectory files organized as follows:

```
{analysis_dir}/
├── med4/
│   ├── 0000/          # reduced condition
│   │   ├── rep1/
│   │   │   ├── solute_fixed.pdb       # Reference structure with chain IDs
│   │   │   ├── solute.tpr             # GROMACS topology
│   │   │   └── solute_fit.xtc         # MD trajectory (fitted)
│   │   ├── rep2/
│   │   │   └── (same files)
│   │   └── rep3/
│   └── 0003/          # oxidized condition
│       ├── rep1/
│       └── (same as 0000)
└── hm2/
    ├── 0000/
    │   └── rep{1,2,3}/
    └── 0003/
        └── rep{1,2,3}/
```

**File descriptions:**
- **solute_fixed.pdb**: PDB reference structure with explicit chain IDs in ATOM records
- **solute.tpr**: GROMACS binary topology file
- **solute_fit.xtc**: GROMACS trajectory file (atoms fitted/aligned to reference)

### Condition Codes

- **0000**: Reduced form (NAD+ bound to PRK in reduced state)
- **0003**: Oxidized form (NADP+ or NAD in oxidized state)

### Protein Components

The script automatically identifies protein components from PDB chain IDs:

| Component | Chains | Count | Notes |
|-----------|--------|-------|-------|
| **CP12** | C, H, M, P | 4 copies | Central regulatory protein (homotetrameric) |
| **PRK** | B, G, J, L | 4 copies | Phosphoribulokinase (paired with CP12 chains) |
| **GAP2** | A, D, E, F, I, K, N, O | 8 copies | Glyceraldehyde-3-phosphate dehydrogenase |
| **NAD** | Q | 1+ copy | Cofactor (optional) |

---

## Usage

### Basic Usage

```bash
python complex_sasa_gmx.py \
  --analysis-dir /path/to/trajectory/root \
  --out-dir /path/to/output \
  --mode dynamic
```

### Example 1: Dynamic Analysis (Full Trajectory)

```bash
python complex_sasa_gmx.py \
  --analysis-dir data_dir/0000_0003_analysis \
  --out-dir results/complex_sasa_bsa_results \
  --strains med4 hm2 \
  --systems 0000 0003 \
  --replicates rep1 rep2 rep3 \
  --mode dynamic \
  --begin-ps 10000 \
  --end-ps 100000 \
  --dt-ps 1000
```

**Parameters:**
- `--begin-ps 10000`: Skip first 10 ns (equilibration)
- `--end-ps 100000`: Analyze up to 100 ns
- `--dt-ps 1000`: Subsample every 1000 ps (sample ~every 20 frames)

### Example 2: Static Analysis (Single Frame)

```bash
python complex_sasa_gmx.py \
  --analysis-dir data_dir/0000_0003_analysis \
  --out-dir results/complex_sasa_bsa_results \
  --mode static \
  --static-source pdb
```

### Example 3: Mixed Static + Dynamic

```bash
python complex_sasa_gmx.py \
  --analysis-dir data_dir/0000_0003_analysis \
  --out-dir results/complex_sasa_bsa_results \
  --mode both
```

### Example 4: With HMM Alignment Mapping

```bash
python complex_sasa_gmx.py \
  --analysis-dir data_dir/0000_0003_analysis \
  --out-dir results/complex_sasa_bsa_results \
  --mode dynamic \
  --cp12-combined-hmm-map CP12_family_alignment.csv \
  --med4-cp12-accession WP_011131859.1 \
  --hm2-cp12-accession YP_004323596.1
```

---

## Command-Line Arguments

### Input/Output

| Argument | Default | Description |
|----------|---------|-------------|
| `--analysis-dir` | `data_dir/0000_0003_analysis` | Root directory containing trajectory data |
| `--out-dir` | `data_dir/0000_0003_analysis/complex_sasa_bsa_results` | Output directory for results |

### Trajectory Selection

| Argument | Default | Description |
|----------|---------|-------------|
| `--strains` | `med4 hm2` | Protein strains/hosts to process |
| `--systems` | `0000 0003` | Condition codes (redox states) |
| `--replicates` | `rep1 rep2 rep3` | Replicate directory names |

### File Names

| Argument | Default | Description |
|----------|---------|-------------|
| `--pdb-name` | `solute_fixed.pdb` | Reference PDB filename |
| `--tpr-name` | `solute.tpr` | GROMACS topology filename |
| `--xtc-name` | `solute_fit.xtc` | Trajectory filename |

### Analysis Mode

| Argument | Choices | Default | Description |
|----------|---------|---------|-------------|
| `--mode` | `static`, `dynamic`, `both` | `dynamic` | Analysis type (single frame vs trajectory) |
| `--static-source` | `pdb`, `trajectory` | `pdb` | Coordinate source for static mode |
| `--static-time-ps` | float | `0.0` | Frame time for `--static-source trajectory` |

### Dynamic (Trajectory) Analysis

| Argument | Default | Description |
|----------|---------|-------------|
| `--begin-ps` | None | Start time in picoseconds (skip equilibration) |
| `--end-ps` | None | End time in picoseconds |
| `--dt-ps` | None | Sampling interval in ps (e.g., 1000 = sample every ~10 frames if dt=100fs) |

### GROMACS Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--probe-radius` | `0.14` | Probe radius in nm (0.14 nm = 1.4 Å) |
| `--gmx-bin` | `gmx` | GROMACS executable path |

### HMM Alignment Mapping (Optional)

| Argument | Default | Description |
|----------|---------|-------------|
| `--med4-cp12-hmm-map` | None | MED4 CP12 residue-to-HMM mapping CSV |
| `--hm2-cp12-hmm-map` | None | HM2 CP12 residue-to-HMM mapping CSV |
| `--cp12-combined-hmm-map` | None | Family-wide CP12 residue-to-HMM mapping (preferred) |
| `--med4-cp12-accession` | `WP_011131859.1` | MED4 CP12 protein accession ID |
| `--hm2-cp12-accession` | `YP_004323596.1` | HM2 CP12 protein accession ID |

### Execution Control

| Argument | Description |
|----------|-------------|
| `--skip-existing` | Reuse existing XVG files instead of rerunning `gmx sasa` |
| `--dry-run` | Write index files and commands, but do not execute `gmx sasa` |
| `--skip-residue-sasa` | Skip atom-level per-residue SASA computation (faster) |
| `--skip-empty-xtc` | Skip missing or zero-byte trajectory files (default: True) |

---

## Output Files

### GROMACS XVG Files

Raw output from `gmx sasa` in GROMACS XVG format:

```
{out_dir}/{strain}_{system}_{replicate}_*_static.xvg    (static analysis)
{out_dir}/{strain}_{system}_{replicate}_*_dynamic.xvg   (dynamic analysis)
```

**Naming convention:** `{metric}` corresponds to different SASA calculations:
- `protein_complex_sasa_A2`: Protein complex (CP12+PRK+GAP2) SASA
- `full_complex_with_nad_sasa_A2`: Full complex (includes NAD) SASA
- `cp12_isolated_sasa_A2`: Isolated CP12 SASA
- `prk_isolated_sasa_A2`: Isolated PRK SASA
- `gap2_isolated_sasa_A2`: Isolated GAP2 SASA
- `cp12_prk_pair_sasa_A2`: CP12+PRK complex SASA
- `cp12_gap2_pair_sasa_A2`: CP12+GAP2 complex SASA
- `cp12_prk_gap2_pair_sasa_A2`: CP12+PRK+GAP2 complex SASA

### Computed Results (CSV)

Parsed and calculated results in CSV format:

```
{out_dir}/
├── sasa_summary.csv                    # Aggregated SASA results
├── bsa_summary.csv                     # Buried surface area calculations
└── (optional) residue_sasa/
    ├── {strain}_{system}_{replicate}_cp12_per_residue_sasa.csv
    ├── {strain}_{system}_{replicate}_prk_per_residue_sasa.csv
    └── {strain}_{system}_{replicate}_gap2_per_residue_sasa.csv
```

### GROMACS Index Files

Chain-based index files (one per replicate):

```
{out_dir}/{strain}_{system}_{replicate}.ndx
```

Contains groups:
- `CP12`, `PRK`, `GAP2`, `NAD` (individual components)
- `ProteinComplex` (CP12+PRK+GAP2)
- `FullComplex` (CP12+PRK+GAP2+NAD)
- `CP12_PRK`, `CP12_GAP2`, `PRK_GAP2` (pair groups)
- `CP12_PRK_GAP2` (all three proteins)

---

## Output CSV Schemas

### `sasa_summary.csv`

Aggregated SASA data across all replicates and analysis modes.

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| `strain` | str | Host strain (med4 or hm2) |
| `system` | str | Condition code (0000=reduced, 0003=oxidized) |
| `condition` | str | Human-readable condition (reduced or oxidized) |
| `replicate` | str | Replicate identifier (rep1, rep2, rep3) |
| `analysis` | str | Analysis type (static or dynamic) |
| `metric` | str | SASA metric (e.g., cp12_isolated_sasa_A2) |
| `time_ps` | float | Time in picoseconds (dynamic only, NaN for static) |
| `sasa_nm2` | float | Solvent-accessible surface area in nm² |

**Example rows:**
```
strain,system,condition,replicate,analysis,metric,time_ps,sasa_nm2
med4,0000,reduced,rep1,dynamic,cp12_isolated_sasa_A2,0.0,45.23
med4,0000,reduced,rep1,dynamic,cp12_isolated_sasa_A2,1000.0,45.15
med4,0000,reduced,rep1,dynamic,cp12_prk_pair_sasa_A2,0.0,78.42
```

### `bsa_summary.csv`

Buried surface area calculated from isolated and complex SASA values.

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| `strain` | str | Host strain |
| `system` | str | Condition code |
| `condition` | str | Human-readable condition |
| `replicate` | str | Replicate identifier |
| `analysis` | str | Analysis type (static or dynamic) |
| `bsa_type` | str | BSA pair type (cp12_prk, cp12_gap2, or cp12_prk_gap2) |
| `time_ps` | float | Time in picoseconds |
| `bsa_nm2` | float | Buried surface area in nm² |

**Example rows:**
```
strain,system,condition,replicate,analysis,bsa_type,time_ps,bsa_nm2
med4,0000,reduced,rep1,dynamic,bsa_cp12_prk_A2,0.0,12.34
med4,0000,reduced,rep1,dynamic,bsa_cp12_prk_A2,1000.0,12.41
hm2,0003,oxidized,rep2,static,bsa_cp12_gap2_A2,,15.67
```

### Per-Residue SASA (`*_per_residue_sasa.csv`)

Atom-level SASA aggregated to per-residue values (generated with `--skip-residue-sasa` not set).

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| `resid` | int | Residue ID |
| `resname` | str | Residue name (3-letter code) |
| `chain` | str | Chain identifier |
| `time_ps` | float | Time in picoseconds (dynamic) |
| `sasa_nm2` | float | Per-residue solvent-accessible surface area in nm² |

---

## Analysis Pipeline

```
Trajectory Files (solute_fixed.pdb, solute.tpr, solute_fit.xtc)
    ↓
Parse PDB → Extract chain IDs
    ↓
Create GROMACS index file (.ndx) with protein groups
    ↓
Run gmx sasa for each configuration:
    • Isolated components (CP12, PRK, GAP2)
    • Pair complexes (CP12+PRK, CP12+GAP2)
    • Full complex
    ↓
Parse XVG output → Calculate average SASA
    ↓
Compute BSA using: BSA(A,B) = SASA(A) + SASA(B) - SASA(A+B)
    ↓
Aggregate across replicates
    ↓
Output CSVs: sasa_summary.csv, bsa_summary.csv
```

---

## Technical Details

### SASA Calculation

- **Probe radius**: 0.14 nm (1.4 Å, standard for protein SASA)
- **Algorithm**: GROMACS uses the Shrake-Rupley algorithm
- **Output units**: Ångströms² (converted to nm² for output)

### Buried Surface Area (BSA)

BSA quantifies the surface area that becomes "buried" (inaccessible to solvent) upon complex formation:

```
BSA(A,B) = SASA(A_isolated) + SASA(B_isolated) - SASA(A+B_complex)
```

**Interpretation:**
- **High BSA**: Significant burial → tight binding interface
- **Low BSA**: Minimal burial → weak or transient interaction
- **Positive BSA**: Expected for stable complexes
- **Negative BSA**: Indicates unusual solvation; may signal dynamic dissociation

### Time Averaging

For dynamic analysis, the script:
1. Reads trajectory and extracts all frames in the specified time range
2. Computes SASA for each frame independently
3. Outputs per-frame SASA values and frame times

External analysis can then compute frame-averaged BSA.

### Static vs. Dynamic

- **Static mode**: 
  - Fast single-point calculation
  - Useful for quick assessments or structure validation
  - `--static-source pdb`: uses solute_fixed.pdb directly
  - `--static-source trajectory`: extracts one frame from solute_fit.xtc (avoids mass guessing warnings)

- **Dynamic mode**:
  - Full trajectory analysis
  - Time-dependent SASA values reveal fluctuations
  - Enables statistical comparison across conditions

---

## Common Use Cases

### 1. Quick Static Comparison

Compare interface burial in oxidized vs. reduced CP12:

```bash
python complex_sasa_gmx.py \
  --analysis-dir data \
  --out-dir results \
  --mode static
```

Then examine `bsa_summary.csv` for system/condition differences.

### 2. Full Trajectory Analysis with Equilibration Skip

Compute time-resolved BSA, skipping first 10 ns:

```bash
python complex_sasa_gmx.py \
  --analysis-dir data \
  --out-dir results \
  --mode dynamic \
  --begin-ps 10000 \
  --dt-ps 1000
```

### 3. Per-Residue Interface Analysis

Identify which CP12 residues contribute most to interface burial:

```bash
python complex_sasa_gmx.py \
  --analysis-dir data \
  --out-dir results \
  --mode dynamic \
  --skip-residue-sasa=False
```

Then analyze `*_per_residue_sasa.csv` files.

### 4. Cross-Strain Residue Mapping

Enable HMM-aligned residue identifiers for comparing MED4 vs HM2:

```bash
python complex_sasa_gmx.py \
  --analysis-dir data \
  --out-dir results \
  --mode dynamic \
  --cp12-combined-hmm-map CP12_family_alignment.csv
```

---

## Troubleshooting

### GROMACS Not Found

**Error:** `FileNotFoundError: gmx: command not found`

**Solution:** Install GROMACS or specify its path:
```bash
python complex_sasa_gmx.py --gmx-bin /path/to/gmx ...
```

Or use a GROMACS module:
```bash
module load gromacs
python complex_sasa_gmx.py ...
```

### Missing Trajectory Files

**Error:** `FileNotFoundError: /path/to/solute.tpr`

**Solution:** Verify directory structure matches expected layout. Check:
- All three files exist: `solute_fixed.pdb`, `solute.tpr`, `solute_fit.xtc`
- Directory names match `--strains`, `--systems`, `--replicates` arguments

### Chain ID Mismatch

**Error:** `ValueError: Missing required chain-derived groups`

**Solution:** Verify `solute_fixed.pdb` has explicit chain IDs. Expected chains:
- CP12: C, H, M, P
- PRK: B, G, J, L
- GAP2: A, D, E, F, I, K, N, O

Check PDB file column 22 (chain ID) contains these values.

### Empty or Missing Index Groups

**Error:** `gmx sasa` fails with "group not found"

**Solution:** Check the generated `.ndx` file (e.g., `med4_0000_rep1.ndx`). Groups should list atom indices. If groups are empty, verify PDB parsing.

### No Output Files

**Use `--dry-run` to debug:**
```bash
python complex_sasa_gmx.py --dry-run --mode dynamic > commands.sh
```

This writes expected commands without executing them. Review for errors.

---

## Performance Tips

1. **Use `--skip-residue-sasa`** to speed up analysis (skips atom-level aggregation)
2. **Use `--dt-ps`** to reduce frame sampling (e.g., `--dt-ps 1000` samples every ~10 frames)
3. **Use `--skip-existing`** to reuse XVG outputs if changing only analysis parameters
4. **Run in parallel**: Process multiple replicates/strains independently and combine results

---

## References

- GROMACS SASA documentation: https://manual.gromacs.org/
- Shrake-Rupley algorithm: Shrake & Rupley (1973), J. Mol. Biol. 79(2):351-364
