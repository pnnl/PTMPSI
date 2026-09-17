# CP12 Protein Family HMM Column Metrics

Computes per-HMM-column sequence-diversity metrics across the CP12 protein family and uses those features to classify sequences by phylogenetic domain (Bacteria vs. Viruses) with logistic regression.

The pipeline proceeds in four phases:

1. **Architecture features** — derive per-sequence cysteine and acidity features from a metadata table.
2. **Sequence grouping** — assign each sequence to a `PrimaryGroup` and `PrimaryGroupCombined` category for stratified counting.
3. **Per-column amino acid counts** — parse a Stockholm HMM alignment and count amino acids at each HMM match column, stratified by group.
4. **Per-column metrics & classification** — compute Shannon entropy, cysteine fraction, acidic fraction, occupancy, and uniqueness per column; integrate interface annotations; train and cross-validate a logistic regression classifier.

---

## Directory Layout

```
CP12_protein_family_HMM_column_metrics/
├── data/
│   ├── CP12_MED4_vs_HM2_alignment_annotation.csv          # pairwise alignment of MED4 vs HM2 reference sequences (input)
│   ├── HMM_seq_align.fasta                                # HMM-aligned FASTA (reference)
│   ├── cp12_seq_groups.txt                                # sequence group assignments (pipeline output)
│   ├── phase5_cp12.step1_architecture_v2.txt              # per-sequence architecture features (pipeline output)
│   ├── figure2_per_column_aa_counts_v2.txt                # per-column AA counts, all groups (pipeline output)
│   ├── figure2_per_column_stats_long_static_rep2_5A.txt   # per-column diversity metrics (pipeline output)
│   ├── figure2_per_column_redo_get_more_per_lineage_stats_long_step2v6.txt        # per-column metrics v6, by lineage (pipeline output)
│   ├── figure2_per_column_redo_get_more_per_lineage_stats_long_step2v6_with_iface_static_rep2_5A_0526_FINAL.txt  # metrics with interface annotation (pipeline output)
│   ├── cp12_hmm_diversity_interface_static_rep2_5A.tsv    # column diversity + interface annotation (pipeline output)
│   ├── combo1_position_fisher_fdr_with_coef_rank_outNEW_0526_FINAL.tsv  # per-position Fisher test + classifier coefficients (pipeline output)
│   ├── combo1_featureclass_fisher_fdr_0526_FINAL.tsv      # per-feature-class Fisher FDR results (pipeline output)
│   ├── combo1_logreg_coefficients.tsv                     # logistic regression coefficients (pipeline output)
│   ├── combo1_cv_predictions.tsv                          # cross-validation per-sequence predictions (pipeline output)
│   └── combo1_confusion_matrix_cv.tsv                     # cross-validation confusion matrix (pipeline output)
└── scripts/
    ├── phase5_step1_cp12_architecture_features_v2.py      # Step 1: derive architecture features
    ├── phase5_step2_make_seq_groups.py                    # Step 2a: assign sequence groups
    ├── phase5_step2_per_column_aa_counts_v2.py            # Step 2b: count AAs per HMM column
    ├── phase5_step2_derive_per_column_metrics.py          # Step 2c: compute diversity metrics
    ├── phase5_step2_per_column_generateOldToNewColumnMapping.py  # utility: map residue index → HMM column
    ├── phase5_cp12_diversity_with_interfaces.py           # Step 3: merge diversity with interface annotation
    ├── redo_step2_cp12_combo1_v6.py                       # Step 4a: build per-accession long table (v6)
    └── combo1_classifier_cp12_v2.py                       # Step 4b: logistic regression classifier
```


---

## Required Input Data Structure

### 1. Stockholm multiple-sequence alignment (`*.sto`)

An HMM-aligned Stockholm file produced by `hmmalign` (HMMER 3) against the CP12 Pfam profile (PF02672). It **must** contain a `#=GC RF` annotation line marking HMM match columns (`x`) vs. insert columns (`.`).

```
# STOCKHOLM 1.0
#=GF ...
<ACCESSION>   <ALIGNED_SEQUENCE>
...
#=GC RF       xxxx....xxxx
//
```

- Sequence names must match the `Accession` column in the metadata table.
- Wrapped (multi-block) alignments are supported.

### 2. Combined CP12 metadata table (`cp12_combined.metadata_v3.1.txt`)

A tab-separated file produced by the `CP12_protein_family_build` pipeline. Required columns:

| Column | Description |
|--------|-------------|
| `Accession` | Unique sequence identifier matching alignment sequence names |
| `SourceType` | `Reference` or `Metagenome` |
| `SequenceStatus` | `full` or `partial` (metagenome sequences only) |
| `NumCys` | Number of cysteine residues in the sequence |
| `CysPositions` | Python list literal of cysteine positions, e.g. `[12, 21]` |
| `CysSpacings` | Python list literal of spacings between consecutive cysteines, e.g. `[9]` |
| `AcidicFraction` | Fraction of residues that are D or E |
| `Lineage_consensus_cleaned` | Colon-delimited lineage string, e.g. `Bacteria:Cyanobacteria` |

Optional columns used by classifiers:

| Column | Description |
|--------|-------------|
| `TaxID` | NCBI taxonomy ID |
| `PrimaryGroup` | Pre-computed group label (computed by pipeline if absent) |
| `PrimaryGroupCombined` | Coarser group label (computed by pipeline if absent) |
| `Lineage_consensus_cleaned` | Used to derive `Domain` (`Bacteria` / `Viruses`) for classification |

### 3. Interface annotation table

A tab-separated file that maps HMM column numbers to interface residue annotations derived from the MED4 vs. HM2 structural comparison. Required columns:

| Column | Description |
|--------|-------------|
| `HMM_col` | 1-based HMM match column number |
| `MED4_<ID>_Interface` | Boolean/flag for MED4 interface residue |
| `MED4_<ID>_Type` | Interaction type for MED4 residue |
| `HM2_<ID>_Interface` | Boolean/flag for HM2 interface residue |
| `HM2_<ID>_Type` | Interaction type for HM2 residue |
| `ChangeClass` | Class of interface change between MED4 and HM2 |

The file `data/CP12_MED4_vs_HM2_alignment_annotation.csv` provides the pairwise residue mapping used to construct this table.

---

## Usage

All scripts accept command-line arguments. Run any script with `--help` for the full option list.

### Step 1 — Architecture features

Derives per-sequence cysteine pair and acidity features.

```bash
python scripts/phase5_step1_cp12_architecture_features_v2.py \
    --meta  cp12_combined.metadata_v3.1.txt \
    --out   data/phase5_cp12.step1_architecture_v2.txt
```

**Output columns:** `Accession`, `PrimaryGroup`, `SourceType`, `SequenceStatus`, `NumCys`, `NumCysPairs`, `MeanCysSpacing`, `AcidicFraction`

---

### Step 2a — Sequence groups

Assigns `PrimaryGroup`, `PrimaryGroupCombined`, and `Lineage_consensus_cleaned` labels used for stratified counting.

```bash
python scripts/phase5_step2_make_seq_groups.py \
    --meta  cp12_combined.metadata_v3.1.txt \
    --out   data/cp12_seq_groups.txt
```

**Output columns:** `Accession`, `PrimaryGroup`, `PrimaryGroupCombined`, `SourceType`, `SequenceStatus`, `Lineage_consensus_cleaned`

Group definitions:

| `PrimaryGroup` | Criteria |
|----------------|----------|
| `Reference` | `SourceType == "Reference"` |
| `Metagenome_full` | `SourceType == "Metagenome"` and `SequenceStatus == "full"` |
| `Metagenome_partial` | `SourceType == "Metagenome"` and `SequenceStatus == "partial"` |

---

### Step 2b — Per-column amino acid counts

Parses the Stockholm alignment and produces a long-format count table stratified by `PrimaryGroup`, `PrimaryGroupCombined`, and `Lineage`.

```bash
python scripts/phase5_step2_per_column_aa_counts_v2.py \
    --sto   cp12_combined_pf02672.sto \
    --meta  data/cp12_seq_groups.txt \
    --out   data/figure2_per_column_aa_counts_v2.txt
```

**Output columns:** `HMM_col`, `GroupType`, `GroupValue`, `NumSeq`, `NumNonGap`, `Gap_count`, `A_count` … `Y_count` (one column per standard amino acid)

---

### Step 2c — Per-column diversity metrics

Computes Shannon entropy, cysteine fraction, acidic fraction, occupancy, and amino acid uniqueness per column per group, and merges interface annotations.

```bash
python scripts/phase5_step2_derive_per_column_metrics.py \
    --counts      data/figure2_per_column_aa_counts_v2.txt \
    --interfaces  <interface_annotation.tsv> \
    --out         data/figure2_per_column_stats_long_static_rep2_5A.txt \
    --cys_thresh  0.7
```

The `--cys_thresh` flag sets the minimum cysteine fraction for a column to be flagged as `IsConservedCys` (default `0.7`).

**Output metrics per column/group:** `C_fraction`, `IsConservedCys`, `NumUnique`, `Shannon`, `AcidicFraction`, `Occupancy`

---

### Step 3 — Diversity with interface annotation

Merges raw alignment-based diversity statistics with interface residue annotations.

```bash
python scripts/phase5_cp12_diversity_with_interfaces.py \
    --sto         cp12_combined_pf02672.sto \
    --interfaces  <interface_annotation.tsv> \
    --out         data/cp12_hmm_diversity_interface_static_rep2_5A.tsv
```

---

### Step 4a — Per-accession long table (v6)

Builds a long-format table with one row per accession per HMM column, including residue identity, chemical class, and metadata labels. This is the primary input to the classifier.

```bash
python scripts/redo_step2_cp12_combo1_v6.py \
    --sto   cp12_combined_pf02672.sto \
    --meta  cp12_combined.metadata_v3.1.txt \
    --out_prefix  combo1 \
    --out_dir     data/
```

**Output columns:** `Accession`, `HMM_col`, `AA_raw`, `AA`, `IsGap`, `ChemClass`, `SourceType`, `SequenceStatus`, `PrimaryGroup`, `PrimaryGroupCombined`, `LineageCleaned`, `TaxID`

Chemical class assignments (`ChemClass`):

| Class | Residues |
|-------|----------|
| `Cysteine` | C |
| `Positive` | K, R, H |
| `Negative` | D, E |
| `Polar` | S, T, N, Q |
| `Hydrophobic` | A, V, I, L, M, F, W, Y |
| `Gly_Pro` | G, P |

---

### Step 4b — Logistic regression classifier

Classifies sequences as Bacteria (0) or Viruses (1) using per-HMM-column chemical class and/or amino acid identity as features. Performs stratified k-fold cross-validation.

```bash
python scripts/combo1_classifier_cp12_v2.py \
    --accession_long  data/<combo1_accession_per_column_step2v6.txt> \
    --out_dir         data/ \
    --feature_set     chem \
    --kfold           5 \
    --seed            1
```

`--feature_set` options:

| Value | Features used |
|-------|--------------|
| `chem` | `ChemClass` per HMM column (default) |
| `aa` | raw amino acid identity per HMM column |
| `both` | `ChemClass` and amino acid identity |

Pass `--use_gaps` to retain gap characters as a feature category (default: gap positions are excluded for each sequence).

**Outputs written to `--out_dir`:**
- `combo1_logreg_coefficients.tsv` — per-feature logistic regression coefficients
- `combo1_cv_predictions.tsv` — per-sequence cross-validation predicted label and probability
- `combo1_confusion_matrix_cv.tsv` — aggregated cross-validation confusion matrix
- `combo1_position_fisher_fdr_with_coef_rank_outNEW_0526_FINAL.tsv` — per-position Fisher exact test with FDR correction and coefficient rank
- `combo1_featureclass_fisher_fdr_0526_FINAL.tsv` — per-feature-class Fisher test results

---

### Utility — Residue-to-HMM-column mapping

Exports a per-sequence mapping of residue index to HMM match column number for downstream analyses (e.g., contact analysis, RMSF mapping).

```bash
python scripts/phase5_step2_per_column_generateOldToNewColumnMapping.py \
    --sto   cp12_combined_pf02672.sto \
    --meta  cp12_combined.metadata_v3.1.txt \
    --out   cp12_residue_to_hmmcol_mapping.tsv
```

---

## Data Flow Summary

```
cp12_combined.metadata_v3.1.txt  ─────────────────────────────────────────────────────────────┐
                                                                                               │
cp12_combined_pf02672.sto ──► phase5_step2_per_column_aa_counts_v2.py                         │
                          ──► phase5_step2_per_column_generateOldToNewColumnMapping.py         │
                          ──► phase5_cp12_diversity_with_interfaces.py                         │
                          ──► redo_step2_cp12_combo1_v6.py ──► combo1_classifier_cp12_v2.py   │
                                                                                               │
interface_annotation.tsv ───► phase5_step2_derive_per_column_metrics.py                       │
                          ──► phase5_cp12_diversity_with_interfaces.py                         │
                                                                                               │
                ┌──────────────────────────────────────────────────────────────────────────────┘
                ▼
        phase5_step1_cp12_architecture_features_v2.py  →  architecture feature table
        phase5_step2_make_seq_groups.py                →  sequence group assignments
        figure2_per_column_aa_counts_v2.txt            →  AA count table
        figure2_per_column_stats_long_*.txt            →  diversity metrics
        combo1_*                                       →  classification results
```

---

## Notes

- All HMM column numbers are **1-based** and correspond to match-state positions in the PF02672 profile.
- Sequences not found in the metadata table are silently skipped during alignment parsing; check script stderr for counts.
- The classifier targets only sequences with `Domain` equal to `Bacteria` or `Viruses` (derived from `Lineage_consensus_cleaned`); other lineages are excluded.
- This analysis depends on outputs from `CP12_protein_family_build/`; run that pipeline first to generate `cp12_combined.metadata_v3.1.txt` and `cp12_combined_pf02672.sto`.
