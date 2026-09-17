# Static Interface HMM Mapping

Maps CP12 protein–partner interface residues from static PDB structures onto columns of the PF02672 CP12 family HMM profile, enabling cross-species comparison of interface conservation.

## Overview

The analysis is split into four sequential phases:

| Phase | Script | Description |
|-------|--------|-------------|
| 1 | `phase1_chain_map.py` | Map PDB chains to named proteins |
| 2 | `phase2_extract_contacts.py` | Extract heavy-atom residue contacts at a distance cutoff |
| 3 | `phase3_collapse_cp12_interfaces.py` | Collapse contacts to per-residue interface annotations |
| 4a | `phase4_cp12_resid_to_hmm.py` | Align CP12 chains to HMM, map residues to HMM columns |
| 4b | `phase4_project_and_collapse_cp12_interfaces.py` | Project interfaces onto HMM columns and collapse across chains |


## Input Data

### Required files

| File | Description |
|------|-------------|
| `<system>_solute_fixed.pdb` | Solute-only PDB of the complex with chain IDs assigned. Hydrogens may be present; they are filtered automatically. Residue records must use standard `ATOM` records (HETATM residues are ignored). |
| `PF02672.hmm` | HMMER profile HMM for the CP12 family (Pfam accession PF02672). Download from [Pfam](https://pfam.xfam.org/) or [InterPro](https://www.ebi.ac.uk/interpro/). |

### PDB chain assignment convention

Each chain in the PDB must belong to exactly one named protein. The mapping is provided on the command line in `Phase 1` using `--protein-map`:

```
CP12:C,H,M,P  PRK:A,B  GAPDH:D,E,F
```

CP12 is the focal protein; all other proteins listed in `--partner-proteins` during Phase 2 are treated as binding partners.

## Usage

Run each phase in order. All paths below assume you are working from the `scripts/` directory or provide full paths to inputs/outputs.

### Phase 1 — Chain map

```bash
python phase1_chain_map.py \
    --pdb      ../data/MED4_0003_rep2_solute_fixed.pdb \
    --system   MED4_0003 \
    --protein-map  CP12:C,H,M,P  PRK:A,B  GAPDH:D,E,F \
    --out      MED4_0003_chain_map.tsv
```

**Output:** `MED4_0003_chain_map.tsv`

| Column | Description |
|--------|-------------|
| System | System name |
| Chain | PDB chain ID |
| Protein | Assigned protein name |
| ResidueCount | Number of standard residues in the chain |

---

### Phase 2 — Extract contacts

```bash
python phase2_extract_contacts.py \
    --pdb             ../data/MED4_0003_rep2_solute_fixed.pdb \
    --chain-map       MED4_0003_chain_map.tsv \
    --system          MED4_0003 \
    --focal-protein   CP12 \
    --partner-proteins  PRK GAPDH \
    --distance-cutoff 5.0 \
    --out             MED4_0003_phase2_contacts.tsv
```

`--distance-cutoff` is in Ångströms (default 5.0 Å). Only heavy atoms are considered.

**Output:** `MED4_0003_phase2_contacts.tsv` — one row per CP12 residue–partner residue contact pair.

---

### Phase 3 — Collapse to interface annotations

```bash
python phase3_collapse_cp12_interfaces.py \
    --contacts       MED4_0003_phase2_contacts.tsv \
    --system         MED4_0003 \
    --focal-protein  CP12 \
    --out-prefix     MED4_0003_phase3_cp12_interface
```

**Output:** `MED4_0003_phase3_cp12_interface.tsv`

| Column | Description |
|--------|-------------|
| System | System name |
| CP12_Chain | CP12 chain ID |
| CP12_ResID | PDB residue number |
| CP12_ResName | Three-letter residue name |
| InterfaceType | Partner(s) in contact (e.g. `PRK_only`, `PRK+GAPDH`) |
| PartnerProteins | Comma-separated partner protein names |
| MinDistance | Minimum heavy-atom distance to any partner residue (Å) |
| NumContacts | Number of partner residues in contact |
| IsInterface | Always `True` (non-interface residues are not written) |

---

### Phase 4a — Map CP12 residues to HMM columns

```bash
python phase4_cp12_resid_to_hmm.py \
    --pdb        ../data/MED4_0003_rep2_solute_fixed.pdb \
    --chains     C,H,M,P \
    --hmm        ../data/PF02672.hmm \
    --system     MED4_0003 \
    --out-prefix MED4_0003_phase4a
```

Requires `hmmalign` on `PATH`. The `--hmmalign-bin` flag can override the executable path.

**Output files:**

| File | Description |
|------|-------------|
| `MED4_0003_phase4a_cp12_chains.faa` | FASTA sequences extracted from each CP12 chain |
| `MED4_0003_phase4a_cp12_hmmalign.sto` | Stockholm-format alignment to PF02672 |
| `MED4_0003_phase4a_cp12_resid_to_hmmcol.tsv` | Per-residue mapping to HMM column |

Columns in `_cp12_resid_to_hmmcol.tsv`:

| Column | Description |
|--------|-------------|
| System | System name |
| CP12_Chain | Chain ID |
| ResID | PDB residue number |
| SeqPos | 1-based position in extracted chain sequence |
| AA | One-letter amino acid code |
| HMM_col | Corresponding PF02672 HMM column number (`NA` if outside the model) |

---

### Phase 4b — Project interfaces to HMM columns and collapse

```bash
python phase4_project_and_collapse_cp12_interfaces.py \
    --interfaces  MED4_0003_phase3_cp12_interface.tsv \
    --resid-map   MED4_0003_phase4a_cp12_resid_to_hmmcol.tsv \
    --system      MED4_0003 \
    --out-prefix  MED4_0003_phase4b
```

**Output files:**

| File | Description |
|------|-------------|
| `MED4_0003_phase4b_cp12_interface_hmm_chain.tsv` | Chain-resolved interface mapped to HMM columns |
| `MED4_0003_phase4b_cp12_interface_hmm_collapsed.tsv` | Per-HMM-column interface, collapsed across chains |

Columns in `_cp12_interface_hmm_collapsed.tsv`:

| Column | Description |
|--------|-------------|
| System | System name |
| HMM_col | PF02672 HMM column |
| Interface_present | 1 if any chain at this HMM column is in contact |
| InterfaceType_collapsed | Interface type; `Mixed` if chains disagree |
| ChainsSupporting | Comma-separated CP12 chains in contact at this HMM column |
| NumChainsSupporting | Count of supporting chains |

---

## Data Directory

```
data/
├── MED4_0003_rep2_solute_fixed.pdb          # Example input PDB (MED4 organism)
├── HM2_0003_rep2_solute_fixed.pdb           # Example input PDB (HM2 organism)
├── PF02672.hmm                              # CP12 family HMM profile
└── MED4_vs_HM2_cp12_interface_hmm_static_rep2_5A.tsv  # Final comparative output
```

### Comparative output format

`MED4_vs_HM2_cp12_interface_hmm_static_rep2_5A.tsv` is produced by running the full pipeline on two systems and joining the collapsed Phase 4b outputs on `HMM_col`. The `ChangeClass` column categorises each HMM position as:

| ChangeClass | Meaning |
|-------------|---------|
| `Conserved` | Interface present in both systems |
| `Lost_in_<system>` | Interface present only in the reference system |
| `Gained_in_<system>` | Interface present only in the comparison system |

## Notes

- PDB files should contain only the solute (no solvent or ions) to avoid spurious contacts.
- Multiple CP12 chains (copies) in the complex are processed independently in Phases 1–4a and collapsed in Phase 4b. `NumChainsSupporting` reflects how many copies of CP12 show the interface at a given HMM position.
- Residues that do not align within the HMM profile (insertions relative to the consensus) are assigned `HMM_col = NA` and excluded from Phase 4b outputs.
- The `_5A` suffix in the example output filename indicates the 5 Å heavy-atom distance cutoff used in Phase 2.
