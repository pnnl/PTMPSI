# CP12 Protein Family Build

Builds a combined CP12 protein family database from two sources:

1. **Metagenome-derived sequences** — environmental metagenomes annotated with KofamScan, Pfam (hmmscan), and CAT taxonomic classification.
2. **Public reference sequences** — NCBI RefSeq and reviewed UniProt entries carrying the PF02672 (CP12) Pfam domain.

The final outputs are a harmonized metadata table (`cp12_combined.metadata_v3.1.txt`) and a Stockholm multiple-sequence alignment (`cp12_combined_pf02672.sto`) used by downstream HMM-based analyses.

---

## Directory Layout

```
CP12_protein_family_build/
├── data/
│   ├── PF02672.hmm                        # Pfam CP12 HMM profile (required input)
│   ├── CP12_motif_and_classification.tsv  # CP12 motif/class lookup (required input)
│   ├── combined_sequences.fasta           # metagenome full + partial sequences (pipeline output)
│   ├── combined_annotation.tsv            # merged metagenome annotation (pipeline output)
│   ├── combined_sequences_derep.fasta     # 100%-identity dereplicated sequences (pipeline output)
│   ├── combined_annotation_derep.tsv      # annotation for dereplicated set (pipeline output)
│   ├── combined_annotation_derep_cp12.tsv # CP12-classified annotation (pipeline output)
│   ├── cp12_combined.metadata_v3.1.txt    # final combined metadata v3.1 (pipeline output)
│   └── cp12_combined_pf02672.sto          # Stockholm alignment against PF02672 (pipeline output)
└── scripts/
    ├── megahit_contigs_filtered_1kb_prot_kofam_cleaned.sb     # SLURM: per-sample KofamScan
    ├── megahit_contigs_filtered_1kb_prot_pfam_cleaned.sb      # SLURM: per-sample Pfam hmmscan
    ├── megahit_1kb_contig_CAT_v2_NRaddNames.sb                # SLURM: per-contig CAT taxonomy (GTDB + NR)
    ├── NR_withName_formatConsistent_gtdb.py                   # convert NR taxonomy to GTDB-style
    ├── batch_super_merge_v3.sh                                # batch: merge all per-sample annotations
    ├── super_merge_v10.py                                     # per-sample annotation merge logic
    ├── batch_filters_function_taxonomy_v4.sh                  # batch: KO/PF + taxonomy filter
    ├── function_filter.py                                     # filter by KO / PF accession
    ├── taxonomy_filter_v7.py                                  # filter by taxonomic category
    ├── batch_extract_fasta_by_new_gene_v6.sh                  # batch: FASTA extraction + master tables
    ├── extract_by_id_split_partial_full_v2.py                 # extract full/partial sequences by ID
    ├── public_seq_fetching_v6.sh                              # fetch CP12 references from NCBI + UniProt
    ├── pipeline_v3.sh                                         # combine, dereplicate, and classify metagenome hits
    ├── phase11_build_cp12_combined.py                         # final merge of metagenome + references
    └── curated_exclude_ids.example.txt                        # example accession exclusion list
```

---

## External Dependencies

| Tool | Purpose |
|------|---------|
| Python ≥ 3.8 | all Python scripts |
| KofamScan (`exec_annotation`) | per-sample KO annotation |
| HMMER 3 (`hmmscan`, `hmmbuild`, `hmmalign`) | Pfam annotation and final alignment |
| CD-HIT (`cd-hit`) | 100%-identity dereplication |
| CAT_pack | contig-level taxonomic classification |
| NCBI Entrez Direct (`esearch`, `efetch`, `xtract`) | public reference fetching |
| `curl` | UniProt REST API calls |
| geNomad *(optional)* | per-contig viral classification |

Set `ENTREZ_EMAIL` in your environment before running `public_seq_fetching_v6.sh`.

---

## Required Input Data Structure

### Metagenome branch

All per-sample files are expected to follow a consistent naming convention derived from the sample ID. The sample ID is the hyphen-delimited token matching the pattern `<DATE>-PNNL-<TYPE>-Seq<N>` (e.g., `20241106-PNNL-NRB-Seq149`).

#### Per-sample protein FASTA
```
renamed_prodigal_out_<SAMPLE>_megahit_contigs_filtered_1kb_prot.fasta[.gz]
```
Amino-acid sequences produced by Prodigal from ≥1 kb assembled contigs.

#### Per-sample KofamScan output
```
<SAMPLE>_megahit_contigs_filtered_1kb_prot_kofamscan_out
```
Produced by:
```bash
exec_annotation -f detail-tsv -E 0.001 --cpu <N> -o <out> <protein.fasta>
```

#### Per-sample Pfam output (tblout)
```
<SAMPLE>_megahit_contigs_filtered_1kb_prot_pfam_out
```
Produced by:
```bash
hmmscan --tblout <out> --cut_tc --cpu <N> Pfam-A.hmm <protein.fasta>
```

#### Per-contig CAT GTDB lineage file
Tab-delimited file with at least a contig ID column and a `lineage` column in GTDB-style format (`d__;p__;c__;o__;f__;g__;s__`), produced by CAT_pack against the GTDB database.

#### Per-contig CAT NR lineage file
Same structure as the GTDB file but from CAT_pack against the NR database. Convert to GTDB-style prefixes before merging using:
```bash
python scripts/NR_withName_formatConsistent_gtdb.py <CAT_NR_*_withNames.txt>
```

#### Per-contig geNomad viral summary *(optional)*
geNomad `*_virus_summary.tsv` files providing per-contig viral taxonomy. Absence of this file causes the `virus:*` columns to be populated with `NA`.

#### Sample list file
Plain text file (one sample ID per line) passed to the batch scripts:
```
20240418-PNNL-NRB-Seq11
20240718-PNNL-NRB-Seq58
20240805-PNNL-MPF-Seq79
...
```

### Public reference branch

No pre-existing files are required. `public_seq_fetching_v6.sh` fetches everything from NCBI and UniProt and writes its outputs to a `cp12_public_refs/` subdirectory. Manual RefSeq accessions and length thresholds are configured at the top of the script.

---

## Pipeline Execution

The pipeline has six ordered phases. Phases 1–5 (metagenome) and Phase 5 (public refs) can run concurrently; Phase 6 requires both to be complete.

### Phase 1 — Per-sample annotation (HPC / SLURM array jobs)

Run on your HPC cluster as SLURM array jobs. Edit the path variables at the top of each `.sb` script to match your installation.

```bash
# KofamScan (one job per sample)
sbatch --array=1-<N> scripts/megahit_contigs_filtered_1kb_prot_kofam_cleaned.sb

# Pfam hmmscan (one job per sample)
sbatch --array=1-<N> scripts/megahit_contigs_filtered_1kb_prot_pfam_cleaned.sb

# CAT contig taxonomy — GTDB and NR (one job per sample)
sbatch --array=1-<N> scripts/megahit_1kb_contig_CAT_v2_NRaddNames.sb

# Convert NR taxonomy to GTDB-style lineage strings
python scripts/NR_withName_formatConsistent_gtdb.py <CAT_NR_*_withNames.txt>
```

### Phase 2 — Per-sample annotation merge

Merges KofamScan, Pfam, GTDB, NR, and (optionally) geNomad annotations into one TSV per sample.

Edit the directory path variables at the top of `batch_super_merge_v3.sh`, then:

```bash
bash scripts/batch_super_merge_v3.sh samples.txt [outer|inner] [outdir]
```

Output per sample: `<SAMPLE>-merged_kofam_pfam_union_gtdb_nr_viral.txt`

### Phase 3 — Function and taxonomy filtering

Filter merged tables to CP12-relevant genes. Use `--pf PF02672` to target the CP12 Pfam domain. Default taxonomy targets are Viruses, Bacteria, Archaea, Fungi, Viridiplantae, Metazoa, and SAR.

```bash
bash scripts/batch_filters_function_taxonomy_v4.sh samples.txt \
    --merged-dir <path/to/merged_dir> \
    --pf PF02672 \
    --ko NA \
    --targets "Viruses,Bacteria,Archaea,Fungi,Viridiplantae,Metazoa,SAR" \
    --taxa-mode lenient \
    --outdir filter_results/
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--merged-dir` | hardcoded path | directory containing per-sample merged TSVs |
| `--ko` | `NA` | KO accession to filter on (`NA` = skip KO filter) |
| `--pf` | `NA` | PF accession to filter on (e.g., `PF02672`) |
| `--targets` | all major domains | comma-separated taxonomy categories |
| `--taxa-mode` | `lenient` | `lenient` or `strict` taxonomy matching |
| `--chain` | `yes` | apply taxonomy filter after function filter |
| `--outdir` | `.` | root output directory |

Outputs written to `<outdir>/func_filter/`, `<outdir>/func_then_taxa/`, etc.

### Phase 4 — FASTA extraction and master table assembly

Extract matched sequences from per-sample protein FASTAs, splitting into full-length and partial hits. Assembles master FASTA and annotation tables across all samples.

```bash
bash scripts/batch_extract_fasta_by_new_gene_v6.sh \
    filter_results/ \
    <path/to/protein_fasta_dir/> \
    "func_filter func_then_taxa"
```

Outputs written to the filter results root:
```
master_full_sequences.fasta
master_partial_sequences.fasta
master_full_annotation.tsv
master_partial_annotation.tsv
```

### Phase 5 — Public reference sequences

Fetches CP12 references from NCBI RefSeq (explicit accessions + title search) and reviewed UniProt entries (PF02672). Requires NCBI Entrez Direct and `ENTREZ_EMAIL` to be set.

```bash
export ENTREZ_EMAIL="your@email.edu"
bash scripts/public_seq_fetching_v6.sh
```

Outputs written to `cp12_public_refs/`:
```
cp12_public_refs.clean.fasta
cp12_public_refs.metadata.tsv
```

### Phase 5 (metagenome) — Combination and dereplication

Run from the directory containing `master_full_sequences.fasta` and `master_partial_sequences.fasta` (output of Phase 4). Combines full and partial metagenome sequences, dereplicates at 100% identity with CD-HIT, and produces a CP12-classified annotation table.

```bash
cd <master_tables_dir>
bash /path/to/scripts/pipeline_v3.sh
```

Outputs:
```
combined_sequences.fasta
combined_annotation.tsv
combined_sequences_derep.fasta
combined_annotation_derep.tsv
combined_annotation_derep_cp12.tsv
```

### Phase 6 — Final combined database

Combines metagenome and public reference sequences, harmonizes metadata to the v3.1 schema, and generates the Stockholm alignment.

```bash
python scripts/phase11_build_cp12_combined.py \
    --metagenome-annotation combined_annotation_derep_cp12.tsv \
    --metagenome-fasta      combined_sequences_derep.fasta \
    --reference-metadata    cp12_public_refs/cp12_public_refs.metadata.tsv \
    --reference-fasta       cp12_public_refs/cp12_public_refs.clean.fasta \
    --hmm                   data/PF02672.hmm \
    --outdir                data/
```

Optional flags:

| Flag | Description |
|------|-------------|
| `--exclude-accessions <file>` | text file with one accession per line to drop (curated exclusions); see `curated_exclude_ids.example.txt` for format |
| `--strict-pfam-cp12` | keep metagenome rows only when `pfamsearch:target_accession` starts with `PF02672` |
| `--metadata-name <name>` | override output metadata filename (default: `cp12_combined.metadata_v3.1.txt`) |
| `--combined-fasta-name <name>` | override combined FASTA filename |
| `--sto-name <name>` | override Stockholm alignment filename |

Final outputs:
```
data/cp12_combined.metadata_v3.1.txt
data/cp12_combined_combined.faa
data/cp12_combined_pf02672.sto
```

---

## Output Metadata Schema (v3.1)

The file `cp12_combined.metadata_v3.1.txt` is tab-delimited with the following columns:

| Column | Description |
|--------|-------------|
| `Accession` | unique protein identifier (`new_gene_name` for metagenome; accession for references) |
| `SourceType` | `Metagenome` or `Reference` |
| `SourceDB` | source database (e.g., `Metagenome`, `NCBI_refseq_CP12`, `UniProt`) |
| `Organism` | organism name (NA for metagenome sequences) |
| `TaxID` | NCBI taxonomy ID (NA for metagenome) |
| `Lineage_ncbi` | NCBI-style lineage string |
| `Lineage_consensus` | consensus lineage from CAT GTDB + NR (metagenome) |
| `pfamsearch:E-value` | Pfam hmmscan E-value |
| `pfamsearch:query_name` | query gene name |
| `pfamsearch:score` | Pfam bit score |
| `pfamsearch:target_accession` | Pfam accession matched (e.g., `PF02672.20`) |
| `pfamsearch:target_name` | Pfam target name |
| `SequenceStatus` | `full` or `partial` (based on contig edge proximity) |
| `Length` | amino acid sequence length |
| `NumCys` | number of cysteine residues |
| `CysPositions` | list of 1-based cysteine positions |
| `CysSpacings` | pairwise distances between cysteines |
| `AcidicFraction` | fraction of D + E residues |
| `CP12_Class` | `Canonical_CP12`, `Divergent_CP12`, or `Likely_false_positive` |
| `CP12_ClassReason` | plain-text justification for the classification |
| `Lineage_consensus_cleaned` | short `Domain:Phylum` label for visualization |

### CP12 Classification Rules

| Class | Criteria |
|-------|---------|
| `Likely_false_positive` | length < 50 or > 250 aa, fewer than 2 cysteines, or acidic fraction < 0.15 |
| `Canonical_CP12` | length 60–150 aa, exactly 2 cysteines, acidic fraction ≥ 0.20, single Cys spacing of 15–80 residues |
| `Divergent_CP12` | passes basic biochemical thresholds but does not meet canonical criteria |

---

## Notes

- Gene IDs in metagenome-derived sequences follow the convention `<SAMPLE>_PNNL_<GENE>`. If the original gene name already contains `PNNL`, it is used as-is.
- CD-HIT dereplication uses 100% identity (`-c 1.0`), preserving all non-redundant sequences.
- The SLURM `.sb` scripts contain hardcoded paths specific to the Tahoma HPC cluster at EMSL. Update `BASE_DIR`, `KOF_DIR`, `PFAM_DIR`, `GTDB_DIR`, `NR_DIR`, `VIRUS_DIR`, and tool binary paths before running on a different system.
- An example curated exclusion list (`curated_exclude_ids.example.txt`) is provided to illustrate the one-accession-per-line format expected by `--exclude-accessions`.
