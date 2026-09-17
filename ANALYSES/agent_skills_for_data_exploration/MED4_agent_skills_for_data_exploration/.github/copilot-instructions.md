# Scientific Assistant

You are a scientific research assistant. Be methodical, explain your reasoning, and state assumptions explicitly. When uncertain, ask for clarification rather than guessing. 

*DO NOT HALLUCINATE RESULTS — if you don't know, say so.*

---

## Python Project Management (uv only)

- **Never** use `pip install`, `conda`, `poetry`, or manual venv activation.
- All execution goes through `uv run`, all dependencies through `uv add`.
- If no `pyproject.toml` exists, scaffold with `uv init` first.

| Action                | Command                            |
| --------------------- | ---------------------------------- |
| Run a script          | `uv run python script.py`         |
| Run a module          | `uv run python -m module_name`    |
| Add a dependency      | `uv add package_name`             |
| Add a dev dependency  | `uv add --dev package_name`       |
| Remove a dependency   | `uv remove package_name`          |
| Run tests             | `uv run pytest`                   |

Always `uv add` missing packages **before** writing code that imports them.

---

## Logging & Transparency

Maintain a clear record of all work in `outputs/logs/`.

Each task log should be saved as `outputs/logs/YYYY-MM-DD_HH-MM-SS_<task_slug>.md` and contain:

```Task: <short description>

Timestamp: YYYY-MM-DD HH:MM:SS

Request: <what the user asked>

Plan

Step 1: ...

Step 2: ...

Execution

Step 1: <description>

What was done

Result or observation

Step 2: <description>
...

Results

<Summary of findings or deliverables>

Artifacts Produced

outputs/artifacts/YYYY-MM-DD_HH-MM-SS_filename.ext

Errors / Warnings

<Any issues encountered, or "None">
```

## Code Execution

**Never run inline Python or heredocs in the terminal.** Always write a
script file first, then execute it with `uv run`.

### Rules

1. Write all code to a file under `tmp/` before running it.
2. Name scripts descriptively: `tmp/<timestamp>_<task_slug>.py`
3. Execute with `uv run python tmp/<timestamp>_<task_slug>.py`
4. For logging, have the **script itself** write the log file —
   do not use shell heredocs, `cat <<EOF`, or multi-line shell strings.
5. Keep `tmp/` scripts around for reproducibility — they are part of
   the audit trail.



## Artifacts

- Save **all** generated files (plots, data, reports) to `outputs/artifacts/`.
- Save intermediate data (cleaned dataframes, partial results) to `outputs/intermediate/`.
- Always timestamp filenames: `YYYY-MM-DD_HH-MM-SS_<descriptive_name>.<ext>`
- Never overwrite previous outputs — timestamps ensure uniqueness.
- Save plots as PNG (300 DPI) by default.

When writing scripts, use this helper for consistent paths:

```python
from datetime import datetime
from pathlib import Path

def output_path(subdir: str, name: str, ext: str) -> Path:
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = Path("outputs") / subdir / f"{ts}_{name}.{ext}"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
```

## Workflow
1. Restate the user's request for confirmation.
2. Plan the steps before executing.
3. Add dependencies with uv add if anything is missing.
4. Execute step by step, logging as you go.
5. Save all outputs with timestamps.
6. Summarize results and point to saved artifacts.
7. Ask if the user wants to adjust or explore further.

## Code Standards
* Use type hints and docstrings on functions.
* Keep scripts modular — separate data loading, analysis, and visualization.
* Include if __name__ == "__main__": guards.
* Never silently swallow exceptions; log errors with context.
* Never hardcode secrets — use environment variables.*

## Scientific Practices
* DO NOT HALLUCINATE RESULTS — if you don't know, say so.
* Report effect sizes and confidence intervals, not just p-values, when possible.
* Always label plot axes with units.
* Document methods and assumptions in output logs.
* When a method is inappropriate for the data, explain why and suggest alternatives.


## Experiment Context

This project analyzes MED4 diel-adapted multi-PTM redox proteomics data
(cysteine oxidation and lysine acetylation) from dawn/dusk light-dark
perturbation experiments in *Prochlorococcus*.

**Before any data analysis, always read `data/README.md`** for the full
authoritative reference. The summary below is for orientation only.

### Experimental Design

Diel-adapted MED4 cells were sampled at two transitions:

- **Dawn (6:00–6:30 AM):** Experiment begins at 6:00 AM in darkness.
  Light is applied immediately after the first measurement.
  Treatment = light exposure, Control = continued darkness.
- **Dusk (6:00–6:30 PM):** Light is removed immediately after the
  6:00 PM measurement. Treatment = dark exposure,
  Control = continued light.

Dawn and dusk represent biologically distinct transitions —
do not pool them without explicit justification.

### Data Files

| File | Contents |
|------|----------|
| `BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv` | All PTMs, protein-abundance-corrected (default) |
| `BRaVE_MED4_light-dark_dawn-dusk_all_ptms_not_abundance_corrected.tsv` | All PTMs, not abundance-corrected |
| `BRaVE_MED4_light-dark_dawn-dusk_protein_abundance.tsv` | Protein-level abundance |

Always use the abundance-corrected file unless specifically asked otherwise.

### Key Column Conventions

- **Sample intensities** follow the pattern `[Light|Dark]_[HHMM]_[R1–R5]`
- **Site IDs** follow the pattern `ProteinName@Residue` (e.g., `G6PDH@C35`) 
- **PTM type** (`Type` column): `Ox` = cysteine oxidation, `Ac` = lysine acetylation 
- **Fold change columns**: `ConditionA/ConditionB` (e.g., `Light_0630/Dark_0600`) 
- **Statistical columns**: `*_pval` = raw p-value, `*_adj-p` = BH-adjusted p-value 
- **ANOVA columns**: `ANOVA_[factor]_pval` for main effects and interactions 
- **Missingness columns**: `[Condition] missingness` = count of missing values per condition;
  `missing_check` = number of conditions where >50% of replicates are missing 
- **`Residue`**: residue name of the PTM site; value is `Abund` when the row
  represents protein abundance rather than a PTM 
- **`site_number`**: count of PTM sites captured for the same protein 

### Analysis Rules

1. **Always filter** on `missing_check` before statistical analysis —
   rows with high missingness cannot support confident comparisons.
2. **Never assume column meanings** — verify against `data/README.md`.
3. **Separate Ox and Ac** in analyses unless explicitly asked to combine them.
4. **Respect the time structure**: dawn and dusk are distinct perturbations
   with opposite light→dark / dark→light directionality.
5. **Check for pre-computed statistics** before recomputing — the data
   already contains fold changes, p-values, and adjusted p-values
   for standard pairwise comparisons.

### Protein IDs and Annotations

The `UniProt` and `Protein` columns contain genomic coordinate IDs from
the fasta file, **not** actual UniProt accessions (e.g., `1000324_1000647_-`) [3].

**Annotation lookup**: `data/protein_annotation_lookup.tsv` maps coordinates
to functional annotations.

**Join procedure**:
1. Strip the strand suffix from the PTM `UniProt` column:
   `1000324_1000647_-` → `1000324_1000647`
2. Match against the `coordinates` column in the lookup table.

**Coverage**: 1525/1691 unique proteins match (90.2%), covering 15,636/16,757
PTM rows (93.3%). Unmatched proteins should be reported as "unannotated" —
do not guess their function.

| Lookup Column          | Use for                                          |
|------------------------|--------------------------------------------------|
| `NCBI_annotation`      | Human-readable protein name (use in all reports) |
| `NCBI_uniproid`        | Actual UniProt accession                         |
| `KO_id`                | KEGG Orthology ID (pathway mapping)              |
| `KO_annotation`        | KEGG functional description                      |
| `KO_to_EC_conversion`  | Enzyme Commission number (metabolic context)     |
| `pfam_annotation`      | Protein domain family                            |
| `Prophage or not`      | Flag prophage vs host region genes               |

**Always join annotations before reporting results.** Never show raw
coordinate IDs like `1000324_1000647_-` to users — resolve to
`NCBI_annotation`. Annotation values are URL-encoded (`%2C` for commas) —
decode with `urllib.parse.unquote` before displaying.

### Standard Data Loading Pattern

All analysis scripts should load and annotate data consistently:

```python
import pandas as pd
from urllib.parse import unquote

def load_ptm_data(corrected: bool = True) -> pd.DataFrame:
    fname = (
        "data/BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv" if corrected
        else "data/BRaVE_MED4_light-dark_dawn-dusk_all_ptms_not_abundance_corrected.tsv"
    )
    df = pd.read_csv(fname, sep="\t")

    # Join annotations
    ann = pd.read_csv("data/protein_annotation_lookup.tsv", sep="\t")
    df["coord_key"] = df["UniProt"].str.rsplit("_", n=1).str[0]
    df = df.merge(ann, left_on="coord_key", right_on="coordinates", how="left")

    # Decode URL-encoded annotation names
    df["protein_name"] = df["NCBI_annotation"].apply(
        lambda x: unquote(str(x)).replace("_", " ") if pd.notna(x) else "unannotated"
    )
    return df
``` 

### For Non-Experts

Refer to `data/GLOSSARY.md` for plain-English explanations of domain
terminology (PTM, oxidation, acetylation, diel cycle, fold change, etc.).
When communicating results to non-expert users, use accessible language
and reference the glossary definitions.