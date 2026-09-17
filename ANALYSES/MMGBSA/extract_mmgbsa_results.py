#!/usr/bin/env python3
"""Extract per-system MM/GBSA binding energies from
<ENS>/postprocess_final/<sys>/mmgbsa/FINAL_RESULTS_MMPBSA.dat.

The gmx_MMPBSA configuration in this run is Generalized-Born (no Poisson).
The receptor mask covers the GAPDH + PRK chains; the ligand mask is the
four CP12 chains. ΔTOTAL therefore measures **CP12 binding energy** to
the GAPDH·PRK supercomplex.

Outputs
-------
PRK_SSG_RERUN/ensemble/mmgbsa_summary_structural_only_90pct/
   mmgbsa_per_system_<ENS>.csv     -- one row per system with Δ block + counts
   mmgbsa_per_system_combined.csv  -- both ensembles concatenated with ensemble col
   extraction_summary.txt
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

SCRATCH = Path("/scratch/path/to")
ENS_BASES = {
    "MED4": SCRATCH / "MED4_NAD" / "postprocess_final",
    "HM2":  SCRATCH / "HM2_NAD"  / "postprocess_final",
}
ENS_TOP5_SOURCE = {
    "MED4": SCRATCH / "MED4_NAD" / "postprocess_final" / "ptm_enrichment" / "ptm_occupancy_binary.csv",
    "HM2":  SCRATCH / "HM2_NAD"  / "postprocess_final" / "ptm_enrichment" / "ptm_occupancy_binary.csv",
}
OUT_DIR = SCRATCH / "PRK_SSG_RERUN" / "ensemble" / "mmgbsa_summary_structural_only_90pct"

# Components to extract from the Δ block. We capture the average, SD, and SEM.
DELTA_COMPONENTS = ("BOND", "ANGLE", "DIHED", "VDWAALS", "EEL",
                    "1-4 VDW", "1-4 EEL", "EGB", "ESURF",
                    "GGAS", "GSOLV", "TOTAL")


def parse_delta_block(path: Path) -> dict[str, float] | None:
    """Return {component: average, ..., component_sd, component_sem, n_frames}.

    Returns None if the file is malformed or the Δ block is missing.
    """
    try:
        text = path.read_text()
    except OSError:
        return None

    # Pull n_frames from the header (line like 'Calculations performed using 201 complex frames')
    m_frames = re.search(r"performed using\s+(\d+)\s+complex frames", text)
    n_frames = int(m_frames.group(1)) if m_frames else None

    # Slice to the Δ block: after the literal "Delta (Complex - Receptor - Ligand):"
    marker = "Delta (Complex - Receptor - Ligand):"
    idx = text.find(marker)
    if idx < 0:
        return None
    # Each Δ component only ever appears in the Δ block, so we can safely
    # scan the rest of the file without worrying about cross-section leakage.
    delta_block = text[idx:]

    out: dict[str, float] = {"n_frames": n_frames}
    for comp in DELTA_COMPONENTS:
        # Match a line like: "ΔVDWAALS               -933.15         24.38      63.15         1.72       4.45"
        # (Allow ASCII or Δ prefix; allow + or - sign; one or more spaces)
        pat = rf"[Δ]?{re.escape(comp)}\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        m = re.search(pat, delta_block)
        if m:
            avg, sd_prop, sd, sem_prop, sem = (float(x) for x in m.groups())
            key = comp.replace(" ", "").replace("-", "_")  # 1-4 VDW -> 14VDW
            out[f"d{key}_mean"] = avg
            out[f"d{key}_sd"]   = sd
            out[f"d{key}_sem"]  = sem
    return out


def extract_ensemble(ens: str) -> pd.DataFrame:
    """Extract top-5% systems plus REF (unmodified baseline).

    Per-system `mmgbsa/` subdirectories from prior PCA rankings exist on disk
    (Jan-2026 batches) and must be excluded; they used a different selection
    and are scientifically invalid against the current PCA. The truth gate
    for the top-5% set is the `in_top5pct=True` column of the current ranked
    occupancy table.

    REF (the unmodified reference trajectory) is included unconditionally so
    that downstream code can compute ΔΔTOTAL = ΔTOTAL_system − ΔTOTAL_REF.
    REF's MM/GBSA run uses the same trajectory and topology across PCA
    selections, so its result is methodologically valid even when it dates
    from a prior batch.
    """
    occ = pd.read_csv(ENS_TOP5_SOURCE[ens], dtype={"system": str})
    top5 = set(occ.loc[occ["in_top5pct"] == True, "system"])

    base = ENS_BASES[ens]
    rows = []
    skipped_parse = []
    skipped_stale = []
    for s in sorted(d for d in os.listdir(base) if (base / d).is_dir()):
        p = base / s / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
        if not p.exists():
            continue
        if s != "REF" and s not in top5:
            skipped_stale.append(s)
            continue
        parsed = parse_delta_block(p)
        if parsed is None:
            skipped_parse.append(s)
            continue
        row = {"system": s, "ensemble": ens,
               "is_ref": (s == "REF"), **parsed}
        rows.append(row)
    df = pd.DataFrame(rows)
    n_top = int((~df["is_ref"]).sum())
    n_ref = int(df["is_ref"].sum())
    print(f"[{ens}] extracted {n_top} top-5% systems + {n_ref} REF; "
          f"excluded {len(skipped_stale)} stale mmgbsa dirs (prior PCA)")
    if skipped_parse:
        print(f"[{ens}] WARNING: {len(skipped_parse)} systems with unparseable mmgbsa: {skipped_parse}")
    return df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frames = []
    for ens in ("MED4", "HM2"):
        df = extract_ensemble(ens)
        df.to_csv(OUT_DIR / f"mmgbsa_per_system_{ens}.csv", index=False)
        print(f"[{ens}] wrote {len(df)} systems "
              f"(median ΔTOTAL = {df['dTOTAL_mean'].median():.1f} kcal/mol, "
              f"IQR {df['dTOTAL_mean'].quantile(0.25):.1f}..{df['dTOTAL_mean'].quantile(0.75):.1f})")
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUT_DIR / "mmgbsa_per_system_combined.csv", index=False)

    summary_lines = [
        f"Extracted from FINAL_RESULTS_MMPBSA.dat under "
        f"{{MED4,HM2}}_NAD/postprocess_final/<sys>/mmgbsa/",
        f"Method                 : Generalized Born (gmx_MMPBSA v1.6.4)",
        f"Receptor               : GAPDH chains + PRK chains (residues per .dat masks)",
        f"Ligand                 : CP12 chains (4 × 74 residues)",
        f"ΔTOTAL semantics       : binding energy of CP12 to the GAPDH·PRK supercomplex (more negative = stronger)",
        "",
    ]
    for ens, df in zip(("MED4", "HM2"), frames):
        summary_lines += [
            f"--- {ens} ---",
            f"  systems with mmgbsa : {len(df)}",
            f"  mean n_frames       : {df['n_frames'].mean():.0f}",
            f"  ΔTOTAL median       : {df['dTOTAL_mean'].median():.1f} kcal/mol",
            f"  ΔTOTAL mean ± SD    : {df['dTOTAL_mean'].mean():.1f} ± {df['dTOTAL_mean'].std():.1f}",
            f"  ΔTOTAL min / max    : {df['dTOTAL_mean'].min():.1f} / {df['dTOTAL_mean'].max():.1f}",
            "",
        ]
    text = "\n".join(summary_lines)
    (OUT_DIR / "extraction_summary.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
