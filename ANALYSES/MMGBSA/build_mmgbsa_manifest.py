#!/usr/bin/env python3
"""Build MM/GBSA downselection and replacement manifests from final PCA ranks."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path("/scratch/path/to/PRK_SSG_RERUN")
OUT_DIR = ROOT / "ensemble" / "mmgbsa_selection"

ENSEMBLES = {
    "HM2": {
        "rank_csv": Path(
            "/scratch/path/to/HM2_NAD/postprocess_final/"
            "suman_pca_structural_only_90pct/pca_distance_to_ref_90pct_ranked.csv"
        ),
        "affected": ROOT / "HM2" / "PRK_SSG_RERUN_INTERACTIVE" / "affected_systems.txt",
        "postprocess_new": Path("/scratch/path/to/HM2_NAD/postprocess_new"),
        "postprocess_final": Path("/scratch/path/to/HM2_NAD/postprocess_final"),
        "rerun_input": ROOT / "HM2" / "PRK_SSG_RERUN_INTERACTIVE",
        "frontier_tuples": Path(
            "/path/to"
            "HM2_NAD/UID/1tuples"
        ),
    },
    "MED4": {
        "rank_csv": Path(
            "/scratch/path/to/MED4_NAD/postprocess_final/"
            "suman_pca_structural_only_90pct/pca_distance_to_ref_90pct_ranked.csv"
        ),
        "affected": ROOT / "MED4" / "PRK_SSG_RERUN_INTERACTIVE" / "affected_systems.txt",
        "postprocess_new": Path("/scratch/path/to/MED4_NAD/postprocess_new"),
        "postprocess_final": Path("/scratch/path/to/MED4_NAD/postprocess_final"),
        "rerun_input": ROOT / "MED4" / "PRK_SSG_RERUN_INTERACTIVE",
        "frontier_tuples": Path(
            "/path/to"
            "MED4_NAD/UID/1tuples"
        ),
    },
}


def read_ids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def read_ranked(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for rank, row in enumerate(csv.DictReader(handle), start=1):
            system = row["system"].strip()
            row["rank"] = str(rank)
            rows[system] = row
    return rows


def complete_mmgbsa_systems(root: Path) -> set[str]:
    systems: set[str] = set()
    for result in root.glob("*/mmgbsa/FINAL_RESULTS_MMPBSA.dat"):
        systems.add(result.parents[1].name)
    return systems


def final_inputs_exist(root: Path, system: str) -> bool:
    system_dir = root / system
    required = [
        "solute_fixed.pdb",
        "solute_fit.xtc",
        "protein.pdb",
        "index.ndx",
    ]
    return system_dir.is_dir() and any((system_dir / name).exists() for name in required[:2])


def topol_source(paths: dict[str, Path], system: str) -> str:
    candidates = [
        paths["rerun_input"] / system / "topol.top",
        paths["postprocess_new"] / system / "mmgbsa" / "topol.top",
        paths["postprocess_new"] / system / "topol.top",
        paths["frontier_tuples"] / system / "topol.top",
        paths["postprocess_final"] / system / "mmgbsa" / "topol.top",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return ""


def decision(system: str, affected: bool, old_complete: bool, final_exists: bool) -> str:
    if system == "REF":
        return "reuse_ref_control_if_mmgbsa_present" if old_complete else "run_ref_control_from_postprocess_final"
    if affected:
        return "rerun_from_postprocess_final_replace_affected_old_result"
    if old_complete:
        return "reuse_postprocess_new_mmgbsa"
    if final_exists:
        return "run_from_postprocess_final_missing_old_result"
    return "missing_corrected_postprocess_inputs"


def write_list(path: Path, ids: list[str]) -> None:
    path.write_text("\n".join(ids) + ("\n" if ids else ""))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    ranked: dict[str, dict[str, dict[str, str]]] = {}
    affected: dict[str, set[str]] = {}
    old_complete: dict[str, set[str]] = {}

    for ensemble, paths in ENSEMBLES.items():
        ranked[ensemble] = read_ranked(paths["rank_csv"])
        affected[ensemble] = read_ids(paths["affected"])
        old_complete[ensemble] = complete_mmgbsa_systems(paths["postprocess_new"])

    non_ref_count = len([s for s in ranked["HM2"] if s != "REF"])
    top_n = math.ceil(non_ref_count * 0.05)
    top = {
        ensemble: {
            row["system"]
            for row in sorted(rows.values(), key=lambda r: int(r["rank"]))[:top_n]
            if row["system"] != "REF"
        }
        for ensemble, rows in ranked.items()
    }

    overlap = sorted(top["HM2"] & top["MED4"])
    union = sorted(top["HM2"] | top["MED4"])

    for ensemble in ENSEMBLES:
        ranked_top = [
            row["system"]
            for row in sorted(ranked[ensemble].values(), key=lambda r: int(r["rank"]))
            if row["system"] != "REF"
        ][:top_n]
        write_list(OUT_DIR / f"{ensemble.lower()}_top5pct_structural90_systems.txt", ranked_top)

    write_list(OUT_DIR / "tierA_top5pct_overlap_systems.txt", overlap)
    write_list(OUT_DIR / "tierB_top5pct_union_systems.txt", union)

    for ensemble in ENSEMBLES:
        invalid_old = sorted(old_complete[ensemble] & affected[ensemble])
        valid_old = sorted((old_complete[ensemble] - affected[ensemble]) - {"REF"})
        write_list(OUT_DIR / f"{ensemble.lower()}_mmgbsa_invalid_affected_old_results.txt", invalid_old)
        write_list(OUT_DIR / f"{ensemble.lower()}_mmgbsa_valid_reusable_old_results.txt", valid_old)

    rows = []
    for system in ["REF"] + union:
        row = {
            "system": system,
            "tier": "REF_control" if system == "REF" else ("A_top5pct_overlap" if system in overlap else "B_top5pct_union_only"),
            "in_hm2_top5pct": str(system in top["HM2"]),
            "in_med4_top5pct": str(system in top["MED4"]),
        }
        for ensemble, paths in ENSEMBLES.items():
            rank_row = ranked[ensemble].get(system, {})
            is_affected = system in affected[ensemble]
            is_old_complete = system in old_complete[ensemble]
            has_final = (paths["postprocess_final"] / system).is_dir()
            has_final_inputs = final_inputs_exist(paths["postprocess_final"], system)
            topol = topol_source(paths, system)
            source = decision(system, is_affected, is_old_complete, has_final)
            needs_new = source.startswith("run_") or source.startswith("rerun_")
            row.update(
                {
                    f"{ensemble.lower()}_rank": rank_row.get("rank", ""),
                    f"{ensemble.lower()}_variance_normalized_distance_to_ref": rank_row.get(
                        "variance_normalized_distance_to_ref", ""
                    ),
                    f"{ensemble.lower()}_affected": str(is_affected),
                    f"{ensemble.lower()}_old_mmgbsa_complete": str(is_old_complete),
                    f"{ensemble.lower()}_postprocess_final_exists": str(has_final),
                    f"{ensemble.lower()}_postprocess_final_inputs_exist": str(has_final_inputs),
                    f"{ensemble.lower()}_topol_source": topol,
                    f"{ensemble.lower()}_mmgbsa_action": source,
                    f"{ensemble.lower()}_needs_new_mmgbsa": str(needs_new),
                    f"{ensemble.lower()}_ready_for_new_mmgbsa": str(needs_new and has_final_inputs and bool(topol)),
                }
            )
        rows.append(row)

    fieldnames = list(rows[0].keys())
    with (OUT_DIR / "mmgbsa_downselection_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    for ensemble in ENSEMBLES:
        needs = sorted(
            row["system"]
            for row in rows
            if row["system"] != "REF" and row[f"{ensemble.lower()}_needs_new_mmgbsa"] == "True"
        )
        ready = sorted(
            row["system"]
            for row in rows
            if row["system"] != "REF" and row[f"{ensemble.lower()}_ready_for_new_mmgbsa"] == "True"
        )
        blocked = sorted(set(needs) - set(ready))
        write_list(OUT_DIR / f"{ensemble.lower()}_systems_needing_new_mmgbsa.txt", needs)
        write_list(OUT_DIR / f"{ensemble.lower()}_systems_ready_for_new_mmgbsa.txt", ready)
        write_list(OUT_DIR / f"{ensemble.lower()}_systems_needing_new_mmgbsa_missing_inputs.txt", blocked)

        top_ready = [
            row["system"]
            for row in rows
            if row["system"] != "REF"
            and row[f"in_{ensemble.lower()}_top5pct"] == "True"
            and row[f"{ensemble.lower()}_ready_for_new_mmgbsa"] == "True"
        ]
        top_needs = [
            row["system"]
            for row in rows
            if row["system"] != "REF"
            and row[f"in_{ensemble.lower()}_top5pct"] == "True"
            and row[f"{ensemble.lower()}_needs_new_mmgbsa"] == "True"
        ]
        write_list(OUT_DIR / f"{ensemble.lower()}_top5pct_ready_for_new_mmgbsa.txt", sorted(top_ready))
        write_list(OUT_DIR / f"{ensemble.lower()}_top5pct_needing_new_mmgbsa.txt", sorted(top_needs))

    summary_rows = [
        ("non_ref_system_count", str(non_ref_count)),
        ("top5pct_n", str(top_n)),
        ("tierA_top5pct_overlap_count", str(len(overlap))),
        ("tierB_top5pct_union_count", str(len(union))),
        ("hm2_old_mmgbsa_complete_count_in_postprocess_new", str(len(old_complete["HM2"]))),
        ("med4_old_mmgbsa_complete_count_in_postprocess_new", str(len(old_complete["MED4"]))),
        ("hm2_invalid_affected_old_mmgbsa_count", str(len(old_complete["HM2"] & affected["HM2"]))),
        ("med4_invalid_affected_old_mmgbsa_count", str(len(old_complete["MED4"] & affected["MED4"]))),
        (
            "hm2_selected_nonref_needing_new_mmgbsa_count",
            str(sum(row["hm2_needs_new_mmgbsa"] == "True" and row["system"] != "REF" for row in rows)),
        ),
        (
            "hm2_selected_nonref_ready_for_new_mmgbsa_count",
            str(sum(row["hm2_ready_for_new_mmgbsa"] == "True" and row["system"] != "REF" for row in rows)),
        ),
        (
            "hm2_selected_nonref_needing_new_mmgbsa_missing_inputs_count",
            str(
                sum(
                    row["hm2_needs_new_mmgbsa"] == "True"
                    and row["hm2_ready_for_new_mmgbsa"] != "True"
                    and row["system"] != "REF"
                    for row in rows
                )
            ),
        ),
        (
            "med4_selected_nonref_needing_new_mmgbsa_count",
            str(sum(row["med4_needs_new_mmgbsa"] == "True" and row["system"] != "REF" for row in rows)),
        ),
        (
            "med4_selected_nonref_ready_for_new_mmgbsa_count",
            str(sum(row["med4_ready_for_new_mmgbsa"] == "True" and row["system"] != "REF" for row in rows)),
        ),
        (
            "med4_selected_nonref_needing_new_mmgbsa_missing_inputs_count",
            str(
                sum(
                    row["med4_needs_new_mmgbsa"] == "True"
                    and row["med4_ready_for_new_mmgbsa"] != "True"
                    and row["system"] != "REF"
                    for row in rows
                )
            ),
        ),
    ]
    with (OUT_DIR / "mmgbsa_downselection_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    readme = f"""# MM/GBSA Downselection Manifest

This directory records the PCA-based MM/GBSA selection set.

Selection rule:
- Use the final 90% structural-only PCA rankings for downselection.
- Select the top 5% of non-REF systems per ensemble by variance-normalized distance to REF.
- With {non_ref_count} ranked non-REF systems, top 5% is {top_n} systems per ensemble.
- Tier A is the HM2/MED4 top-5% overlap.
- Tier B is the HM2/MED4 top-5% union.
- REF is kept as a separate control and is not counted in the top-5% denominator.

Replacement rule:
- Existing MM/GBSA results under `postprocess_new` are reusable only when the system is not in the affected/rerun list.
- Existing MM/GBSA results for affected systems are marked invalid and should be replaced from corrected `postprocess_final` trajectories.
- Missing selected non-affected systems should be run from `postprocess_final` so each selected system has HM2 and MED4 MM/GBSA coverage.

Primary files:
- `mmgbsa_downselection_manifest.csv`: per-system HM2/MED4 ranks, affected status, old-result status, and action.
- `hm2_top5pct_structural90_systems.txt` and `med4_top5pct_structural90_systems.txt`: direct top-5% lists from the primary 90% structural-only PCA rankings.
- `tierA_top5pct_overlap_systems.txt`: focused cross-ensemble overlap set.
- `tierB_top5pct_union_systems.txt`: inclusive selected set.
- `hm2_systems_needing_new_mmgbsa.txt` and `med4_systems_needing_new_mmgbsa.txt`: corrected/new runs needed by ensemble.
- `hm2_systems_ready_for_new_mmgbsa.txt` and `med4_systems_ready_for_new_mmgbsa.txt`: needed runs with both corrected trajectories and a located `topol.top`.
- `hm2_top5pct_ready_for_new_mmgbsa.txt` and `med4_top5pct_ready_for_new_mmgbsa.txt`: corrected/new runs for the original top-5%-per-ensemble MM/GBSA plan.
- `*_missing_inputs.txt`: selected systems that still need input/topology reconstruction before MM/GBSA can run.
"""
    (OUT_DIR / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
