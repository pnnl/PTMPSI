#!/usr/bin/env python3
"""Stage MM/GBSA input directories under postprocess_final for selected systems."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path("/scratch/path/to/PRK_SSG_RERUN")
SELECT_ROOT = ROOT / "ensemble" / "mmgbsa_selection"

ENSEMBLES = {
    "HM2": {
        "base": Path("/scratch/path/to/HM2_NAD/postprocess_final"),
        "old_base": Path("/scratch/path/to/HM2_NAD/postprocess_new"),
        "rerun_input": ROOT / "HM2" / "PRK_SSG_RERUN_INTERACTIVE",
        "frontier_tuples": Path(
            "/path/to"
            "HM2_NAD/UID/1tuples"
        ),
        "template": Path("/scratch/path/to/HM2_NAD/postprocess_new/REF/mmgbsa"),
        "list": SELECT_ROOT / "hm2_systems_ready_for_new_mmgbsa.txt",
    },
    "MED4": {
        "base": Path("/scratch/path/to/MED4_NAD/postprocess_final"),
        "old_base": Path("/scratch/path/to/MED4_NAD/postprocess_new"),
        "rerun_input": ROOT / "MED4" / "PRK_SSG_RERUN_INTERACTIVE",
        "frontier_tuples": Path(
            "/path/to"
            "MED4_NAD/UID/1tuples"
        ),
        "template": Path("/scratch/path/to/MED4_NAD/postprocess_new/REF/mmgbsa"),
        "helper_template": Path("/scratch/path/to/MED4_NAD/postprocess_new/0327/mmgbsa"),
        "list": SELECT_ROOT / "med4_systems_ready_for_new_mmgbsa.txt",
    },
}


def read_ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def topol_source(paths: dict[str, Path], system: str) -> Path | None:
    candidates = [
        paths["rerun_input"] / system / "topol.top",
        paths["old_base"] / system / "mmgbsa" / "topol.top",
        paths["old_base"] / system / "topol.top",
        paths["frontier_tuples"] / system / "topol.top",
        paths["base"] / system / "mmgbsa" / "topol.top",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def forcefield_source(paths: dict[str, Path]) -> Path:
    source = paths["rerun_input"] / "amber99sb.ff"
    if not source.is_dir():
        raise FileNotFoundError(f"missing custom force field directory: {source}")
    return source


def symlink_force(target: str, link: Path) -> None:
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(target)


def stage(ensemble: str, paths: dict[str, Path], list_override: Path | None = None) -> tuple[int, list[str]]:
    staged = 0
    errors: list[str] = []
    for system in read_ids(list_override or paths["list"]):
        source_dir = paths["base"] / system
        workdir = source_dir / "mmgbsa"
        topol = topol_source(paths, system)
        required = ["solute_fixed.pdb", "solute_fit.xtc", "solute.tpr"]

        missing = [name for name in required if not (source_dir / name).exists()]
        if topol is None:
            missing.append("topol.top")
        if missing:
            errors.append(f"{ensemble} {system}: missing {','.join(missing)}")
            continue

        workdir.mkdir(exist_ok=True)
        shutil.copy2(paths["template"] / "mmpbsa.in", workdir / "mmpbsa.in")
        helper_dir = paths.get("helper_template", paths["template"])
        shutil.copy2(helper_dir / "make_index_groups_withNAD.py", workdir / "make_index_groups_withNAD.py")
        if topol.resolve() != (workdir / "topol.top").resolve():
            shutil.copy2(topol, workdir / "topol.top")
        ff_link = workdir / "amber99sb.ff"
        if ff_link.is_symlink() or ff_link.exists():
            if ff_link.is_dir() and not ff_link.is_symlink():
                shutil.rmtree(ff_link)
            else:
                ff_link.unlink()
        ff_link.symlink_to(forcefield_source(paths))
        symlink_force("../solute_fixed.pdb", workdir / "solute_fixed.pdb")
        symlink_force("../solute_fit.xtc", workdir / "solute_fit.xtc")
        symlink_force("../solute.tpr", workdir / "solute.tpr")
        staged += 1

    return staged, errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hm2-list", type=Path)
    parser.add_argument("--med4-list", type=Path)
    args = parser.parse_args()
    overrides = {"HM2": args.hm2_list, "MED4": args.med4_list}

    all_errors: list[str] = []
    for ensemble, paths in ENSEMBLES.items():
        staged, errors = stage(ensemble, paths, overrides[ensemble])
        print(f"{ensemble}: staged {staged} systems")
        all_errors.extend(errors)

    if all_errors:
        print("Errors:")
        for error in all_errors:
            print(f"  {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
