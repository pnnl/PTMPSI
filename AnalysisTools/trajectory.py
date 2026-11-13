"""Trajectory preprocessing and structural analysis for PTM-Psi GROMACS runs.

Standardizes post-processing around two tools:
1. GROMACS (`gmx_mpi` on HPC / `gmx` locally) for periodic-boundary-condition
   unwrapping (`trjconv -pbc whole` -> `-pbc nojump`) with PTM-Psi's custom
   `residuetypes.dat` staged so modified residues remain in `Protein`.
2. `MDAnalysis` loading `.tpr` + unwrapped `.xtc` so exact atomic masses, bonds,
   and nonstandard PTM residue names are preserved across RMSD, mass-weighted
   radius of gyration (RoG), iteratively aligned RMSF, minimum-image contact
   maps, and mass-weighted Cartesian PCA (single-system and joint WT-vs-PTM).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, Mapping, Sequence

import MDAnalysis as mda
from MDAnalysis.analysis import align, pca as mda_pca, rms
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib import distances as mda_distances
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent.parent
FORCEFIELD_DIR = REPO_ROOT / "ptmpsi" / "gromacs" / "forcefield"
RESIDUETYPES_DAT = FORCEFIELD_DIR / "residuetypes.dat"
RESIDUES_INIT_PY = REPO_ROOT / "ptmpsi" / "residues" / "__init__.py"

# Cofactors listed as Protein in GROMACS residuetypes.dat for group convenience,
# excluded from polypeptide backbone/protein selections by default.
COFACTOR_RESNAMES = frozenset({"NAD", "NAI", "NADP", "NDP", "FAD", "FMN", "HEM", "HEME"})

# Fallback PTM codes in case repository forcefield files are relocated.
_FALLBACK_PTM_RESNAMES = frozenset(
    {
        "SNC", "CSO", "CSD", "OCS", "CSS", "CGL", "SMC", "QCS", "XCN", "IYY",
        "ALY", "PTM", "CYZ", "CCYZ", "EAC", "ABA", "CEAC", "NEAC",
    }
)


def _load_ptmpsi_residue_dicts() -> tuple[set[str], set[str]]:
    """Read standard and PTM residue names from ptmpsi without importing rdkit/meeko."""
    if "ptmpsi.residues" in sys.modules:
        mod = sys.modules["ptmpsi.residues"]
        std = set(getattr(mod, "three2one", {}).keys())
        ptm_vals = set(getattr(mod, "ptm2nonstandard", {}).values())
        res_keys = set(getattr(mod, "resdict", {}).keys())
        return std, (ptm_vals | res_keys) - std

    std: set[str] = set()
    ptm_candidates: set[str] = set()
    if RESIDUES_INIT_PY.is_file():
        tree = ast.parse(RESIDUES_INIT_PY.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
                names = {t.id for t in node.targets if isinstance(t, ast.Name)}
                if "three2one" in names:
                    for k in node.value.keys:
                        if isinstance(k, ast.Constant) and isinstance(k.value, str):
                            std.add(k.value)
                elif "ptm2nonstandard" in names:
                    for v in node.value.values:
                        if isinstance(v, ast.Constant) and isinstance(v.value, str):
                            ptm_candidates.add(v.value)
                elif "resdict" in names:
                    for k in node.value.keys:
                        if isinstance(k, ast.Constant) and isinstance(k.value, str):
                            ptm_candidates.add(k.value)
    return std, ptm_candidates - std


def get_ptmpsi_residue_sets(include_cofactors: bool = False) -> dict[str, frozenset[str]]:
    """Return authoritative PTM-Psi residue-name sets (`protein`, `ptm`, `water`, `ion`, `nucleic`)."""
    std_from_code, ptm_from_code = _load_ptmpsi_residue_dicts()

    by_kind: dict[str, set[str]] = {
        "Protein": set(std_from_code),
        "Water": {"SOL", "WAT", "HOH", "TIP3", "TP3M"},
        "Ion": {"NA", "CL", "K", "MG", "CA", "ZN", "NA+", "CL-"},
        "DNA": set(),
        "RNA": set(),
    }

    if RESIDUETYPES_DAT.is_file():
        for raw in RESIDUETYPES_DAT.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                resname, rtype = parts[0], parts[1]
                by_kind.setdefault(rtype, set()).add(resname)

    ptm_names = set(_FALLBACK_PTM_RESNAMES) | ptm_from_code
    # Ensure canonical unmodified cysteines/amino acids are never marked as PTM-modified
    ptm_names -= std_from_code | {"CYS", "CYX", "CYM", "HIS", "HID", "HIE", "HIP"}

    protein_names = (by_kind["Protein"] | ptm_names | {"CCYX", "CCYZ"})
    if not include_cofactors:
        protein_names = protein_names - set(COFACTOR_RESNAMES)

    return {
        "protein": frozenset(protein_names),
        "ptm": frozenset(ptm_names),
        "standard_protein": frozenset(protein_names - ptm_names),
        "water": frozenset(by_kind["Water"]),
        "ion": frozenset(by_kind["Ion"]),
        "nucleic": frozenset(by_kind["DNA"] | by_kind["RNA"]),
    }


def protein_selection_str(
    kind: str = "all",
    *,
    include_ptms: bool = True,
    include_cofactors: bool = False,
    extra_resnames: Iterable[str] = (),
) -> str:
    """Build an MDAnalysis selection string that includes PTM-Psi modified residues.

    Parameters
    ----------
    kind : {"all", "backbone", "ca", "heavy"}
        Atom subset within the protein + PTM residue selection.
    include_ptms : bool, default True
        Include PTM-Psi modified residue names (`SNC`, `CGL`, `IYY`, `ALY`, etc.).
    include_cofactors : bool, default False
        Include cofactors such as `NAD`/`NAI`.
    extra_resnames : iterable of str
        Additional custom residue names to treat as protein.
    """
    sets = get_ptmpsi_residue_sets(include_cofactors=include_cofactors)
    resnames = set(sets["protein"] if include_ptms else sets["standard_protein"])
    resnames.update(extra_resnames)
    res_clause = "protein"
    if resnames:
        sorted_names = " ".join(sorted(resnames))
        if include_ptms or extra_resnames or include_cofactors:
            res_clause = f"(protein or resname {sorted_names})"
        else:
            res_clause = f"(protein and resname {sorted_names})"
    if not include_cofactors:
        cof_names = " ".join(sorted(COFACTOR_RESNAMES - set(extra_resnames)))
        if cof_names:
            res_clause = f"({res_clause} and not resname {cof_names})"

    kind_norm = kind.lower().strip()
    if kind_norm == "all":
        return res_clause
    if kind_norm == "backbone":
        # Explicit atom names instead of MDAnalysis 'backbone', which drops nonstandard PTM residues!
        return f"({res_clause} and name N CA C O)"
    if kind_norm == "ca":
        return f"({res_clause} and name CA)"
    if kind_norm == "heavy":
        return f"({res_clause} and not name H* [123]H* D* and not type H D DU)"
    raise ValueError(f"Unsupported selection kind {kind!r}; expected 'all', 'backbone', 'ca', or 'heavy'")


def select_protein(
    universe_or_ag: mda.Universe | mda.AtomGroup,
    kind: str = "all",
    *,
    include_ptms: bool = True,
    include_cofactors: bool = False,
    extra_resnames: Iterable[str] = (),
) -> mda.AtomGroup:
    """Select protein atoms while preserving PTM-Psi modified residues."""
    sel = protein_selection_str(
        kind=kind,
        include_ptms=include_ptms,
        include_cofactors=include_cofactors,
        extra_resnames=extra_resnames,
    )
    return universe_or_ag.select_atoms(sel)


def find_ptm_residues(universe_or_ag: mda.Universe | mda.AtomGroup) -> mda.ResidueGroup:
    """Return all residues in `universe_or_ag` carrying a PTM-Psi modified residue name."""
    sets = get_ptmpsi_residue_sets()
    non_ptm = sets["standard_protein"] | sets["water"] | sets["ion"] | sets["nucleic"] | COFACTOR_RESNAMES
    residues = universe_or_ag.residues
    mask = [
        (r.resname in sets["ptm"]) or (r.resname not in non_ptm and "CA" in set(r.atoms.names))
        for r in residues
    ]
    if not any(mask):
        return residues[[]]
    return residues[np.flatnonzero(mask)]


def find_gmx(gmx: str | None = None) -> str:
    """Locate `gmx_mpi` (Deception / HPC) or `gmx` (local desktop)."""
    if gmx:
        resolved = shutil.which(gmx) or (gmx if Path(gmx).is_file() else None)
        if not resolved:
            raise FileNotFoundError(f"Requested GROMACS executable {gmx!r} was not found on PATH")
        return resolved
    for candidate in ("gmx_mpi", "gmx"):
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise FileNotFoundError(
        "Neither 'gmx_mpi' nor 'gmx' was found on PATH. Load the GROMACS module on a compute node or install GROMACS locally."
    )


def ensure_residuetypes_dat(directory: str | Path) -> Path | None:
    """Copy PTM-Psi's `residuetypes.dat` into `directory` if absent; returns path if created."""
    target = Path(directory) / "residuetypes.dat"
    if target.exists() or not RESIDUETYPES_DAT.is_file():
        return None
    shutil.copy2(RESIDUETYPES_DAT, target)
    return target


def load_universe(
    topology: str | Path,
    trajectory: str | Path | Sequence[str | Path] | None = None,
    *,
    in_memory: bool = False,
    unwrap: bool = False,
) -> mda.Universe:
    """Load an `MDAnalysis.Universe` (prefer `.tpr` + `.xtc`) with optional on-the-fly PBC unwrap."""
    top_str = str(topology)
    if trajectory is None:
        u = mda.Universe(top_str, in_memory=in_memory)
    elif isinstance(trajectory, (str, Path)):
        u = mda.Universe(top_str, str(trajectory), in_memory=in_memory)
    else:
        u = mda.Universe(top_str, [str(p) for p in trajectory], in_memory=in_memory)

    if unwrap:
        from MDAnalysis import transformations as trans

        prot = select_protein(u, kind="all")
        solvent = u.atoms - prot
        ops = [trans.unwrap(u.atoms)]
        if len(prot) > 0:
            ops.append(trans.center_in_box(prot, center="mass"))
        if len(solvent) > 0:
            ops.append(trans.wrap(solvent, compound="residues"))
        u.trajectory.add_transformations(*ops)
    return u


def check_pbc_jumps(
    universe: mda.Universe,
    selection: str | None = None,
    *,
    max_step_angstrom: float = 15.0,
    max_ca_span_fraction: float = 0.85,
) -> dict[str, object]:
    """Detect likely un-unwrapped PBC jumps across trajectory frames.

    Checks two signatures of broken PBC wrapping:
    1. Consecutive-frame C-alpha (or selected atom) COM/coordinate jumps exceeding
       `max_step_angstrom`.
    2. Adjacent-residue C-alpha distances > 6.0 Å within the same chain/segid
       (indicating a chain split across a periodic boundary face).
    """
    ag = universe.select_atoms(selection) if selection else select_protein(universe, kind="ca")
    if len(ag) == 0:
        ag = select_protein(universe, kind="all")
    if len(ag) == 0:
        raise ValueError("Selection matched 0 atoms for PBC jump check")

    prev_pos: np.ndarray | None = None
    max_jump = 0.0
    jump_frames: list[int] = []
    split_chain_frames: list[int] = []

    same_chain_pairs: list[tuple[int, int]] = []
    if all(name == "CA" for name in ag.names) and len(ag) > 1:
        for idx in range(len(ag) - 1):
            a, b = ag[idx], ag[idx + 1]
            if a.segid == b.segid and (b.resid - a.resid) == 1:
                same_chain_pairs.append((idx, idx + 1))

    for ts in universe.trajectory:
        pos = ag.positions.copy()
        if prev_pos is not None:
            disp = np.linalg.norm(pos - prev_pos, axis=1)
            step_max = float(disp.max())
            if step_max > max_jump:
                max_jump = step_max
            if step_max > max_step_angstrom:
                jump_frames.append(int(ts.frame))
        if same_chain_pairs:
            idx_a = [p[0] for p in same_chain_pairs]
            idx_b = [p[1] for p in same_chain_pairs]
            ca_dists = np.linalg.norm(pos[idx_b] - pos[idx_a], axis=1)
            if float(ca_dists.max()) > 6.0:
                split_chain_frames.append(int(ts.frame))
        elif ts.dimensions is not None and np.all(ts.dimensions[:3] > 0):
            span = pos.max(axis=0) - pos.min(axis=0)
            if np.any(span > max_ca_span_fraction * ts.dimensions[:3]):
                pass
        prev_pos = pos

    has_jumps = bool(jump_frames or split_chain_frames)
    return {
        "has_pbc_jumps": has_jumps,
        "max_step_angstrom": max_jump,
        "jump_frames": jump_frames,
        "split_chain_frames": split_chain_frames,
        "n_frames": len(universe.trajectory),
    }


def unwrap_trajectory(
    topology: str | Path,
    trajectory: str | Path,
    output_xtc: str | Path,
    *,
    gmx: str | None = None,
    output_group: str = "System",
    center_group: str | None = "Protein",
    fit_group: str | None = None,
    index: str | Path | None = None,
) -> Path:
    """Unwrap PBC jumps in a GROMACS trajectory via `trjconv -pbc whole` -> `-pbc nojump`.

    Automatically stages PTM-Psi's `residuetypes.dat` in a temporary working
    directory so `Protein` / `Backbone` groups include PTM-Psi modified residues.
    """
    gmx_bin = find_gmx(gmx)
    top_path = Path(topology).resolve()
    traj_path = Path(trajectory).resolve()
    out_path = Path(output_xtc).resolve()
    idx_args = ["-n", str(Path(index).resolve())] if index is not None else []

    if not top_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {top_path}")
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="ptmpsi_trjconv_") as tmpdir:
        tmp = Path(tmpdir)
        ensure_residuetypes_dat(tmp)

        whole_xtc = tmp / "step1_whole.xtc"
        cmd_whole = [
            gmx_bin, "trjconv",
            "-s", str(top_path),
            "-f", str(traj_path),
            "-o", str(whole_xtc),
            "-pbc", "whole",
            *idx_args,
        ]
        if center_group:
            cmd_whole.extend(["-center"])
            stdin_whole = f"{center_group}\n{output_group}\n"
        else:
            stdin_whole = f"{output_group}\n"

        res1 = subprocess.run(cmd_whole, cwd=tmp, input=stdin_whole, text=True, capture_output=True)
        if res1.returncode != 0:
            raise RuntimeError(f"gmx trjconv -pbc whole failed:\n{res1.stderr}")

        nojump_target = (tmp / "step2_nojump.xtc") if fit_group else out_path
        cmd_nojump = [
            gmx_bin, "trjconv",
            "-s", str(top_path),
            "-f", str(whole_xtc),
            "-o", str(nojump_target),
            "-pbc", "nojump",
            *idx_args,
        ]
        res2 = subprocess.run(cmd_nojump, cwd=tmp, input=f"{output_group}\n", text=True, capture_output=True)
        if res2.returncode != 0:
            raise RuntimeError(f"gmx trjconv -pbc nojump failed:\n{res2.stderr}")

        if fit_group:
            cmd_fit = [
                gmx_bin, "trjconv",
                "-s", str(top_path),
                "-f", str(nojump_target),
                "-o", str(out_path),
                "-fit", "rot+trans",
                *idx_args,
            ]
            res3 = subprocess.run(
                cmd_fit,
                cwd=tmp,
                input=f"{fit_group}\n{output_group}\n",
                text=True,
                capture_output=True,
            )
            if res3.returncode != 0:
                raise RuntimeError(f"gmx trjconv -fit rot+trans failed:\n{res3.stderr}")

    return out_path


def _ensure_positive_masses(ag: mda.AtomGroup, context: str) -> np.ndarray:
    if not hasattr(ag, "masses"):
        raise ValueError(
            f"{context} requires atomic masses, which are missing from this topology. "
            "Load a GROMACS .tpr topology so exact force-field and PTM masses are present."
        )
    masses = np.asarray(ag.masses, dtype=float)
    # Dummy TI atoms (DU) legitimately have 0 mass or small mass; require positive total mass
    # and warn/raise if non-dummy heavy atoms all have zero mass.
    non_dummy = np.array([t != "DU" and not n.startswith("DU") for t, n in zip(
        getattr(ag, "types", [""] * len(ag)), ag.names
    )], dtype=bool)
    check_masses = masses[non_dummy] if np.any(non_dummy) else masses
    if check_masses.size == 0 or np.all(check_masses <= 0.0) or np.any(np.isnan(masses)):
        raise ValueError(
            f"{context} requires valid atomic masses, but found missing/zero masses. "
            "Pass a GROMACS .tpr file as topology rather than a mass-less .pdb/.xyz."
        )
    if np.any(check_masses <= 0.0):
        bad = [f"{a.resname}{a.resid}:{a.name}" for a, ok in zip(ag, masses > 0.0) if not ok and a.name != "DU"]
        if bad:
            raise ValueError(
                f"{context} found zero mass on atoms {bad[:5]}. Load a .tpr topology so PTM atoms carry exact masses."
            )
    return masses


def compute_rmsd(
    universe: mda.Universe,
    *,
    reference: mda.Universe | None = None,
    ref_frame: int = 0,
    select: str | None = None,
    groupselections: Sequence[str] | None = None,
    mass_weighted: bool = False,
) -> pd.DataFrame:
    """Compute superposed RMSD (in Å) across trajectory frames using `MDAnalysis.analysis.rms.RMSD`.

    Defaults to PTM-aware backbone (`N CA C O`) so modified residues are included.
    """
    sel = select or protein_selection_str("backbone")
    ref = reference if reference is not None else universe
    ag = universe.select_atoms(sel)
    if len(ag) == 0:
        raise ValueError(f"RMSD selection matched 0 atoms: {sel}")
    weights = "mass" if mass_weighted else None
    if mass_weighted:
        _ensure_positive_masses(ag, "Mass-weighted RMSD")

    rmsd_calc = rms.RMSD(
        universe,
        ref,
        select=sel,
        groupselections=list(groupselections) if groupselections else None,
        weights=weights,
        ref_frame=ref_frame,
    ).run()

    arr = rmsd_calc.results.rmsd
    cols = ["frame", "time_ps", "rmsd_angstrom"]
    if groupselections:
        cols.extend([f"group_{i}_rmsd_angstrom" for i in range(len(groupselections))])
    df = pd.DataFrame(arr, columns=cols)
    df["frame"] = df["frame"].astype(int)
    return df


def compute_rog(
    universe: mda.Universe,
    *,
    select: str | None = None,
    include_components: bool = False,
) -> pd.DataFrame:
    """Compute mass-weighted radius of gyration (in Å) across trajectory frames."""
    sel = select or protein_selection_str("all")
    ag = universe.select_atoms(sel)
    if len(ag) == 0:
        raise ValueError(f"RoG selection matched 0 atoms: {sel}")
    _ensure_positive_masses(ag, "Mass-weighted radius of gyration")

    rows: list[dict[str, float | int]] = []
    for ts in universe.trajectory:
        row: dict[str, float | int] = {
            "frame": int(ts.frame),
            "time_ps": float(ts.time),
            "rg_angstrom": float(ag.radius_of_gyration()),
        }
        if include_components:
            moments = ag.moment_of_inertia()
            w, _ = np.linalg.eigh(moments)
            total_m = float(ag.total_mass())
            radii = np.sqrt(np.maximum(w / total_m, 0.0))
            row["rg_p1_angstrom"], row["rg_p2_angstrom"], row["rg_p3_angstrom"] = (
                float(radii[0]), float(radii[1]), float(radii[2])
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _aligned_in_memory_copy(
    universe: mda.Universe,
    *,
    align_select: str,
    reference: mda.Universe | None = None,
    iterative_average: bool = True,
    mass_weighted: bool = True,
    ref_frame: int = 0,
) -> tuple[mda.Universe, mda.Universe]:
    """Return an aligned in-memory copy of `universe` without mutating the caller's trajectory."""
    coords = np.empty((len(universe.trajectory), universe.atoms.n_atoms, 3), dtype=np.float64)
    times = np.empty(len(universe.trajectory), dtype=np.float64)
    dt = float(getattr(universe.trajectory, "dt", 1.0) or 1.0)
    for i, ts in enumerate(universe.trajectory):
        coords[i] = ts.positions
        times[i] = float(ts.time)

    u_copy = mda.Merge(universe.atoms)
    u_copy.load_new(coords, format=MemoryReader, dt=dt)
    for i, ts in enumerate(u_copy.trajectory):
        ts.time = times[i]

    weights = "mass" if mass_weighted else None
    if mass_weighted:
        _ensure_positive_masses(u_copy.select_atoms(align_select), "Trajectory alignment")

    if reference is not None:
        ref_u = reference
    elif iterative_average and len(u_copy.trajectory) > 1:
        avg = align.AverageStructure(
            u_copy,
            u_copy,
            select=align_select,
            ref_frame=ref_frame,
            weights=weights,
        ).run()
        ref_u = avg.results.universe
    else:
        ref_u = mda.Merge(universe.atoms)
        universe.trajectory[ref_frame]
        ref_u.load_new(universe.atoms.positions[None, :, :].copy(), format=MemoryReader)

    align.AlignTraj(
        u_copy,
        ref_u,
        select=align_select,
        weights=weights,
        in_memory=True,
    ).run()
    return u_copy, ref_u


def compute_rmsf(
    universe: mda.Universe,
    *,
    select: str | None = None,
    align_select: str | None = None,
    reference: mda.Universe | None = None,
    by_residue: bool = True,
    mass_weighted_align: bool = True,
) -> pd.DataFrame:
    """Compute iteratively aligned RMSF (in Å) using `AverageStructure` + `AlignTraj` + `RMSF`.

    When `by_residue=True` (default), aggregates atomic fluctuations per residue using
    mass weights ($\\sqrt{\\sum_i m_i \\,\\mathrm{RMSF}_i^2 / \\sum_i m_i}$), matching `gmx rmsf -res`.
    """
    if len(universe.trajectory) < 2:
        raise ValueError("RMSF requires at least 2 trajectory frames")

    fit_sel = align_select or protein_selection_str("backbone")
    calc_sel = select or (protein_selection_str("heavy") if by_residue else protein_selection_str("ca"))

    u_aligned, _ = _aligned_in_memory_copy(
        universe,
        align_select=fit_sel,
        reference=reference,
        iterative_average=(reference is None),
        mass_weighted=mass_weighted_align,
    )
    ag = u_aligned.select_atoms(calc_sel)
    if len(ag) == 0:
        raise ValueError(f"RMSF selection matched 0 atoms: {calc_sel}")

    rmsf_vals = rms.RMSF(ag).run().results.rmsf
    ptm_set = get_ptmpsi_residue_sets()["ptm"]

    if not by_residue:
        return pd.DataFrame(
            {
                "atom_index": ag.indices,
                "atom_name": ag.names,
                "resid": ag.resids,
                "resname": ag.resnames,
                "segid": ag.segids,
                "is_ptm": [r in ptm_set for r in ag.resnames],
                "rmsf_angstrom": rmsf_vals,
            }
        )

    masses = _ensure_positive_masses(ag, "Mass-weighted residue RMSF")
    rows: list[dict[str, object]] = []
    for res in ag.residues:
        mask = ag.resindices == res.resindex
        m_r = masses[mask]
        f_r = rmsf_vals[mask]
        positive = m_r > 0
        if np.any(positive):
            res_rmsf = float(np.sqrt(np.sum(m_r[positive] * (f_r[positive] ** 2)) / np.sum(m_r[positive])))
        else:
            res_rmsf = float(np.sqrt(np.mean(f_r ** 2)))
        rows.append(
            {
                "resid": int(res.resid),
                "resname": str(res.resname),
                "segid": str(res.segid),
                "label": f"{res.segid}:{res.resname}{res.resid}" if res.segid else f"{res.resname}{res.resid}",
                "is_ptm": res.resname in ptm_set,
                "rmsf_angstrom": res_rmsf,
            }
        )
    return pd.DataFrame(rows)


@dataclass
class ContactMapResult:
    """Residue-residue contact probability and mean minimum-distance matrices."""

    residues: pd.DataFrame
    frequency: np.ndarray
    mean_distance_angstrom: np.ndarray
    cutoff_angstrom: float
    representation: str
    n_frames: int


def compute_contact_map(
    universe: mda.Universe,
    *,
    select: str | None = None,
    representation: str = "ca",
    cutoff_angstrom: float | None = None,
    pbc: bool = True,
    stride: int = 1,
) -> ContactMapResult:
    """Compute time-averaged residue-residue contact frequency matrix in a single trajectory pass.

    Parameters
    ----------
    representation : {"ca", "heavy"}
        `"ca"` uses C-alpha distances (default cutoff 8.0 Å).
        `"heavy"` uses minimum heavy-atom distance between every residue pair (default cutoff 4.5 Å).
    pbc : bool, default True
        Use minimum-image periodic boundary conditions when unit-cell dimensions are present.
    """
    rep = representation.lower().strip()
    if rep not in ("ca", "heavy"):
        raise ValueError(f"Unsupported contact representation {representation!r}; expected 'ca' or 'heavy'")
    if cutoff_angstrom is None:
        cutoff_angstrom = 8.0 if rep == "ca" else 4.5

    sel = select or protein_selection_str("ca" if rep == "ca" else "heavy")
    ag = universe.select_atoms(sel)
    if len(ag) == 0:
        raise ValueError(f"Contact-map selection matched 0 atoms: {sel}")

    residues = ag.residues
    n_res = len(residues)
    ptm_set = get_ptmpsi_residue_sets()["ptm"]
    res_df = pd.DataFrame(
        {
            "res_pos": np.arange(n_res, dtype=int),
            "resid": [int(r.resid) for r in residues],
            "resname": [str(r.resname) for r in residues],
            "segid": [str(r.segid) for r in residues],
            "label": [f"{r.segid}:{r.resname}{r.resid}" if r.segid else f"{r.resname}{r.resid}" for r in residues],
            "is_ptm": [r.resname in ptm_set for r in residues],
        }
    )

    counts = np.zeros((n_res, n_res), dtype=np.float64)
    dist_sum = np.zeros((n_res, n_res), dtype=np.float64)
    n_frames = 0

    # Map each selected atom to its local 0..n_res-1 index
    resindex_to_local = {int(r.resindex): idx for idx, r in enumerate(residues)}
    atom_to_res = np.array([resindex_to_local[int(ri)] for ri in ag.resindices], dtype=np.intp)

    for ts in universe.trajectory[::stride]:
        box = ts.dimensions if (pbc and ts.dimensions is not None and np.all(ts.dimensions[:3] > 0)) else None
        pos = ag.positions
        if rep == "ca" and len(ag) == n_res:
            dmat = mda_distances.distance_array(pos, pos, box=box)
        else:
            atom_dmat = mda_distances.distance_array(pos, pos, box=box)
            # Reduce atom-atom distances to minimum per residue pair via vectorized reduceat
            order = np.argsort(atom_to_res, kind="stable")
            sorted_res = atom_to_res[order]
            splits = np.flatnonzero(np.r_[True, sorted_res[1:] != sorted_res[:-1]])
            d_sorted = atom_dmat[np.ix_(order, order)]
            row_min = np.minimum.reduceat(d_sorted, splits, axis=0)
            col_min = np.minimum.reduceat(row_min, splits, axis=1)
            present_res = sorted_res[splits]
            if len(present_res) == n_res:
                dmat = col_min
            else:
                dmat = np.full((n_res, n_res), np.inf, dtype=np.float64)
                dmat[np.ix_(present_res, present_res)] = col_min

        np.fill_diagonal(dmat, 0.0)
        counts += (dmat < cutoff_angstrom).astype(np.float64)
        dist_sum += dmat
        n_frames += 1

    if n_frames == 0:
        raise ValueError("Trajectory slice yielded 0 frames for contact map")

    return ContactMapResult(
        residues=res_df,
        frequency=counts / n_frames,
        mean_distance_angstrom=dist_sum / n_frames,
        cutoff_angstrom=float(cutoff_angstrom),
        representation=rep,
        n_frames=n_frames,
    )


def contact_difference_map(
    target: ContactMapResult,
    reference: ContactMapResult,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute `(target - reference)` contact-probability difference aligned by `(segid, resid)`.

    Aligning by `(segid, resid)` rather than `resname` allows direct comparison
    between Wild-Type (`CYS54`) and PTM-modified (`SNC54` / `CGL54`) runs.
    """
    ref_keys = {(row.segid, int(row.resid)): int(row.res_pos) for row in reference.residues.itertuples()}
    tgt_keys = {(row.segid, int(row.resid)): int(row.res_pos) for row in target.residues.itertuples()}
    shared = [k for k in tgt_keys if k in ref_keys]
    if not shared:
        raise ValueError("Target and reference contact maps share no (segid, resid) residue positions")

    tgt_idx = [tgt_keys[k] for k in shared]
    ref_idx = [ref_keys[k] for k in shared]
    diff = target.frequency[np.ix_(tgt_idx, tgt_idx)] - reference.frequency[np.ix_(ref_idx, ref_idx)]

    meta = target.residues.iloc[tgt_idx].reset_index(drop=True).copy()
    meta["ref_resname"] = reference.residues.iloc[ref_idx]["resname"].to_numpy()
    meta["res_pos"] = np.arange(len(shared), dtype=int)
    return meta, diff


@dataclass
class CartesianPCAResult:
    """Results of Cartesian coordinate PCA (single-trajectory or joint multi-system)."""

    eigenvalues: np.ndarray
    explained_variance_ratio: np.ndarray
    cumulative_variance_ratio: np.ndarray
    components: np.ndarray
    mean_coordinates: np.ndarray
    atom_weights: np.ndarray
    atoms: pd.DataFrame
    projections: pd.DataFrame
    residue_loadings: pd.DataFrame
    mass_weighted: bool


def _apply_deterministic_pc_signs(
    components: np.ndarray,
    projections: np.ndarray,
    ref_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Make eigenvector signs reproducible across runs and platforms."""
    comps = components.copy()
    projs = projections.copy()
    for k in range(comps.shape[1]):
        flip = False
        if ref_mask is not None and np.any(ref_mask) and not np.all(ref_mask):
            ref_mean = float(projs[ref_mask, k].mean())
            other_mean = float(projs[~ref_mask, k].mean())
            if abs(other_mean - ref_mean) > 1e-9:
                # Orient axis from reference (< 0 or smaller) toward modified (> ref_mean)
                flip = ref_mean > other_mean
            else:
                max_idx = int(np.argmax(np.abs(comps[:, k])))
                flip = comps[max_idx, k] < 0.0
        else:
            max_idx = int(np.argmax(np.abs(comps[:, k])))
            flip = comps[max_idx, k] < 0.0
        if flip:
            comps[:, k] *= -1.0
            projs[:, k] *= -1.0
    return comps, projs


def _run_coordinate_pca_from_arrays(
    coords_by_system: Sequence[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
    atom_meta: pd.DataFrame,
    masses: np.ndarray,
    *,
    n_components: int = 5,
    mass_weighted: bool = True,
    ref_system: str | None = None,
) -> CartesianPCAResult:
    """Run Cartesian PCA via `MDAnalysis.analysis.pca.PCA` on aligned coordinates.

    When `mass_weighted=True`, coordinates are weighted by $\\sqrt{m_i}$ after
    mass-weighted Kabsch superposition, yielding the exact mass-weighted covariance
    eigenproblem $M^{1/2}\\,\\mathrm{Cov}(\\mathbf{r})\\,M^{1/2}$ (identical to `gmx covar -mwa`,
    with eigenvalues in $\\mathrm{amu}\\cdot\\mathrm{\\AA}^2$ and projections in $\\sqrt{\\mathrm{amu}}\\,\\mathrm{\\AA}$).
    """
    all_coords = np.concatenate([item[1] for item in coords_by_system], axis=0)
    n_total_frames, n_atoms, _ = all_coords.shape
    if n_total_frames < 2:
        raise ValueError("Cartesian PCA requires at least 2 total frames")

    weights = np.sqrt(masses) if mass_weighted else np.ones(n_atoms, dtype=np.float64)
    mean_coords = all_coords.mean(axis=0)
    scaled_coords = all_coords * weights[None, :, None]

    n_dim = n_atoms * 3
    k = min(int(n_components), n_dim, max(n_total_frames - 1, 1))

    # Exact real-symmetric covariance / SVD solver: MDAnalysis.analysis.pca.PCA calls
    # non-symmetric np.linalg.eig, which emits ComplexWarning on the 6 zero rigid-body modes
    # created by rotational/translational superposition. Symmetric eigh (small 3N) and SVD
    # (large 3N) match MDAnalysis's exact sample covariance (ddof=1) without complex roundoff.
    centered = scaled_coords.reshape(n_total_frames, n_dim) - (mean_coords * weights[:, None]).ravel()
    if n_dim <= 1500 and n_total_frames > n_dim:
        cov = (centered.T @ centered) / (n_total_frames - 1)
        evals_asc, evecs_asc = np.linalg.eigh(cov)
        all_evals = np.maximum(evals_asc[::-1], 0.0)
        total_var = float(all_evals.sum())
        eigenvalues = all_evals[:k]
        ev_ratio = eigenvalues / total_var if total_var > 0 else np.zeros(k)
        cum_var = np.cumsum(ev_ratio)
        components = evecs_asc[:, ::-1][:, :k]
        projections = centered @ components
    else:
        _, svals, vt = np.linalg.svd(centered, full_matrices=False)
        all_evals = (svals ** 2) / (n_total_frames - 1)
        total_var = float(all_evals.sum())
        eigenvalues = all_evals[:k]
        ev_ratio = eigenvalues / total_var if total_var > 0 else np.zeros(k)
        cum_var = np.cumsum(ev_ratio)
        components = vt[:k].T
        projections = centered @ components

    sys_labels: list[str] = []
    frame_ids: list[int] = []
    times_ps: list[float] = []
    for sys_name, _, f_arr, t_arr in coords_by_system:
        sys_labels.extend([sys_name] * len(f_arr))
        frame_ids.extend([int(x) for x in f_arr])
        times_ps.extend([float(x) for x in t_arr])

    ref_mask = np.array([s == ref_system for s in sys_labels], dtype=bool) if ref_system else None
    components, projections = _apply_deterministic_pc_signs(components, projections, ref_mask=ref_mask)

    proj_data: dict[str, object] = {
        "system": sys_labels,
        "frame": frame_ids,
        "time_ps": times_ps,
    }
    for idx in range(k):
        proj_data[f"PC{idx + 1}"] = projections[:, idx]
    proj_df = pd.DataFrame(proj_data)

    # Per-residue root-mean-square eigenvector loading for each PC
    atom_comp_3d = components.reshape(n_atoms, 3, k)
    atom_loading_norm = np.linalg.norm(atom_comp_3d, axis=1)  # (n_atoms, k)
    res_rows: list[dict[str, object]] = []
    for (segid, resid, resname, is_ptm), grp in atom_meta.groupby(
        ["segid", "resid", "resname", "is_ptm"], sort=False
    ):
        a_idx = grp["atom_pos"].to_numpy(dtype=int)
        row: dict[str, object] = {
            "segid": segid,
            "resid": int(resid),
            "resname": resname,
            "label": f"{segid}:{resname}{resid}" if segid else f"{resname}{resid}",
            "is_ptm": bool(is_ptm),
        }
        for idx in range(k):
            row[f"PC{idx + 1}_loading"] = float(np.sqrt(np.mean(atom_loading_norm[a_idx, idx] ** 2)))
        res_rows.append(row)

    return CartesianPCAResult(
        eigenvalues=eigenvalues,
        explained_variance_ratio=ev_ratio,
        cumulative_variance_ratio=cum_var,
        components=components,
        mean_coordinates=mean_coords,
        atom_weights=masses if mass_weighted else np.ones_like(masses),
        atoms=atom_meta,
        projections=proj_df,
        residue_loadings=pd.DataFrame(res_rows),
        mass_weighted=mass_weighted,
    )


def compute_cartesian_pca(
    universe: mda.Universe,
    *,
    select: str | None = None,
    align_select: str | None = None,
    n_components: int = 5,
    mass_weighted: bool = True,
    system_name: str = "system",
) -> CartesianPCAResult:
    """Perform aligned (mass-weighted by default) Cartesian PCA on a single trajectory."""
    calc_sel = select or protein_selection_str("backbone")
    fit_sel = align_select or calc_sel

    u_aligned, _ = _aligned_in_memory_copy(
        universe,
        align_select=fit_sel,
        iterative_average=True,
        mass_weighted=mass_weighted,
    )
    ag = u_aligned.select_atoms(calc_sel)
    if len(ag) == 0:
        raise ValueError(f"PCA selection matched 0 atoms: {calc_sel}")
    masses = _ensure_positive_masses(ag, "Mass-weighted Cartesian PCA") if mass_weighted else np.ones(len(ag))

    coords = np.empty((len(u_aligned.trajectory), len(ag), 3), dtype=np.float64)
    frames = np.empty(len(u_aligned.trajectory), dtype=int)
    times = np.empty(len(u_aligned.trajectory), dtype=np.float64)
    for i, ts in enumerate(u_aligned.trajectory):
        coords[i] = ag.positions
        frames[i] = int(ts.frame)
        times[i] = float(ts.time)

    ptm_set = get_ptmpsi_residue_sets()["ptm"]
    atom_meta = pd.DataFrame(
        {
            "atom_pos": np.arange(len(ag), dtype=int),
            "atom_name": list(ag.names),
            "resid": [int(x) for x in ag.resids],
            "resname": list(ag.resnames),
            "segid": list(ag.segids),
            "mass": masses,
            "is_ptm": [r in ptm_set for r in ag.resnames],
        }
    )

    return _run_coordinate_pca_from_arrays(
        [(system_name, coords, frames, times)],
        atom_meta,
        masses,
        n_components=n_components,
        mass_weighted=mass_weighted,
    )


def compute_joint_cartesian_pca(
    systems: Mapping[str, mda.Universe | tuple[str | Path, str | Path]],
    *,
    ref_system: str | None = None,
    atom_kind: str = "backbone",
    n_components: int = 5,
    mass_weighted: bool = True,
) -> CartesianPCAResult:
    """Perform joint mass-weighted Cartesian PCA across Wild-Type and PTM trajectories.

    Matches the shared polypeptide backbone (`backbone` or `ca`) across all systems
    by `(segid, resid, atom_name)` so side-chain atom-count differences (e.g. `CYS`
    vs `SNC` vs `CGL`) never mismatch feature dimensions or drop modified sites.
    """
    if len(systems) < 2:
        raise ValueError("compute_joint_cartesian_pca expects at least 2 systems (e.g., {'WT': ..., 'PTM': ...})")
    if atom_kind not in ("backbone", "ca"):
        raise ValueError("Joint WT-vs-PTM Cartesian PCA requires atom_kind='backbone' or 'ca' so atom counts match")

    universes: dict[str, mda.Universe] = {}
    for name, item in systems.items():
        if isinstance(item, mda.Universe):
            universes[name] = item
        else:
            top, traj = item
            universes[name] = load_universe(top, traj)

    ref_key = ref_system or next(iter(universes))
    if ref_key not in universes:
        raise ValueError(f"ref_system {ref_key!r} not in systems {list(universes)}")

    # Identify common (segid, resid, atom_name) keys across all systems, including PTM sites
    key_maps: dict[str, dict[tuple[str, int, str], int]] = {}
    ordered_keys: list[tuple[str, int, str]] = []
    for idx, (name, u) in enumerate(universes.items()):
        ag = select_protein(u, kind=atom_kind, include_ptms=True)
        if len(ag) == 0:
            raise ValueError(f"System {name!r} matched 0 {atom_kind} atoms")
        kmap: dict[tuple[str, int, str], int] = {}
        for atom in ag:
            key = (str(atom.segid), int(atom.resid), str(atom.name))
            kmap[key] = int(atom.index)
            if idx == 0:
                ordered_keys.append(key)
        key_maps[name] = kmap

    common_keys = [k for k in ordered_keys if all(k in key_maps[name] for name in universes)]
    if not common_keys:
        raise ValueError("Systems share no common (segid, resid, atom_name) backbone atoms")

    ref_u = universes[ref_key]
    ref_indices = [key_maps[ref_key][k] for k in common_keys]
    ref_ag = ref_u.atoms[ref_indices]
    masses = _ensure_positive_masses(ref_ag, "Joint mass-weighted Cartesian PCA") if mass_weighted else np.ones(len(ref_ag))

    # Compute iterative reference structure on ref_system's common atoms
    ref_u.trajectory[0]
    ref_positions = ref_ag.positions.copy()
    ref_com = np.average(ref_positions, axis=0, weights=masses)
    ref_centered = ref_positions - ref_com

    coords_by_system: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []
    ptm_set = get_ptmpsi_residue_sets()["ptm"]
    any_ptm_at_key = {k: False for k in common_keys}

    for name, u in universes.items():
        idxs = [key_maps[name][k] for k in common_keys]
        ag = u.atoms[idxs]
        for k, atom in zip(common_keys, ag):
            if atom.resname in ptm_set:
                any_ptm_at_key[k] = True

        n_frames = len(u.trajectory)
        coords = np.empty((n_frames, len(common_keys), 3), dtype=np.float64)
        frames = np.empty(n_frames, dtype=int)
        times = np.empty(n_frames, dtype=np.float64)

        for i, ts in enumerate(u.trajectory):
            pos = ag.positions.copy()
            mob_com = np.average(pos, axis=0, weights=masses)
            mob_centered = pos - mob_com
            rot, _ = align.rotation_matrix(mob_centered, ref_centered, weights=masses)
            coords[i] = mob_centered @ rot.T + ref_com
            frames[i] = int(ts.frame)
            times[i] = float(ts.time)

        coords_by_system.append((name, coords, frames, times))

    atom_meta = pd.DataFrame(
        {
            "atom_pos": np.arange(len(common_keys), dtype=int),
            "atom_name": [k[2] for k in common_keys],
            "resid": [k[1] for k in common_keys],
            "resname": list(ref_ag.resnames),
            "segid": [k[0] for k in common_keys],
            "mass": masses,
            "is_ptm": [any_ptm_at_key[k] for k in common_keys],
        }
    )

    return _run_coordinate_pca_from_arrays(
        coords_by_system,
        atom_meta,
        masses,
        n_components=n_components,
        mass_weighted=mass_weighted,
        ref_system=ref_key,
    )


def save_standard_report(
    topology: str | Path,
    trajectory: str | Path,
    output_dir: str | Path,
    *,
    prefix: str = "analysis",
    contact_rep: str = "ca",
    n_pca_components: int = 3,
    make_plots: bool = True,
) -> dict[str, Path]:
    """Run all five starter analyses on `(topology, trajectory)` and write CSV/NPY (+ PNG) outputs."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    u = load_universe(topology, trajectory)

    rmsd_df = compute_rmsd(u)
    rog_df = compute_rog(u)
    rmsf_df = compute_rmsf(u, by_residue=True)
    cmap = compute_contact_map(u, representation=contact_rep)
    pca_res = compute_cartesian_pca(u, n_components=n_pca_components)

    written: dict[str, Path] = {}
    for key, df in (
        ("rmsd", rmsd_df),
        ("rog", rog_df),
        ("rmsf", rmsf_df),
        ("contact_residues", cmap.residues),
        ("pca_projections", pca_res.projections),
        ("pca_residue_loadings", pca_res.residue_loadings),
    ):
        p = out / f"{prefix}_{key}.csv"
        df.to_csv(p, index=False)
        written[key] = p

    cmap_path = out / f"{prefix}_contact_frequency.npy"
    np.save(cmap_path, cmap.frequency)
    written["contact_frequency"] = cmap_path

    pca_var_df = pd.DataFrame(
        {
            "component": [f"PC{i + 1}" for i in range(len(pca_res.eigenvalues))],
            "eigenvalue": pca_res.eigenvalues,
            "explained_variance_ratio": pca_res.explained_variance_ratio,
            "cumulative_variance_ratio": pca_res.cumulative_variance_ratio,
        }
    )
    pca_var_path = out / f"{prefix}_pca_variance.csv"
    pca_var_df.to_csv(pca_var_path, index=False)
    written["pca_variance"] = pca_var_path

    if make_plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
        axes[0, 0].plot(rmsd_df["time_ps"], rmsd_df["rmsd_angstrom"], lw=1.4)
        axes[0, 0].set_xlabel("Time (ps)")
        axes[0, 0].set_ylabel("Backbone RMSD (Å)")
        axes[0, 0].set_title("Aligned Backbone RMSD")

        axes[0, 1].plot(rog_df["time_ps"], rog_df["rg_angstrom"], lw=1.4, color="tab:orange")
        axes[0, 1].set_xlabel("Time (ps)")
        axes[0, 1].set_ylabel("Radius of Gyration (Å)")
        axes[0, 1].set_title("Mass-Weighted Radius of Gyration")

        axes[1, 0].plot(rmsf_df["resid"], rmsf_df["rmsf_angstrom"], lw=1.4, color="tab:green")
        ptm_rows = rmsf_df[rmsf_df["is_ptm"]]
        if not ptm_rows.empty:
            axes[1, 0].scatter(ptm_rows["resid"], ptm_rows["rmsf_angstrom"], color="crimson", zorder=4, label="PTM")
            axes[1, 0].legend(frameon=False)
        axes[1, 0].set_xlabel("Residue ID")
        axes[1, 0].set_ylabel("Residue RMSF (Å)")
        axes[1, 0].set_title("Iteratively Aligned Per-Residue RMSF")

        im = axes[1, 1].imshow(cmap.frequency, origin="lower", vmin=0.0, vmax=1.0, cmap="viridis")
        fig.colorbar(im, ax=axes[1, 1], fraction=0.046, pad=0.04, label="Contact probability")
        axes[1, 1].set_xlabel("Residue index")
        axes[1, 1].set_ylabel("Residue index")
        axes[1, 1].set_title(f"Contact Map ({cmap.representation.upper()}, < {cmap.cutoff_angstrom:g} Å)")

        fig_path = out / f"{prefix}_summary.png"
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        written["summary_plot"] = fig_path

    return written
