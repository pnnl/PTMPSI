from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import pytest

from AnalysisTools import (
    RESIDUETYPES_DAT,
    check_pbc_jumps,
    compute_cartesian_pca,
    compute_contact_map,
    compute_joint_cartesian_pca,
    compute_rmsd,
    compute_rmsf,
    compute_rog,
    contact_difference_map,
    find_ptm_residues,
    get_ptmpsi_residue_sets,
    load_universe,
    save_standard_report,
    select_protein,
    unwrap_trajectory,
)


def _make_peptide_universe(
    resnames: list[str],
    coords: np.ndarray,
    *,
    extra_sidechain_on_res: dict[int, list[tuple[str, float, np.ndarray]]] | None = None,
    box: np.ndarray | None = None,
) -> mda.Universe:
    """Build an in-memory multi-frame peptide Universe with backbone + optional PTM atoms."""
    backbone_defs = [("N", 14.007), ("CA", 12.011), ("C", 12.011), ("O", 15.999)]
    extra = extra_sidechain_on_res or {}

    atom_names: list[str] = []
    atom_masses: list[float] = []
    atom_resindex: list[int] = []
    n_frames = coords.shape[0]
    per_atom_traj: list[np.ndarray] = []

    bb_cursor = 0
    for r_idx in range(len(resnames)):
        ca_traj = coords[:, bb_cursor + 1, :]
        for a_name, a_mass in backbone_defs:
            atom_names.append(a_name)
            atom_masses.append(a_mass)
            atom_resindex.append(r_idx)
            per_atom_traj.append(coords[:, bb_cursor, :])
            bb_cursor += 1
        for sc_name, sc_mass, offset in extra.get(r_idx, []):
            atom_names.append(sc_name)
            atom_masses.append(sc_mass)
            atom_resindex.append(r_idx)
            per_atom_traj.append(ca_traj + offset[None, :])

    full_coords = np.stack(per_atom_traj, axis=1).astype(np.float32)
    n_atoms = len(atom_names)
    n_res = len(resnames)

    u = mda.Universe.empty(
        n_atoms,
        n_residues=n_res,
        atom_resindex=atom_resindex,
        residue_segindex=[0] * n_res,
        trajectory=True,
    )
    u.add_TopologyAttr("name", atom_names)
    u.add_TopologyAttr("resname", resnames)
    u.add_TopologyAttr("resid", list(range(1, n_res + 1)))
    u.add_TopologyAttr("segid", ["A"])
    u.add_TopologyAttr("mass", atom_masses)

    u.load_new(full_coords, format=MemoryReader, dt=2.0)
    if box is not None:
        for ts in u.trajectory:
            ts.dimensions = box
    return u


def _ideal_backbone_frame(n_res: int) -> np.ndarray:
    """Generate simple realistic peptide backbone coordinates (4 atoms per residue)."""
    pts = []
    for r in range(n_res):
        base_x = r * 3.8
        pts.extend(
            [
                [base_x - 1.2, 0.0, 0.0],   # N
                [base_x, 0.4, 0.0],         # CA
                [base_x + 1.3, 0.0, 0.0],   # C
                [base_x + 1.3, -1.2, 0.0],  # O
            ]
        )
    return np.asarray(pts, dtype=np.float64)


def test_ptmpsi_residue_sets_and_ptm_selections():
    sets = get_ptmpsi_residue_sets()
    for ptm_code in ("SNC", "CSO", "CSD", "OCS", "CSS", "CGL", "SMC", "QCS", "XCN", "IYY", "ALY", "CYZ"):
        assert ptm_code in sets["ptm"]
        assert ptm_code in sets["protein"]
    for std_code in ("ALA", "CYS", "CYX", "LYS"):
        assert std_code not in sets["ptm"]
        assert std_code in sets["protein"]
    assert "NAD" not in sets["protein"]
    assert "NAD" in get_ptmpsi_residue_sets(include_cofactors=True)["protein"]

    base = _ideal_backbone_frame(4)[None, :, :]
    u = _make_peptide_universe(
        ["ALA", "SNC", "CGL", "NAD"],
        base,
        extra_sidechain_on_res={1: [("SG", 32.06, np.array([0.0, 1.8, 0.0])), ("NO", 30.01, np.array([0.0, 3.1, 0.0]))]},
    )

    # MDAnalysis built-in 'protein' / 'backbone' drops SNC and CGL; select_protein keeps them and excludes NAD!
    assert set(u.select_atoms("protein").resnames) == {"ALA"}
    prot_all = select_protein(u, kind="all")
    assert list(prot_all.residues.resnames) == ["ALA", "SNC", "CGL"]
    prot_bb = select_protein(u, kind="backbone")
    assert len(prot_bb) == 3 * 4
    prot_ca = select_protein(u, kind="ca")
    assert len(prot_ca) == 3

    ptm_res = find_ptm_residues(u)
    assert list(ptm_res.resnames) == ["SNC", "CGL"]


def test_pbc_jump_check_and_rmsd_rog_rmsf():
    n_res = 6
    n_frames = 12
    f0 = _ideal_backbone_frame(n_res)
    rng = np.random.default_rng(42)

    # Rigid body rotations + translations on residues 0..4, genuine motion only on PTM residue 2 (SNC)
    coords = np.zeros((n_frames, n_res * 4, 3), dtype=np.float64)
    for t in range(n_frames):
        frame = f0.copy()
        # Oscillate SNC (residue index 2, atoms 8..11) along z
        frame[8:12, 2] += 1.5 * np.sin(2.0 * np.pi * t / n_frames)
        # Apply a large rigid rotation + translation that alignment must remove
        theta = 0.25 * t
        rot = np.array(
            [[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]]
        )
        coords[t] = frame @ rot.T + np.array([3.0 * t, -2.0 * t, 0.5 * t])

    u = _make_peptide_universe(
        ["ALA", "GLY", "SNC", "VAL", "LEU", "ALA"],
        coords,
        extra_sidechain_on_res={2: [("SG", 80.0, np.array([0.0, 2.5, 0.0]))]},
    )

    # 1. RMSD: frame 0 vs frame 0 is 0; rigid translation/rotation removed
    rmsd_df = compute_rmsd(u)
    assert rmsd_df.shape == (n_frames, 3)
    assert rmsd_df["rmsd_angstrom"].iloc[0] == pytest.approx(0.0, abs=1e-5)
    assert rmsd_df["rmsd_angstrom"].max() < 1.0

    # 2. Mass-weighted RoG: verify exact match with analytical mass-weighted formula and difference from unweighted Rg
    rog_df = compute_rog(u, include_components=True)
    u.trajectory[0]
    ag = select_protein(u, "all")
    masses = ag.masses
    com = np.average(ag.positions, axis=0, weights=masses)
    expected_mw_rg = np.sqrt(np.sum(masses * np.sum((ag.positions - com) ** 2, axis=1)) / np.sum(masses))
    unweighted_cog = ag.positions.mean(axis=0)
    unweighted_rg = np.sqrt(np.mean(np.sum((ag.positions - unweighted_cog) ** 2, axis=1)))
    assert rog_df["rg_angstrom"].iloc[0] == pytest.approx(expected_mw_rg, rel=1e-6)
    assert abs(expected_mw_rg - unweighted_rg) > 0.05

    # Missing/zero masses must fail loudly instead of falling back to unweighted geometry
    u_zero = _make_peptide_universe(["ALA", "SNC"], coords[:3, :8, :])
    u_zero.atoms.masses = 0.0
    with pytest.raises(ValueError, match="masses"):
        compute_rog(u_zero)

    # 3. Iteratively aligned RMSF: rigid residues ~0, mobile PTM residue 3 ('SNC') clearly elevated
    rmsf_df = compute_rmsf(u, by_residue=True)
    assert len(rmsf_df) == n_res
    assert bool(rmsf_df.loc[rmsf_df["resname"] == "SNC", "is_ptm"].iloc[0]) is True
    snc_rmsf = float(rmsf_df.loc[rmsf_df["resname"] == "SNC", "rmsf_angstrom"].iloc[0])
    rigid_max = float(rmsf_df.loc[rmsf_df["resname"] != "SNC", "rmsf_angstrom"].max())
    assert snc_rmsf > 0.75
    assert snc_rmsf > 3.0 * rigid_max
    assert rigid_max < 0.30

    # 4. PBC jump detector catches an artificial 50 Å box jump
    wrapped_coords = coords.copy()
    wrapped_coords[5, 12:, 0] += 50.0
    u_wrapped = _make_peptide_universe(["ALA", "GLY", "SNC", "VAL", "LEU", "ALA"], wrapped_coords)
    report = check_pbc_jumps(u_wrapped)
    assert report["has_pbc_jumps"] is True
    assert 5 in report["split_chain_frames"]


def test_minimum_image_contact_map_and_wt_ptm_difference():
    # Place residue 1 near x=1.0 Å and residue 4 near x=49.0 Å in a 50 Å orthorhombic box:
    # raw distance is 48 Å, minimum-image PBC distance is 2.0 Å (< 8.0 Å cutoff)!
    f = _ideal_backbone_frame(4)
    f[12:16, 0] = f[0:4, 0] + 48.0
    coords_wt = np.repeat(f[None, :, :], 4, axis=0)
    box = np.array([50.0, 50.0, 50.0, 90.0, 90.0, 90.0], dtype=np.float32)
    u_wt = _make_peptide_universe(["ALA", "CYS", "GLY", "VAL"], coords_wt, box=box)

    cmap_pbc = compute_contact_map(u_wt, representation="ca", cutoff_angstrom=8.0, pbc=True)
    cmap_nopbc = compute_contact_map(u_wt, representation="ca", cutoff_angstrom=8.0, pbc=False)
    assert cmap_pbc.frequency[0, 3] == pytest.approx(1.0)
    assert cmap_nopbc.frequency[0, 3] == pytest.approx(0.0)

    # PTM system moves residue 4 away so contact (0, 3) breaks, and CYS2 -> CGL2
    coords_ptm = coords_wt.copy()
    coords_ptm[:, 12:16, 1] += 20.0
    u_ptm = _make_peptide_universe(
        ["ALA", "CGL", "GLY", "VAL"],
        coords_ptm,
        extra_sidechain_on_res={1: [("SG", 32.06, np.array([0.0, 1.8, 0.0])), ("G1", 12.01, np.array([0.0, 3.2, 0.0]))]},
        box=box,
    )
    cmap_ptm = compute_contact_map(u_ptm, representation="ca", cutoff_angstrom=8.0, pbc=True)
    meta, diff = contact_difference_map(cmap_ptm, cmap_pbc)
    assert list(meta["resname"]) == ["ALA", "CGL", "GLY", "VAL"]
    assert list(meta["ref_resname"]) == ["ALA", "CYS", "GLY", "VAL"]
    assert diff[0, 3] == pytest.approx(-1.0)


def test_mass_weighted_and_joint_wt_ptm_cartesian_pca():
    n_res = 5
    n_frames = 20
    f0 = _ideal_backbone_frame(n_res)
    t_grid = np.linspace(0.0, 2.0 * np.pi, n_frames, endpoint=False)

    # Mode 1: hinge motion along z on residue 3; WT centered at z=-1.5 Å, PTM shifted to z=+1.5 Å
    wt_coords = np.zeros((n_frames, n_res * 4, 3), dtype=np.float64)
    ptm_coords = np.zeros((n_frames, n_res * 4, 3), dtype=np.float64)
    for i, t in enumerate(t_grid):
        w = f0.copy()
        p = f0.copy()
        w[8:12, 2] += -1.5 + 0.6 * np.sin(t)
        p[8:12, 2] += +1.5 + 0.6 * np.cos(t)
        wt_coords[i] = w
        ptm_coords[i] = p

    u_wt = _make_peptide_universe(["ALA", "GLY", "CYS", "VAL", "ALA"], wt_coords)
    # CGL has 4 extra bulky sidechain atoms vs CYS -> different total atom counts!
    u_ptm = _make_peptide_universe(
        ["ALA", "GLY", "CGL", "VAL", "ALA"],
        ptm_coords,
        extra_sidechain_on_res={
            2: [
                ("SG", 32.06, np.array([0.0, 1.8, 0.0])),
                ("CB2", 12.011, np.array([0.5, 2.9, 0.0])),
                ("N2", 14.007, np.array([1.2, 3.8, 0.0])),
                ("O2", 15.999, np.array([-0.8, 3.6, 0.0])),
            ]
        },
    )
    assert u_ptm.atoms.n_atoms > u_wt.atoms.n_atoms

    single_pca = compute_cartesian_pca(u_ptm, n_components=3, mass_weighted=True)
    assert single_pca.explained_variance_ratio[0] > 0.95
    # Dominant loading is on residue 3 (CGL)
    top_res = single_pca.residue_loadings.sort_values("PC1_loading", ascending=False).iloc[0]
    assert top_res["resname"] == "CGL"
    assert bool(top_res["is_ptm"]) is True

    joint = compute_joint_cartesian_pca(
        {"WT": u_wt, "CGL3": u_ptm},
        ref_system="WT",
        atom_kind="backbone",
        n_components=3,
        mass_weighted=True,
    )
    assert len(joint.atoms) == n_res * 4
    wt_pc1 = joint.projections.loc[joint.projections["system"] == "WT", "PC1"].mean()
    ptm_pc1 = joint.projections.loc[joint.projections["system"] == "CGL3", "PC1"].mean()
    # Deterministic sign convention orients WT < PTM along separating PC1
    assert wt_pc1 < ptm_pc1
    assert joint.explained_variance_ratio[0] > 0.70
    assert joint.cumulative_variance_ratio[1] > 0.99


@pytest.mark.skipif(shutil.which("gmx") is None and shutil.which("gmx_mpi") is None, reason="GROMACS CLI not installed")
def test_end_to_end_gromacs_snc_unwrap_pca_and_cli(tmp_path: Path):
    """Build a real GROMACS .tpr + .xtc with a PTM-Psi SNC residue and validate unwrap + analyses + CLI."""
    gmx = shutil.which("gmx_mpi") or shutil.which("gmx")
    assert gmx is not None

    # Copy PTM-Psi amber99sb.ff + residuetypes.dat + specbond.dat into tmp_path
    ff_root = RESIDUETYPES_DAT.parent
    shutil.copytree(ff_root / "amber99sb.ff", tmp_path / "amber99sb.ff")
    shutil.copy2(RESIDUETYPES_DAT, tmp_path / "residuetypes.dat")
    shutil.copy2(ff_root / "specbond.dat", tmp_path / "specbond.dat")

    # Build ACE-SNC-ALA-NME using PTM-Psi Protein builder if available, otherwise write minimal heavy-atom PDB
    pdb_in = tmp_path / "snc_pep.pdb"
    try:
        from ptmpsi.protein import Protein
        from ptmpsi.gromacs.utils import amber_to_gromacs_names

        prot = Protein()
        prot.prepend("A", "NME")
        prot.prepend("A", "ALA", phi="alpha")
        prot.prepend("A", "CYS", phi="alpha")
        prot.prepend("A", "ACE")
        prot.modify("A:CYS2", "nitrosylation")
        amber_to_gromacs_names(prot)
        prot.write_pdb(str(pdb_in))
    except Exception:
        # Fallback heavy-atom ACE-SNC-ALA-NME PDB matching amber99sb.ff/aminoacids.rtp atom names
        pdb_in.write_text(
            "ATOM      1  CH3 ACE A   1      -1.200   0.400   0.000  1.00  0.00           C\n"
            "ATOM      2  C   ACE A   1       0.100   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  O   ACE A   1       0.400  -1.200   0.000  1.00  0.00           O\n"
            "ATOM      4  N   SNC A   2       1.100   0.900   0.000  1.00  0.00           N\n"
            "ATOM      5  CA  SNC A   2       2.450   0.450   0.000  1.00  0.00           C\n"
            "ATOM      6  CB  SNC A   2       2.800   0.950   1.400  1.00  0.00           C\n"
            "ATOM      7  SG  SNC A   2       2.700   2.750   1.550  1.00  0.00           S\n"
            "ATOM      8  ND  SNC A   2       3.950   3.250   2.600  1.00  0.00           N\n"
            "ATOM      9  OE  SNC A   2       4.050   4.450   2.800  1.00  0.00           O\n"
            "ATOM     10  C   SNC A   2       3.450   0.950  -1.020  1.00  0.00           C\n"
            "ATOM     11  O   SNC A   2       3.200   2.000  -1.620  1.00  0.00           O\n"
            "ATOM     12  N   ALA A   3       4.600   0.300  -1.200  1.00  0.00           N\n"
            "ATOM     13  CA  ALA A   3       5.650   0.700  -2.120  1.00  0.00           C\n"
            "ATOM     14  CB  ALA A   3       6.850  -0.220  -2.050  1.00  0.00           C\n"
            "ATOM     15  C   ALA A   3       6.100   2.130  -1.850  1.00  0.00           C\n"
            "ATOM     16  O   ALA A   3       6.300   2.920  -2.770  1.00  0.00           O\n"
            "ATOM     17  N   NME A   4       6.250   2.450  -0.570  1.00  0.00           N\n"
            "ATOM     18  CH3 NME A   4       6.700   3.780  -0.180  1.00  0.00           C\n"
            "TER\nEND\n",
            encoding="utf-8",
        )

    subprocess.run(
        [gmx, "pdb2gmx", "-f", str(pdb_in), "-o", "processed.gro", "-ff", "amber99sb", "-water", "tip3p", "-ignh"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [gmx, "editconf", "-f", "processed.gro", "-o", "box.gro", "-c", "-d", "1.2", "-bt", "dodecahedron"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    mdp = tmp_path / "md.mdp"
    mdp.write_text(
        "integrator = steep\nnsteps = 5\nemtol = 1000.0\nemstep = 0.01\n"
        "cutoff-scheme = Verlet\ncoulombtype = cutoff\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n",
        encoding="utf-8",
    )
    tpr_path = tmp_path / "md.tpr"
    subprocess.run(
        [gmx, "grompp", "-f", str(mdp), "-c", "box.gro", "-p", "topol.top", "-o", str(tpr_path), "-maxwarn", "2"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )

    # Create a 10-frame .xtc with smooth thermal motion PLUS a periodic-box wrap jump on frames 5..9
    u_gro = mda.Universe(str(tpr_path), str(tmp_path / "box.gro"))
    base_pos = u_gro.atoms.positions.copy()
    box_dims = u_gro.dimensions.copy()
    rng = np.random.default_rng(7)
    raw_xtc = tmp_path / "raw_wrapped.xtc"
    with mda.Writer(str(raw_xtc), u_gro.atoms.n_atoms) as writer:
        for frame_idx in range(10):
            disp = 0.25 * np.sin(0.5 * frame_idx) + rng.normal(0.0, 0.08, size=base_pos.shape)
            pos = base_pos + disp.astype(np.float32)
            if frame_idx >= 5:
                # Shift whole peptide by one periodic box vector along x to mimic a box crossing jump
                pos[:, 0] += box_dims[0]
            u_gro.atoms.positions = pos
            u_gro.trajectory.ts.time = float(frame_idx * 2.0)
            u_gro.trajectory.ts.dimensions = box_dims
            writer.write(u_gro.atoms)

    u_raw = load_universe(tpr_path, raw_xtc)
    assert check_pbc_jumps(u_raw)["has_pbc_jumps"] is True

    clean_xtc = unwrap_trajectory(tpr_path, raw_xtc, tmp_path / "md_clean.xtc", center_group=None)
    u_clean = load_universe(tpr_path, clean_xtc)
    assert check_pbc_jumps(u_clean)["has_pbc_jumps"] is False
    assert [r.resname for r in find_ptm_residues(u_clean)] == ["SNC"]

    # Compare mass-weighted Cartesian PCA eigenvalues against GROMACS `gmx covar -mwa` on Backbone (group 4)!
    # GROMACS 'Backbone' is N, CA, C (8 atoms across ACE-SNC-ALA-NME when residuetypes.dat is staged).
    subprocess.run(
        [gmx, "covar", "-s", str(tpr_path), "-f", str(clean_xtc), "-o", "eigenval_fit.xvg", "-mwa"],
        cwd=tmp_path,
        input="4\n4\n",
        check=True,
        text=True,
        capture_output=True,
    )
    xvg_vals_nm2 = [
        float(line.split()[1])
        for line in (tmp_path / "eigenval_fit.xvg").read_text().splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    ]
    # 1 nm^2 = 100 Å^2; GROMACS `gmx covar` normalizes by 1/N while MDAnalysis/PCA uses 1/(N-1)
    gmx_trace_ang2 = sum(xvg_vals_nm2) * 100.0
    bb_n_ca_c = compute_cartesian_pca(
        u_clean,
        select="name N CA C",
        align_select="name N CA C",
        n_components=24,
        mass_weighted=True,
    )
    n_frames = len(u_clean.trajectory)
    our_pop_trace_ang2 = float(bb_n_ca_c.eigenvalues.sum()) * (n_frames - 1) / n_frames
    assert our_pop_trace_ang2 == pytest.approx(gmx_trace_ang2, rel=1e-3)

    # Test CLI analyze subcommand and legacy ptmpsi.gromacs.analysis wrappers
    out_dir = tmp_path / "report"
    cli = subprocess.run(
        [sys.executable, "-m", "AnalysisTools", "analyze", "-s", str(tpr_path), "-f", str(clean_xtc), "-o", str(out_dir)],
        capture_output=True,
        text=True,
    )
    assert cli.returncode == 0, cli.stderr
    assert (out_dir / "analysis_summary.png").is_file()
    assert (out_dir / "analysis_pca_variance.csv").is_file()

