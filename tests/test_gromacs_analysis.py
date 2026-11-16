"""ptmpsi.gromacs.analysis against a real GROMACS ACE-SNC-ALA-NME system."""

from __future__ import annotations

import os
import shutil

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import pytest

from AnalysisTools import check_pbc_jumps, load_universe
from ptmpsi.gromacs import analysis as an


def _xvg(path):
    return np.array(
        [[float(x) for x in line.split()] for line in open(path) if line.strip() and not line.startswith(("#", "@"))]
    )


@pytest.fixture(scope="module")
def fitted(snc_system, tmp_path_factory):
    out = tmp_path_factory.mktemp("analysis")
    xtc = an.post_process(
        snc_system["raw_xtc"], snc_system["tpr"], prefix=out / "pp", group="Protein", fit_group=4, gmx=snc_system["gmx"]
    )[0]
    return {"xtc": xtc, "tpr": str(snc_system["tpr"]), "gmx": snc_system["gmx"], "out": out}


def test_find_gmx_prefers_gmx_mpi(tmp_path, monkeypatch):
    for name in ("gmx", "gmx_mpi"):
        exe = tmp_path / name
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert an.find_gmx() == str(tmp_path / "gmx_mpi")
    assert an.find_gmx("gmx") == str(tmp_path / "gmx")
    (tmp_path / "gmx_mpi").unlink()
    assert an.find_gmx() == str(tmp_path / "gmx")
    (tmp_path / "gmx").unlink()
    with pytest.raises(FileNotFoundError):
        an.find_gmx()


def test_interactive_false_is_not_silently_ignored(tmp_path):
    with pytest.raises(NotImplementedError, match="interactive=True"):
        an.rmsd("a.xtc", "a.tpr", tmp_path, interactive=False)


def test_add_residue_file_uses_packaged_copy(tmp_path):
    target = an.add_residue_file(None, destination=tmp_path)
    assert os.path.basename(target) == "residuetypes.dat"
    assert "SNC" in open(target).read()


def test_process_dssp_names_outputs_after_the_dat_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # "spade" starts and ends with characters that str.strip("_dssp.dat") used to remove
    (tmp_path / "spade_dssp.dat").write_text("~HHHH~\n~HHH~~\n")
    an.process_dssp(["spade_dssp.dat"], ["spade.xtc"])
    assert (tmp_path / "spade_summary.dat").is_file()
    assert (tmp_path / "spade_result.out").is_file()
    assert not (tmp_path / "e_summary.dat").exists()


@pytest.mark.skipif(shutil.which("gmx") is None and shutil.which("gmx_mpi") is None, reason="GROMACS CLI not installed")
def test_post_process_removes_box_jumps(snc_system, fitted):
    assert check_pbc_jumps(load_universe(snc_system["tpr"], snc_system["raw_xtc"]))["has_pbc_jumps"]
    assert not check_pbc_jumps(load_universe(snc_system["tpr"], fitted["xtc"]))["has_pbc_jumps"]
    assert os.path.isfile(str(fitted["xtc"]).replace("_whole_nojump_fit.xtc", ".gro"))


def test_rmsd_matches_mdanalysis_against_tpr_structure(fitted):
    (xvg,) = an.rmsd(fitted["xtc"], fitted["tpr"], fitted["out"] / "rmsd", gmx=fitted["gmx"])
    gmx_rmsd = _xvg(xvg)[:, 1] * 10.0
    u = mda.Universe(fitted["tpr"], fitted["xtc"])
    ref = mda.Universe(fitted["tpr"])
    mine = rms.RMSD(u, ref, select="name N CA C", weights="mass").run().results.rmsd[:, 2]
    assert gmx_rmsd == pytest.approx(mine, abs=1e-4)
    assert (fitted["out"] / "rmsd" / "DELSPEC_rmsd.out").is_file()


def test_gyration_matches_mdanalysis(fitted):
    (xvg,) = an.gyration(fitted["xtc"], fitted["tpr"], fitted["out"] / "gyr", gmx=fitted["gmx"])
    gmx_rg = _xvg(xvg)[:, 1] * 10.0
    u = mda.Universe(fitted["tpr"], fitted["xtc"])
    bb = u.select_atoms("name N CA C")
    mine = np.array([bb.radius_of_gyration() for _ in u.trajectory])
    assert gmx_rg == pytest.approx(mine, abs=1e-3)


def test_rmsf_per_residue_matches_numpy(fitted):
    ((xvg, pdb),) = an.rmsf(fitted["xtc"], fitted["tpr"], fitted["out"] / "rmsf", gmx=fitted["gmx"])
    assert os.path.isfile(pdb)
    gmx_rmsf = _xvg(xvg)[:, 1] * 10.0
    u = mda.Universe(fitted["tpr"], fitted["xtc"])
    bb = u.select_atoms("name N CA C")
    coords = np.array([bb.positions.copy() for _ in u.trajectory])
    per_atom = np.sqrt(((coords - coords.mean(0)) ** 2).sum(-1).mean(0))
    per_res = [np.sqrt((per_atom[bb.resindices == r] ** 2).mean()) for r in sorted(set(bb.resindices))]
    assert len(gmx_rmsf) == len(per_res)
    # gmx fits on the atoms it averages over, so allow a few percent
    assert gmx_rmsf == pytest.approx(per_res, rel=0.05)


def test_dssp_and_contact_map(fitted):
    (dat,) = an.dssp(fitted["xtc"], fitted["tpr"], fitted["out"] / "dssp", gmx=fitted["gmx"])
    assert os.path.isfile(dat)
    (matrix,) = an.contact_map(
        fitted["xtc"], fitted["tpr"], None, fitted["out"] / "cm", dt=None, selection='group "Protein"', gmx=fitted["gmx"]
    )
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 1.0)


def test_process_contact_uses_fraction_of_frames_below_cutoff(tmp_path):
    xvg = tmp_path / "x_dist.xvg"
    # time, then a 2x2 pair matrix: the (0, 1) distance is 0.4 nm in 3 of 4 frames
    xvg.write_text(
        "# comment\n@ title\n"
        "0 0.0 0.4 0.4 0.0\n2 0.0 0.4 0.4 0.0\n4 0.0 0.4 0.4 0.0\n6 0.0 0.9 0.9 0.0\n"
    )
    (matrix,) = an.process_contact([str(xvg)], contact_cutoff=0.5)
    assert matrix == pytest.approx(np.array([[1.0, 0.75], [0.75, 1.0]]))
    assert (tmp_path / "x_dist_fraction.dat").is_file()


def test_unknown_keyword_arguments_are_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(an, "find_gmx", lambda gmx=None: "gmx")
    for name in ("rmsd", "rmsf", "gyration", "dssp"):
        with pytest.raises(TypeError, match="bogus"):
            getattr(an, name)(__file__, __file__, tmp_path, bogus=1)
    with pytest.raises(NotImplementedError, match="selection"):
        an.contact_map(__file__, __file__, None, tmp_path, start_residue=3)
