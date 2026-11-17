"""ptmpsi.gromacs.analysis against a real GROMACS ACE-SNC-ALA-NME system."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys

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
        snc_system["raw_xtc"], snc_system["tpr"], prefix=out / "pp", group="Protein", fit_group=4, gmx=snc_system["gmx"], interactive=True
    )[0]
    return {"xtc": xtc, "tpr": str(snc_system["tpr"]), "gmx": snc_system["gmx"], "out": out}


def test_find_gmx_prefers_gmx_and_falls_back_to_gmx_mpi(tmp_path, monkeypatch):
    for name in ("gmx", "gmx_mpi"):
        exe = tmp_path / name
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert an.find_gmx() == str(tmp_path / "gmx")
    assert an.find_gmx("gmx_mpi") == str(tmp_path / "gmx_mpi")
    (tmp_path / "gmx").unlink()
    assert an.find_gmx() == str(tmp_path / "gmx_mpi")
    (tmp_path / "gmx_mpi").unlink()
    with pytest.raises(FileNotFoundError):
        an.find_gmx()


def test_interactive_false_still_does_nothing_but_says_so(tmp_path, capsys):
    assert an.rmsd("a.xtc", "a.tpr", tmp_path) is None
    assert "interactive=False" in capsys.readouterr().out
    assert not list(tmp_path.iterdir())


def test_add_residue_file_uses_packaged_copy(tmp_path):
    target = an.add_residue_file(None, destination=tmp_path)
    assert os.path.basename(target) == "residuetypes.dat"
    assert "SNC" in open(target).read()


@pytest.mark.parametrize("failure", ["exit", "missing"])
def test_post_process_propagates_command_and_missing_output_failures(tmp_path, monkeypatch, failure):
    trajectory = tmp_path / "md.xtc"
    structure = tmp_path / "md.tpr"
    trajectory.touch()
    structure.touch()
    monkeypatch.setattr(an, "find_gmx", lambda gmx: "gmx")

    def run(command, **kwargs):
        if failure == "exit":
            return subprocess.CompletedProcess(command, 7, "", "trjconv failed")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(an.subprocess, "run", run)
    match = "exit code 7" if failure == "exit" else "expected output"
    with pytest.raises(RuntimeError, match=match):
        an.post_process(trajectory, structure, prefix=tmp_path / "clean", interactive=True)


@pytest.mark.parametrize("missing_suffix", [".gro", "_whole.xtc", "_whole_nojump.xtc", "_whole_nojump_fit.xtc", "_Protein.tpr"])
def test_post_process_checks_every_output(tmp_path, monkeypatch, missing_suffix):
    trajectory, structure = tmp_path / "md.xtc", tmp_path / "md.tpr"
    trajectory.touch()
    structure.touch()
    monkeypatch.setattr(an, "find_gmx", lambda gmx: "gmx")

    def run(command, cwd, stdin):
        output = Path(command[command.index("-o") + 1])
        if output.name != f"md{missing_suffix}":
            output.touch()

    monkeypatch.setattr(an, "_run", run)
    with pytest.raises(RuntimeError, match="expected output"):
        an.post_process(trajectory, structure, prefix=tmp_path / "clean", group="Protein", interactive=True)


def test_unwrap_cli_exits_nonzero_on_failure(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "AnalysisTools", "unwrap", "-s", str(tmp_path / "md.tpr"),
         "-f", str(tmp_path / "md.xtc"), "-o", str(tmp_path / "clean"), "--gmx", str(tmp_path / "missing-gmx")],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "Unwrap failed:" in result.stderr
    assert "fitted_trajectory:" not in result.stdout


@pytest.mark.parametrize("label", ["0001", "spade"])
def test_process_dssp_names_outputs_after_the_dat_file(tmp_path, monkeypatch, label):
    # str.strip("_dssp.dat") removes any of those characters from both ends: "0001" is unharmed, "spade" is not
    monkeypatch.chdir(tmp_path)
    (tmp_path / f"{label}_dssp.dat").write_text("~HHHH~\n~HHH~~\n")
    an.process_dssp([f"{label}_dssp.dat"], [f"{label}.xtc"])
    assert (tmp_path / f"{label}_summary.dat").is_file()
    assert (tmp_path / f"{label}_result.out").is_file()


@pytest.mark.skipif(shutil.which("gmx") is None and shutil.which("gmx_mpi") is None, reason="GROMACS CLI not installed")
def test_post_process_removes_box_jumps(snc_system, fitted):
    assert check_pbc_jumps(load_universe(snc_system["tpr"], snc_system["raw_xtc"]))["has_pbc_jumps"]
    assert not check_pbc_jumps(load_universe(snc_system["tpr"], fitted["xtc"]))["has_pbc_jumps"]
    assert os.path.isfile(str(fitted["xtc"]).replace("_whole_nojump_fit.xtc", ".gro"))


def test_rmsd_matches_mdanalysis_against_tpr_structure(fitted):
    (xvg,) = an.rmsd(fitted["xtc"], fitted["tpr"], fitted["out"] / "rmsd", gmx=fitted["gmx"], interactive=True)
    gmx_rmsd = _xvg(xvg)[:, 1] * 10.0
    u = mda.Universe(fitted["tpr"], fitted["xtc"])
    ref = mda.Universe(fitted["tpr"])
    mine = rms.RMSD(u, ref, select="name N CA C", weights="mass").run().results.rmsd[:, 2]
    assert gmx_rmsd == pytest.approx(mine, abs=1e-4)
    assert (fitted["out"] / "rmsd" / "DELSPEC_rmsd.out").is_file()


def test_gyration_matches_mdanalysis(fitted):
    (xvg,) = an.gyration(fitted["xtc"], fitted["tpr"], fitted["out"] / "gyr", gmx=fitted["gmx"], interactive=True)
    gmx_rg = _xvg(xvg)[:, 1] * 10.0
    u = mda.Universe(fitted["tpr"], fitted["xtc"])
    bb = u.select_atoms("name N CA C")
    mine = np.array([bb.radius_of_gyration() for _ in u.trajectory])
    assert gmx_rg == pytest.approx(mine, abs=1e-3)


def test_rmsf_per_residue_matches_numpy(fitted):
    ((xvg, pdb),) = an.rmsf(fitted["xtc"], fitted["tpr"], fitted["out"] / "rmsf", gmx=fitted["gmx"], interactive=True)
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


def test_dssp_and_contact_map(fitted, capsys):
    (dat,) = an.dssp(fitted["xtc"], fitted["tpr"], fitted["out"] / "dssp", gmx=fitted["gmx"], interactive=True)
    assert os.path.isfile(dat)
    # group 1 is Protein; start_residue/end_residue/bogus are accepted and ignored, as before
    matrix = an.contact_map(
        fitted["xtc"], fitted["tpr"], None, fitted["out"] / "cm", dt=2, selection="1",
        gmx=fitted["gmx"], interactive=True, start_residue=0, end_residue=3, bogus=1,
    )
    assert matrix.shape == (4, 4)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 1.0)
    assert (fitted["out"] / "cm" / "contact_fraction.png").is_file()


def test_process_contact_averages_fraction_of_frames_below_cutoff(tmp_path):
    first = tmp_path / "a_dist_0.xvg"
    second = tmp_path / "a_dist_1000.xvg"
    # time, then a 2x2 pair matrix: the (0, 1) distance is 0.4 nm in 3 of 4 frames, then 0.9 nm in all of them
    first.write_text("# comment\n@ title\n0 0.0 0.4 0.4 0.0\n2 0.0 0.4 0.4 0.0\n4 0.0 0.4 0.4 0.0\n6 0.0 0.9 0.9 0.0\n")
    second.write_text("0 0.0 0.9 0.9 0.0\n2 0.0 0.9 0.9 0.0\n")
    matrix = an.process_contact([str(first), str(second)], contact_cutoff=0.5)
    assert matrix == pytest.approx(np.array([[1.0, 0.375], [0.375, 1.0]]))
    assert (tmp_path / "contact_fraction.dat").is_file()
