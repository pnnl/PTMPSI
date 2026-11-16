"""Shared fixtures: a real GROMACS ACE-SNC-ALA-NME system with a box-crossing raw trajectory."""

from __future__ import annotations

import os

os.environ["MPLBACKEND"] = "Agg"   # plt.show() in process_contact must not open a window or need a GUI backend

import shutil
import subprocess

import MDAnalysis as mda
import numpy as np
import pytest

from AnalysisTools import RESIDUETYPES_DAT

GMX = shutil.which("gmx_mpi") or shutil.which("gmx")


@pytest.fixture(scope="session")
def snc_system(tmp_path_factory):
    """Return a dict with the GROMACS binary, a .tpr, and a wrapped .xtc whose frames 5..9 jump one box vector."""
    if GMX is None:
        pytest.skip("GROMACS CLI not installed")
    gmx = GMX
    tmp_path = tmp_path_factory.mktemp("snc_system")

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

    return {"gmx": gmx, "dir": tmp_path, "tpr": tpr_path, "raw_xtc": raw_xtc, "box": box_dims}
