---
name: ptmpsi
description: Use the PTM-Psi Python package (import ptmpsi) to build proteins and polymers, apply post-translational modifications (PTMs), and prepare simulations. Covers installing ptmpsi, loading or building structures (PDB, AlphaFold DB, local file, from scratch), AlphaFold job scripts, nylon surfaces, point mutations, AMBER99 parameters via NWChem, PTMs, docking, QM/MM input, and generating GROMACS .mdp/.itp/.tpr run directories for PTM combinations with thermodynamic integration. Encodes tutorials/Tutorial.ipynb as the worked example. Use when the user mentions PTM-Psi, ptmpsi, PTMPSI, Protein.modify, get_ptm_combinations, gen_ptm_files, get_qm_data, or asks to set up GROMACS or thermodynamic integration runs for modified residues.
---

# PTM-Psi

PTM-Psi is a Python package (`import ptmpsi`) for building proteins with post-translational modifications and preparing simulations. This skill covers the package. The steps below follow `tutorials/Tutorial.ipynb`, which is the worked example of every main feature.

Read this file first. Open the references only when a step tells you to:

- `references/api.md`: every function the tutorial calls, with file and line.
- `references/templates.md`: the GROMACS `.mdp`, `.itp`, `.tpr` and other template files. Read it before changing or checking any generated run directory.

## Ground rules

1. Run every Python snippet from a scratch directory, not from the repo root. `Protein(pdbid=...)` writes `<pdbid>.pdb` into the current directory, and `gen_ptm_files` writes a UUID-named folder there.
2. Use the exact strings below. PTM names are lowercase and must be keys of `ptmdict` (`ptmpsi/residues/__init__.py`, line 178). The tutorial text spells one of them wrong ("sulfynilation"). The valid key is `sulfinylation`.
3. Residue selectors use `"<chain>:<RES><number>"`, for example `"A:ALA130"`. The residue name may be left out (`"A:130"`). The number is the 1-based position in the chain's residue list, not the PDB number.
4. Do not invent keyword arguments. If a keyword is not listed in `references/api.md`, read the source function before using it.
5. Needs internet: `pdbid=`, `uniprotid=`.
6. Run `scripts/check_env.sh` before anything else. It reports what is missing.
7. Any script that calls `gen_ptm_files` must put its code under `if __name__ == "__main__":`. On macOS and Windows, Python starts worker processes by re-importing the script, and without the guard the run dies with `RuntimeError: An attempt has been made to start a new process before the current process has finished its bootstrapping phase`, followed by `BrokenProcessPool`. `scripts/ptm_combinations.py` shows the pattern.
8. On macOS, read "Platform notes: macOS (Apple Silicon)" below before running any GROMACS command.

## Step 0: environment

```bash
# From the repository root. Creates the venv (Python 3.9 to 3.12), installs ptmpsi, downloads Vina, runs check_env.sh
bash .agents/skills/ptmpsi/scripts/setup_env.sh ~/ptmpsi-venv --install-system-tools
source ~/ptmpsi-venv/bin/activate     # puts pdb2pqr30 and vina on PATH; do this in every new shell
```

`--install-system-tools` runs `brew install gromacs nwchem gnu-sed` for whichever of those is missing (macOS only). Without the flag the script lists what is missing and prints the command. On Linux, install GROMACS and NWChem with the system package manager; that path was not tested.

Manual equivalent: `python3.12 -m venv ~/ptmpsi-venv`, `source ~/ptmpsi-venv/bin/activate`, `pip install -e <repo root>`, `brew install gromacs nwchem gnu-sed`, then `bash .agents/skills/ptmpsi/scripts/check_env.sh`.

`Protein.protonate` (`ptmpsi/protein/__init__.py`, line 725) calls `which("pdb2pqr30")` and raises `KeyError("Cannot find PDB2PQR")` if it is missing. `gen_ptm_files` calls `protonate` by default, so the venv must be activated, or its `bin` directory exported on `PATH`, in the shell that runs Python.

External programs by step:

| Step | Needs |
|---|---|
| 1-4, 6 (AlphaFold scripts, structures, chains, mutation, PTMs) | `ptmpsi` only |
| 5 get_qm_data | `ptmpsi`; the generated NWChem inputs are run later with NWChem |
| 7 docking | AutoDock Vina (`vina`) |
| 8 QM/MM | `pdb2pqr30`, NWChem input generation needs `gsed` on macOS (`brew install gnu-sed`); running the input needs NWChem |
| 9 PTM combinations | `pdb2pqr30`, internet for the example PDB; GROMACS only to run the output |

## Where things are

All paths are relative to the repository root. Do not search for these; open them directly.

| To do this | Open |
|---|---|
| Load, build, mutate, modify, dock, protonate | `ptmpsi/protein/__init__.py` (`Protein`) |
| Combination logic | `ptmpsi/protein/tools.py` (`ptm_combination`, line 138) |
| Point mutation and PTM geometry | `ptmpsi/protein/mutate.py`, `ptmpsi/residues/ptms.py` |
| Valid residue and PTM names | `ptmpsi/residues/__init__.py` (`resdict`, `ptmdict`, `ptm2nonstandard`) |
| GROMACS run folders and `.mdp` files | `ptmpsi/gromacs/__init__.py`, `ptmpsi/gromacs/templates.py` |
| Force fields, `specbond.dat` | `ptmpsi/gromacs/forcefield/` |
| Cluster and partition definitions | `ptmpsi/slurm/__init__.py` |
| NWChem inputs (`get_qm_data`, `qmmm_optimize`) | `ptmpsi/nwchem/__init__.py`, `ptmpsi/nwchem/templates.py`, `ptmpsi/nwchem/qmmm.py` |
| Docking | `ptmpsi/docking/__init__.py` |
| Nylon surfaces | `ptmpsi/polymers/nylon6.py`, `nylon4.py` |
| AlphaFold scripts | `ptmpsi/alphafold/__init__.py` |
| The tutorial | `tutorials/Tutorial.ipynb` (`Tutorial_lcf.ipynb`, `Tutorial_frontier.ipynb` for those clusters) |

Machines and partitions (names are what `machine=` and `partition=` accept):

| Machine object | Partitions (CPUs / GPUs per node, max nodes) |
|---|---|
| `ptmpsi.slurm.aqe_ldrd` | `basic` (120/0, 200), `standard` (120/0, 10), `premium` (40/8, 8), `prm96c4g` (96/4, 3) |
| `ptmpsi.slurm.aqe_h100` | `premium` / `prm96c8g` (96/8, 28) |
| `ptmpsi.slurm.tahoma` | `normal` (36/0, 10, 72 h), `analysis` (36/1, 1, 3 h) |
| `ptmpsi.slurm.perlmutter` | `regular` (64/0, 10, 48 h), `preempt`, `shared` |
| `ptmpsi.slurm.frontier` | `batch` (56/8, 12 h), `extended` (56/8, 24 h, 64 nodes) |
| `ptmpsi.slurm.polaris` | `debug` (32/4, 1 h), `prod` (32/4, 24 h, 476 nodes), `preemptable` (32/4, 72 h) |
| `ptmpsi.slurm.deception` | `slurm` (64/0, 200 nodes, 168 h; default) |

The default machine is `tahoma`. Requests beyond a partition's limits raise `KeyError`.

## Platform notes: macOS (Apple Silicon)

Checked on macOS 26, arm64, Python 3.12.14, Homebrew GROMACS 2026.4 (CPU only, `thread_mpi`).

| Issue | What happens | What to do |
|---|---|---|
| Worker start method | `spawn` is the default, so unguarded scripts crash (ground rule 7) | Use `if __name__ == "__main__":` |
| BSD `sed` | The generated scripts used `sed -i '...\n...'`, which macOS `sed` rejects (`undefined label`). The restraint `#ifdef` blocks never reach `topol.top` and the later `grompp` calls fail | Fixed in this branch: `ptmpsi/gromacs/__init__.py` now writes a `perl -0pi -e` line. Folders generated by older code still contain `sed -i`; `scripts/localize_sbatch.py` rewrites it |
| Cluster launch lines | Generated `*_slurm.sbatch` files call `mpirun -np N --rankfile ... apptainer exec --nv ... gmx_mpi` and pass `-gpu_id`, `-nb gpu`, `-pme gpu`, `-update gpu`. None of that exists on a Mac | `scripts/localize_sbatch.py` strips them and emits a plain `gmx` script |
| `pdb2pqr30` | Not on `PATH` unless the venv is active; `protonate()` then raises `KeyError: Cannot find PDB2PQR` | `source ~/ptmpsi-venv/bin/activate`, or export `<venv>/bin` on `PATH` |
| Receptor hydrogens | `Protein.dock` fails for a receptor without hydrogens | Call `protonate()` before `write_pdb` (see step 7) |
| AutoDock Vina | Step 7 needs `vina` on `PATH`. `brew` has no `autodock-vina` formula and `pip download vina` fails to build on this Mac | Download the official arm64 release: `curl -L -o "$VIRTUAL_ENV/bin/vina" https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64 && chmod +x "$VIRTUAL_ENV/bin/vina"`, then `vina --version` (expect `v1.2.7`). Use `vina_1.2.7_mac_x86_64` on an Intel Mac (not tested) |
| Python version | `requires-python` is `>=3.6,<3.13` | Use 3.12, not the newest Homebrew Python |

Running the generated GROMACS workflow on a Mac (cluster steps and time limits do not apply):

```bash
cd <uid>/1tuples/0000
python <skill>/scripts/localize_sbatch.py 0000_slurm.sbatch --smoke 10 > smoke.sh   # 10 steps per mdrun
bash smoke.sh                                                                        # 148 s on a 15-core M-series Mac
```

Drop `--smoke 10` for the full-length stages. The localized script ends at `grompp ... -o md.tpr`; run production yourself with `gmx mdrun -deffnm md`.

## Step 1: AlphaFold prediction (writes job files, does not run AlphaFold)

```python
from ptmpsi.alphafold import prediction
prediction("GGGGGGGGGGG", multimer=False)
```

Writes `alphafold.sbatch`, `run_singularity.py`, `temp1.fasta` in the current directory. Database paths default to the Tahoma cluster. A list of sequences gives one prediction each. For multimer pass `multimer=True` and a path to a FASTA file with two or more `>name` records. Run the output on a cluster, not locally.

## Step 2: load a structure

```python
from ptmpsi.protein import Protein
local   = Protein(filename="/path/to/local.pdb")
fetched = Protein(pdbid="6tht")          # downloads from RCSB, saves 6tht.pdb
af      = Protein(uniprotid="G9BY57")    # downloads the current AlphaFold DB model
```

`delwat=True` and `delhet=True` are the defaults, so waters and hetero atoms are dropped on load.

## Step 3: build a chain from scratch

```python
chain = Protein()
chain.prepend(chain="A", residue="EAC")
chain.append(chain="A", residue="EAC", psi=120.0)
chain.append(chain="A", residue="EAC", psi=120.0)
chain.prepend(chain="A", residue="ACE", phi=-120.0)
chain.append(chain="A", residue="NME", psi=120.0)
```

- Give `phi` and `psi` as numbers or as `"alpha"` / `"beta"`.
- Defaults: `psi=-40.0` for `append`, `phi=-60.0` for `prepend`. Both accept `"alpha"` or `"beta"` as well as numbers (`append` needed a one-word fix for the string form; it is in this branch).
- Valid residue names are the keys of `resdict` (`ptmpsi/residues/__init__.py`, line 132). A one-letter code also works.
- Write results with `write_pdb("x.pdb")` or `write_xyz("x.xyz")`.

## Step 3b: nylon surfaces

```python
import ptmpsi
ptmpsi.polymers.nylon6.gen_surface(rangea=3, rangeb=3, rangec=3, model="Malta1979")
ptmpsi.polymers.nylon4.gen_surface(rangea=3, rangeb=3, rangec=3)
```

Nylon-6 models: `"Malta1979"`, `"Holmes"`, `"Arimoto"`. Nylon-4 has one model. Two PDB files are written: one with all strands in one MODEL, one with a MODEL per strand.

## Step 4: point mutation

```python
cutinase = Protein(pdbid="6tht")
cutinase.mutate("A:ALA130", "SER")
```

## Step 5: AMBER99 parameters for a residue (NWChem inputs)

```python
from ptmpsi.nwchem import get_qm_data
eac = Protein()
eac.prepend(chain="A", residue="EAC")
get_qm_data(eac)
```

Generates NWChem inputs for RESP charges and bonded force constants. Two ACE/NME-capped conformers are made. RESP is fitted on both together, so there is one set of charges. Bond and angle constants come out per conformer and must be averaged by hand. Defaults, all overridable as keyword arguments (`ptmpsi/nwchem/__init__.py`, line 91): `mult=1`, `charge=0`, `memory=2000`, `aobasis="def2-tzvp"`, `tdbasis="def2-svp"`, `cdbasis="def2-universal-jfit"`, `xcfun="r2scan"`, `grid="lebedev 120 14"`, `tdgrid="lebedev 100 14"`, `nscf=100`, `nopt=60`, `disp="disp vdw 4"`, `delta=0.0189`, `lshift=0.1`.

## Step 6: post-translational modifications

```python
snc = Protein()
snc.prepend(chain="A", residue="CYS")
snc.modify("A:CYS1", "nitrosylation")
```

Cysteine PTMs: nitrosylation, sulfenylation, sulfinylation, sulfonylation, sulfhydration, glutathionylation, cysteinylation, methylation, carbamoylation, cyanylation. The mapping to residue names is `ptm2nonstandard` (`ptmpsi/residues/__init__.py`, line 206): nitrosylation to SNC, sulfhydration to CSS, sulfenylation to CSO, sulfinylation to CSD, sulfonylation to OCS, glutathionylation to CGL, methylation to SMC, carbamoylation to QCS, cyanylation to XCN, cysteinylation to IYY, reduction to CYS, oxidation to CYX.

Other residues:

| PTM | Residues |
|---|---|
| acetylation | LYS, LYN |
| methylation | GLU, GLH, LYS, LYN, ARG, HIS, HIP, HID, HIE (GLU and GLH are listed by `check_ptm` but fail with `KeyError` in `modify`: they have no entry in `ptmsite`, `ptmpsi/residues/ptms.py` line 84) |
| dimethylation | LYS, LYN, ARG |
| trimethylation | LYS, LYN |
| symmetric / asymmetric dimethylation | ARG |
| phosphorylation | SER, THR, TYR, ARG, HIS, HIP, HID, HIE, LYS, LYN, ASP, ASH, CYS, CYM |

The check lives in `check_ptm` (`ptmpsi/residues/ptms.py`, line 176). Anything else raises `MyDockingError`.

Histidine methylation and phosphorylation stop and ask `Selection:` on stdin (1 = pros, 2 = tele). A non-interactive agent must supply the answer (for example by replacing `builtins.input`), or the call blocks forever.

## Step 7: docking (needs `vina`)

```python
cutinase = Protein(pdbid="6tht")
cutinase.mutate("A:ALA130", "SER")
cutinase.protonate()                     # required: Meeko needs hydrogens; sets cutinase.charge
cutinase.write_pdb("cutinase.pdb")

nylon6 = Protein()
nylon6.prepend("A", "EAC"); nylon6.prepend("A", "ACE")
nylon6.append("A", "EAC");  nylon6.append("A", "NME")
nylon6.write_xyz("nylon6.xyz")

cutinase.dock(ligand="nylon6.xyz", receptor="cutinase.pdb",
              boxcenter="A:SER130", boxsize=20.0)
```

`engine` is `"vina"` (default) or `"adgpu"`. The output file is named `<receptor>_<ligand>_<engine>.pdbqt`. `ligand` and `receptor` need a 3-letter extension, because the code strips the last four characters.

## Step 8: NWChem QM/MM input

Same setup as step 7 (which already protonates), then:

```python
ptmpsi.docking.write_pdb(cutinase.docking, filename="test2.pdb", pose=0, template=nylon6)
ptmpsi.nwchem.qmmm_optimize(filename="test2.pdb",
    qmres=[130, 207, 259, 260, 261, 262], counter=True,
    totcharge=cutinase.charge, xcfun="pbe96", grid="medium", aobasis="def2-svp")
```

`cutinase.charge` is only set by `protonate()`.

## Step 9: PTM combinations and GROMACS run directories

This is the part that produces `.mdp`, `.itp` and `.tpr` inputs. Read `references/templates.md` before touching its output. A runnable version is `scripts/ptm_combinations.py`.

```python
import ptmpsi
from ptmpsi.protein import Protein

protein = Protein(pdbid="6gve")          # GAP2/CP12/PRK complex
protein.findresidue("CYS")               # prints "<chain>:CYS<n>" lines

c1 = protein.get_ptm_combinations(
    [("CYS", "nitrosylation"), ("CYS", "sulfhydration")], ntuple=1)

exclude = [f"{ch}:{r}" for ch in ["A","D","E","F","I","K","N","O"] for r in ["CYS19","CYS77"]]
c2 = protein.get_ptm_combinations(
    [("CYS", "nitrosylation"), ("CYS", "sulfhydration")], exclude=exclude, ntuple=2)

protein.findresidue("CYX")               # disulfide cysteines
c3 = protein.get_ptm_combinations([("M:CYX", "reduction")], ntuple=4)
# c3 is a list; c3[k] holds the (k+1)-tuples. Lengths for this input: [4, 6, 4, 1]

uid = protein.gen_ptm_files(
    [c3[3]],                                         # note the extra list
    machine=ptmpsi.slurm.aqe_ldrd,
    **ptmpsi.slurm.aqe_ldrd.partitions["premium"].options["gromacs"])
```

Facts to rely on:

- `get_ptm_combinations` returns a list indexed by tuple size. Each entry is a list of combinations. Each combination is a tuple of `(residue_string, ptm_name)`.
- `gen_ptm_files` takes a list of lists of tuples. Passing `c3[3]` alone triggers a warning and an automatic wrap; wrap it yourself with `[c3[3]]`.
- It returns the UID string. The run directory is `<cwd>/<uid>/` (or `path=...`).
- The 4-tuple `[("M:CYX","reduction")]` case reduces both disulfides of chain M. 2-tuples that split a disulfide pair are chemically meaningless. CYX19 bonds to CYX29, CYX61 bonds to CYX70.
- A modified `specbond.dat` is copied into the run directory so that GROMACS makes a disulfide bond between nearby CYX residues and not between nearby CYS residues.
- Default simulation: minimization, heating, NPT equilibration, 100 ns production (`lenmd=100.0` in `generate_mdp`), plus 13 lambda windows of thermodynamic integration (TI) when `thermo=True` (default).
- `machine=` and the spread `options["gromacs"]` describe the PNNL AQE-LDRD cluster (apptainer image, rankfile, `gmx_mpi`). They only shape the generated sbatch scripts. Generating the files ran on a Mac with no cluster access.
- `pdb2pqr=False` skips protonation and writes the input structure as `protonated.pdb`. It failed with `UnboundLocalError` before this branch's fix.
- Other machines are defined in `ptmpsi/slurm/__init__.py`: `deception` (line 32), `aqe_ldrd` (106), `aqe_h100`, `tahoma` (207), `perlmutter` (234), `frontier` (250), `polaris` (293).
- `ff` is `"amber99sb"` (default), `"amber99zn"` or `"amber14sb"`. An invalid name does not raise (the `KeyError` in `gen_ptm_files` is built but never raised), and the later copy fails instead.

## Checking your work

After step 9, list `<uid>/` and compare it with the layout in `references/templates.md`. A complete combination folder has 14 `.mdp` files, `<NNNN>_slurm.sbatch`, `md.sbatch`, `submit_lambdas.sh` and a `dualti/` folder with `lam-00` to `lam-12`. If something is missing, do not edit generated files by hand. Fix the input and regenerate.

`assets/` holds real output from this workflow so you can compare without running it: `assets/mdp/` (rendered `.mdp` files, including two TI windows), `assets/example_run/` (the run `README`, `specbond.dat`, and the `amber99sb` `forcefield.itp`).

## Known gaps

- The tutorial's last cell is empty.
- `Tutorial_lcf.ipynb` and `Tutorial_frontier.ipynb` in the same folder are the same workflow for the LCF and Frontier machines. Use them for those clusters.
- Run on macOS arm64 (outputs inspected): step 1 (`prediction` writes `alphafold.sbatch`, `run_singularity.py`, `temp1.fasta`); step 2 (`pdbid=`, `filename=`, `uniprotid=`); step 3 (the trimer with `phi`/`psi`, and `"alpha"`/`"beta"` strings); step 3b (both `gen_surface` calls wrote `*_333.pdb` and `*_models_333.pdb`); step 5 (`get_qm_data` wrote `conf0.nw`, `conf1.nw`, the two `*_hess.nw`, `espfit.py`, `modseminario.py`, `torsions.dat`, `get_qm_data.sbatch`); step 6 (every CYS PTM, plus acetylation, methylation, dimethylation, trimethylation of LYS/LYN, ARG methylation and dimethylations, HIS/HIP/HID/HIE methylation, and phosphorylation of SER, THR, TYR, ARG, HIS, LYS, ASP, ASH, CYS, CYM, plus CYS reduction and oxidation); step 7 (docking, 9 poses); step 8 (`dock_cutinase.py --qmmm` exit 0; `prepare.nw` and `qmmm_optimization.sbatch` written; the cleaned PDB has 1174 lines with the renamed-hydrogen pattern and none with `HB2`/`HB3`); step 9 (combination counts 4, 6, 4, 1; `gen_ptm_files` for the 4-tuple and for a 1-combination run with `pdb2pqr=False`; the localized GROMACS chain through `grompp ... md.tpr` with 10-step `mdrun` calls, exit 0).
- Failing as shipped: `modify` methylation of GLU and GLH (`KeyError`).
- Not run: executing any generated NWChem input, dimethylation and trimethylation on residues other than the ones above, and `Protein.dock` with `flexible=` or `engine="adgpu"`.
- Not tested: production MD, the 13 TI windows, and the generated scripts on a cluster. The smoke run only shows that every `grompp` accepts the files and every `mdrun` takes steps.
- The local GROMACS is 2026.4; the cluster image is 2023.4.
