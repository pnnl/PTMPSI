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

## Step 0: environment

```bash
python3.12 -m venv ~/ptmpsi-venv      # requires-python is >=3.6,<3.13
source ~/ptmpsi-venv/bin/activate     # puts pdb2pqr30 on PATH
pip install -e /path/to/PTMPSI        # repo root, where pyproject.toml is
bash .agents/skills/ptmpsi/scripts/check_env.sh
```

`Protein.protonate` (`ptmpsi/protein/__init__.py`, line 722) calls `which("pdb2pqr30")` and raises `KeyError("Cannot find PDB2PQR")` if it is missing. `gen_ptm_files` calls `protonate` by default, so the venv must be activated, or its `bin` directory exported on `PATH`, in the shell that runs Python.

External programs by step:

| Step | Needs |
|---|---|
| 1-4, 6 (AlphaFold scripts, structures, chains, mutation, PTMs) | `ptmpsi` only |
| 5 get_qm_data | `ptmpsi`; the generated NWChem inputs are run later with NWChem |
| 7 docking | AutoDock Vina (`vina`) |
| 8 QM/MM | `pdb2pqr30`, NWChem |
| 9 PTM combinations | `pdb2pqr30`, internet for the example PDB; GROMACS only to run the output |

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
af      = Protein(uniprotid="G9BY57")    # downloads AlphaFold DB model v4
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

- Give `phi` and `psi` as numbers. Passing the strings `"alpha"` or `"beta"` to `append` fails (line 215 area has a bug: it reads an undefined `phi`). `prepend` accepts the strings.
- Defaults: `psi=-40.0` for `append`, `phi=-60.0` for `prepend`.
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
| methylation | GLU, GLH, LYS, LYN, ARG, HIS, HIP, HID, HIE |
| dimethylation | LYS, LYN, ARG |
| trimethylation | LYS, LYN |
| symmetric / asymmetric dimethylation | ARG |
| phosphorylation | SER, THR, TYR, ARG, HIS, HIP, HID, HIE, LYS, LYN, ASP, ASH, CYS, CYM |

The check lives in `check_ptm` (`ptmpsi/residues/ptms.py`, line 176). Anything else raises `MyDockingError`.

## Step 7: docking (needs `vina`)

```python
cutinase = Protein(pdbid="6tht")
cutinase.mutate("A:ALA130", "SER")
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

Same setup as step 7, plus `cutinase.protonate()` before writing the receptor, then:

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
- `machine=` and the spread `options["gromacs"]` describe the PNNL AQE-LDRD cluster (apptainer image, rankfile, `gmx_mpi`). They only shape the generated sbatch scripts. Whether generating files needs the cluster is not yet verified here; see Known gaps.
- Other machines are defined in `ptmpsi/slurm/__init__.py`: `deception` (line 32), `aqe_ldrd` (106), `aqe_h100`, `tahoma` (207), `perlmutter` (234), `frontier` (250), `polaris` (293).
- `ff` is `"amber99sb"` (default), `"amber99zn"` or `"amber14sb"`. An invalid name does not raise (the `KeyError` in `gen_ptm_files` is built but never raised), and the later copy fails instead.

## Checking your work

After step 9, list `<uid>/` and compare it with the layout in `references/templates.md`. If something is missing, do not edit generated files by hand. Fix the input and regenerate.

## Known gaps

- The tutorial's last cell is empty.
- `Tutorial_lcf.ipynb` and `Tutorial_frontier.ipynb` in the same folder are the same workflow for the LCF and Frontier machines. Use them for those clusters.
- Nothing in this skill has been run end to end yet. Steps 1-8 follow the source code. Step 9 listed the combinations correctly (counts 4, 6, 4, 1) but `gen_ptm_files` has not completed here; a test run stopped at `Cannot find PDB2PQR` because `pdb2pqr30` was not on `PATH`. A rerun with the venv's interpreter is pending.
- Steps 7 and 8 need `vina` and NWChem, which are not confirmed installed.
