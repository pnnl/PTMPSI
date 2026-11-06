# Template files in PTM-Psi

"Template" in this project means two groups of files. The first group matters most for the tutorial's GROMACS output.

## A. GROMACS inputs: `.mdp`, `.itp`, `.tpr`

| Type | Where it comes from | Role |
|---|---|---|
| `.mdp` | Python string templates in `ptmpsi/gromacs/templates.py`, filled by `generate_mdp` (`ptmpsi/gromacs/__init__.py`, line 39) | Run parameters for `gmx grompp` |
| `.itp` | Force field folders under `ptmpsi/gromacs/forcefield/`; restraint files `<force>p.itp` made at run time by `gmx genrestr` | Topology includes: force field terms, water models, ions, restraints |
| `.tpr` | Never stored. Made by `gmx grompp` inside the generated sbatch scripts | Binary run input |

### `.mdp` templates (`ptmpsi/gromacs/templates.py`)

| Variable | Line | Written as | Placeholders |
|---|---|---|---|
| `ions` | 1 | `ions.mdp` | none |
| `minim` | 18 | `minim<force>.mdp` (one per position restraint) and `minim.mdp` | `{restraint}` |
| `minimcg` | 44 | conjugate-gradient variant | `{restraint}` |
| `heating` | 69 | `heating.mdp`, `fheating.mdp` | `{timestep}`, `{nsteps}` |
| `npt` | 135 | `npt.mdp`, `fnpt.mdp` | `{timestep}`, `{nsteps}` |
| `md` | 200 | `md.mdp` | `{timestep}`, `{nsteps}` |
| `qlambdas` | 267 | `dualti/lam-NN/01-q/grompp.mdp` | `{lambda_state}`, `{nsteps}`, `{timestep}`, `{temp}` |
| `vdwlambdas` | 337 | `dualti/lam-NN/02-vdw/grompp.mdp` | `{lambda_state}`, `{nsteps}`, `{timestep}`, `{temp}` |

`generate_mdp` options (all via `kwargs`, line 39 onward): `lennvt=500.0` and `lennpt=500.0` (picoseconds), `lenmd=100.0` (nanoseconds), `timestep=2.0` (femtoseconds), `nstxout_compressed=50000`, `posres=[1000,500,100,50,10,5,1]`, `temp=300`.

The TI templates fill `init_lambda_state` with the window number 0-12 (13 windows). Do not hard-code lambda arrays elsewhere.

### `.itp` and force field files

Folders (each is a complete GROMACS force field directory):

- `ptmpsi/gromacs/forcefield/amber99sb.ff/` (default)
- `ptmpsi/gromacs/forcefield/amber99zn.ff/`
- `ptmpsi/gromacs/forcefield/amber14sb.ff/`

Key files inside: `forcefield.itp` (entry point), `ffnonbonded.itp`, `ffbonded.itp`, `ions.itp`, water models (`tip3p.itp`, `spce.itp`, `spc.itp`, `tip4p.itp`, `tip4pew.itp`, `tip5p.itp`), `aminoacids.rtp` (residue topologies, where the modified residues live), `aminoacids.hdb`, `*.tdb`, `atomtypes.atp`, `watermodels.dat`.

Sibling files copied next to the force field folder in each run directory:

- `ptmpsi/gromacs/forcefield/specbond.dat`: the modified disulfide bond rule (CYX bonds, CYS does not).
- `ptmpsi/gromacs/forcefield/residuetypes.dat`: residue type table.

Where each is placed: `Protein.gen_ptm_files` copies the force field folder, `residuetypes.dat` and `specbond.dat` to `<uid>/` (line 459 onward), then `process_combination` (line 352) makes relative symlinks to them in each `<N>tuples/<NNNN>/` folder and in every TI window folder.

### Per-combination topology edits

`ptmpsi/gromacs/templates.py` also holds the residue parameter blocks and the script that rewrites topologies for the modified cysteines: `CYZ`, `SNC`, `CSO`, `CGL`, `IYY` strings (lines 842-926) and `update_topology` (line 951). Treat them as templates too. Change them only when adding or fixing a modification.

### Run directory layout (what `gen_ptm_files` writes)

```
<uid>/
  README                      list of combination -> folder
  submit.sh                   submits everything
  <ff>.ff/  residuetypes.dat  specbond.dat
  protonated.pdb  protonated.pqr
  <N>tuples/<NNNN>/
      <NNNN>.pdb  <NNNN>_slurm.sbatch  md.sbatch  submit_lambdas.sh  process_combination.log
      ions.mdp  minim.mdp  minim{1000,500,100,50,10,5,1}.mdp  heating.mdp  npt.mdp  fheating.mdp  fnpt.mdp  md.mdp
      <ff>.ff  residuetypes.dat  specbond.dat      symlinks up to <uid>/
      dualti/topol.top  dualti/update_topology.py
      dualti/lam-00 ... lam-12/
          <NNNN>_lamNN_slurm.sbatch
          01-q/   grompp.mdp  + symlinks (topol.top, index.ndx, md.gro, md.cpt)
          02-vdw/ grompp.mdp  + symlinks (topol.top, index.ndx, md.gro)
      (topol.top, index.ndx, md.gro, md.cpt do not exist until the run; the symlinks dangle until then)
```

`.tpr` files appear only after the sbatch scripts run `gmx grompp` (`ptmpsi/gromacs/__init__.py`, lines 206-345).

## B. Python-side templates

| File | What it is |
|---|---|
| `ptmpsi/residues/template.py` | `Template` class: name, coordinates, elements, backbone and attachment vectors for one residue. `find()`/`find_coord()` look up atoms. |
| `ptmpsi/residues/aminoacids.py`, `nonstandard.py` | One `Template()` instance per residue. `nonstandard.py` has the PTM residues (SNC, CSO, CSD, OCS, CSS, CGL, QCS, SMC, XCN, IYY) plus EAC and ABA. |
| `ptmpsi/residues/ptms.py` | Attachment geometry for methyl, acetyl and phosphoryl groups. |
| `ptmpsi/nwchem/templates.py` | NWChem input strings and Slurm headers (`hfresp`, `respnw`, `hessnw`, `torsnw`, `forcebalance_input`, and others). |
| `ptmpsi/alphafold/templates.py` | `slurm_header` and `slurm_body` per machine. |
| `ptmpsi/gromacs/templates.py` | `slurm_header`, `flux_header`, `parsl_*` job script heads for each machine, plus the `.mdp` strings above. |

## C. Chimera plugin templates

`plugins/v4.0/CysPTM/src/templates/*.idatmres` (and a `v3.0` copy) are UCSF Chimera atom-type templates for CYS, SNC, CSO, CSD, OCS, CSS, CGL, SMC, QCS, XCN. They are separate from the Python workflow.

## Files in the repo that every run uses

- `generate_mdp()` in `ptmpsi/gromacs/__init__.py` writes the non-TI `.mdp` files. Run it in an empty folder to see them.
- `ptmpsi/gromacs/templates.py` holds the TI window templates. Window `01-q` sets `couple-lambda0 = vdw-q` (line 314). Window `02-vdw` sets `couple-lambda0 = vdw` (line 384), and an inline comment there says `vdw-q` was wrong. Keep `vdw`.
- `ptmpsi/gromacs/forcefield/specbond.dat` defines SG-SG pairs for `CYX-CYX` (0.2 and 0.25 nm), `CYM-CYM` and `CYZ-CYZ`. There is no `CYS-CYS` line, which is what keeps nearby CYS residues from bonding.
- `ptmpsi/gromacs/forcefield/amber99sb.ff/forcefield.itp` is the entry point of the force field. The other `.itp` files are in the same folder.
- `gen_ptm_files` writes a run index `README` that lists each combination and its folder.

Every stage runs the same way in the macOS check: `ions.tpr`, seven `minim*.tpr`, `heating.tpr`, `npt.tpr`, `minim.tpr`, `fheating.tpr`, `fnpt.tpr` and `md.tpr` were all produced by `gmx grompp` 2026.4.

## Rules

- Do not edit a template string to fix one run. Fix the input arguments.
- A new modification needs four matching edits: the residue `Template` in `nonstandard.py`, the entries in `resdict`, `ptmdict` and `ptm2nonstandard` (`ptmpsi/residues/__init__.py`), and the force field entry in each `aminoacids.rtp` it should support.
