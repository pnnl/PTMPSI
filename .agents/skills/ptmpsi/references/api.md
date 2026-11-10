# API pointers

Paths are relative to the repo root. Line numbers are for the `skill` branch (base commit `f5223119` plus the fixes in this branch). Re-check with `grep -n` if the file changed.

## Protein (`ptmpsi/protein/__init__.py`)

| Method | Line | Notes |
|---|---|---|
| `Protein(filename=None, pdbid=None, uniprotid=None, interactive=False, delwat=True, delhet=True)` | 35 | Exactly one source, or none for an empty chain "A" |
| `write_pdb(pdbfile)`, `write_xyz(xyzfile)`, `write_fasta(fastafile)` | 97-103 | |
| `mutate(original, new)` | 213 | Calls `point_mutation` in `ptmpsi/protein/mutate.py`, line 14 |
| `append(chain, residue, psi=None)` | 218 | Numeric `psi`, or `"alpha"` / `"beta"` |
| `prepend(chain, residue, phi=None)` | 259 | `phi` numeric or `"alpha"` / `"beta"` |
| `delchain(chain)` | 303 | |
| `findresidue(name)` | 313 | Prints `chain:NAMEnumber` |
| `modify(original, modification, logger=None)` | 321 | Calls `post_translational_modification`, `mutate.py` line 84 |
| `get_ptm_combinations(ptms, exclude=None, ntuple=-1)` | 325 | Calls `ptm_combination`, `ptmpsi/protein/tools.py` line 138 |
| `check_combinations_format(combinations)` | 328 | Needs list of lists of tuples |
| `process_combination(...)` | 355 | One combination to one folder; used by `gen_ptm_files` |
| `gen_ptm_files(combinations, path=None, prefix=None, ff='amber99sb', **kwargs)` | 462 | Returns the UID |
| `dock(ligand, receptor, boxcenter, boxsize=10, output=None, flexible=None, engine=None, exhaustiveness=None, spacing=None)` | 666 | |
| `protonate(pdbin=None, pdb=None, pqr=None, ph=7, path=".")` | 725 | Runs `pdb2pqr30` with AMBER and propka; sets `self.charge` |

`gen_ptm_files` keyword arguments read in the source (lines 459-540): `thermo=True`, `pdb2pqr=True`, `cofactor=None`, `lenlambda=2.0`, `timestep=2.0`, `temp=300.0`, `nsplit=1`, `gpu_id=""`, `max_script_workers=1`, `machine`, `auto_submit_ti=True`, `checkpointing=True`, `bundling=False`, `time=12`, `combination_logging=True`. Remaining keywords go to `ptmpsi.gromacs.generate` and then `generate_slurm` (`mpirun`, `container`, `gmx`, `ntasks`, `nthreads`, `nstlist`, `partition`, `account`, `nnodes`).

## Residues and PTMs (`ptmpsi/residues/`)

| Name | File | Line |
|---|---|---|
| `resdict` (valid residue names) | `__init__.py` | 132 |
| `ptmdict` (valid PTM names) | `__init__.py` | 178 |
| `ptm2nonstandard` | `__init__.py` | 206 |
| `check_ptm(residue, ptm)` | `ptms.py` | 179 |
| `Template` class | `template.py` | 4 |

Valid PTM names (keys of `ptmdict`): acetylation, citrullination, cysteinylation, glutathionylation, glycosylation, hydroxylation, methylation, myristoylation, nitration, nitrosylation, oxidation, palmitoylation, phosphorylation, prenylation, reduction, sulfenylation, sulfhydration, sulfinylation, sulfonylation, dimethylation, trimethylation, symmetric dimethylation, asymmetric dimethylation, cyanylation, carbamoylation. Being a key does not mean it is implemented: only the ones in the `SKILL.md` tables work.

## GROMACS generation (`ptmpsi/gromacs/__init__.py`)

| Function | Line |
|---|---|
| `generate_mdp(temp=300, posres=[...], **kwargs)` | 39 |
| `generate_slurm(infile, posres=[...], ...)` | 95 |
| `sltcap(protein, bt="dodecahedron", size=None, conc=0.154, charge=0)` | 355 |
| `generate(protein, filename=None, size=None, conc=None, charge=0, **kwargs)` | 376 |

## Slurm machines (`ptmpsi/slurm/__init__.py`)

`Machine` (1), `Partition` (21), `deception` (38; partitions `shared`, `slurm`, `short`), `aqe_ldrd` (77; partitions `basic`, `standard`, `premium`, `prm96c4g`), `get_machine_name` (142), `aqe_h100` (150), `tahoma` (178), `perlmutter` (190), `frontier` (206), `polaris` (249), `Slurm` (306). `Slurm(caller, **kwargs)` takes `partition`, `ngpus`, `ntasks`, `nthreads`, `time` (hours), `account`, `nnodes`, `scratch`, `jobname`. It raises `KeyError` when a request exceeds the partition.

## NWChem (`ptmpsi/nwchem/__init__.py`)

`get_qm_data(residue, ligand=False, metal=False, ff="AMBER99", dohfresp=True, path=None, **kwargs)` at 91. `qmmm_optimize(filename, qmres=None, counter=False, center=False, orient=False, totcharge=0, path=None, **kwargs)` at 673. `forcebalance` (418), `torsionscan` (583), `bonds_and_angles` (619).

## Docking (`ptmpsi/docking/__init__.py`)

`dock_ligand(cls)` at 27 and `write_pdb(cls, filename=None, pose=0, template=None)` at 115.

## Polymers (`ptmpsi/polymers/`)

`nylon4.gen_surface(rangea=1, rangeb=1, rangec=1, model="Fredericks")` at `nylon4.py` 130. `nylon6.gen_surface(rangea=1, rangeb=1, rangec=1, model="Malta1979")` at `nylon6.py` 223.

## AlphaFold (`ptmpsi/alphafold/__init__.py`)

`AlphaFoldOptions` takes `model="monomer"`, `version="2.3.2"`, `dbs="full_dbs"`, `machine`, `data_dir`, and others. `prediction(fasta, **kwargs)` is at line 147. Frontier and Polaris force version 2.3.2 and `full_dbs`.
