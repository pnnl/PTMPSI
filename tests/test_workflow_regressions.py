"""Dependency-free workflow regressions: python3 -m unittest -v test_workflow_regressions.

Load the selected production definitions without importing scientific dependencies.
All generated inputs and mock programs remain in temporary directories. No MD runs.
"""
import ast
import contextlib
from datetime import datetime
import io
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / '.agents/skills/ptmpsi/scripts'


def definitions(path, names, **namespace):
    tree = ast.parse(path.read_text(), filename=str(path))
    tree.body = [node for node in tree.body
                 if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names]
    exec(compile(tree, str(path), 'exec'), namespace)
    return SimpleNamespace(**namespace)


@contextlib.contextmanager
def scratch():
    original = Path.cwd()
    with tempfile.TemporaryDirectory() as directory:
        os.chdir(directory)
        try:
            yield Path(directory)
        finally:
            os.chdir(original)


class ResidueSelectors(unittest.TestCase):
    def setUp(self):
        class Residue:
            name = 'ALA'
        self.residues = [Residue() for _ in range(130)]
        self.other = [Residue() for _ in range(130)]
        self.protein = SimpleNamespace(chains=[
            SimpleNamespace(name='A', residues=self.residues),
            SimpleNamespace(name='B', residues=self.other)])
        self.get_residue = definitions(
            ROOT / 'ptmpsi/protein/tools.py', {'get_residue'},
            Residue=Residue, one2three={'A': 'ALA', 'G': 'GLY'},
            MyDockingError=RuntimeError).get_residue

    def test_numeric_and_named_positions(self):
        for position in (1, 2, 30, 130):
            for selector in (position, str(position), f'A:{position}',
                             f'ALA{position}', f'A:ALA{position}',
                             f'A{position}', f'A:A{position}', f'A:ala{position}'):
                with self.subTest(selector=selector):
                    self.assertIs(self.get_residue(self.protein, selector),
                                  self.residues[position - 1])
            self.assertIs(self.get_residue(self.protein, f'B:{position}'),
                          self.other[position - 1])

    def test_residue_instance(self):
        self.assertIs(self.get_residue(self.protein, self.other[0]), self.other[0])

    def test_invalid_selectors(self):
        for selector in ('', 'A:', ':2', 'A', 'ALA', 'AL2', 'ALA-2', 'A:2x',
                         'A:0', '0', 0, -1, 'A:-1', 'A:1:2', 'A:²', 'Z2'):
            with self.subTest(selector=selector):
                with self.assertRaises(KeyError):
                    self.get_residue(self.protein, selector)
        with self.assertRaises(RuntimeError):
            self.get_residue(self.protein, 'Z:2')
        with self.assertRaises(IndexError):
            self.get_residue(self.protein, 'A:131')
        for selector in ('A:GLY2', 'A:G2'):
            with self.assertRaises(KeyError):
                self.get_residue(self.protein, selector)


class AlphaFoldInputs(unittest.TestCase):
    def setUp(self):
        slurm = runpy.run_path(str(ROOT / 'ptmpsi/slurm/__init__.py'))
        templates = runpy.run_path(str(ROOT / 'ptmpsi/alphafold/templates.py'))
        self.templates = SimpleNamespace(**templates)
        self.af = definitions(
            ROOT / 'ptmpsi/alphafold/__init__.py',
            {'AlphaFoldOptions', 'gen_fasta_paths', 'gen_pull', 'prediction'},
            datetime=datetime, isfile=os.path.isfile,
            Slurm=slurm['Slurm'], slurm_body=templates['slurm_body'],
            singularity_232='oras://example/alphafold:2.3.2',
            tahoma_datasets_231='/datasets', gen_script=lambda options: None)

    def test_grouping_and_independent_paths(self):
        with scratch():
            Path('target.fasta').write_text('>target\nGGG\n')
            paths = self.af.gen_fasta_paths(['AAA', 'CCC', 'target.fasta', 'DDD'])
            self.assertEqual(paths, ['./temp1.fasta', './target.fasta', './temp2.fasta'])
            self.assertEqual(Path('temp1.fasta').read_text(), '>sequence_1\nAAA\n>sequence_2\nCCC\n')
            self.assertEqual(Path('temp2.fasta').read_text(), '>sequence_1\nDDD\n')
            Path('second.fasta').write_text('>second\nCCC\n')
            self.assertEqual(self.af.gen_fasta_paths(['target.fasta', 'second.fasta']),
                             ['./target.fasta', './second.fasta'])

    def test_model_presets_in_generated_jobs(self):
        for model in ('monomer', 'multimer'):
            with self.subTest(model=model), scratch(), contextlib.redirect_stdout(io.StringIO()), \
                    patch.dict(sys.modules, {'ptmpsi.alphafold.templates': self.templates}):
                self.af.prediction(['AAA', 'CCC'], model=model)
                self.assertIn(f'--model_preset={model}', Path('alphafold.sbatch').read_text())
                self.assertEqual(Path('temp1.fasta').read_text().count('>'), 2)
        self.assertEqual(self.af.AlphaFoldOptions([]).model, 'monomer')

    def test_notebook_guidance(self):
        for name in ('Tutorial', 'Tutorial_lcf', 'Tutorial_frontier'):
            notebook = json.loads((ROOT / f'tutorials/{name}.ipynb').read_text())
            source = '\n'.join(''.join(cell['source']) for cell in notebook['cells'])
            self.assertNotIn('multimer=True', source)
            self.assertNotIn('multimer=False', source)
            self.assertIn('model="multimer"', source)
            self.assertIn('separate single-record FASTA files', source)


class TopologyUpdater(unittest.TestCase):
    def setUp(self):
        self.templates = runpy.run_path(str(ROOT / 'ptmpsi/gromacs/templates.py'))
        self.script = self.templates['update_topology']

    def test_exact_generated_script_compiles(self):
        compile(self.script, 'update_topology.py', 'exec')

    def test_supported_cyz_topology_execution(self):
        with scratch():
            names = ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'SG', 'DU', 'C', 'O']
            parameters = self.templates['CYZ'].splitlines()
            original = '[ moleculetype ]\nProtein 3\n[ atoms ]\n'
            for residue in range(2):
                original += f'; residue {residue + 1} CYZ\n'
                for offset, (name, parameter) in enumerate(zip(names, parameters), 1):
                    atom = residue * 11 + offset
                    atomtype, charge, mass = parameter.split()
                    original += f'{atom} {atomtype} {residue+1} CYZ {name} {atom} {charge} {mass}\n'
            original += '[ bonds ]\n8 19 1\n[ angles ]\n5 8 19 1\n[ dihedrals ]\n5 8 19 16 9\n[ system ]\nTwo CYZ updater fixture\n[ molecules ]\nProtein 1\n'
            Path('topol.top').write_text(original)
            Path('update_topology.py').write_text(self.script)
            result = subprocess.run([sys.executable, 'update_topology.py'],
                                    text=True, capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(Path('topol.top').read_text(), original)
            output = Path('TItop.top').read_text()
            atoms = output.split('[ atoms ]\n', 1)[1].split('[ bonds ]', 1)[0]
            atom_lines = [line.split() for line in atoms.splitlines() if line and not line.startswith(';')]
            self.assertEqual(len(atom_lines), 22)
            self.assertTrue(all(len(fields) == 11 for fields in atom_lines))
            self.assertEqual(atom_lines[7][-3:], ['SH', '-0.3119', '32.06'])
            self.assertIn('0.20380   138908.8   0.45000   0.000000', output)
            self.assertIn('103.700    569.024   103.700   0.00000', output)
            self.assertIn('14.64400  2  0.0   0.00000', output)
            self.assertIn('[ molecules ]\nProtein 1\n', output)


class LocalScripts(unittest.TestCase):
    def test_generated_default_and_launcher_free_stages(self):
        slurm = runpy.run_path(str(ROOT / 'ptmpsi/slurm/__init__.py'))
        generator = definitions(ROOT / 'ptmpsi/gromacs/__init__.py', {'generate_slurm'},
                                Slurm=slurm['Slurm']).generate_slurm
        templates = SimpleNamespace(**runpy.run_path(str(ROOT / 'ptmpsi/gromacs/templates.py')))
        for launcher in (None, ''):
            with self.subTest(launcher=launcher), scratch(), \
                    patch.dict(sys.modules, {'ptmpsi.gromacs.templates': templates}):
                kwargs = {} if launcher is None else {'mpirun': launcher}
                generator('0000.pdb', conc=0.154, thermo=False, checkpointing=False, **kwargs)
                result = subprocess.run([sys.executable, str(SCRIPTS / 'localize_sbatch.py'),
                                         '0000_slurm.sbatch', '--smoke', '10'],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                lines = result.stdout.splitlines()
                runs = [line for line in lines if line.startswith('gmx mdrun ')]
                self.assertEqual(len(runs), 12)
                self.assertTrue(all(line.endswith('-nsteps 10 -ntmpi 1') for line in runs))
                self.assertEqual(sum(line.startswith('gmx grompp ') for line in lines), 14)
                self.assertTrue(lines[-1].endswith('-o md.tpr'))
                self.assertNotIn('-deffnm md ', result.stdout)
                Path('local.sh').write_text(result.stdout)
                syntax = subprocess.run(['/bin/bash', '-n', 'local.sh'], capture_output=True)
                self.assertEqual(syntax.returncode, 0, syntax.stderr)

    def test_deception_scripts_use_modules_and_cpu_launcher(self):
        slurm = runpy.run_path(str(ROOT / 'ptmpsi/slurm/__init__.py'))
        generator = definitions(ROOT / 'ptmpsi/gromacs/__init__.py', {'generate_slurm'},
                                Slurm=slurm['Slurm']).generate_slurm
        templates = SimpleNamespace(**runpy.run_path(str(ROOT / 'ptmpsi/gromacs/templates.py')))
        machine = slurm['deception']
        options = machine.partitions['shared'].options['gromacs']
        with scratch(), patch.dict(sys.modules, {'ptmpsi.gromacs.templates': templates}):
            generator('0000.pdb', conc=0.154, thermo=False, checkpointing=False,
                      machine=machine, account='nwbrave', time=1, **options)
            script = Path('0000_slurm.sbatch').read_text()
        self.assertIn('#SBATCH --account=nwbrave\n#SBATCH --partition=shared', script)
        self.assertIn('module load gcc/14.2.0 openmpi/5.0.7 gromacs/2025.2', script)
        self.assertIn('mpirun -np 8 --bind-to none  gmx_mpi mdrun -deffnm npt', script)
        for absent in ('apptainer', '-nb gpu', '-gpu_id'):
            self.assertNotIn(absent, script)
        self.assertEqual(machine.default_partition().name, 'shared')
        self.assertIn('set -eo pipefail', script)
        with scratch(), patch.dict(sys.modules, {'ptmpsi.gromacs.templates': templates}):
            with self.assertRaises(ValueError):
                generator('0000.pdb', conc=0.154, thermo=False, machine=machine, account='nwbrave',
                          time=1, partition='a100_shared', ngpus=1)

    def test_deception_nwchem_and_alphafold_templates(self):
        nwchem = runpy.run_path(str(ROOT / 'ptmpsi/nwchem/templates.py'))
        header = nwchem['slurm_header']['DECEPTION']
        self.assertIn('module load gcc/14.2.0 openmpi/5.0.7 nwchem/7.2.3', header)
        self.assertIn('NWCHEM_BASIS_LIBRARY=/share/apps/nwchem/7.2.3/share/libraries/', header)
        self.assertIn('amber_s', nwchem['qmmm_slurm']['DECEPTION'])
        self.assertIn('mpirun -np $SLURM_NTASKS', nwchem['runsingularity']['DECEPTION'])
        alphafold = runpy.run_path(str(ROOT / 'ptmpsi/alphafold/templates.py'))
        self.assertIn('{gpuline}', alphafold['slurm_header']['DECEPTION'])
        body = alphafold['slurm_body']['Deception'].format(
            header='', data_dir='/d', container='/c.sif', nv='--nv ', dbs='full_dbs', model='monomer',
            run_relax=True, gpu_relax=True, fasta_paths='./a.fasta', date='2022-01-01')
        self.assertIn('apptainer run --nv', body)
        self.assertIn('--data_dir=$Data_Loc', body)

    def test_deception_nwchem_inputs_omit_cando_flags(self):
        try:
            from ptmpsi.protein import Protein
            from ptmpsi.nwchem import get_qm_data
            from ptmpsi.slurm import deception, tahoma
        except ImportError as error:
            self.skipTest(str(error))
        for machine, expected in ((deception, 0), (tahoma, 2)):
            with scratch() as directory:
                residue = Protein()
                residue.prepend(chain='A', residue='EAC')
                get_qm_data(residue, machine=machine, account='nwbrave')
                for name in ('conf0.nw', 'conf0_hess.nw', 'conf1.nw'):
                    text = (directory / name).read_text()
                    self.assertEqual(text.count('set int:cando_'), expected, (machine.name, name))

    def test_bare_piped_and_launched_commands(self):
        with scratch():
            commands = ['# comment -o md.tpr', 'echo SOL | gmx_mpi genion -s ions.tpr',
                        '  gmx mdrun -deffnm minim']
            for launcher in ('srun -n 2', 'mpirun -np 2', 'mpiexec -n 2', 'flux run -n 2'):
                commands.append(f'{launcher} apptainer exec --nv image gmx_mpi mdrun -gpu_id 0 -nb gpu -pme gpu -npme 1 -bonded gpu -update gpu -deffnm npt')
            commands += ["sed -i '/; Include water topology/i #ifdef P1000' topol.top",
                         'gmx grompp -f md.mdp -o   md.tpr', 'gmx mdrun -deffnm md']
            Path('input.sbatch').write_text('\n'.join(commands))
            result = subprocess.run([sys.executable, str(SCRIPTS / 'localize_sbatch.py'),
                                     'input.sbatch'], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('echo SOL | gmx genion', result.stdout)
            self.assertEqual(result.stdout.count('gmx mdrun'), 5)
            self.assertIn('perl -0pi', result.stdout)
            for removed in ('apptainer', 'gpu', 'mpirun', 'srun', '-nsteps', '-deffnm md'):
                self.assertNotIn(removed, result.stdout)
            Path('input.sbatch').write_text('gmx mdrun -deffnm minim\n')
            failure = subprocess.run([sys.executable, str(SCRIPTS / 'localize_sbatch.py'),
                                      'input.sbatch'], text=True, capture_output=True)
            self.assertNotEqual(failure.returncode, 0)
            self.assertEqual(failure.stdout, '')

    def test_environment_failures(self):
        script = SCRIPTS / 'check_env.sh'
        with tempfile.TemporaryDirectory() as directory:
            tools = Path(directory)
            for name in ('pdb2pqr30', 'mk_prepare_receptor.py', 'vina', 'gmx', 'nwchem', 'perl', 'sed', 'grep'):
                tool = tools / name
                tool.write_text('#!/bin/sh\nexit 0\n')
                tool.chmod(0o755)
            uname = tools / 'uname'
            uname.write_text('#!/bin/sh\necho Linux\n')
            uname.chmod(0o755)
            for version, importable, expected in ((12, True, 0), (13, True, 1), (5, True, 1), (12, False, 1)):
                with self.subTest(version=version, importable=importable):
                    python = tools / 'python3'
                    python.write_text(f'#!{sys.executable}\n'
                                      'import sys, types\nfrom collections import namedtuple\n'
                                      f'sys.version_info = namedtuple("Version", "major minor micro")(3, {version}, 0)\n'
                                      f'sys.modules["ptmpsi"] = types.SimpleNamespace(__file__="mock") if {importable!r} else None\n'
                                      'exec(compile(sys.stdin.read(), "<check_env>", "exec"))\n')
                    python.chmod(0o755)
                    result = subprocess.run(['/bin/bash', str(script)],
                                            env=dict(os.environ, PATH=directory), text=True, capture_output=True)
                    self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
                    self.assertIn('Programs on PATH:', result.stdout)
                    if not importable:
                        self.assertIn('ptmpsi NOT importable', result.stdout)
            python.write_text(f'#!{sys.executable}\n'
                              'import sys, types\nfrom collections import namedtuple\n'
                              'sys.version_info = namedtuple("Version", "major minor micro")(3, 12, 0)\n'
                              'sys.modules["ptmpsi"] = types.SimpleNamespace(__file__="mock")\n'
                              'exec(compile(sys.stdin.read(), "<check_env>", "exec"))\n')
            (tools / 'vina').unlink()
            external = subprocess.run(['/bin/bash', str(script)], env=dict(os.environ, PATH=directory),
                                      text=True, capture_output=True)
            self.assertEqual(external.returncode, 1)
            self.assertIn('MISSING  vina', external.stdout)
            python.unlink()
            missing = subprocess.run(['/bin/bash', str(script)], env=dict(os.environ, PATH=directory),
                                     text=True, capture_output=True)
            self.assertEqual(missing.returncode, 1)
            self.assertIn('python3', missing.stderr)
            self.assertIn('Programs on PATH:', missing.stdout)


if __name__ == '__main__':
    unittest.main()
