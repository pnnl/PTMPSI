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

ROOT = Path(__file__).resolve().parent
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
