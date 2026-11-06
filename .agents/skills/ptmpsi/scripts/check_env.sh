#!/usr/bin/env bash
# Report which PTM-Psi steps can run in this shell.
# Usage: bash check_env.sh
ok=0
check() { # name, command, steps-that-need-it
  if command -v "$2" >/dev/null 2>&1; then
    printf '  found    %-22s %s\n' "$1" "$(command -v "$2")"
  else
    printf '  MISSING  %-22s needed for: %s\n' "$1" "$3"; ok=1
  fi
}
echo "Platform: $(uname -s) $(uname -m)"
echo "Python:"
if ! python3 - <<'PY'
import sys, multiprocessing as mp
v = sys.version_info
ok = (3, 6) <= v[:2] < (3, 13)
print(f"  {v.major}.{v.minor}.{v.micro}", "(ok)" if (3, 6) <= v[:2] < (3, 13) else "(ptmpsi needs >=3.6,<3.13)")
print("  multiprocessing start method:", mp.get_start_method(),
      "(scripts calling gen_ptm_files need an if __name__ == '__main__' guard)" if mp.get_start_method() != "fork" else "")
try:
    import ptmpsi
    print("  ptmpsi importable from", ptmpsi.__file__)
except Exception as e:
    print("  ptmpsi NOT importable:", e)
    ok = False
sys.exit(0 if ok else 1)
PY
then
  ok=1
fi
echo "Programs on PATH:"
check pdb2pqr30 pdb2pqr30 "protonate(), qmmm_optimize(), gen_ptm_files() (activate the venv)"
check mk_prepare_receptor mk_prepare_receptor.py "Protein.dock() (installed with Meeko)"
check vina vina "Protein.dock()"
check gmx gmx "running generated GROMACS inputs locally"
check nwchem nwchem "running generated NWChem inputs"
check perl perl "topology edits in generated GROMACS scripts"
if [ "$(uname -s)" = "Darwin" ] && ! sed --version >/dev/null 2>&1; then
  echo "Note: sed is BSD sed. Generated scripts from before the perl fix contain 'sed -i' and will fail; run scripts/localize_sbatch.py on them."
fi
if command -v gmx >/dev/null 2>&1; then
  gmx --version 2>&1 | grep -E 'GROMACS version|MPI library|GPU support' | sed 's/^/  /'
fi
exit $ok
