#!/usr/bin/env bash
# Report which PTM-Psi tutorial steps can run in this shell.
# Usage: bash check_env.sh
ok=0
check() { # name, command, steps-that-need-it
  if command -v "$2" >/dev/null 2>&1; then
    printf '  found    %-10s %s\n' "$1" "$(command -v "$2")"
  else
    printf '  MISSING  %-10s needed for: %s\n' "$1" "$3"; ok=1
  fi
}
echo "Python:"
python3 - <<'PY'
import sys
v = sys.version_info
print(f"  {v.major}.{v.minor}.{v.micro}", "(ok)" if (3, 6) <= v[:2] < (3, 13) else "(ptmpsi needs >=3.6,<3.13)")
try:
    import ptmpsi
    print("  ptmpsi importable from", ptmpsi.__file__)
except Exception as e:
    print("  ptmpsi NOT importable:", e)
PY
echo "Programs on PATH:"
check pdb2pqr30 pdb2pqr30 "protonate(), qmmm_optimize(), gen_ptm_files() (activate the venv)"
check gmx gmx "running the generated GROMACS inputs"
check nwchem nwchem "running generated NWChem inputs"
check vina vina "Protein.dock()"
exit $ok
