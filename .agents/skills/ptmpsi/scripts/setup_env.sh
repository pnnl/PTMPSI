#!/usr/bin/env bash
# Create a Python environment for PTM-Psi and install everything the skill's scripts need
# that pip can provide, plus the AutoDock Vina binary.
#
#   bash setup_env.sh [VENV_DIR] [--install-system-tools]
#
# Default VENV_DIR: $HOME/ptmpsi-venv. Run it from anywhere; it finds the repository from its
# own location (<repo>/.agents/skills/ptmpsi/scripts).
#
# Without --install-system-tools it only reports missing system programs (GROMACS, NWChem and,
# on macOS, GNU sed) and prints the Homebrew command. With the flag it runs that command.
#
# Tested only on macOS arm64. The Vina asset names for the other platforms come from the
# AutoDock-Vina v1.2.7 release page and were not run.
set -euo pipefail

INSTALL_SYSTEM=0
ARGS=()
for a in "$@"; do
  if [ "$a" = "--install-system-tools" ]; then INSTALL_SYSTEM=1; else ARGS+=("$a"); fi
done
VENV=${ARGS[0]:-$HOME/ptmpsi-venv}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
[ -f "$REPO/pyproject.toml" ] || { echo "ERROR: $REPO is not the PTMPSI repository root (no pyproject.toml)"; exit 1; }

# pyproject.toml says requires-python >=3.6,<3.13, so the newest Homebrew Python is too new.
PY=""
for c in python3.12 python3.11 python3.10 python3.9; do
  if command -v "$c" >/dev/null 2>&1; then PY=$(command -v "$c"); break; fi
done
[ -n "$PY" ] || { echo "ERROR: need python3.9 to python3.12 on PATH (ptmpsi requires <3.13)"; exit 1; }
echo "Using $PY"

[ -d "$VENV" ] || "$PY" -m venv "$VENV"
"$VENV/bin/python" -m pip install --quiet --upgrade pip
echo "Installing ptmpsi from $REPO (needs git and network; several minutes the first time)"
"$VENV/bin/python" -m pip install --quiet -e "$REPO"

# AutoDock Vina: not on pip for this use and not in Homebrew.
if [ ! -x "$VENV/bin/vina" ] && ! command -v vina >/dev/null 2>&1; then
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64)   asset=vina_1.2.7_mac_aarch64 ;;
    Darwin-x86_64)  asset=vina_1.2.7_mac_x86_64 ;;
    Linux-x86_64)   asset=vina_1.2.7_linux_x86_64 ;;
    Linux-aarch64)  asset=vina_1.2.7_linux_aarch64 ;;
    *) echo "WARNING: no Vina asset known for $(uname -s)-$(uname -m); install vina yourself"; asset="" ;;
  esac
  if [ -n "$asset" ]; then
    curl -fL --silent --show-error -o "$VENV/bin/vina" \
      "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/$asset"
    chmod +x "$VENV/bin/vina"
  fi
fi
"$VENV/bin/vina" --version 2>/dev/null | head -1 || true

# System programs that pip cannot provide.
if [ "$(uname -s)" = "Darwin" ]; then
  need=()
  command -v gmx    >/dev/null 2>&1 || need+=(gromacs)
  command -v nwchem >/dev/null 2>&1 || need+=(nwchem)
  command -v gsed   >/dev/null 2>&1 || need+=(gnu-sed)   # qmmm_optimize swaps sed for gsed on macOS
  if [ ${#need[@]} -gt 0 ]; then
    if [ "$INSTALL_SYSTEM" = 1 ]; then
      command -v brew >/dev/null 2>&1 || { echo "ERROR: Homebrew not found; install ${need[*]} another way"; exit 1; }
      brew install "${need[@]}"
    else
      echo "Missing system programs: ${need[*]}"
      echo "Install them with:  brew install ${need[*]}     (or rerun this script with --install-system-tools)"
    fi
  fi
fi

echo
echo "Environment ready: $VENV"
echo "Activate it in every shell that runs ptmpsi (puts pdb2pqr30 and vina on PATH):"
echo "  source $VENV/bin/activate"
echo
PATH="$VENV/bin:$PATH" bash "$HERE/check_env.sh" || echo "(some programs above are missing; install them for the steps that need them)"
