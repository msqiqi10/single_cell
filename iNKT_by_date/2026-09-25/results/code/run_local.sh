#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export MPLCONFIGDIR="$PWD/cache/matplotlib"
export NUMBA_CACHE_DIR="$PWD/cache/numba"
export XDG_CACHE_HOME="$PWD/cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
case "${1:-verify}" in
 verify) .venv/bin/python src/verify_local.py ;;
 go-pilot) .venv/bin/python src/go_pilot.py ;;
 followup) .venv/bin/python src/run_followup.py ;;
 verify-followup) .venv/bin/python src/validate_followup.py ;;
 *) echo 'Usage: ./run_local.sh verify | go-pilot | followup | verify-followup' >&2; exit 2 ;;
esac
