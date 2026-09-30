#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH=/home/zzz0054/bio3/notebooks/scripts/iNKT
out=output/iNKT_discovery_20260905
if .venv/bin/python -m unittest discover -s notebooks/scripts/iNKT -p 'test_*.py' -v > "$out/logs/verification_tests.log" 2>&1; then
  echo 0 > "$out/logs/verification_tests.exit"
else
  echo 1 > "$out/logs/verification_tests.exit"
  exit 1
fi
.venv/bin/python notebooks/scripts/iNKT/verify_inkt_discovery_outputs.py
