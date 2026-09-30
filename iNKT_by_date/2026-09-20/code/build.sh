#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=iNKT_by_date/2026-09-20
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 'echo $? > "$out/logs/build.exit"' EXIT
.venv/bin/python -u "$out/package/code/rebuild_decks.py" > "$out/logs/build.log" 2>&1
.venv/bin/python -u "$out/code/verify_package.py" > "$out/logs/verify.log" 2>&1
