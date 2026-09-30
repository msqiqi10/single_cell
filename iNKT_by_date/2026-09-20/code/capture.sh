#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=iNKT_by_date/2026-09-20
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 'echo $? > "$out/logs/capture.exit"' EXIT
.venv/bin/python -u "$out/code/capture_sources.py" > "$out/logs/capture.log" 2>&1
