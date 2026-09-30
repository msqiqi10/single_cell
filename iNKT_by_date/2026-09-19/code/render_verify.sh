#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=iNKT_by_date/2026-09-19
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
rm -f "$out/results/logs/render_verify.exit"
trap 'echo $? > "$out/results/logs/render_verify.exit"' EXIT
.venv/bin/python -u "$out/code/test_analysis.py" > "$out/results/logs/tests.log" 2>&1
.venv/bin/python -u "$out/code/summarize.py" > "$out/results/logs/summary.log" 2>&1
.venv/bin/python -u "$out/code/presentation.py" > "$out/results/logs/presentation.log" 2>&1
.venv/bin/python -u "$out/code/verify.py" > "$out/results/logs/verification.log" 2>&1
