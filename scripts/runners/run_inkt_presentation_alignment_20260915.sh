#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
out=output/iNKT_presentation_alignment_20260915
mkdir -p "$out/logs"
trap 'echo $? > "$out/logs/build.exit"' EXIT
.venv/bin/python -u notebooks/scripts/iNKT/build_inkt_reference_aligned_presentation.py > "$out/logs/build.log" 2>&1
