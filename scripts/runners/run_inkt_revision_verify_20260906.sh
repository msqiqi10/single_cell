#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=output/iNKT_revision_responses_20260906
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 'echo $? > "$out/logs/verify.exit"' EXIT
.venv/bin/python notebooks/scripts/iNKT/verify_inkt_revision_responses.py
