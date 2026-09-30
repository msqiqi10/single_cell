#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python notebooks/scripts/iNKT/build_inkt_discovery_brief.py
