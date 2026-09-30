#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out_dir=output/iNKT_meeting_followup_20260905
stage=${1:-scores}
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 's=$?; echo "$s" > "$out_dir/logs/$stage.exit"' EXIT
.venv/bin/python notebooks/scripts/iNKT/run_inkt_meeting_followup.py "$stage"
