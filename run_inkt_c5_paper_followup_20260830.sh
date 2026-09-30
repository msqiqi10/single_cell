#!/usr/bin/env bash
set -euo pipefail

cd /home/zzz0054/bio3
run_output_dir="output/iNKT_reproduction_deck/20260830_C5_paper_Fig3DEF_followup"
kegg_gmt="$run_output_dir/tables/20260830_KEGG_2019_Mouse.gmt"
if [[ ! -s "$kegg_gmt" ]]; then
  echo "Missing frozen KEGG GMT: $kegg_gmt" >&2
  exit 1
fi
mkdir -p "$run_output_dir"

.venv/bin/python notebooks/scripts/iNKT/run_inkt_c5_paper_followup.py \
  --out-dir "$run_output_dir" \
  --kegg-gmt "$kegg_gmt" \
  "$@" \
  2>&1 | tee "$run_output_dir/20260830_run.log"
