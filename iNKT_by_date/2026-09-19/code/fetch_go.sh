#!/usr/bin/env bash
set -euo pipefail
cd /home/zzz0054/bio3
out=iNKT_by_date/2026-09-19
trap 'echo $? > "$out/results/logs/fetch.exit"' EXIT
.venv/bin/python -u "$out/code/fetch_go.py" > "$out/results/logs/fetch.log" 2>&1
