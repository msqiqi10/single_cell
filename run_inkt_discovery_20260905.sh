set -eu
cd /home/zzz0054/bio3
out=output/iNKT_discovery_20260905
mkdir -p "$out/logs"
trap 'echo $? > "$out/logs/run.exit"' EXIT
.venv/bin/python -u notebooks/scripts/iNKT/run_inkt_discovery.py all > "$out/logs/run.log" 2>&1
