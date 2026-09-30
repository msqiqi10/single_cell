set -eu
cd /home/zzz0054/bio3
out=output/iNKT_meeting_followup_20260905
trap 'echo $? > "$out/logs/report.exit"' EXIT
while [ ! -f "$out/logs/comparisons.completed.json" ]; do
 if ! tmux has-session -t bio3_meeting_pathways_20260905 2>/dev/null; then
  echo 'Pathway/comparison job ended without completion receipt' > "$out/logs/report.log"
  exit 1
 fi
 sleep 5
done
.venv/bin/python -u notebooks/scripts/iNKT/build_inkt_meeting_report.py > "$out/logs/report.log" 2>&1
