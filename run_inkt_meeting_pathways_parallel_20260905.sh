set -eu
cd /home/zzz0054/bio3
out=output/iNKT_meeting_followup_20260905
trap 'echo $? > "$out/logs/pathways.exit"' EXIT
pids=()
for worker in 0 1 2 3; do
 MEETING_WORKER=$worker MEETING_WORKERS=4 .venv/bin/python -u notebooks/scripts/iNKT/run_inkt_meeting_pathways.py enrichment > "$out/logs/pathways_worker${worker}.log" 2>&1 &
 pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
if [ "$failed" -ne 0 ]; then exit 1; fi
.venv/bin/python - <<'PY'
from pathlib import Path
import json,pandas as pd
p=Path('output/iNKT_meeting_followup_20260905');pd.concat([pd.read_csv(p/f'tables/enrichment_ranking_audit_worker{i}.csv') for i in range(4)]).to_csv(p/'tables/enrichment_ranking_audit.csv',index=False)
j=[json.loads((p/f'logs/pathways_worker{i}.completed.json').read_text()) for i in range(4)];result=j[0];result['n_original_min20_contrasts']=sum(r['n_original_min20_contrasts'] for r in j);result['workers']=4;(p/'logs/pathways.completed.json').write_text(json.dumps(result,indent=2))
PY
.venv/bin/python -u notebooks/scripts/iNKT/run_inkt_meeting_pathways.py comparisons > "$out/logs/comparisons.log" 2>&1
