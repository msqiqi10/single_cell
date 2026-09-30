"""Explicit offline reproduction entry point for the frozen meeting-followup analysis."""
from pathlib import Path
import subprocess,sys,json,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup'
commands=[]
def run(args):
    cmd=[sys.executable,*args];start=time.time()
    print('Running',args[0],flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True,timeout=1800)
    commands.append(dict(command=cmd,elapsed_seconds=time.time()-start,exit_code=0))
for name in ['build_function_inputs.py','rescue_markers.py','run_fusion.py','prepare_gnn_adapter.py']:
    run(['src/'+name])
for seed in [42,43,44]:
    run(['src/gnn_mouse_holdout.py','--edge_file',str(OUT/'network/mouse_GO_edges.tsv'),
         '--feature_npz',str(OUT/'network/semantic_embeddings.npz'),'--weight_column','JACCARD','--edge_weight_transform','none',
         '--hidden_dim','256','--num_layers','2','--epochs','120','--device','cpu','--seed',str(seed),'--split_seed',str(seed),
         '--output_csv',str(OUT/'models'/f'gnn_seed{seed}.csv'),'--output_npz',str(OUT/'models'/f'gnn_seed{seed}.npz')])
for name in ['summarize_network.py','targeted_sensitivity.py','integrate_existing.py','ontology_context.py','plot_followup.py','validate_followup.py']:
    run(['src/'+name])
(OUT/'rerun_commands.json').write_text(json.dumps(commands,indent=2))
print('Frozen-snapshot analysis reproduced. Narrative reports must be reviewed if inputs or methods change.',flush=True)
