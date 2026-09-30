"""Keep ontology ancestry and expert-prior bridges distinct from overlap edges."""
from pathlib import Path
import json
import networkx as nx
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup';SRC=ROOT/'data/bio3/iNKT_by_date/2026-09-19/sources'
graph=nx.DiGraph();current=None
for line in (SRC/'go-basic.obo').read_text().splitlines():
    if line=='[Term]':current=''
    elif line.startswith('['):current=None
    elif current is not None:
        if line.startswith('id: '):current=line[4:];graph.add_node(current)
        elif line.startswith('is_a: '):graph.add_edge(current,line.split()[1],relation='is_a')
        elif line.startswith('relationship: part_of '):graph.add_edge(current,line.split()[2],relation='part_of')
assert nx.is_directed_acyclic_graph(graph)
nodes=pd.read_csv(OUT/'network/annotated_nodes.csv').set_index('GOID');rows=[]
for child in nodes.index:
    for parent in sorted(nx.descendants(graph,child)&set(nodes.index)):
        path=nx.shortest_path(graph,child,parent)
        relations=[graph.edges[a,b]['relation'] for a,b in zip(path,path[1:])]
        rows.append(dict(child_GO=child,parent_GO=parent,child_module=nodes.loc[child,'gnn_module'],parent_module=nodes.loc[parent,'gnn_module'],
                         path_length=len(path)-1,path=';'.join(path),relations=';'.join(relations),
                         interpretation='ontology ancestry; child to parent; not a causal or regulatory edge'))
pd.DataFrame(rows).to_csv(OUT/'network/ontology_ancestry_paths.csv',index=False)
sets={}
for line in (SRC/'GO_mouse_measured.gmt').read_text().splitlines():
    gid,name,*genes=line.split('\t');sets[gid]=set(genes)
priors={'AP1_transcription_factor_complex':('GO:0035976',sets['GO:0035976']),
        'MAPK_cascade':('GO:0000165',sets['GO:0000165']),
        'iNKT17_identity_prior':('expert_marker_panel',{'Rorc','Il23r','Ccr6'}),
        'cytotoxic_core_prior':('existing_PAGER_core',{'Prf1','Gzma','Gzmb'})}
evidence=pd.read_csv(OUT/'tables/module_gene_evidence.csv.gz');evidence=evidence[evidence.method=='gnn_module']
go=pd.read_csv(OUT/'tables/all_original_condition_GO.csv.gz',keep_default_na=False);bridges=[]
for (module,uid,direction),group in evidence.groupby(['module','unit_id','direction']):
    for name,(source,members) in priors.items():
        overlap=set(group.gene)&members
        q=go[(go.unit_id==uid)&(go.direction==direction)&(go.go_id==source)]
        bridges.append(dict(module=module,unit_id=uid,direction=direction,prior=name,prior_source=source,
                            shared_significant_genes=';'.join(sorted(overlap)),n_shared=len(overlap),
                            prior_term_q_family=q.q_family.iloc[0] if len(q) else None,
                            interpretation='shared observed genes only; prior GO significance is separate; marker panels are not GO enrichment tests'))
pd.DataFrame(bridges).to_csv(OUT/'tables/expert_prior_to_module_bridges.csv',index=False)
(OUT/'ontology_context.json').write_text(json.dumps(dict(ontology_paths=len(rows),prior_bridges=len(bridges),relations=['is_a','part_of'],graph_acyclic=True,
    scope='Separate explanatory layers, not added to or used to retrain the Jaccard GNN'),indent=2))
print('Ontology paths',len(rows),'expert bridges',len(bridges))
