"""Freeze full candidate coverage, GO nodes, and explicit mouse gene-overlap edges."""
from pathlib import Path
import json, hashlib, re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/bio3'
SRC = DATA/'iNKT_by_date/2026-09-19'
OUT = ROOT/'runs/meeting_followup'
for sub in ['tables', 'gene_lists', 'network', 'figures', 'models']:
    (OUT/sub).mkdir(parents=True, exist_ok=True)
def save(frame, name):
    frame.to_csv(OUT/name, index=False)
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
go = pd.read_csv(SRC/'results/tables/GO_ORA_all.csv.gz', keep_default_na=False, low_memory=False)
go = go[(go.analysis=='condition') & go.definition.isin(['original','tissue'])].copy()
catalog = pd.read_csv(SRC/'results/tables/GO_term_catalog.csv').set_index('go_id')
sets = {}
for line in (SRC/'sources/GO_mouse_measured.gmt').read_text().splitlines():
    term, name, *genes = line.split('\t'); sets[term] = set(genes)
save(go, 'tables/all_original_condition_GO.csv.gz')
status = pd.read_csv(DATA/'output/iNKT_meeting_followup_20260905/tables/de_status_original.csv')
save(status, 'tables/all_comparison_status.csv')
full, coverage, membership = [], [], []
focus = 'Dusp1 Fos Jun Junb Jund Fosb Fosl1 Fosl2 Nr4a1 Il1r1 Il6ra Cd8a Ciita Rorc Il23r Ccr6 Il17a Il17f Tbx21 Gata3 Zbtb16 Nkg7 Prf1 Gzma Gzmb Ccl5'.split()
focus_rows=[]
for uid, family in go.groupby('unit_id', sort=True):
    de = pd.read_csv(SRC/'results/de'/f'{uid}.csv.gz')
    assert len(de)==10670 and de.gene.is_unique
    de['unit_id']=uid; de['tissue']=family.tissue.iloc[0]; de['cluster']=family.cluster.iloc[0]
    full.append(de)
    expr=set(de.loc[(de.pct_expressing_tumor>0)|(de.pct_expressing_control>0),'gene'])
    focus_rows.append(de[de.gene.isin(focus)])
    for direction in ['T2_up','Ctrl_up']:
        mask = de.logfoldchanges>=.25 if direction=='T2_up' else de.logfoldchanges<=-.25
        selected=set(de.loc[(de.pvals_adj<=.05)&mask,'gene'])
        (OUT/'gene_lists'/f'{uid}__{direction}.txt').write_text('\n'.join(sorted(selected))+'\n')
        for ns in ['BP','MF','CC']:
            f=family[(family.direction==direction)&(family.namespace==ns)]
            annotated=expr & set().union(*(sets[g] for g in catalog.index[catalog.namespace==ns]))
            tested=set().union(*(sets[g]&expr for g in f.go_id))
            significant=f[f.q_family<=.05]
            hits={g:[] for g in selected}
            for row in significant.itertuples():
                for gene in row.genes.split(';'):
                    if gene: hits[gene].append(row.go_id)
            covered={g for g, terms in hits.items() if terms}
            globally=set(g for value in f.loc[f.q_analysis_global<=.05,'genes'] for g in value.split(';') if g)
            coverage.append(dict(unit_id=uid,tissue=family.tissue.iloc[0],cluster=family.cluster.iloc[0],direction=direction,namespace=ns,
                n_selected=len(selected),n_annotated=len(selected&annotated),n_in_tested_terms=len(selected&tested),
                n_in_significant_terms=len(covered),n_in_global_significant_terms=len(globally),
                n_without_branch_annotation=len(selected-annotated),n_not_in_significant_terms=len(selected-covered),
                n_significant_terms=len(significant),n_tested_terms=len(f)))
            for gene in sorted(selected):
                membership.append(dict(unit_id=uid,direction=direction,namespace=ns,gene=gene,
                    annotated=gene in annotated,in_tested_terms=gene in tested,
                    significant_go_ids=';'.join(sorted(hits[gene])),n_significant_terms=len(hits[gene])))
    print(uid, 'complete', flush=True)
save(pd.concat(full), 'tables/all_20_condition_DE.csv.gz')
save(pd.concat(focus_rows), 'tables/Rob_focus_all_20_comparisons.csv')
save(pd.DataFrame(coverage), 'tables/full_candidate_GO_coverage.csv')
save(pd.DataFrame(membership), 'tables/candidate_to_GO_membership.csv.gz')
# Hypothesis-name retrieval is an audit, not an independently significant gene set.
pattern = r'MAPK|mitogen-activated|AP-1|interleukin-17|T-helper 17|cytotoxic|cell killing|interferon|ribosom|translation|oxidative phosphorylation'
audit=go[go.term.str.contains(pattern,case=False,regex=True)].copy()
save(audit,'tables/hypothesis_GO_all_results.csv.gz')
# Nodes are all BP terms significant in at least one original condition comparison.
sig=go[(go.namespace=='BP')&(go.q_family<=.05)]
ids=sorted(sig.go_id.unique())
obo={}; current=None
for line in (SRC/'sources/go-basic.obo').read_text().splitlines():
    if line=='[Term]': current={}
    elif line.startswith('['): current=None
    elif current is not None:
        if line.startswith('id: '): current['id']=line[4:]; obo[line[4:]]=current
        elif line.startswith('def: '): current['definition']=line.split('"')[1]
        elif line.startswith('name: '): current['name']=line[6:]
nodes=[]
for gid in ids:
    rows=sig[sig.go_id==gid]
    nodes.append(dict(GOID=gid,name=catalog.loc[gid,'term'],DESCRIPTION=obo[gid].get('definition',catalog.loc[gid,'term']),
        measured_members=';'.join(sorted(sets[gid])),n_measured_members=len(sets[gid]),
        significant_comparisons=rows.unit_id.nunique(),min_q_family=rows.q_family.min()))
save(pd.DataFrame(nodes),'network/mouse_GO_nodes.csv')
save(sig,'network/node_condition_evidence.csv')
matrix=np.eye(len(ids)); edges=[]
for i,a in enumerate(ids):
    for j in range(i+1,len(ids)):
        b=ids[j]; shared=sets[a]&sets[b]; jac=len(shared)/len(sets[a]|sets[b])
        matrix[i,j]=matrix[j,i]=jac
        if jac>=.25 and len(shared)>=3:
            edges.append(dict(GS_A_ID=a,GS_B_ID=b,JACCARD=jac,n_shared=len(shared),shared_measured_genes=';'.join(sorted(shared))))
save(pd.DataFrame(edges),'network/mouse_GO_edges.tsv')
# GNN loader requires a tab-delimited edge file.
pd.DataFrame(edges).to_csv(OUT/'network/mouse_GO_edges.tsv',sep='\t',index=False)
adj=np.zeros_like(matrix);node_index={gid:i for i,gid in enumerate(ids)}
for edge in edges:
    i,j=node_index[edge['GS_A_ID']],node_index[edge['GS_B_ID']]
    adj[i,j]=adj[j,i]=edge['JACCARD']
np.savez_compressed(OUT/'network/overlap_matrices.npz',ID=np.array(ids),jaccard=matrix,adjacency=adj)
assert all((sets[r['GS_A_ID']]&sets[r['GS_B_ID']])==set(r['shared_measured_genes'].split(';')) for r in edges)
clusters={}
tree=linkage(squareform(1-matrix,checks=True),method='complete')
for minimum in [.25,.5,.75]:
    labels=fcluster(tree,t=1-minimum,criterion='distance')
    clusters[f'jaccard_{minimum}']=labels
save(pd.DataFrame({'GOID':ids,**clusters}),'network/overlap_only_complete_linkage.csv')
manifest=dict(stage='full candidate audit and mouse GO network',comparisons=go.unit_id.nunique(),nodes=len(ids),edges=len(edges),
    candidate_rule='gene BH FDR<=0.05 and abs(log2FC)>=0.25; separate directions',
    node_rule='BP q_family<=0.05 in at least one of 20 original/tissue condition comparisons; stable and identity comparisons excluded from node selection',
    edge_rule='Jaccard of complete measured gene members >=0.25 and at least 3 shared genes; not a regulatory or causal edge',
    clustering_baseline='complete linkage on 1-Jaccard, thresholds 0.25,0.5,0.75; exploratory redundancy grouping',
    membership='complete frozen measured gene sets, not query-only overlap',
    input_hashes={str(p.relative_to(DATA)):sha(p) for p in [SRC/'sources/GO_mouse_measured.gmt',SRC/'sources/go-basic.obo',SRC/'results/tables/GO_ORA_all.csv.gz']})
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest,indent=2))
