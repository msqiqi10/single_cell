"""Compare learned graph representations to simple baselines and export traceability."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import roc_auc_score, average_precision_score, silhouette_score, adjusted_rand_score
from sklearn.preprocessing import normalize
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup'
nodes=pd.read_csv(OUT/'network/fusion_nodes.csv');ids=nodes.GOID.to_numpy()
semantic=np.load(OUT/'network/semantic_embeddings.npz')['embeddings']
metrics=[];ks=[];solutions={};representations={}
for seed in [42,43,44]:
    prefix=OUT/'models'/f'gnn_seed{seed}'
    result=json.loads(prefix.with_suffix('.summary.json').read_text())
    base=json.loads(Path(str(prefix)+'.csv.baseline.json').read_text())
    splits=np.load(str(prefix)+'.csv.splits.npz');assert np.array_equal(splits['ID'],ids)
    graph=nx.Graph();graph.add_nodes_from(range(len(ids)));graph.add_edges_from(splits['train_positive'])
    pairs=np.vstack([splits['test_positive'],splits['test_negative']]);truth=np.r_[np.ones(len(splits['test_positive'])),np.zeros(len(splits['test_negative']))]
    # The baseline can use training edges only, just as the encoder does.
    common=np.array([len(set(graph[u])&set(graph[v])) for u,v in pairs])
    aa=np.array([sum(1/np.log(graph.degree(w)) for w in set(graph[u])&set(graph[v])) for u,v in pairs])
    metrics.append(dict(seed=seed,n_test_positive=len(splits['test_positive']),n_test_negative=len(splits['test_negative']),
        GNN_AUC=result['final_at_best']['test_auc'],GNN_AP=result['final_at_best']['test_ap'],
        semantic_AUC=base['semantic_cosine_test_auc'],semantic_AP=base['semantic_cosine_test_ap'],
        common_neighbors_AUC=roc_auc_score(truth,common),common_neighbors_AP=average_precision_score(truth,common),
        Adamic_Adar_AUC=roc_auc_score(truth,aa),Adamic_Adar_AP=average_precision_score(truth,aa),best_epoch=result['best']['epoch']))
    z=np.load(prefix.with_suffix('.npz'));assert np.array_equal(z['ID'],ids)
    features=np.concatenate([.5*semantic,.5*normalize(z['embeddings'])],axis=1)
    representations[seed]=features
    for k in range(2,13):
        labels=AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(features)
        solutions[seed,k]=labels
        ks.append(dict(seed=seed,k=k,silhouette=silhouette_score(features,labels),min_size=int(np.bincount(labels).min())))
pd.DataFrame(metrics).to_csv(OUT/'network/gnn_baseline_comparison.csv',index=False)
ks=pd.DataFrame(ks);ks.to_csv(OUT/'network/gnn_fusion_k_sensitivity.csv',index=False)
# k is chosen from seed 42 without choosing the best-looking seed or inspecting biological labels.
k=int(ks[ks.seed==42].sort_values(['silhouette','k'],ascending=[False,True]).iloc[0].k)
ref=solutions[42,k];order=sorted(set(ref),key=lambda x:min(ids[ref==x]));mapping={c:f'G{i+1:02d}' for i,c in enumerate(order)}
nodes['gnn_module']=[mapping[x] for x in ref]
nodes.to_csv(OUT/'network/annotated_nodes.csv',index=False)
pd.DataFrame([dict(seed=seed,k=k,ARI_to_seed42=adjusted_rand_score(ref,solutions[seed,k])) for seed in [42,43,44]]).to_csv(OUT/'network/gnn_partition_stability.csv',index=False)
evidence=pd.read_csv(OUT/'network/node_condition_evidence.csv',keep_default_na=False)
evidence=evidence.merge(nodes[['GOID','fusion_module','gnn_module']],left_on='go_id',right_on='GOID',validate='many_to_one')
evidence.to_csv(OUT/'tables/module_to_GO_evidence.csv',index=False)
summaries=[];gene_rows=[]
de=pd.read_csv(OUT/'tables/all_20_condition_DE.csv.gz').set_index(['unit_id','gene'])
for method in ['fusion_module','gnn_module']:
    for (module,uid,direction),group in evidence.groupby([method,'unit_id','direction']):
        genes=sorted({gene for text in group.genes for gene in text.split(';') if gene})
        summaries.append(dict(method=method,module=module,unit_id=uid,tissue=group.tissue.iloc[0],cluster=group.cluster.iloc[0],direction=direction,
                              n_GO=len(group),n_unique_drivers=len(genes),min_term_q=group.q_family.min(),
                              genes=';'.join(genes),GO_ids=';'.join(sorted(group.go_id)),
                              representative_terms=' | '.join(group.sort_values('q_family').term.head(4))))
        for gene in genes:
            row=de.loc[uid,gene]
            hit=group[group.genes.map(lambda text:gene in text.split(';'))]
            gene_rows.append(dict(method=method,module=module,unit_id=uid,direction=direction,gene=gene,
                                 log2FC=row.logfoldchanges,gene_FDR=row.pvals_adj,
                                 pct_T2=row.pct_expressing_tumor,pct_Ctrl=row.pct_expressing_control,
                                 GO_ids=';'.join(sorted(hit.go_id))))
pd.DataFrame(summaries).to_csv(OUT/'tables/module_by_condition.csv',index=False)
pd.DataFrame(gene_rows).to_csv(OUT/'tables/module_gene_evidence.csv.gz',index=False)
edges=pd.read_csv(OUT/'network/mouse_GO_edges.tsv',sep='\t')
edges=edges.merge(nodes[['GOID','gnn_module']],left_on='GS_A_ID',right_on='GOID').rename(columns={'gnn_module':'module_A'}).drop(columns='GOID')
edges=edges.merge(nodes[['GOID','gnn_module']],left_on='GS_B_ID',right_on='GOID').rename(columns={'gnn_module':'module_B'}).drop(columns='GOID')
module_edges=[]
for a in sorted(nodes.gnn_module.unique()):
    for b in sorted(nodes.gnn_module.unique()):
        if a>=b:continue
        part=edges[((edges.module_A==a)&(edges.module_B==b))|((edges.module_A==b)&(edges.module_B==a))]
        if len(part):module_edges.append(dict(module_A=a,module_B=b,n_GO_edges=len(part),sum_jaccard=part.JACCARD.sum(),
            shared_measured_genes=';'.join(sorted({g for text in part.shared_measured_genes for g in text.split(';')}))))
pd.DataFrame(module_edges).to_csv(OUT/'network/module_edges.csv',index=False)
validation=dict(passed=True,gnn_k=k,seeds=[42,43,44],node_alignment=True,
    external_validation=False,scope='Internal transductive reconstruction of a small GO overlap graph, balanced held-out edge/non-edge samples; not biological validation or a method benchmark',
    gnn_fusion='concat(0.5 * L2 semantic, 0.5 * L2 GNN); exploratory extension; reference seed42, k maximizes silhouette2..12',
    metrics=metrics)
(OUT/'network_validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2))
