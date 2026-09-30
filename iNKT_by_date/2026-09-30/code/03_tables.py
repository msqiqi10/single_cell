"""Step 3: anchor x method tables, group tables, k-sensitivity of anchor grouping, driver genes per group x comparison."""
import sys; sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from common import *
import numpy as np, pandas as pd, itertools
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize
T=RES/'tables'
METHODS=['overlap_louvain','semantic_only','golden_fusion','gnn_semantic']
res=pd.read_csv(T/'ORA_all_terms_both_versions.csv.gz',keep_default_na=False,na_values=[''])
def canon(labels,q):
    labels=np.asarray(labels); gs=sorted(set(labels),key=lambda g:(q[labels==g].min(),g)); m={g:i+1 for i,g in enumerate(gs)}; return np.array([m[x] for x in labels])
for version in ['all','excl']:
    vd=RES/'network'/version; n=pd.read_csv(vd/'nodes.csv',keep_default_na=False,na_values=['']); n['best_q']=pd.to_numeric(n.best_q)
    # group tables
    rows=[]
    for m in METHODS+['gnn_only_seed42']:
        for g,d in n.groupby(m):
            d=d.sort_values('best_q',na_position='last'); rep=d.iloc[0]
            rows.append(dict(version=version,method=m,group=g,n_terms=len(d),representative_term=rep.term,representative_id=rep.term_id,representative_q=rep.best_q,
                top_terms=' | '.join(f'{t} (q={q:.2g})' if q==q else f'{t} (anchor, ns)' for t,q in zip(d.term.head(5),d.best_q.head(5))),
                anchors_in_group=' | '.join(d[d.is_anchor].term),n_sources=';'.join(f'{k}:{v}' for k,v in d.source.value_counts().items())))
    G=pd.DataFrame(rows); G.to_csv(T/f'group_summary_{version}.csv',index=False)
    # anchor x method
    A=n[n.is_anchor].copy()
    tab=A[['term_id','term','source','anchor_family','significant','best_q','n_measured_genes']+METHODS+['gnn_only_seed42']].copy()
    for m in METHODS:
        rep=G[G.method==m].set_index('group').representative_term
        tab[m+'_group_rep']=tab[m].map(rep)
    tab.to_csv(T/f'anchor_by_method_{version}.csv',index=False)
    # same-group summary per method across the 4 families + TNF
    fams={'AP-1':['GO:0035976'],'MAPK':['KEGG:MAPK signaling pathway','GO:0000165','GO:0070371','GO:0038066','GO:0007254'],
          'IL-17/Th17':['KEGG:IL-17 signaling pathway','KEGG:Th17 cell differentiation','GO:0072538','GO:0032620'],'TCR':['KEGG:T cell receptor signaling pathway']}
    srows=[]
    for m in METHODS+['gnn_only_seed42']:
        gm=n.set_index('term_id')[m]
        def fg(f): return set(gm[t] for t in fams[f])
        for a,b in itertools.combinations(fams,2):
            srows.append(dict(version=version,method=m,family_A=a,family_B=b,shares_a_group=bool(fg(a)&fg(b)),groups_A=';'.join(sorted(fg(a))),groups_B=';'.join(sorted(fg(b)))))
        srows.append(dict(version=version,method=m,family_A='ALL 4 families',family_B='',shares_a_group=len(set().union(*[fg(f) for f in fams]))==1,groups_A=';'.join(sorted(set().union(*[fg(f) for f in fams]))),groups_B=''))
    pd.DataFrame(srows).to_csv(T/f'anchor_family_cogrouping_{version}.csv',index=False)
    # k sensitivity of anchor grouping (semantic, GOLDEN fusion, GNN+semantic seed42) k=2..15
    ids=n.term_id.to_numpy(); sem=np.load(vd/'semantic_embeddings.npz')['embeddings']; Aj=np.load(vd/'adjacency.npz')['adjacency']
    z=normalize(np.load(vd/'models/gnn_seed42.npz')['embeddings'])
    feats={'semantic_only':sem,'golden_fusion':np.concatenate([.5*sem,.5*Aj],1),'gnn_semantic':np.concatenate([.5*sem,.5*z],1)}
    q=n.best_q.fillna(9).to_numpy(); krows=[]; anc=n.set_index('term_id').loc[[a for a in A.term_id]].term
    for m,f in feats.items():
        for k in range(2,16):
            lab=canon(AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(f),q); s=silhouette_score(f,lab)
            gm=dict(zip(ids,lab)); r=dict(version=version,method=m,k=k,silhouette=s,largest_group=int(np.bincount(lab).max()),
                n_groups_containing_anchors=len({gm[t] for t in A.term_id}),
                AP1_with_MAPK=gm['GO:0035976'] in {gm[t] for t in fams['MAPK']},AP1_with_IL17Th17=gm['GO:0035976'] in {gm[t] for t in fams['IL-17/Th17']},
                AP1_with_TCR=gm['GO:0035976']==gm['KEGG:T cell receptor signaling pathway'],MAPK_with_IL17Th17=bool({gm[t] for t in fams['MAPK']}&{gm[t] for t in fams['IL-17/Th17']}),
                MAPK_with_TCR=gm['KEGG:T cell receptor signaling pathway'] in {gm[t] for t in fams['MAPK']},IL17Th17_with_TCR=gm['KEGG:T cell receptor signaling pathway'] in {gm[t] for t in fams['IL-17/Th17']})
            for t,nm in zip(A.term_id,A.term): r['grp_'+nm]=int(gm[t])
            krows.append(r)
    pd.DataFrame(krows).to_csv(T/f'anchor_k_sensitivity_{version}.csv',index=False)
    # node attributes export
    n.drop(columns=['measured_genes']).to_csv(T/f'node_attributes_{version}.csv',index=False)
    # driver genes per group x comparison x direction (representative + merged members, sig rows only, k>=3)
    R=res[(res.version==version)&(res.q<=.05)&(res.k>=3)]
    mem=pd.read_csv(T/f'redundancy_merges_{version}.csv')
    members={t:[t]+list(mem[mem.representative==t].merged_term_id) for t in ids}
    drows=[]
    for m in METHODS:
        for g,d in n.groupby(m):
            ids_g=[x for t in d.term_id for x in members[t]]; r=R[R.term_id.isin(ids_g)]
            for (u,dr),rr in r.groupby(['unit_id','direction']):
                gs=sorted({x for s in rr.genes for x in str(s).split(';') if x})
                drows.append(dict(version=version,method=m,group=g,unit_id=u,direction=dr,n_driver_genes=len(gs),n_terms=rr.term_id.nunique(),driver_genes=';'.join(gs)))
    pd.DataFrame(drows).to_csv(T/f'group_driver_genes_{version}.csv',index=False)
    edges=pd.read_csv(vd/'edges.tsv',sep='\t'); print(version,len(edges))
