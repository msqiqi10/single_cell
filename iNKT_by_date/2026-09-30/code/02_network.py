"""Step 2: node set (significant + anchors), Jaccard>=0.75 collapse, edges, embeddings, GNN runs, 4 grouping methods, ablation metrics."""
import sys,json,subprocess,os; sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from common import *
import numpy as np, pandas as pd, networkx as nx, collections, itertools
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, adjusted_rand_score, roc_auc_score, average_precision_score
from sklearn.preprocessing import normalize
res=pd.read_csv(RES/'tables/ORA_all_terms_both_versions.csv.gz',keep_default_na=False,na_values=[''])
de=pd.read_csv(P25/'results/tables/all_20_condition_DE.csv.gz',usecols=['gene']); measured0=set(de.gene)
cat=pd.read_csv(P19/'results/tables/GO_term_catalog.csv').set_index('go_id')
gmt={}
for l in (P19/'sources/GO_mouse_measured.gmt').read_text().splitlines():
    t,n,*g=l.split('\t'); gmt[t]=set(g)
up=collections.defaultdict(set)
for g in measured0: up[g.upper()].add(g)
names={t:cat.loc[t,'term'] for t in gmt}; sets=dict(gmt); srcof={t:cat.loc[t,'namespace'] for t in gmt}
for l in KEGG.read_text().splitlines():
    n,_,*g=l.split('\t'); sets['KEGG:'+n]=set().union(*(up.get(x,set()) for x in g)); names['KEGG:'+n]=n; srcof['KEGG:'+n]='KEGG'
# GO definitions and names for terms absent from measured catalogue (interleukin-17 production)
obo={};cur=None
for line in (P19/'sources/go-basic.obo').read_text().splitlines():
    if line=='[Term]': cur={}
    elif line.startswith('['): cur=None
    elif cur is not None:
        if line.startswith('id: '): cur['id']=line[4:]; obo[line[4:]]=cur
        elif line.startswith('def: '): cur['def']=line.split('"')[1]
        elif line.startswith('name: '): cur['name']=line[6:]
for t,v in obo.items():
    names.setdefault(t,v.get('name',t))
go_name2id={}
for t in gmt: go_name2id.setdefault(names[t],t)
go_name2id['interleukin-17 production']='GO:0032620'; srcof['GO:0032620']='BP'; sets.setdefault('GO:0032620',set())
anchors={}
for n in ANCHOR_KEGG: anchors['KEGG:'+n]=('KEGG')
for a in ANCHOR_GO:
    t=a if a.startswith('GO:') else go_name2id[a]; anchors[t]=srcof[t]
fam={'GO:0035976':'AP-1','KEGG:MAPK signaling pathway':'MAPK','GO:0000165':'MAPK','GO:0070371':'MAPK','GO:0038066':'MAPK','GO:0007254':'MAPK',
     'KEGG:IL-17 signaling pathway':'IL-17/Th17','KEGG:Th17 cell differentiation':'IL-17/Th17','GO:0072538':'IL-17/Th17','GO:0032620':'IL-17/Th17',
     'KEGG:T cell receptor signaling pathway':'TCR','KEGG:TNF signaling pathway':'TNF'}
def srclabel(t,is_anchor_cc=False):
    s=srcof[t]; return {'BP':'GO BP','MF':'GO MF','CC':'GO CC anchor','KEGG':'KEGG'}[s]
def semtext(t):
    return names[t] if t.startswith('KEGG:') else names[t]+'. '+obo.get(t,{}).get('def',names[t])
def canon(labels,order_key):
    """relabel groups by best (smallest) q of members: order_key gives per-node sortable value"""
    labels=np.asarray(labels); groups=sorted(set(labels),key=lambda g:(min(order_key[labels==g]),g)); m={g:i+1 for i,g in enumerate(groups)}
    return np.array([m[x] for x in labels])
def pick_k(features,kmax=15):
    rows=[];sol={}
    for k in range(2,min(kmax,len(features)-1)+1):
        lab=AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(features)
        rows.append(dict(k=k,silhouette=silhouette_score(features,lab),min_size=int(np.bincount(lab).min()))); sol[k]=lab
    d=pd.DataFrame(rows); k=int(d.sort_values(['silhouette','k'],ascending=[False,True]).iloc[0].k); return k,d,sol
def modularity(G,labels,ids):
    comm=[set(ids[labels==g]) for g in set(labels)]
    return nx.community.modularity(G,comm,weight='weight')
from sentence_transformers import SentenceTransformer
model=SentenceTransformer(str(MINILM),device='cpu',local_files_only=True)
summary={}
for version in ['all','excl']:
    vd=OUT/'results'/'network'/version; vd.mkdir(parents=True,exist_ok=True)
    keep=(lambda g:True) if version=='all' else (lambda g:not EXCL.match(g))
    R=res[res.version==version]
    S={t:{g for g in gs if keep(g)} for t,gs in sets.items()}
    sig=R[(R.q<=.05)&(R.k>=3)&R.source.isin(['BP','MF','KEGG'])]
    node_ids=set(sig.term_id)|set(anchors)
    rows=[]
    for t in node_ids:
        allr=R[R.term_id==t]; s=allr[(allr.q<=.05)&(allr.k>=3)]
        is_sig=len(s)>0
        best=(s if is_sig else allr).q.min() if len(allr) else np.nan
        genes=sorted({g for x in s.genes for g in str(x).split(';') if g})
        rows.append(dict(term_id=t,term=names[t],source=srclabel(t),is_anchor=t in anchors,anchor_family=fam.get(t,''),significant=is_sig,
            best_q=best,best_neglog10q=-np.log10(max(best,1e-300)) if best==best else np.nan,
            sig_comparisons=';'.join(sorted(f'{u}|{d}' for u,d in zip(s.unit_id,s.direction))),n_sig_comparisons=len(s),
            hit_genes_union=';'.join(genes),n_measured_genes=len(S[t]),measured_genes=';'.join(sorted(S[t]))))
    nodes=pd.DataFrame(rows).sort_values(['best_q','term_id'],na_position='last').reset_index(drop=True)
    n_pre=len(nodes)
    # --- redundancy collapse (Jaccard>=.75 on measured sets). Anchors are protected: never absorbed into another node, and may absorb others.
    def jac(a,b):
        u=len(S[a]|S[b]); return len(S[a]&S[b])/u if u else 0.
    order=list(nodes[nodes.is_anchor].term_id)+list(nodes[~nodes.is_anchor].term_id)
    order=sorted(order,key=lambda t:(not (t in anchors),nodes.set_index('term_id').best_q.get(t,np.nan) if nodes.set_index('term_id').best_q.get(t,np.nan)==nodes.set_index('term_id').best_q.get(t,np.nan) else 9))
    assigned={};members=collections.defaultdict(list)
    qmap=nodes.set_index('term_id').best_q.fillna(9).to_dict()
    order=sorted(node_ids,key=lambda t:(t not in anchors,qmap[t],t))
    for t in order:
        if t in assigned: continue
        assigned[t]=t
        for u in order:
            if u in assigned or u in anchors: continue
            if jac(t,u)>=.75: assigned[u]=t; members[t].append(u)
    reps=[t for t in order if assigned[t]==t]
    merge=[dict(representative=r,representative_term=names[r],merged_term_id=u,merged_term=names[u],jaccard=jac(r,u),merged_source=srclabel(u)) for r in reps for u in members[r]]
    pd.DataFrame(merge).to_csv(OUT/'results/tables'/f'redundancy_merges_{version}.csv',index=False)
    nodes=nodes[nodes.term_id.isin(reps)].copy()
    nodes['merged_members']=nodes.term_id.map(lambda r:' | '.join(f'{u} {names[u]}' for u in members[r]))
    nodes['n_merged']=nodes.term_id.map(lambda r:len(members[r]))
    nodes=nodes.sort_values(['best_q','term_id'],na_position='last').reset_index(drop=True)
    ids=nodes.term_id.to_numpy(); N=len(ids); qkey=nodes.best_q.fillna(9).to_numpy()
    # --- edges
    edges=[];A=np.zeros((N,N))
    for i,j in itertools.combinations(range(N),2):
        sh=S[ids[i]]&S[ids[j]]; jj=jac(ids[i],ids[j])
        if jj>=.25 and len(sh)>=3:
            edges.append(dict(GS_A_ID=ids[i],GS_B_ID=ids[j],JACCARD=jj,n_shared=len(sh),shared_measured_genes=';'.join(sorted(sh)))); A[i,j]=A[j,i]=jj
    E=pd.DataFrame(edges); E.to_csv(vd/'edges.tsv',sep='\t',index=False)
    G=nx.Graph(); G.add_nodes_from(range(N)); G.add_weighted_edges_from([(i,j,A[i,j]) for i,j in zip(*np.nonzero(np.triu(A)))])
    # --- semantic embeddings
    sem=model.encode([semtext(t) for t in ids],normalize_embeddings=True,show_progress_bar=False)
    np.savez_compressed(vd/'semantic_embeddings.npz',ID=ids.astype(str),embeddings=sem); np.savez_compressed(vd/'adjacency.npz',ID=ids.astype(str),adjacency=A)
    # --- (i) Louvain
    comm=nx.community.louvain_communities(G,weight='weight',seed=42,resolution=1.0)
    lou=np.zeros(N,int)
    for c,ms in enumerate(comm): lou[list(ms)]=c
    # --- (ii) semantic, (iii) GOLDEN fusion
    k_sem,d_sem,sol_sem=pick_k(sem); feats=np.concatenate([.5*sem,.5*A],axis=1); k_fus,d_fus,sol_fus=pick_k(feats)
    d_sem.assign(method='semantic').append if False else None
    # --- (iv) GNN
    seeds=[42,43,44]; gnn_feats={}; metr=[]; gnn_sol={}
    (vd/'models').mkdir(exist_ok=True)
    for sd in seeds:
        pre=vd/'models'/f'gnn_seed{sd}'
        subprocess.run([VENV_PY,str(OUT/'code/gnn_mouse_holdout.py'),'--edge_file',str(vd/'edges.tsv'),'--feature_npz',str(vd/'semantic_embeddings.npz'),'--weight_column','JACCARD','--edge_weight_transform','none','--hidden_dim','256','--num_layers','2','--epochs','120','--device','cpu','--seed',str(sd),'--split_seed',str(sd),'--output_csv',str(pre)+'.csv','--output_npz',str(pre)+'.npz'],check=True,capture_output=True,env={**os.environ,'CUDA_VISIBLE_DEVICES':''})
        r=json.loads(open(str(pre)+'.csv'.replace('.csv','')+'.summary.json').read()) if os.path.exists(str(pre)+'.summary.json') else json.loads(open(str(pre)+'.csv.summary.json').read())
        base=json.loads(open(str(pre)+'.csv.baseline.json').read()); sp=np.load(str(pre)+'.csv.splits.npz'); assert np.array_equal(sp['ID'],ids.astype(str))
        g=nx.Graph(); g.add_nodes_from(range(N)); g.add_edges_from(sp['train_positive'])
        pairs=np.vstack([sp['test_positive'],sp['test_negative']]); y=np.r_[np.ones(len(sp['test_positive'])),np.zeros(len(sp['test_negative']))]
        cn=np.array([len(set(g[u])&set(g[v])) for u,v in pairs]); aa=np.array([sum(1/np.log(g.degree(w)) for w in set(g[u])&set(g[v])) for u,v in pairs])
        fa=r['final_at_best']
        metr.append(dict(version=version,seed=seed if False else sd,n_test_pos=len(sp['test_positive']),n_test_neg=len(sp['test_negative']),GNN_AUC=fa['test_auc'],GNN_AP=fa['test_ap'],
            semantic_cosine_AUC=base['semantic_cosine_test_auc'],common_neighbours_AUC=roc_auc_score(y,cn),Adamic_Adar_AUC=roc_auc_score(y,aa),
            GNN_AP_=fa['test_ap'],semantic_cosine_AP=base['semantic_cosine_test_ap'],common_neighbours_AP=average_precision_score(y,cn),Adamic_Adar_AP=average_precision_score(y,aa),best_epoch=r['best']['epoch']))
        z=np.load(str(pre)+'.npz'); assert np.array_equal(z['ID'],ids.astype(str))
        gnn_feats[sd]=np.concatenate([.5*sem,.5*normalize(z['embeddings'])],axis=1)
    pd.DataFrame(metr).drop(columns='GNN_AP_').to_csv(OUT/'results/tables'/f'gnn_vs_baselines_{version}.csv',index=False)
    gk={sd:pick_k(gnn_feats[sd]) for sd in seeds}
    k_gnn=gk[42][0]  # k from reference seed 42 as in 2026-09-25
    lab={'overlap_louvain':lou,'semantic_only':sol_sem[k_sem],'golden_fusion':sol_fus[k_fus],'gnn_semantic':gk[42][2][k_gnn]}
    lab={m:canon(v,qkey) for m,v in lab.items()}
    seed_lab={sd:canon(AgglomerativeClustering(n_clusters=k_gnn,linkage='ward').fit_predict(gnn_feats[sd]),qkey) for sd in seeds}
    ks=[]
    for nm,d in [('semantic_only',d_sem),('golden_fusion',d_fus)]+[(f'gnn_semantic_seed{sd}',gk[sd][1]) for sd in seeds]: ks.append(d.assign(method=nm))
    pd.concat(ks).to_csv(OUT/'results/tables'/f'k_selection_{version}.csv',index=False)
    # add-on: additional GNN-only (no semantic) partition to show what the graph embedding alone gives
    z42=normalize(np.load(str(vd/'models/gnn_seed42.npz'))['embeddings']); kz,_,solz=pick_k(z42); lab['gnn_only_seed42']=canon(solz[kz],qkey)
    for m,v in lab.items(): nodes[m]=[f'{m[:3].upper() if False else {"overlap_louvain":"L","semantic_only":"S","golden_fusion":"F","gnn_semantic":"G","gnn_only_seed42":"Z"}[m]}{x:02d}' for x in v]
    for sd in seeds: nodes[f'gnn_semantic_seed{sd}']=[f'G{x:02d}' for x in seed_lab[sd]]
    nodes.to_csv(vd/'nodes.csv',index=False)
    # ARIs
    meths=['overlap_louvain','semantic_only','golden_fusion','gnn_semantic']
    ari=pd.DataFrame([[adjusted_rand_score(lab[a],lab[b]) for b in meths] for a in meths],index=meths,columns=meths); ari.to_csv(OUT/'results/tables'/f'ARI_between_methods_{version}.csv')
    sa=pd.DataFrame([[adjusted_rand_score(seed_lab[a],seed_lab[b]) for b in seeds] for a in seeds],index=seeds,columns=seeds); sa.to_csv(OUT/'results/tables'/f'ARI_across_gnn_seeds_{version}.csv')
    # GNN vs simpler: seed-wise ARI of each method's partition vs gnn seed labels, plus per-seed own-k
    quality=[]
    semmat=sem@sem.T
    for m,v in list(lab.items()):
        ids_arr=np.asarray(v)
        within=np.mean([semmat[np.ix_(ids_arr==g,ids_arr==g)][np.triu_indices((ids_arr==g).sum(),1)].mean() for g in set(ids_arr) if (ids_arr==g).sum()>1])
        quality.append(dict(version=version,method=m,n_groups=len(set(v)),n_singletons=int((np.bincount(v)==1).sum()),modularity_overlap_graph=modularity(G,np.asarray(v),np.arange(N)),mean_within_group_semantic_cosine=within,
            k_rule='Louvain resolution 1, seed 42' if m=='overlap_louvain' else 'Ward, max silhouette k=2..15'))
    pd.DataFrame(quality).to_csv(OUT/'results/tables'/f'partition_quality_{version}.csv',index=False)
    pd.DataFrame({'k_own_seed42':[gk[42][0]],'k_own_seed43':[gk[43][0]],'k_own_seed44':[gk[44][0]]}).to_csv(vd/'gnn_own_k.csv',index=False)
    summary[version]=dict(n_significant_terms_before_collapse=int(len(sig.term_id.unique())),n_nodes_before_collapse=n_pre,n_nodes=N,n_merged_away=n_pre-N,n_edges=len(E),n_isolates=int(sum(1 for i in range(N) if G.degree(i)==0)),
        n_anchor_nodes=int(nodes.is_anchor.sum()),n_anchor_significant=int(nodes[nodes.is_anchor].significant.sum()),
        groups={m:int(len(set(v))) for m,v in lab.items()},k_gnn=k_gnn,k_sem=k_sem,k_fusion=k_fus)
    print(version,json.dumps(summary[version]))
json.dump(summary,open(OUT/'results/network/network_summary.json','w'),indent=2)
