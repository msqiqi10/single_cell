"""Step 8: build single-file interactive explorer (vis-network via CDN): pathway network (a/b) + gene network, linked by tabs/hash."""
from common_dd import *
import pandas as pd, numpy as np
aux=json.load(open(TAB/'_gene_net_aux.json')); mem=pd.read_csv(TAB/'pathway_gene_membership.csv')
E=pd.read_csv(TAB/'gene_comembership_edges.csv')
nb=pd.read_csv(RES/'network/excl/nodes.csv',keep_default_na=False,na_values=[''])
tname=dict(zip(nb.term_id,nb.term)); sigset=set(nb.term_id[nb.significant.astype(str)=='True'])
full_by_id={pid:sorted(mem[mem.pathway_id==pid].gene) for pid,_ in DRILL}
def R(x,d=3): return None if x!=x else round(float(x),d)
def terms(version):
    n=pd.read_csv(RES/'network'/version/'nodes.csv',keep_default_na=False,na_values=['']); e=pd.read_csv(RES/'network'/version/'edges.tsv',sep='\t')
    idx={t:i for i,t in enumerate(n.term_id)}; nodes=[]
    for r in n.itertuples():
        genes=full_by_id.get(r.term_id) or (r.measured_genes.split(';') if isinstance(r.measured_genes,str) else [])
        nodes.append(dict(id=r.term_id,term=r.term,src=r.source,anchor=str(r.is_anchor)=='True',sig=str(r.significant)=='True',
            q=None if pd.isna(r.best_q) else float(r.best_q),nl=None if pd.isna(r.best_neglog10q) else round(float(r.best_neglog10q),2),
            comps=(r.sig_comparisons if isinstance(r.sig_comparisons,str) else '').replace('cluster_tissue__','').replace('tissue__','tissue ').replace('_',' '),
            hits=r.hit_genes_union if isinstance(r.hit_genes_union,str) else '',grp=r.gnn_semantic,genes=genes))
    ed=[[idx[a],idx[b],round(float(j),3),int(k)] for a,b,j,k in zip(e.GS_A_ID,e.GS_B_ID,e.JACCARD,e.n_shared)]
    gs=pd.read_csv(RES/'tables'/f'group_summary_{version}.csv'); gs=gs[gs.method=='gnn_semantic']
    return dict(nodes=nodes,edges=ed,groups={r.group:f'{r.group}: {r.representative_term} (n={r.n_terms})' for r in gs.itertuples()})
TV={'b':terms('excl'),'a':terms('all')}
gm=mem.drop_duplicates('gene').set_index('gene')
genes=[]
for g in aux['universe']:
    r=gm.loc[g]; det=r.detected_in_data=='yes'
    genes.append(dict(g=g,det=det,lfc=[R(r[f'{s} log2FC']) if det else None for _,s in KEY],fdr=[None if not det or r[f'{s} FDR']!=r[f'{s} FDR'] else float(f'{r[f"{s} FDR"]:.2g}') for _,s in KEY],
        deg=[r[f'{s} DEG']=='yes' for _,s in KEY],dp=aux['drill_of'].get(g,[]),nsig=int(r.n_significant_nodes_containing),
        tm=[tname[t] for t in aux['gene_terms'].get(g,[]) if t in sigset][:10]))
gi={x['g']:i for i,x in enumerate(genes)}
ed=[[gi[a],gi[b],int(w)] for a,b,w in zip(E.gene_a,E.gene_b,E.n_shared_terms)]
DATA=dict(TV=TV,genes=genes,gedges=ed,comps=[s for _,s in KEY],drill=[dict(id=i,name=n,genes=full_by_id[i]) for i,n in DRILL],caveats=CAVEATS)
tpl=open(Path(__file__).parent/'explorer_template.html').read()
html=tpl.replace('/*DATA*/null',json.dumps(DATA,separators=(',',':')))
(HTMLD/'inkt_network_explorer.html').write_text(html)
print(len(html)/1e6,'MB')
