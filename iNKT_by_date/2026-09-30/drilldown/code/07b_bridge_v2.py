"""Step 7: per-pathway static gene figures, bridge heatmap, and v2 of the version (b) full network."""
from common_dd import *
import numpy as np, pandas as pd, networkx as nx, matplotlib, textwrap
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','font.family':'DejaVu Sans'})
def save(fig,name,d=FIGD):
    for ext in ['png','svg']: fig.savefig(d/f'{name}.{ext}',dpi=170 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)
mem=pd.read_csv(TAB/'pathway_gene_membership.csv'); E=pd.read_csv(TAB/'gene_comembership_edges.csv')
G_all=nx.Graph(); G_all.add_weighted_edges_from(E.itertuples(index=False,name=None))
FOOT='Exploratory, cell-level; each tissue x condition = 3 pooled mice (n=1 pool). Annotation/enrichment != activation. Edges = co-membership in >=2 significant terms/anchors (annotation, not PPI; no STRING data).'
CS=[('BM C0','BM C0'),('Spleen C3','Spleen C3')]
VM=1.5
DRILLNAME=DRILLNAME
# ---------- bridge heatmap: genes in >=2 drilled KEGG/AP-1 pathways that are DEG somewhere, plus AP-1 members
core=[k for k in DRILLID if k!='ENERGY']
mm=mem[mem.pathway_id.isin(core)]
g_sets=mm.groupby('gene').agg(n=('pathway_id','nunique'),deg=('n_key_comparisons_DEG','max'),det=('detected_in_data','first'))
sel=g_sets[(g_sets.det=='yes')&(((g_sets.n>=2)&(g_sets.deg>=1))|(g_sets.index.isin(mm[mm.pathway_id=='GO:0035976'].gene)))].index
W=mm.drop_duplicates('gene').set_index('gene').loc[sel]
W=W.assign(n=g_sets.loc[sel,'n']).sort_values(['n','n_key_comparisons_DEG'],ascending=False)
cols=[s for _,s in KEY]; L=W[[f'{s} log2FC' for s in cols]].values; D=np.array([[W.at[g,f'{s} DEG']=='yes' for s in cols] for g in W.index])
fig,(a1,a2)=plt.subplots(1,2,figsize=(11,0.3*len(W)+2.5),gridspec_kw=dict(width_ratios=[1,1.2],wspace=.05),sharey=True)
memb=np.array([[g in set(mm[mm.pathway_id==k].gene) for k in core] for g in W.index])
a1.imshow(memb,cmap='Greys',aspect='auto',vmin=0,vmax=1.6); a1.set_xticks(range(len(core))); a1.set_xticklabels([DRILLNAME[k] for k in core],rotation=60,ha='right'); a1.set_yticks(range(len(W))); a1.set_yticklabels(W.index,fontsize=8); a1.set_title('Annotation membership')
im2=a2.imshow(np.clip(L,-VM,VM),cmap='RdBu_r',aspect='auto',vmin=-VM,vmax=VM); a2.set_xticks(range(len(cols))); a2.set_xticklabels(cols,rotation=60,ha='right')
for i,j in zip(*np.where(D)): a2.text(j,i,'*',ha='center',va='center',fontsize=11)
a2.set_title('log2FC (* = DEG)')
fig.suptitle('Genes shared by >=2 anchor pathways (DEG in >=1 key comparison) + AP-1 members',fontsize=11,y=1.0)
fig.subplots_adjust(bottom=0.24,right=0.86)
cax=fig.add_axes([0.88,0.35,0.015,0.35]); cbb=fig.colorbar(im2,cax=cax); cbb.set_label('log2FC (T2 vs Ctrl), clipped +/-%g'%VM,fontsize=8); cbb.ax.tick_params(labelsize=8)
fig.text(.01,.01,FOOT,fontsize=6.5,color='0.3',wrap=True)
save(fig,'bridge_genes_heatmap_v2'); print('heatmap done',len(W))
