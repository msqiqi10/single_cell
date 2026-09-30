from pathlib import Path
import textwrap,json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup';FIG=OUT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
short={'bone_marrow':'BM','spleen':'Spleen','thymus':'Thymus'}
def label(t,c):return f'{short[t]} / {c}'
def save(fig,name):
    for ext in ['png','svg','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)
cov=pd.read_csv(OUT/'tables/full_candidate_GO_coverage.csv');cov=cov[(cov.namespace=='BP')&(cov.direction=='T2_up')&(cov.cluster!='ALL')]
fig,ax=plt.subplots(figsize=(12,8));y=np.arange(len(cov))
v=[cov.n_in_significant_terms,cov.n_in_tested_terms-cov.n_in_significant_terms,cov.n_annotated-cov.n_in_tested_terms,cov.n_selected-cov.n_annotated]
left=np.zeros(len(cov));names=['In significant BP terms','Tested BP annotation, no significant hit','BP annotation outside term-size filter','No BP annotation in frozen background'];colors=['#137C8B','#8EBBC4','#D9B96E','#D7DCE2']
for values,name,color in zip(v,names,colors):ax.barh(y,values,left=left,label=name,color=color);left+=values
for i,r in enumerate(cov.itertuples()):ax.text(r.n_selected+3,i,f'{r.n_in_significant_terms}/{r.n_selected}',va='center',fontsize=9)
ax.set_yticks(y,[label(r.tissue,r.cluster) for r in cov.itertuples()]);ax.invert_yaxis();ax.set_xlim(0,max(cov.n_selected)*1.18);ax.set_xlabel('Unique T2-up genes (gene BH FDR ≤ 0.05; log2FC ≥ 0.25)')
ax.set_title('Full candidate coverage across all 17 eligible tissue / cell-cluster comparisons',loc='left',pad=18)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.1),ncol=2,frameon=False,fontsize=9)
fig.subplots_adjust(bottom=.25)
fig.text(.12,.01,'Numbers show significant-BP-covered genes / all selected genes. Zero selected genes does not establish biological equivalence.',fontsize=9)
save(fig,'01_full_candidate_coverage')
nodes=pd.read_csv(OUT/'network/annotated_nodes.csv');mods=sorted(nodes.gnn_module.unique())
modnames={'G01':'Mixed annotations / stress; inspect drivers','G02':'Translation / ribosome biogenesis','G03':'Shared ribosomal annotations; not muscle identity','G04':'Respiration / ATP production','G05':'Respiratory-chain assembly','G06':'Nucleotide metabolism / overlapping energy genes','G07':'Ion / proton transport'}
module=pd.read_csv(OUT/'tables/module_by_condition.csv');module=module[(module.method=='gnn_module')&(module.direction=='T2_up')&(module.cluster!='ALL')]
uids=cov.unit_id.tolist();matrix=module.pivot(index='unit_id',columns='module',values='n_unique_drivers').reindex(index=uids,columns=mods).fillna(0)
fig,ax=plt.subplots(figsize=(12,8));im=ax.imshow(matrix,cmap='YlGnBu',aspect='auto',vmin=0)
ax.set_yticks(range(len(uids)),[label(r.tissue,r.cluster) for r in cov.itertuples()]);ax.set_xticks(range(len(mods)),mods)
for i in range(len(uids)):
    for j in range(len(mods)):
        value=int(matrix.iloc[i,j]);ax.text(j,i,str(value) if value else '—',ha='center',va='center',color='white' if value>25 else '#223344',fontsize=9)
fig.colorbar(im,ax=ax,label='Unique T2-up drivers in significant BP terms',shrink=.7,pad=.02)
ax.set_title('Functional groups mapped back to tissue and cell clusters',loc='left',pad=16)
fig.subplots_adjust(bottom=.28,right=.84)
fig.text(.13,.04,'\n'.join(f'{m}  {modnames[m]}' for m in mods),fontsize=9,linespacing=1.3)
fig.text(.13,.235,'Genes can contribute to several groups: columns must not be summed. A dash means no significant BP-term evidence.',fontsize=9)
save(fig,'02_function_groups_by_cluster')
edges=pd.read_csv(OUT/'network/mouse_GO_edges.tsv',sep='\t');graph=nx.Graph();graph.add_nodes_from(nodes.GOID)
for r in edges.itertuples():graph.add_edge(r.GS_A_ID,r.GS_B_ID,weight=r.JACCARD)
palette=dict(zip(mods,plt.get_cmap('tab10').colors));which=nodes.set_index('GOID').gnn_module.to_dict()
fig,(ax,bx)=plt.subplots(1,2,figsize=(14,7),gridspec_kw={'width_ratios':[1.35,1]})
pos=nx.spring_layout(graph,seed=42,iterations=200,weight='weight',k=.3)
nx.draw_networkx_edges(graph,pos,ax=ax,alpha=.14,width=.7)
nx.draw_networkx_nodes(graph,pos,ax=ax,node_size=45,node_color=[palette[which[g]] for g in graph],edgecolors='white',linewidths=.3)
ax.set_title('97 GO BP nodes; 577 shared-gene edges',loc='left');ax.axis('off')
mg=nx.Graph();mg.add_nodes_from(mods);me=pd.read_csv(OUT/'network/module_edges.csv')
for r in me.itertuples():mg.add_edge(r.module_A,r.module_B,weight=r.n_GO_edges)
mp=nx.circular_layout(mg)
nx.draw_networkx_edges(mg,mp,ax=bx,width=[1+np.log1p(mg.edges[e]['weight']) for e in mg.edges],alpha=.25)
nx.draw_networkx_nodes(mg,mp,ax=bx,node_size=[450+40*(nodes.gnn_module==m).sum() for m in mods],node_color=[palette[m] for m in mods])
nx.draw_networkx_labels(mg,mp,ax=bx,font_size=11,font_color='white',font_weight='bold')
bx.set_title('Seven exploratory GNN + semantic groups',loc='left');bx.axis('off');bx.margins(.2)
fig.subplots_adjust(bottom=.25)
fig.text(.08,.03,'\n'.join(f'{m}: {modnames[m]}' for m in mods),fontsize=9,linespacing=1.2)
fig.text(.08,.205,'Edges represent overlapping measured gene sets, not regulation or causality. Group G01 is heterogeneous.',fontsize=10)
save(fig,'03_function_network')
de=pd.read_csv(OUT/'tables/Rob_focus_all_20_comparisons.csv');panel='Fos Jun Junb Jund Fosb Fosl2 Dusp1 Nr4a1 Il1r1 Il6ra Cd8a Ciita'.split()
d=de.pivot(index='unit_id',columns='gene',values='logfoldchanges').reindex(index=uids,columns=panel);q=de.pivot(index='unit_id',columns='gene',values='pvals_adj').reindex(index=uids,columns=panel)
fig,ax=plt.subplots(figsize=(13,8));im=ax.imshow(d,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
ax.set_xticks(range(len(panel)),panel,rotation=45,ha='right');ax.set_yticks(range(len(uids)),[label(r.tissue,r.cluster) for r in cov.itertuples()])
for i in range(len(uids)):
    for j in range(len(panel)):
        if q.iloc[i,j]<=.05 and abs(d.iloc[i,j])>=.25:ax.text(j,i,'*',ha='center',va='center',fontsize=14,color='black')
fig.colorbar(im,ax=ax,label='log2FC, T2 / Ctrl (color clipped at ±2)',shrink=.75)
ax.set_title('Rob\'s candidate panel: gene effects across all eligible cell clusters',loc='left',pad=18)
fig.text(.12,.005,'* Existing gene BH FDR ≤ 0.05 and |log2FC| ≥ 0.25. Full effects and detection fractions are retained in the evidence table.',fontsize=9)
save(fig,'04_Rob_gene_panel')
raw=pd.read_csv(OUT/'tables/full_feature_marker_by_group.csv');raw=raw[raw.cluster=='C5-2'];tissues=['bone_marrow','spleen','thymus']
fig,axes=plt.subplots(1,3,figsize=(13,5),sharey=True)
for ax,gene in zip(axes,['Rorc','Il23r','Ccr6']):
    for idx,cond in enumerate(['Ctrl','T2']):
        r=raw[(raw.gene==gene)&(raw.condition==cond)].set_index('tissue').loc[tissues]
        x=np.arange(3)+(idx-.5)*.34;ax.bar(x,r.pct_detected,width=.32,label=cond,color=['#728394','#C96A46'][idx])
        for xx,rr in zip(x,r.itertuples()):ax.text(xx,rr.pct_detected+.9,f'{rr.n_detected}/{rr.n_cells}',ha='center',va='bottom',fontsize=8,rotation=60)
    ax.set_xticks(range(3),['BM','Spleen','Thymus']);ax.set_title(gene);ax.set_ylim(0,100)
axes[0].set_ylabel('C5-2 cells with detected transcript (%)');axes[-1].legend(frameon=False)
fig.suptitle('Full-feature audit: C5-2 marker identity versus condition-specific detection',x=.08,ha='left',fontsize=12)
fig.subplots_adjust(bottom=.23,top=.82)
fig.text(.08,.04,'Labels: detected cells / all C5-2 cells. Thymus Ctrl has only 11 cells.\nIl17a: 0 / 15,532 retained cells; Il17f: 2 / 15,532 (both Ctrl spleen, outside C5-2).\nDetection is not cytokine secretion; these percentages have no independent-animal inference.',fontsize=10)
save(fig,'05_full_feature_iNKT17')
print('Saved five figures in PNG, SVG and PDF',flush=True)
