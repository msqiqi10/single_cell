"""Step 4: figures (PNG+SVG). Colour = GNN+semantic group (seed 42)."""
import sys,textwrap; sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from common import *
import numpy as np, pandas as pd, networkx as nx, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':10,'svg.fonttype':'none','font.family':'DejaVu Sans'})
FIG=RES/'figures'; T=RES/'tables'
MARK={'GO BP':'o','GO MF':'s','KEGG':'D','GO CC anchor':'^'}
KEY=[('cluster_tissue__C0__bone_marrow','BM C0'),('cluster_tissue__C3__spleen','Spleen C3'),('cluster_tissue__C4__bone_marrow','BM C4'),
     ('cluster_tissue__C5-2__bone_marrow','BM C5-2'),('cluster_tissue__C5-2__spleen','Spleen C5-2'),('cluster_tissue__C6__thymus','Thymus C6')]
def save(fig,name):
    for ext in ['png','svg']: fig.savefig(FIG/f'{name}.{ext}',dpi=200 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)
def wrap(s,w=28): return '\n'.join(textwrap.wrap(s,w))
for version,vtag in [('excl','b_primary_excl_ribo_mito_hsp'),('all','a_all_genes')]:
    n=pd.read_csv(RES/'network'/version/'nodes.csv',keep_default_na=False,na_values=['']); n['best_q']=pd.to_numeric(n.best_q)
    E=pd.read_csv(RES/'network'/version/'edges.tsv',sep='\t'); idx={t:i for i,t in enumerate(n.term_id)}
    G=nx.Graph(); G.add_nodes_from(range(len(n)))
    for r in E.itertuples(): G.add_edge(idx[r.GS_A_ID],idx[r.GS_B_ID],weight=r.JACCARD)
    groups=sorted(n.gnn_semantic.unique()); cmap=plt.get_cmap('tab20'); col={g:cmap(i%20) for i,g in enumerate(groups)}
    gs=pd.read_csv(T/f'group_summary_{version}.csv'); gs=gs[gs.method=='gnn_semantic'].set_index('group')
    conn=[i for i in G if G.degree(i)>0]; iso=[i for i in G if G.degree(i)==0]
    comps=sorted([list(c) for c in nx.connected_components(G.subgraph(conn))],key=len,reverse=True)
    pos={}; xcur=0; ycur=0; rowh=0; W=14
    for c in comps:
        H=G.subgraph(c); sc=max(1.2,0.55*np.sqrt(len(c))*2.2)
        p=nx.spring_layout(H,weight='weight',seed=42,k=2.5/np.sqrt(len(c)),iterations=800) if len(c)>1 else {c[0]:np.zeros(2)}
        P=np.array(list(p.values())); P=(P-P.mean(0))/(np.abs(P-P.mean(0)).max()+1e-9)*sc/2
        if xcur+sc>W: xcur=0; ycur+=rowh+1.2; rowh=0
        for (i,_),q in zip(p.items(),P): pos[i]=np.array([xcur+sc/2+q[0],ycur+sc/2+q[1]])
        xcur+=sc+1.2; rowh=max(rowh,sc)
    top=ycur+rowh
    iso_y0=-2.0
    ycur=0
    cols=int(np.ceil(len(iso)/2)) if iso else 1
    for j,i in enumerate(iso): pos[i]=np.array([W*(j%cols)/max(cols-1,1),-2.0-1.4*(j//cols)])
    fig,ax=plt.subplots(figsize=(17,14))
    for u,v,d in G.edges(data=True): ax.plot(*zip(pos[u],pos[v]),color='0.6',lw=0.4+2.5*d['weight'],alpha=.5,zorder=1)
    size=lambda q:40+28*q if q==q else 60
    for src,mk in MARK.items():
        d=n[n.source==src]
        for solid in [True,False]:
            dd=d[(~d.is_anchor)|(d.significant==solid)] if solid else d[d.is_anchor&(~d.significant)]
            if solid: dd=d[(d.significant)]
            if not len(dd): continue
            xy=np.array([pos[i] for i in dd.index]); 
            if solid: ax.scatter(xy[:,0],xy[:,1],s=[size(q) for q in -np.log10(dd.best_q.clip(lower=1e-300))],c=[col[g] for g in dd.gnn_semantic],marker=mk,edgecolors=['k' if a else 'w' for a in dd.is_anchor],linewidths=[1.8 if a else .6 for a in dd.is_anchor],zorder=3)
            else: ax.scatter(xy[:,0],xy[:,1],s=110,facecolors='none',c=None,marker=mk,edgecolors=[col[g] for g in dd.gnn_semantic],linewidths=2.0,linestyle='--',zorder=3)
    texts=[]
    for g in groups:
        r=n[n.gnn_semantic==g].sort_values('best_q',na_position='last').iloc[0]; 
        if r.is_anchor and r.best_q!=r.best_q: r=n[n.gnn_semantic==g].iloc[0]
        texts.append((r.name,f'[{g}] '+r.term,'k',9,False))
    for i in n.index[n.is_anchor]: texts.append((i,'ANCHOR: '+n.term[i],'darkred',8.5,True))
    seen={}
    for i,t,c,fs,anc in texts:
        x,y=pos[i]; k=seen.get(i,0); seen[i]=k+1
        ax.annotate(wrap(t,30),(x,y),xytext=(10+k*6,10+k*14+(hash(t)%4)*10),textcoords='offset points',fontsize=fs,color=c,fontweight='bold' if not anc else 'normal',
                    bbox=dict(boxstyle='round,pad=0.2',fc='white',ec=c,lw=.6,alpha=.85),arrowprops=dict(arrowstyle='-',lw=.5,color=c),zorder=6)
    ax.text(0,-1.2,'Isolated terms (no edge with Jaccard>=0.25 and >=3 shared genes)',fontsize=9,style='italic')
    ax.set_axis_off()
    h=[Line2D([],[],marker=m,ls='',mfc='0.8',mec='k',ms=9,label=s) for s,m in MARK.items()]+[Line2D([],[],marker='o',ls='',mfc='none',mec='0.3',ms=10,label='non-significant anchor (hollow, dashed)'),
       Line2D([],[],marker='o',ls='',mfc='0.8',mec='k',mew=2,ms=9,label='anchor (black rim)')]
    h+=[Line2D([],[],marker='o',ls='',mfc=col[g],mec='w',ms=9,label=f'{g}: {wrap(gs.loc[g,"representative_term"],26).splitlines()[0]}... (n={gs.loc[g,"n_terms"]})') for g in groups]
    ax.legend(handles=h,loc='center left',bbox_to_anchor=(1.0,.5),fontsize=9,frameon=False,title='Shape = source; colour = GNN+semantic group;\nsize ∝ best -log10(q)')
    ver='b PRIMARY (excl. Rpl/Rps/Mrpl/Mrps/mt-/Hsp/Dnaj)' if version=='excl' else 'a (all genes)'
    ax.set_title(f'iNKT term network, version {ver}: {len(n)} nodes, {len(E)} edges (exploratory; enrichment != activation)',fontsize=12)
    save(fig,f'network_full_{vtag}')
    # ---- super-PAG panels
    ng=len(groups); nc=4; nr=int(np.ceil(ng/nc)); fig,axs=plt.subplots(nr,nc,figsize=(5.2*nc,4.6*nr),squeeze=False)
    for a in axs.ravel(): a.set_axis_off()
    for a,g in zip(axs.ravel(),groups):
        m=list(n.index[n.gnn_semantic==g]); H=G.subgraph(m); p=nx.spring_layout(H,weight='weight',seed=1,k=None) if len(m)>1 else {m[0]:np.zeros(2)}
        for u,v,d in H.edges(data=True): a.plot(*zip(p[u],p[v]),color='0.6',lw=.5+3*d['weight'],alpha=.6,zorder=1)
        for i in m:
            r=n.loc[i]; sig=r.significant
            a.scatter(*p[i],s=size(-np.log10(r.best_q)) if sig else 110,marker=MARK[r.source],c=[col[g]] if sig else 'none',edgecolors='k' if r.is_anchor else ('w' if sig else col[g]),linewidths=1.6 if r.is_anchor else (.6 if sig else 2),zorder=3)
        lab=n.loc[m].sort_values('best_q',na_position='last'); show=list(lab.index[:3])+[i for i in m if n.is_anchor[i] and i not in lab.index[:3]]
        for j,i in enumerate(show):
            a.annotate(wrap(n.term[i],24),p[i],xytext=(5,5+7*(j%3)),textcoords='offset points',fontsize=7.5,color='darkred' if n.is_anchor[i] else 'k',zorder=6,bbox=dict(boxstyle='round,pad=.1',fc='white',ec='none',alpha=.75))
        rep=gs.loc[g,'representative_term']; a.set_title(f'{g}: {wrap(rep,40)}\n(n={len(m)})',fontsize=10,color='k',loc='left')
        a.margins(.3)
    fig.suptitle(f'Super-PAG panels, GNN+semantic groups, version {ver}. Hollow = non-significant anchor; dark-red label = anchor',fontsize=13,y=1.0)
    fig.tight_layout(); save(fig,f'superPAG_panels_{vtag}')
    # ---- heatmap
    dg=pd.read_csv(T/f'group_driver_genes_{version}.csv'); dg=dg[dg.method=='gnn_semantic']
    fig,axs=plt.subplots(1,2,figsize=(15,0.5*len(groups)+2.6),sharey=True)
    for ax,dr in zip(axs,['T2_up','Ctrl_up']):
        M=np.zeros((len(groups),len(KEY)),int)
        for gi,g in enumerate(groups):
            for ki,(u,_) in enumerate(KEY):
                s=dg[(dg.group==g)&(dg.unit_id==u)&(dg.direction==dr)]; M[gi,ki]=int(s.n_driver_genes.iloc[0]) if len(s) else 0
        im=ax.imshow(M,cmap='Blues',aspect='auto',vmin=0,vmax=max(M.max(),1))
        ax.set_xticks(range(len(KEY)),[k[1] for k in KEY],rotation=35,ha='right'); ax.set_title(f'{dr}: unique driver genes in group terms (q<=0.05, >=3 hits)',fontsize=10)
        for (i,j),v in np.ndenumerate(M): ax.text(j,i,v if v else '',ha='center',va='center',fontsize=9,color='w' if v>M.max()*.6 else 'k')
    axs[0].set_yticks(range(len(groups)),[f'{g}: {gs.loc[g,"representative_term"][:38]}' for g in groups],fontsize=9)
    fig.suptitle(f'Driver genes per GNN+semantic group x key comparison, version {ver}',fontsize=12); fig.tight_layout(); save(fig,f'heatmap_group_driver_genes_{vtag}')
print('done')
