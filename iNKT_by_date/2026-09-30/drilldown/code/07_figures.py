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
def layout(H,members):
    iso=[g for g in members if g not in H or H.degree(g)==0]; conn=[g for g in members if g in H and H.degree(g)>0]
    pos={}
    if conn:
        S=H.subgraph(conn); n=len(conn)
        p=nx.spring_layout(S,weight='weight',seed=7,k=2.2/np.sqrt(n),iterations=400)
        P=np.array(list(p.values())); P=(P-P.mean(0)); P/=np.abs(P).max()
        for g,q in zip(p,P): pos[g]=q*10
    top=-11.5; nc=max(8,int(np.ceil(np.sqrt(len(iso)*3))))
    for j,g in enumerate(sorted(iso)): pos[g]=np.array([-10+20*(j%nc)/max(nc-1,1),top-1.1*(j//nc)])
    return pos,conn,iso
cmap=plt.get_cmap('RdBu_r'); VM=1.5
for pid,pname in DRILL:
    m=mem[mem.pathway_id==pid].set_index('gene'); members=list(m.index)
    H=G_all.subgraph(members).copy(); big=len(members)>120
    if big: H.remove_edges_from([(u,v) for u,v,d in list(H.edges(data=True)) if d['weight']<4])
    pos,conn,iso=layout(H,members)
    fig,axs=plt.subplots(1,2,figsize=(19,10.5 if not iso else 10.5+0.35*np.ceil(len(iso)/12)))
    for ax,(cs,cl) in zip(axs,CS):
        deg={g for g in members if m.at[g,f'{cs} DEG']=='yes'}
        near=set(deg)
        for g in deg:
            if g in H: near|=set(H.neighbors(g))
        for u,v,d in H.edges(data=True):
            hi=(u in near and v in near) or not big
            ax.plot(*zip(pos[u],pos[v]),color='0.55',lw=.3+.25*d['weight'],alpha=.35 if hi else .06,zorder=1)
        for g in members:
            x,y=pos[g]; det=m.at[g,'detected_in_data']=='yes'; fade=big and g not in near
            if not det: ax.scatter(x,y,s=55,facecolors='none',edgecolors='0.55',linewidths=1.2,zorder=3,alpha=.5 if fade else 1); continue
            v=m.at[g,f'{cs} log2FC']; d=g in deg
            ax.scatter(x,y,s=95 if d else 55,c=[cmap((np.clip(v,-VM,VM)+VM)/(2*VM))],edgecolors='red' if d else '0.3',linewidths=2.2 if d else .5,alpha=.2 if fade else 1,zorder=4 if d else 3)
        lab=[g for g in members if g in deg]
        for g in lab:
            ax.annotate('★'+g,pos[g],xytext=(5,5),textcoords='offset points',fontsize=8.5,color='darkred',zorder=6,bbox=dict(boxstyle='round,pad=.1',fc='white',ec='none',alpha=.7))
        if not big:
            for g in members:
                if g not in deg and m.at[g,'detected_in_data']=='yes': ax.annotate(g,pos[g],xytext=(4,-9),textcoords='offset points',fontsize=6.5,color='0.35',zorder=5)
        ax.set_title(f'{cl}: {len(deg)} DEG of {(m.detected_in_data=="yes").sum()} detected',fontsize=12); ax.set_axis_off(); ax.margins(.06)
        if iso: ax.text(-10,-11.0,f'{len(iso)} members with no co-membership edge (grid)',fontsize=8,style='italic')
    sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(-VM,VM)); cb=fig.colorbar(sm,ax=axs,shrink=.35,pad=.01,aspect=25); cb.set_label('log2FC (T2 vs Ctrl), clipped +/-1.5')
    h=[Line2D([],[],marker='o',ls='',mfc='0.85',mec='red',mew=2.2,ms=9,label='DEG (FDR<=0.05, |log2FC|>=0.25) + star'),Line2D([],[],marker='o',ls='',mfc='none',mec='0.55',ms=9,label='not detected in data')]
    axs[0].legend(handles=h,loc='lower left',fontsize=8.5,frameon=False)
    nm=f'{pname} ({len(members)} members, {(m.detected_in_data=="yes").sum()} detected)'+('; non-DEG genes outside DEG neighbourhoods faded' if big else '')
    if pid=='ENERGY': nm=nm.replace(pname,'Energy module (G01 core: gene in >=3 of 15 terms)')
    fig.suptitle(nm,fontsize=14,y=.99); fig.text(.01,.005,FOOT,fontsize=7.5,color='0.3')
    save(fig,'gene_network_'+slug(pname))
print('per-pathway figs done')
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
a2.imshow(np.clip(L,-VM,VM),cmap='RdBu_r',aspect='auto',vmin=-VM,vmax=VM); a2.set_xticks(range(len(cols))); a2.set_xticklabels(cols,rotation=60,ha='right')
for i,j in zip(*np.where(D)): a2.text(j,i,'*',ha='center',va='center',fontsize=11)
a2.set_title('log2FC (* = DEG)')
fig.suptitle('Genes shared by >=2 anchor pathways (DEG in >=1 key comparison) + AP-1 members',fontsize=11,y=1.0); fig.text(.01,-.01,FOOT[:170]+'...',fontsize=6.5)
save(fig,'bridge_genes_heatmap'); print('heatmap done',len(W))
# ---------- v2 network (version b)
n=pd.read_csv(RES/'network/excl/nodes.csv',keep_default_na=False,na_values=['']); n['best_q']=pd.to_numeric(n.best_q)
Ed=pd.read_csv(RES/'network/excl/edges.tsv',sep='\t'); idx={t:i for i,t in enumerate(n.term_id)}
G=nx.Graph(); G.add_nodes_from(range(len(n)))
for r in Ed.itertuples(): G.add_edge(idx[r.GS_A_ID],idx[r.GS_B_ID],weight=r.JACCARD)
groups=sorted(n.gnn_semantic.unique()); cm=plt.get_cmap('tab20'); col={g:cm(i%20) for i,g in enumerate(groups)}
gs=pd.read_csv(RES/'tables/group_summary_excl.csv'); gs=gs[gs.method=='gnn_semantic'].set_index('group')
MARK={'GO BP':'o','GO MF':'s','KEGG':'D','GO CC anchor':'^'}
conn=[i for i in G if G.degree(i)>0]; iso=[i for i in G if G.degree(i)==0]
comps=sorted([list(c) for c in nx.connected_components(G.subgraph(conn))],key=len,reverse=True)
pos={}; xcur=0; ycur=0; rowh=0; Wd=14
for c in comps:
    H=G.subgraph(c); sc=max(1.2,0.55*np.sqrt(len(c))*2.2)
    p=nx.spring_layout(H,weight='weight',seed=42,k=2.5/np.sqrt(len(c)),iterations=800) if len(c)>1 else {c[0]:np.zeros(2)}
    P=np.array(list(p.values())); P=(P-P.mean(0))/(np.abs(P-P.mean(0)).max()+1e-9)*sc/2
    if xcur+sc>Wd: xcur=0; ycur+=rowh+1.2; rowh=0
    for (i,_),q in zip(p.items(),P): pos[i]=np.array([xcur+sc/2+q[0],ycur+sc/2+q[1]])
    xcur+=sc+1.2; rowh=max(rowh,sc)
cols_=int(np.ceil(len(iso)/2)) if iso else 1
for j,i in enumerate(iso): pos[i]=np.array([Wd*(j%cols_)/max(cols_-1,1),-2.0-1.4*(j//cols_)])
fig,ax=plt.subplots(figsize=(18,11))
for u,v,d in G.edges(data=True): ax.plot(*zip(pos[u],pos[v]),color='0.6',lw=0.4+2.5*d['weight'],alpha=.5,zorder=1)
size=lambda q:40+28*q if q==q else 60
for src,mk in MARK.items():
    d=n[n.source==src]
    s_=d[d.significant.astype(str)=='True']
    if len(s_):
        xy=np.array([pos[i] for i in s_.index])
        ax.scatter(xy[:,0],xy[:,1],s=[size(q) for q in -np.log10(s_.best_q.clip(lower=1e-300))],c=[col[g] for g in s_.gnn_semantic],marker=mk,edgecolors=['k' if a else 'w' for a in s_.is_anchor.astype(str)=='True'],linewidths=[1.8 if a else .6 for a in s_.is_anchor.astype(str)=='True'],zorder=3)
    h_=d[(d.is_anchor.astype(str)=='True')&(d.significant.astype(str)!='True')]
    if len(h_):
        xy=np.array([pos[i] for i in h_.index]); ax.scatter(xy[:,0],xy[:,1],s=110,facecolors='none',marker=mk,edgecolors=[col[g] for g in h_.gnn_semantic],linewidths=2.0,linestyle='--',zorder=3)
# labels: group representative + anchors, greedy placement avoiding overlap, title zone and axes edges
wrap=lambda s,w=26:'\n'.join(textwrap.wrap(s,w))
items=[]
for g in groups:
    r=n[n.gnn_semantic==g].sort_values('best_q',na_position='last').iloc[0]; items.append((r.name,f'[{g}] '+r.term,'k',9,True))
for i in n.index[n.is_anchor.astype(str)=='True']: items.append((i,n.term[i],'darkred',8,False))
xs=np.array(list(pos.values())); pad=1.5
ax.set_xlim(xs[:,0].min()-pad,xs[:,0].max()+pad); ax.set_ylim(xs[:,1].min()-pad,xs[:,1].max()+pad+0.5); ax.set_axis_off()
ax.text(0,-1.0,'Isolated terms (no edge with Jaccard>=0.25 and >=3 shared genes)',fontsize=9,style='italic')
fig.canvas.draw(); rend=fig.canvas.get_renderer(); axbb=ax.get_window_extent(rend)
placed=[]
# node bboxes as obstacles (pixels)
pts=ax.transData.transform(np.array([pos[i] for i in n.index]))
node_boxes=[(x-9,y-9,x+9,y+9) for x,y in pts]
def ov(a,b): return max(0,min(a[2],b[2])-max(a[0],b[0]))*max(0,min(a[3],b[3])-max(a[1],b[1]))
cands=[(dx,dy) for r_ in (14,34,58,84,112) for dx,dy in [(r_,r_*.5),(-r_,r_*.5),(r_,-r_*.5),(-r_,-r_*.5),(0,r_),(0,-r_),(r_*1.3,0),(-r_*1.3,0)]]
for i,t,c,fs,bold in sorted(items,key=lambda z:-len(z[1])):
    x,y=ax.transData.transform(pos[i]); tt=ax.text(0,0,wrap(t),fontsize=fs,color=c,fontweight='bold' if bold else 'normal',ha='left',va='bottom',bbox=dict(boxstyle='round,pad=0.2',fc='white',ec=c,lw=.6,alpha=.9),zorder=6)
    bb=tt.get_window_extent(rend); w,h=bb.width,bb.height; best=None
    for dx,dy in cands:
        x0=x+dx-(w if dx<0 else 0); y0=y+dy-(h if dy<0 else 0); box=(x0,y0,x0+w,y0+h)
        if box[0]<axbb.x0+5 or box[2]>axbb.x1-5 or box[1]<axbb.y0+5 or box[3]>axbb.y1-5: continue
        cost=sum(ov(box,b) for b in placed)*5+sum(ov(box,nb) for j,nb in enumerate(node_boxes) if j!=i)*.6+(abs(dx)+abs(dy))*.02
        if best is None or cost<best[0]: best=(cost,box)
    if best is None: best=(0,(x+10,y+10,x+10+w,y+10+h))
    box=best[1]; placed.append(box); tt.remove()
    inv=ax.transData.inverted(); (xa,ya),(xb,yb)=inv.transform((box[0],box[1])),inv.transform((box[2],box[3]))
    ax.text(xa,ya,wrap(t),fontsize=fs,color=c,fontweight='bold' if bold else 'normal',ha='left',va='bottom',bbox=dict(boxstyle='round,pad=0.2',fc='white',ec=c,lw=.6,alpha=.9),zorder=6)
    cx,cy=(box[0]+box[2])/2,(box[1]+box[3])/2; px,py=x,y
    # leader line from node to nearest box point
    qx=min(max(px,box[0]),box[2]); qy=min(max(py,box[1]),box[3]); (lx,ly)=inv.transform((qx,qy)); ax.plot([pos[i][0],lx],[pos[i][1],ly],color=c,lw=.5,zorder=5)
h=[Line2D([],[],marker=m,ls='',mfc='0.8',mec='k',ms=9,label=s) for s,m in MARK.items()]+[Line2D([],[],marker='o',ls='',mfc='none',mec='0.3',ms=10,label='non-significant anchor (hollow)'),Line2D([],[],marker='o',ls='',mfc='0.8',mec='k',mew=2,ms=9,label='anchor (black rim; dark-red label)')]
h+=[Line2D([],[],marker='o',ls='',mfc=col[g],mec='w',ms=9,label=f'{g}: {wrap(gs.loc[g,"representative_term"],30).splitlines()[0]}... (n={gs.loc[g,"n_terms"]})') for g in groups]
ax.legend(handles=h,loc='center left',bbox_to_anchor=(1.0,.5),fontsize=9,frameon=False,title='Shape = source; colour = GNN+semantic group;\nsize ~ best -log10(q)')
fig.suptitle(f'iNKT term network, version b PRIMARY (excl. Rpl/Rps/Mrpl/Mrps/mt-/Hsp/Dnaj): {len(n)} nodes, {len(Ed)} edges',fontsize=13,y=.955)
fig.text(.5,.925,'Exploratory; enrichment != activation. G03 is a catch-all (mostly isolates / 2-node fragments), not evidence that anchors are related.',ha='center',fontsize=9,color='0.3')
save(fig,'network_full_b_primary_excl_ribo_mito_hsp_v2',d=FIGD)
print('v2 done')
