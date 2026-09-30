"""Render paper-style figures, PDF and matching PPTX from the dated numeric outputs."""
from pathlib import Path
import json, textwrap, hashlib
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from scipy.signal import savgol_filter
from pptx import Presentation
from pptx.util import Inches
from PIL import Image, ImageOps, ImageDraw
from analyze import OUT,ROOT,RES,SRC,PREV,KEY,ORDER,TISSUES,MODULES,sha
PRES=OUT/'presentation'; FIG=RES/'figures'; PAGES=[]
NAMES={'bone_marrow':'Bone marrow','spleen':'Spleen','thymus':'Thymus'}
SHORT={'bone_marrow':'BM','spleen':'Spl','thymus':'Thy'}
COLORS={'bone_marrow':{'Ctrl':'#FFFF00','T2':'#FFA500'},'spleen':{'Ctrl':'#FF0000','T2':'#8B0000'},'thymus':{'Ctrl':'#0000FF','T2':'#00008B'}}
C_COLORS=dict(zip(ORDER,['#31a8d3','#ee9920','#087a63','#9467bd','#e5b53c','#4485ab','#d56b73','#67a547','#787878']))
WARN='Exploratory transcriptional results; one sample label per tissue x condition; no animal-level inference.'
PAGER='Huang et al. 2024, PAGER-scFGA, Fig. 5; doi:10.3389/fbinf.2024.1336135'
BORRA='Borra et al. 2026, GAFA, Sec. 3.5 / Fig. 4c; doi:10.3390/biology15070588'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def read(name):return pd.read_csv(RES/'tables'/name)
def new(title,subtitle,reference,note=WARN):
 f=plt.figure(figsize=(16,9),facecolor='white');f.text(.035,.955,title,weight='bold',fontsize=23,va='top');f.text(.035,.900,subtitle,fontsize=12,va='top',color='#444')
 f.text(.035,.065,note,fontsize=9.5,va='top');f.text(.035,.024,reference,fontsize=8.5,color='#555');f.text(.965,.024,str(len(PAGES)+1),ha='right',fontsize=9)
 f._info={'title':title,'subtitle':subtitle,'reference':reference,'note':note};return f

def finish(f,pdf,notes=''):
 i=len(PAGES)+1;stem=f'{i:02}';f.savefig(FIG/(stem+'.png'),dpi=130);f.savefig(FIG/(stem+'.svg'));pdf.savefig(f);PAGES.append({'page':i,**f._info,'notes_zh':notes,'figure':str((FIG/(stem+'.png')).relative_to(OUT))});plt.close(f);print(stem,f._info['title'],flush=True)
def wrapped(s,n=42):return '\n'.join(textwrap.wrap(s,n))
def no_data(ax,text):ax.axis('off');ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes)
def heat(ax,d,title,cmap='RdBu_r',center=True):
 vals=d.to_numpy(float);lim=max(.001,np.nanmax(abs(vals))) if np.isfinite(vals).any() else 1
 im=ax.imshow(np.ma.masked_invalid(vals),aspect='auto',cmap=cmap,vmin=-lim if center else 0,vmax=lim)
 ax.set_xticks(range(len(d.columns)),d.columns,rotation=35,ha='right',fontsize=10);ax.set_yticks(range(len(d)),d.index,fontsize=10);ax.set_title(title,loc='left',weight='bold');plt.colorbar(im,ax=ax,fraction=.035,pad=.025);return im

def choose_terms(d,max_terms=16):
 sig=d[(d.q_family<=.05)&(d.k>=2)].copy();chosen=[];sets=[]
 # Round-robin best terms for every cluster prevents a large cluster taking every row.
 candidates=[]
 for _,g in sig.groupby('cluster',sort=True):candidates.append(g.sort_values(['q_family','enrichment_ratio'],ascending=[True,False]).drop_duplicates('go_id').head(8))
 if not candidates:return []
 joined=pd.concat(candidates).sort_values('q_family').drop_duplicates('go_id')
 for r in joined.itertuples():
  gs=set(str(r.genes).split(';'))
  if any(len(gs&s)/max(1,len(gs|s))>=.75 for s in sets):continue
  chosen.append(r.go_id);sets.append(gs)
  if len(chosen)>=max_terms:break
 return chosen

def go_bubbles(ax,d,terms,condition=False):
 if not terms:no_data(ax,'No GO terms passed BH q <= 0.05');return
 clusters=[c for c in ORDER if c in set(d.cluster)];cols=[c+'/'+dr for c in clusters for dr in ['T2','Ctrl']] if condition else clusters
 lookup=d.drop_duplicates('go_id').set_index('go_id').term.to_dict();ax.set_facecolor('#f7f7f7')
 ax.set_xlim(-.6,len(cols)-.4);ax.set_ylim(len(terms)-.5,-.5)
 ax.set_xticks(range(len(cols)),cols,rotation=55 if condition else 0,ha='right' if condition else 'center',fontsize=10)
 ax.set_yticks(range(len(terms)),[wrapped(lookup[t],43) for t in terms],fontsize=9)
 sub=d[d.go_id.isin(terms)&(d.q_family<=.05)]
 for r in sub.itertuples():
  col=r.cluster+'/'+('T2' if r.direction=='T2_up' else 'Ctrl') if condition else r.cluster
  if col not in cols:continue
  val=min(10,-np.log10(max(r.q_family,1e-300)));ax.scatter(cols.index(col),terms.index(r.go_id),s=25+260*r.gene_ratio,c=[val],vmin=0,vmax=10,cmap='viridis',edgecolors='#333',linewidths=.3)
 plt.colorbar(ScalarMappable(norm=Normalize(0,10),cmap='viridis'),ax=ax,fraction=.024,pad=.03,label='-log10(BH q), capped at 10')
 for ratio in [.05,.2,.5]:ax.scatter([],[],s=25+260*ratio,color='gray',label=f'{ratio:.0%}')
 ax.legend(title='Gene ratio',bbox_to_anchor=(1.14,1),loc='upper left',frameon=False,fontsize=9)
 ax.grid(axis='x',alpha=.13);ax.set_axisbelow(True)

def main():
 a=ad.read_h5ad(PREV/'objects/scored_base.h5ad');obs=read('cell_metadata_scores.csv.gz').set_index('cell');obs=obs.loc[a.obs_names];xy=a.obsm['X_umap']
 contrasts=read('cytotoxicity_contrasts.csv');sens=read('cytotoxicity_QC_sensitivity.csv');coverage=read('module_gene_coverage.csv');dot=read('marker_dotplot_values.csv');go=read('GO_ORA_significant.csv');go_all=read('GO_contrast_coverage.csv');root=read('trajectory_root_audit.csv');curves=read('trajectory_binned_curves.csv');tc=read('trajectory_correlations.csv')
 selected=[]
 with PdfPages(PRES/'iNKT_cytotoxicity_GO_20260919.pdf') as pdf:
  f=new('iNKT cytotoxicity and Gene Ontology','Paper-guided analysis on the latest local data | 19 September 2026',PAGER+'; '+BORRA)
  ax=f.add_axes([.06,.20,.48,.61]);
  for c in ORDER:
   m=obs[KEY].astype(str).eq(c).to_numpy();ax.scatter(xy[m,0],xy[m,1],s=3,c=C_COLORS[c],label=c,rasterized=True,linewidths=0)
  ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(ncol=3,markerscale=3,fontsize=10,frameon=False)
  f.text(.60,.74,'15,532 cells | 10,670 measured genes',fontsize=19,weight='bold');f.text(.60,.63,'3 tissues | Ctrl and T2 | 9 refined clusters',fontsize=15)
  f.text(.60,.50,'PAGER: gene modules and temporal profiles\n\nGAFA: cluster-specific GO enrichment\n\nPrimary contrasts: tissue x original cluster\nSensitivity: stable tissue reclustering + QC matching',fontsize=14,linespacing=1.6)
  finish(f,pdf,'本轮新增的是细胞毒性与GO计算。输入来自9月5日最新分析对象；9月15日仅是展示更新。')
  f=new('Analysis design and reference mapping','Methods adapted to mouse iNKT; biological conclusions are derived from these data.',PAGER+'; '+BORRA)
  rows=[('Cytotoxicity','PAGER Fig. 5C/D','Core: Prf1, Gzma, Gzmb; E1: +Ctla2a','Mean log1p expression; matched-control score sensitivity'),('GO cluster identity','GAFA Fig. 4c','Cluster vs rest within each tissue','Wilcoxon, gene BH <= .05 and log2FC >= .25'),('GO treatment response','GAFA functional comparison','T2-up and Ctrl-up genes tested separately','Hypergeometric; full eligible GO BP/MF/CC families'),('Local trajectories','PAGER Fig. 5D','Recompute diffusion map on tissue-specific graphs','Precursor-rooted DPT; 20 shared bins; 3-root sensitivity'),('Robustness','Adaptation to this dataset','Stable clusters; five QC-matched resamples','Effects and direction consistency; no animal replication')]
  ax=f.add_axes([.04,.24,.92,.55]);ax.axis('off');table=ax.table(cellText=rows,colLabels=['Analysis','Reference','Input / definition','Implementation'],colWidths=[.15,.15,.33,.37],cellLoc='left',colLoc='left',bbox=[0,0,1,1]);table.auto_set_font_size(False);table.set_fontsize(10)
  for (r,c),cell in table.get_celld().items():cell.set_edgecolor('#ddd');cell.set_facecolor('#eef1f4' if r==0 else 'white');cell.get_text().set_wrap(True)
  f.text(.05,.15,'GO uses official mouse annotations. An enrichment of upregulated genes does not prove pathway activation.',fontsize=12)
  finish(f,pdf,'采用论文中的方法和图形组织方式。没有运行SCORPION、随机森林或声称重现PAGER专有网络；GO采用官方MGI小鼠注释。')
  f=new('Cytotoxic gene coverage and scoring','Primary endpoint: mean log-normalized expression of Prf1 / Gzma / Gzmb.',PAGER)
  c=coverage.drop_duplicates('gene');ax=f.add_axes([.08,.23,.40,.56]);ax.bar(c.gene,c.n_detected/len(obs)*100,color='#176b8d');ax.set_ylabel('Cells with detected transcript (%)');ax.set_title('A  Measured marker coverage',loc='left')
  ax=f.add_axes([.61,.23,.31,.56]);values=obs.groupby(['tissue','condition'],observed=True)[['PAGER_core_mean','PAGER_E1_mean','prior_effector_mean']].mean();values.index=[SHORT[x]+' / '+y for x,y in values.index];values.columns=['Core 3','E1 4','Prior 4'];heat(ax,values,'B  Sample-level mean expression',cmap='viridis',center=False)
  finish(f,pdf,'分数是转录表达而不是杀伤实验读数。不同基因集的绝对分数不可视为相同量纲的功能强弱。')
  vmax=np.quantile(obs.PAGER_core_mean,.99)
  for tissue in TISSUES:
   f=new(NAMES[tissue]+': cytotoxic expression and condition effects','Original refined clusters; T2 - Ctrl within the same tissue and cluster.',PAGER)
   m=obs.tissue.eq(tissue).to_numpy();ax=f.add_axes([.04,.48,.28,.32]);im=ax.scatter(xy[m,0],xy[m,1],s=4,c=obs.loc[m,'PAGER_core_mean'],cmap='magma',vmin=0,vmax=vmax,rasterized=True,linewidths=0);ax.set_xticks([]);ax.set_yticks([]);ax.set_title('A  Cytotoxic core expression',loc='left');plt.colorbar(im,ax=ax,fraction=.035)
   dd=dot[(dot.tissue==tissue)&dot.gene.isin(['Prf1','Gzma','Gzmb','Ctla2a','Nkg7'])];clusters=[c for c in ORDER if c in set(dd.cluster)];labels=[c+'/'+cond for c in clusters for cond in ['Ctrl','T2']];genes=['Prf1','Gzma','Gzmb','Ctla2a','Nkg7'];ax=f.add_axes([.41,.22,.22,.58]);vmax_gene=max(dd.mean_log1p.max(),.01)
   for r in dd.itertuples():ax.scatter(genes.index(r.gene),labels.index(r.cluster+'/'+r.condition),s=3+110*r.pct_positive,c=[r.mean_log1p],cmap='Reds',vmin=0,vmax=vmax_gene,edgecolors='#777',linewidths=.2)
   ax.set_xticks(range(5),genes,rotation=40,ha='right');ax.set_yticks(range(len(labels)),labels,fontsize=8);ax.invert_yaxis();ax.set_title('B  Gene expression',loc='left');plt.colorbar(ScalarMappable(norm=Normalize(0,vmax_gene),cmap='Reds'),ax=ax,fraction=.05,pad=.04,label='Mean log1p')
   f.text(.405,.12,'Dot area: fraction of cells expressing gene',fontsize=9)
   d=contrasts[(contrasts.tissue==tissue)&(contrasts.definition=='original')&(contrasts.module=='PAGER_core')].set_index('cluster').reindex(clusters)
   ax=f.add_axes([.74,.25,.21,.52]);colors=['#b2182b' if v>0 else '#2166ac' for v in d.delta_mean.fillna(0)];ax.barh(clusters,d.delta_mean.where(d.status=='tested',0),color=colors);ax.axvline(0,c='#333',lw=.8);ax.invert_yaxis();ax.set_title('C  Core expression difference',loc='left');ax.set_xlabel('Mean(T2) - mean(Ctrl)')
   for j,(c,r) in enumerate(d.iterrows()):
    text='NA' if r.status!='tested' else ('*' if r.q_all<=.05 else '')
    ax.text(r.delta_mean if r.status=='tested' and np.isfinite(r.delta_mean) else 0,j,text,va='center',fontsize=12)
   f.text(.045,.29,'Core = mean log1p(Prf1, Gzma, Gzmb)\n* Global BH q <= 0.05 across tested\n  modules, tissues and cluster definitions\nNA: fewer than 20 cells in either condition',fontsize=10,linespacing=1.5)
   finish(f,pdf,'按组织和原分群展示细胞毒性；红色为T2高，蓝色为Ctrl高。星号是细胞层面探索性检验，不能代表动物重复验证。')
  for analysis,title in [('cluster_identity','GO biological processes: cluster identity'),('condition','GO biological processes: condition response')]:
   for tissue in TISSUES:
    d=go[(go.tissue==tissue)&(go.analysis==analysis)&(go.definition=='original')&(go.namespace=='BP')];terms=choose_terms(d)
    f=new(NAMES[tissue]+': '+title,'Within-tissue comparisons | dot area: query gene ratio | color: enrichment evidence',BORRA,note=WARN+' Blank = no displayed significant term, not proof of absence.')
    ax=f.add_axes([.35,.22,.50,.59]);go_bubbles(ax,d,terms,condition=analysis=='condition')
    for t in terms:selected.append({'page':len(PAGES)+1,'analysis':analysis,'tissue':tissue,'go_id':t})
    f.text(.04,.12,'All BP/MF/CC tests retained. Only clusters with significant BP terms shown; q <= .05; >=2 hits; overlap redundancy reduced.',fontsize=10)
    finish(f,pdf,'GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。')
  for tissue in TISSUES:
   audit=root[root.tissue==tissue].iloc[0];b=ad.read_h5ad(RES/'objects'/f'{tissue}_local_trajectory.h5ad');xyb=b.obsm['X_umap'];f=new(NAMES[tissue]+': local cytotoxicity along pseudotime','PAGER-style functional curves on a newly computed within-tissue diffusion map.',PAGER,note=WARN+' Pseudotime is an assumed ordering, not RNA velocity or lineage validation.')
   ax=f.add_axes([.05,.27,.34,.50]);im=ax.scatter(xyb[:,0],xyb[:,1],s=4,c=b.obs.local_dpt,cmap='viridis',rasterized=True,linewidths=0);ax.set_xticks([]);ax.set_yticks([]);ax.set_title('A  Local DPT on tissue UMAP',loc='left');plt.colorbar(im,ax=ax,fraction=.04)
   ax=f.add_axes([.51,.30,.43,.45]);
   for cond in ['Ctrl','T2']:
    d=curves[(curves.tissue==tissue)&(curves.condition==cond)&(curves.module=='PAGER_core')&(curves.n_cells>=20)].sort_values('bin')
    ax.plot(d.pseudotime,d.mean_expression,'o',mec='#333',mew=.5,color=COLORS[tissue][cond],label=cond+' bin means',ms=5)
    # Smooth only contiguous supported bins; never bridge gaps with low cell coverage.
    groups=(d.bin.diff().fillna(1)>1).cumsum()
    for _,segment in d.groupby(groups):
     window=min(5,len(segment) if len(segment)%2 else len(segment)-1)
     if window>=3:ax.plot(segment.pseudotime,savgol_filter(segment.mean_expression,window,min(2,window-1)),color=COLORS[tissue][cond],lw=2.5)
    rho=tc[(tc.tissue==tissue)&(tc.condition==cond)&(tc.module=='PAGER_core')].iloc[0].rho
    ax.text(.04,.93 if cond=='Ctrl' else .84,f'{cond}: Spearman rho = {rho:+.2f}',transform=ax.transAxes,fontsize=10)
   ax.set_xlabel('Local diffusion pseudotime');ax.set_ylabel('Core mean log1p expression');ax.set_title('B  Common 20-bin functional profiles',loc='left');ax.legend(loc='best',fontsize=9)
   f.text(.06,.15,f'Root state: {audit.root_cluster}; measured root markers: {audit.root_genes.replace(';', ' / ')} (Cd24a unavailable). Root sensitivity rho: {audit.root_rho_1:.2f}, {audit.root_rho_2:.2f}.\nAnalyzed {audit.n_analyzed:,}/{audit.n_cells:,} cells; only bins with >=20 cells shown. Smooth lines are descriptive.',fontsize=11)
   finish(f,pdf,'重新计算了组织内diffusion map，避免复用全局轨迹。用前体marker选择根群，三个根细胞做敏感性。平滑曲线不是置信区间，条件曲线差别不等于分化速率差别。')
  f=new('Core signature changes do not imply a full killing pathway','Two reproducible local core signals; the broader E1 signature is not significant.',PAGER+'; '+BORRA)
  focus=contrasts[(contrasts.definition=='original')&(((contrasts.tissue=='bone_marrow')&(contrasts.cluster=='C4'))|((contrasts.tissue=='spleen')&(contrasts.cluster=='C3')))].copy()
  focus['contrast']=focus.tissue.map(SHORT)+' / '+focus.cluster
  ax=f.add_axes([.10,.41,.35,.35]);mat=focus.pivot(index='contrast',columns='module',values='delta_mean').reindex(columns=['PAGER_core','PAGER_E1','prior_effector']);mat.columns=['Core 3','E1 4','Prior 4'];heat(ax,mat,'A  Signature sensitivity')
  for y,label in enumerate(mat.index):
   for x,module in enumerate(['PAGER_core','PAGER_E1','prior_effector']):
    r=focus[(focus.contrast==label)&(focus.module==module)].iloc[0];ax.text(x,y,f'{r.delta_mean:+.3f}'+(' *' if r.q_all<=.05 else ''),ha='center',va='center',fontsize=11,color='white' if abs(r.delta_mean)>.06 else 'black')
  kill=read('GO_cytotoxicity_targeted_audit.csv');cross=read('focus_stable_cluster_crosswalk.csv')
  text='B  Complementary evidence from the same dataset\n\nBone marrow C4 core: +0.086; q = 0.025\nSpleen C3 core: -0.023; q = 0.0015\nBoth: same direction in all 5 QC resamples\n\nStable-cluster effects: +0.111 and -0.022\nOriginal-cell overlap: 98.7% and 99.5%\n\nE1 (+Ctla2a): q = 0.93 and 0.53\nKilling-related GO terms: no q <= 0.05'
  f.text(.55,.77,text,fontsize=14,va='top',linespacing=1.35)
  f.text(.08,.20,'The evidence supports a localized core-effector transcriptional shift.\nIt does not establish increased/decreased killing activity or coordinated activation of a full cytotoxic pathway.',fontsize=13,linespacing=1.5)
  finish(f,pdf,'骨髓C4和脾脏C3的核心3基因模块及原4基因效应模块有一致局部变化；加入Ctla2a后的论文E1代表模块不显著。杀伤相关GO也未显著，不能推断整个细胞毒性通路增强或减弱。')
  f=new('Cellular components of cytotoxic module genes','GO-CC membership supports localization; this is not a protein-interaction network.',PAGER+'; official mouse GO annotations')
  cc=read('cytotoxic_gene_CC_membership.csv');top=cc.groupby(['go_id','term']).gene.nunique().sort_values(ascending=False).head(14).reset_index();genes=sorted(cc.gene.unique());ax=f.add_axes([.43,.21,.46,.58]);matrix=np.array([[int(((cc.go_id==t)&(cc.gene==g)).any()) for g in genes] for t in top.go_id]);im=ax.imshow(matrix,cmap='Blues',aspect='auto',vmin=0,vmax=1);ax.set_yticks(range(len(top)),[wrapped(t,42) for t in top.term],fontsize=10);ax.set_xticks(range(len(genes)),genes);ax.set_title('Measured genes x propagated GO-CC membership',loc='left')
  finish(f,pdf,'这里只画真实GO细胞组分注释矩阵，不把共同定位误画为PPI相互作用。')
  f=new('Robustness: QC matching and alternative clustering','Five condition-balanced QC resamples within each tissue x cluster.',PAGER+'; sensitivity analyses added for this dataset')
  d=contrasts[(contrasts.module=='PAGER_core')&(contrasts.definition=='original')&(contrasts.status=='tested')].copy();d['label']=d.tissue.map(SHORT)+'/'+d.cluster
  s=sens[(sens.module=='PAGER_core')&(sens.definition=='original')].groupby(['tissue','cluster']).agg(min_delta=('delta_mean','min'),max_delta=('delta_mean','max'),same=('same_direction','sum'),runs=('seed','size'),max_smd=('max_abs_QC_SMD','max')).reset_index();d=d.merge(s,on=['tissue','cluster'],how='left');ax=f.add_axes([.09,.29,.84,.50]);xx=np.arange(len(d));ax.scatter(xx,d.delta_mean,color='black',label='All cells',s=30)
  for x,r in zip(xx,d.itertuples()):
   if np.isfinite(r.min_delta):ax.plot([x,x],[r.min_delta,r.max_delta],lw=5,color='#37a6a0')
  ax.axhline(0,color='#555',lw=.8);ax.set_xticks(xx,d.label,rotation=55,ha='right',fontsize=9);ax.set_ylabel('Core mean difference: T2 - Ctrl');ax.legend(frameon=False);f.text(.075,.14,'Teal range: five QC-matched effect estimates (not a confidence interval).\nStable tissue clusters are a separate sensitivity analysis; their labels are not interchangeable with original C labels.',fontsize=11)
  finish(f,pdf,'绿色范围表示五次QC匹配结果范围，不是置信区间。完整稳定分群分析保存在表里，不把S群直接等同原来的C群。')
  f=new('Evidence summary and next interpretation','Cytotoxic module results, GO support and sensitivity are reported separately.',PAGER+'; '+BORRA)
  d['all5']=d.same.eq(5)&d.runs.eq(5);d['score_agrees']=np.sign(d.delta_mean)==np.sign(d.delta_control_score);d['rank']=d.delta_mean.abs();rank=d.sort_values('rank',ascending=False).head(8)
  ax=f.add_axes([.045,.35,.91,.43]);ax.axis('off');rows=[[r.label,f'{r.delta_mean:+.3f}',f'{r.q_all:.2g}',str(int(r.same))+'/5' if np.isfinite(r.same) else 'NA','yes' if r.score_agrees else 'no',f'{r.max_smd:.2f}' if np.isfinite(r.max_smd) else 'NA'] for r in rank.itertuples()];tab=ax.table(cellText=rows,colLabels=['Contrast','Core delta','BH q (cell level)','QC same direction','Control-score agrees','Max QC |SMD|'],bbox=[0,0,1,1],cellLoc='left',colLoc='left');tab.auto_set_font_size(False);tab.set_fontsize(11)
  f.text(.05,.23,'Ranked by absolute expression difference, not selected for significance.\nA cytotoxic transcript signature is not a killing assay; GO terms can share the same driver genes.\nIndependent biological replicates and functional assays are required to confirm treatment effects.',fontsize=13,linespacing=1.5)
  finish(f,pdf,'总结按表达差异绝对值排序。优先看基因组成、QC匹配、背景校正评分是否同向，再讨论GO机制；不要把转录结果写成已验证杀伤增强或减弱。')
 pd.DataFrame(selected).to_csv(RES/'tables/displayed_GO_terms.csv',index=False);pd.DataFrame(PAGES).to_csv(PRES/'slide_index.csv',index=False)
 prs=Presentation();prs.slide_width=Inches(16);prs.slide_height=Inches(9)
 for p in PAGES:
  s=prs.slides.add_slide(prs.slide_layouts[6]);s.shapes.add_picture(str(OUT/p['figure']),0,0,width=Inches(16),height=Inches(9));s.notes_slide.notes_text_frame.text=p['notes_zh']+'\n'+p['reference']+'\nSource: '+p['figure']
 prs.save(PRES/'iNKT_cytotoxicity_GO_20260919.pptx')
 thumbs=[]
 for p in PAGES:
  im=Image.open(OUT/p['figure']).convert('RGB');im.thumbnail((480,270));tile=Image.new('RGB',(490,300),'#ddd');tile.paste(im,(5,5));ImageDraw.Draw(tile).text((8,280),str(p['page'])+' '+p['title'][:65],fill='black');thumbs.append(tile)
 sheet=Image.new('RGB',(490*3,300*((len(thumbs)+2)//3)),'white')
 for i,tile in enumerate(thumbs):sheet.paste(tile,((i%3)*490,(i//3)*300))
 sheet.save(PRES/'contact_sheet.jpg')
 d.to_csv(RES/'tables/cytotoxicity_integrated_summary.csv',index=False)
 # Chinese report with exact calculated numbers and slide-linked guidance.
 ntest=int((contrasts.status=='tested').sum());qcount=int(((contrasts.status=='tested')&(contrasts.q_all<=.05)).sum());go_prov=json.loads((SRC/'GO_provenance.json').read_text())
 lines=['# 2026-09-19｜iNKT 细胞毒性与 GO 分析','', '[结果 PDF](presentation/iNKT_cytotoxicity_GO_20260919.pdf) · [PPTX（含中文讲稿）](presentation/iNKT_cytotoxicity_GO_20260919.pptx) · [图总览](presentation/contact_sheet.jpg)','', '## 本轮实际完成','',f'- 最新分析输入：9月5日对象，15,532细胞 × 10,670基因；9月15日为展示更新，没有更新表达数据。保留原9个refined cluster，并补充三个组织的稳定分群敏感性。',f'- 新计算：PAGER核心与E1模块、背景校正评分、组织×群条件比较、五次QC匹配、组织内cluster-vs-rest差异表达、GO BP/MF/CC全库ORA、局部DPT与根细胞敏感性。',f'- 条件DE复用9月5日完整结果，核对基因集及各组细胞数，复制到本目录并记录SHA256。所有GO富集为本轮新计算。',f'- 输出 {len(PAGES)} 页论文式图表与PPTX；PPTX采用图像页，保留SVG和绘图源码可修改。代码和结果是真实文件，不是旧目录符号链接。','', '## 当前主要数值','', '**优先关注：骨髓C4杀伤核心模块T2升高（差值+0.086，BH q=0.025）；脾脏C3小幅降低（−0.023，q=0.0015）。两者五次QC匹配均同向，并在主要重叠稳定分群保持方向。**', '', '**签名依赖性必须一并汇报：加Ctla2a后的E1代表模块在这两群均不显著（q=0.93、0.53）。杀伤相关GO条目没有通过q≤0.05，最小q=0.632。因此结论限于局部核心效应转录信号，不是整个细胞毒性通路或杀伤能力的实验证据。**', '', '条件GO较强信号包括骨髓C0与脾脏C3的T2上调基因富集氧化磷酸化；多群另见核糖体/翻译相关条目。诸如“synapse translation”的条目可由核糖体基因驱动，不能据名称推断形成神经突触。', '', '| 组织/原分群 | T2−Ctrl核心均值 | 细胞层面BH q | QC同向次数 | 背景校正评分同向 |','|---|---:|---:|---:|---|']
 for r in rank.itertuples():lines.append(f'| {r.label} | {r.delta_mean:+.3f} | {r.q_all:.3g} | {str(int(r.same))+'/5' if np.isfinite(r.same) else 'NA（匹配后细胞不足）'} | {"是" if r.score_agrees else "否"} |')
 lines+=['','按绝对效应排序，并非独立生物学重复验证。NA表示QC匹配后每组不足20细胞，不是0次同向。模块均值是log1p归一化表达，不是杀伤率，也不是log2FC。',f'总共 {ntest} 条合格模块对照，{qcount} 条通过统一BH q≤0.05；包含原分群、稳定分群和三个重叠基因集，不能当作独立发现数量。','','## GO 方法与结果入口','',f'- 官方GO版本：`{go_prov["ontology_version"]}`。MGI GAF头和下载哈希见 sources/GO_provenance.json 与 download_manifest.json。', '- 仅小鼠注释；排除NOT和ND，保留IEA；只沿is_a/part_of向上继承，不沿regulates继承。精确symbol与明确Ensembl→MGI别名匹配。','- 每个对照的背景为实际检出的基因∩该GO分支已注释基因；集合大小5–500；上调/下调分开。基因筛选：gene BH≤0.05，|log2FC|≥0.25。','- 每个对照×方向×GO分支在全部可检验term中BH校正，包含零命中term；另给出每种分析跨所有对照的q_analysis_global。','- 杀伤相关GO补充核查按term名称匹配cytotox、cell killing、granzyme、exocytosis、lymphocyte mediated；全部显著及不显著条目保留在GO_cytotoxicity_targeted_audit.csv。', '- 图中富集不表示通路激活；不同GO条目可能由同一组基因驱动。主图取显著、有≥2命中的代表条目，Jaccard≥0.75去冗余，完整结果保留。','', '[完整GO结果](results/tables/GO_ORA_all.csv.gz) · [显著GO结果](results/tables/GO_ORA_significant.csv) · [背景与查询覆盖](results/tables/GO_contrast_coverage.csv) · [绘图GO条目](results/tables/displayed_GO_terms.csv)','', '## 拟时序方法与限制','', '- 每组织在9月5日的局部邻居图重新计算diffusion map，避免继承全局坐标。取最大连通分量；未纳入细胞数在root audit中逐项列出。','- 原前体面板Cd24a/Egr2/Hivep3中，Cd24a不在当前10,670基因集，实际以Egr2/Hivep3两基因评分最高的局部群选根。根的前体身份支持有限；三个根细胞做敏感性，不能把此排序当谱系。','- 两条件共享同一组织DPT的20个区间，每bin至少20细胞；仅连续有支持的区间平滑。胸腺细胞集中在不连续的DPT区间，没有可靠连续曲线；骨髓/脾脏根敏感性最低相关约0.76。跨组织不能比较速度。','- 这是转录状态排序，不是RNA velocity、分化速率或已验证谱系。','', '## 可复现与目录','', '- `code/analyze.py`：评分、GO、轨迹；`code/presentation.py`：图与报告；`code/verify.py`：数值及文件验证；`code/test_analysis.py`：统计单元测试。','- 从仓库根目录：`bash iNKT_by_date/2026-09-19/code/run.sh`。长任务在tmux运行，CPU-only，GPU全部禁用。','- `sources/` 保存GO原始注释、版本、输入哈希；`results/de/`含新marker DE和复用的条件DE；`results/figures/`含PNG/SVG；`presentation/`含PDF/PPTX和讲稿索引。','- 重跑需仓库原始输入及现有`.venv`，本日期目录未重复拷贝大体积表达矩阵。原分析文件未改写。','','## 解释边界','', '- 每组织×条件只有1个sample标签，没有可识别动物重复。细胞层面P值和QC匹配不能解决生物重复缺失。','- 这是iNKT免疫杀伤相关转录分析，与图中药物肝/肾毒性项目不同。','- 采用两篇文章的方法与presentation，不声称完整重现PAGER平台、PPI、GAFA随机森林/SCORPION。','', '## 方法来源','', '- [PAGER-scFGA](https://doi.org/10.3389/fbinf.2024.1336135)：Figure 5功能模块、GO细胞组分和拟时序功能曲线。','- [Borra et al. 2026](https://doi.org/10.3390/biology15070588)：§3.5/Figure 4c亚群功能富集展示。','- [GO官方本体](https://geneontology.org/docs/download-ontology/)与[官方注释](https://geneontology.org/docs/download-go-annotations/)。']
 (OUT/'README.md').write_text('\n'.join(lines)+'\n')
 notes=['# 逐页中文讲稿','']
 for p in PAGES:notes.extend([f'## {p["page"]}. {p["title"]}',p['notes_zh'],''])
 (OUT/'notes/presentation_guide.zh.md').write_text('\n'.join(notes))
 (OUT/'notes/method_adaptation.md').write_text('论文只提供方法参考，未将PDF中的指令当作执行要求。\n\nPAGER：保留模块均值与拟时序曲线概念。核心3基因是透明简化；E1四基因是文中描述的代表基因，不冒充完整作者补充表基因集。\nGAFA：复用群特异上调基因→功能富集→气泡矩阵的展示逻辑；未运行预测模型或SCORPION。\n新加官方鼠GO BP/MF/CC、NOT/ND过滤、祖先传播、对照特异背景、跨对照校正列、五次QC匹配。\n原文固定601平滑窗不适合本数据；采用20个公共DPT区间、最小20细胞、连续有效区间内最多5-bin二阶Savitzky–Golay。\n')
 print('Completed presentation:',len(PAGES),'pages',flush=True)
if __name__=='__main__':main()
