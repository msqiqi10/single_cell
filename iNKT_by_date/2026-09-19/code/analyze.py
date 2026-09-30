"""Paper-guided iNKT cytotoxicity/GO analysis; CPU-only, frozen inputs, no source edits."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import argparse, collections, functools, gzip, hashlib, json, time, shutil, importlib.metadata
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.stats import hypergeom, mannwhitneyu, spearmanr
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
RES=OUT/'results'; SRC=OUT/'sources'; PREV=ROOT/'output/iNKT_meeting_followup_20260905'
KEY='cluster_c5_split_20260830'; TISSUES=['bone_marrow','spleen','thymus']
ORDER=['C0','C1','C2','C3','C4','C5-1','C5-2','C6','C7']; SEED=20260919
MODULES={'PAGER_core':['Prf1','Gzma','Gzmb'],'PAGER_E1':['Prf1','Ctla2a','Gzma','Gzmb'],'prior_effector':['Prf1','Gzma','Gzmb','Nkg7']}
NAMESPACES={'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'}
def log(s):print(time.strftime('%H:%M:%S'),s,flush=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def js(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)))
def csv(x,name):x.to_csv(RES/'tables'/name,index=False)
def bh(p):
 p=np.asarray(p,float); order=np.argsort(p); q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q

def read_obo(path):
 terms={}; alt={}; current=None; version=''
 for line in Path(path).read_text().splitlines():
  if line.startswith('data-version:'):version=line.split(': ',1)[1]
  if line=='[Term]':current={'parents':[],'alt':[]};continue
  if line.startswith('['):current=None;continue
  if current is None:continue
  if line.startswith('id: '):current['id']=line[4:];terms[line[4:]]=current
  elif line.startswith('name: '):current['name']=line[6:]
  elif line.startswith('namespace: '):current['namespace']=line[11:]
  elif line.startswith('alt_id: '):current['alt'].append(line[8:])
  elif line.startswith('is_a: '):current['parents'].append(line[6:].split()[0])
  elif line.startswith('relationship: part_of '):current['parents'].append(line.split()[2])
  elif line=='is_obsolete: true':current['obsolete']=True
 terms={k:v for k,v in terms.items() if not v.get('obsolete',False)}
 for k,v in terms.items():
  for x in v['alt']:alt[x]=k
 return terms,alt,version

def go_library(genes):
 terms,alt,version=read_obo(SRC/'go-basic.obo')
 @functools.lru_cache(None)
 def ancestors(t):
  result={t}
  for p in terms[t]['parents']:
   if p in terms:result.update(ancestors(p))
  return frozenset(result)
 # Exact symbol mapping first; Ensembl-to-current-MGI aliases only when unambiguous.
 base=ad.read_h5ad(PREV/'objects/scored_base.h5ad',backed='r')
 ens=pd.read_csv(PREV/'sources/MGI_MRK_ENSEMBL.rpt',sep='\t',header=None,dtype=str).fillna('')
 ensmap=collections.defaultdict(set)
 for r in ens.itertuples(index=False,name=None):ensmap[r[5]].add(r[1])
 symbol_map=collections.defaultdict(set)
 for old,eid in zip(base.var_names,base.var.gene_ids):
  symbol_map[str(old)].add(str(old))
  for sym in ensmap[str(eid)]:symbol_map[sym].add(str(old))
 base.file.close()
 annotations=collections.defaultdict(set); evidence=collections.Counter(); headers=[]; mapping_rows=[]; skipped=collections.Counter()
 with gzip.open(SRC/'mgi.gaf.gz','rt') as f:
  for line in f:
   if line.startswith('!'):headers.append(line.rstrip());continue
   x=line.rstrip('\n').split('\t')
   if len(x)<15:continue
   if 'NOT' in x[3].split('|'):skipped['NOT']+=1;continue
   if x[6]=='ND':skipped['ND']+=1;continue
   if x[12].split('|')[0]!='taxon:10090':skipped['non_mouse']+=1;continue
   t=alt.get(x[4],x[4]);candidates=symbol_map.get(x[2],set())
   if t not in terms or len(candidates)!=1:continue
   gene=next(iter(candidates))
   if gene not in genes:continue
   evidence[x[6]]+=1
   for parent in ancestors(t):annotations[parent].add(gene)
   mapping_rows.append((gene,x[2],t,x[6],x[3]))
 meta=[]
 for t,gs in annotations.items():
  v=terms[t]; meta.append({'go_id':t,'term':v['name'],'namespace':NAMESPACES[v['namespace']],'n_measured_genes':len(gs)})
 csv(pd.DataFrame(meta),'GO_term_catalog.csv')
 csv(pd.DataFrame(mapping_rows,columns=['gene','MGI_symbol','direct_go_id','evidence','qualifier']).drop_duplicates(),'GO_direct_annotation_audit.csv.gz')
 libraries={ns:{t:gs for t,gs in annotations.items() if NAMESPACES[terms[t]['namespace']]==ns} for ns in ['BP','MF','CC']}
 with (SRC/'GO_mouse_measured.gmt').open('w') as f:
  for t in sorted(annotations):f.write('\t'.join([t,terms[t]['name'],*sorted(annotations[t])])+'\n')
 js(SRC/'GO_provenance.json',{'ontology_version':version,'gaf_headers':headers,'evidence_counts':dict(evidence),'excluded':dict(skipped),'mapping':'exact symbol plus unambiguous Ensembl-to-MGI alias from frozen 20260905 report','propagation':['is_a','part_of'],'excluded_propagation':['regulates','positively_regulates','negatively_regulates'],'term_size_bounds':[5,500],'evidence_policy':'all non-NOT, non-ND mouse annotations, including IEA','tested_universe':'per contrast expressed genes intersect namespace-annotated genes'})
 return terms,libraries

def ora(query,universe,sets,terms):
 query=set(query)&set(universe);M=len(universe);n=len(query);rows=[]
 for term,genes in sets.items():
  measured=genes&universe;K=len(measured)
  if 5<=K<=500:
   hit=query&measured;k=len(hit)
   rows.append({'go_id':term,'term':terms[term]['name'],'M':M,'K':K,'n':n,'k':k,'gene_ratio':k/n if n else 0,'enrichment_ratio':k/(n*K/M) if n else 0,'genes':';'.join(sorted(hit))})
 d=pd.DataFrame(rows)
 if len(d):
  d['pvalue']=hypergeom.sf(d.k.to_numpy()-1,M,d.K.to_numpy(),n) if n else 1.
  d['q_family']=bh(d.pvalue)
 return d

def qc_match(obs,seed):
 bins=[]
 for x in ['total_counts','n_genes_by_counts','pct_counts_mt']:
  v=obs[x].astype(float);bins.append(pd.qcut(v,4,labels=False,duplicates='drop').fillna(0).to_numpy(int))
 codes=np.ravel_multi_index(tuple(bins),(4,4,4));rng=np.random.default_rng(seed);cond=obs.condition.astype(str).to_numpy();ids=[]
 for b in np.unique(codes):
  c=np.flatnonzero((codes==b)&(cond=='Ctrl'));t=np.flatnonzero((codes==b)&(cond=='T2'));n=min(len(c),len(t))
  if n:ids.extend(rng.choice(c,n,replace=False));ids.extend(rng.choice(t,n,replace=False))
 return np.asarray(ids,int)

def score_stage(a):
 coverage=[]
 for module,gs in MODULES.items():
  keep=[g for g in gs if g in a.raw.var_names]
  for g in gs:
   ix=list(a.raw.var_names).index(g) if g in a.raw.var_names else None
   coverage.append({'module':module,'gene':g,'measured':ix is not None,'n_detected':int((a.raw.X[:,ix]>0).sum()) if ix is not None else 0})
  if len(keep)<2:raise ValueError('Insufficient module coverage '+module)
  a.obs[module+'_mean']=np.asarray(a.raw[:,keep].X.mean(axis=1)).ravel()
  sc.tl.score_genes(a,keep,score_name=module+'_score',random_state=SEED,use_raw=True,ctrl_size=50,n_bins=25)
 csv(pd.DataFrame(coverage),'module_gene_coverage.csv')
 rows=[];sample_rows=[];sens=[]
 for definition,key in [('original',KEY),('stable','stable_cluster')]:
  for (t,c),o in a.obs.groupby(['tissue',key],observed=True):
   for module in MODULES:
    col=module+'_mean';ctrl=o.loc[o.condition=='Ctrl',col].to_numpy();tum=o.loc[o.condition=='T2',col].to_numpy();ok=min(len(ctrl),len(tum))>=20
    row={'definition':definition,'tissue':str(t),'cluster':str(c),'module':module,'n_Ctrl':len(ctrl),'n_T2':len(tum),'status':'tested' if ok else 'insufficient_cells','delta_mean':np.mean(tum)-np.mean(ctrl) if len(ctrl)*len(tum) else np.nan,'mean_Ctrl':np.mean(ctrl) if len(ctrl) else np.nan,'mean_T2':np.mean(tum) if len(tum) else np.nan,'pvalue':np.nan,'rank_biserial':np.nan}
    if ok:
     u=mannwhitneyu(tum,ctrl,alternative='two-sided');row.update(pvalue=u.pvalue,rank_biserial=2*u.statistic/(len(tum)*len(ctrl))-1)
     row['delta_control_score']=o.loc[o.condition=='T2',module+'_score'].mean()-o.loc[o.condition=='Ctrl',module+'_score'].mean()
     for seed in [0,17,42,123,20260919]:
      m=o.iloc[qc_match(o,seed)];nc=int((m.condition=='Ctrl').sum())
      if nc>=20:
       delta=m.loc[m.condition=='T2',col].mean()-m.loc[m.condition=='Ctrl',col].mean()
       smds={}
       for metric in ['total_counts','n_genes_by_counts','pct_counts_mt']:
        x=m.loc[m.condition=='T2',metric].to_numpy(float);y=m.loc[m.condition=='Ctrl',metric].to_numpy(float)
        if metric!='pct_counts_mt':x=np.log1p(x);y=np.log1p(y)
        sd=np.sqrt((x.var()+y.var())/2);smds[metric]=(x.mean()-y.mean())/sd if sd else 0
       sens.append({**{k:row[k] for k in ['definition','tissue','cluster','module']},'seed':seed,'n_each':nc,'delta_mean':delta,'same_direction':np.sign(delta)==np.sign(row['delta_mean']),'max_abs_QC_SMD':max(abs(v) for v in smds.values())})
    rows.append(row)
   for sample,g in o.groupby('sample',observed=True):
    for module in MODULES:sample_rows.append({'definition':definition,'tissue':str(t),'cluster':str(c),'sample':sample,'condition':str(g.condition.iloc[0]),'n_cells':len(g),'module':module,'mean_expression':g[module+'_mean'].mean(),'control_score':g[module+'_score'].mean()})
 d=pd.DataFrame(rows);d['q_all']=np.nan;ok=d.pvalue.notna();d.loc[ok,'q_all']=bh(d.loc[ok,'pvalue']);csv(d,'cytotoxicity_contrasts.csv');csv(pd.DataFrame(sens),'cytotoxicity_QC_sensitivity.csv');csv(pd.DataFrame(sample_rows),'cytotoxicity_sample_summary.csv')
 a.obs.reset_index(names='cell').to_csv(RES/'tables/cell_metadata_scores.csv.gz',index=False)
 genes=sorted(set(sum(MODULES.values(),[])+['Fasl','Ctsw','Ifng','Tbx21','Il4','Rorc','Cd24a','Egr2','Hivep3']))
 dot=[]
 for (t,c,cond),o in a.obs.groupby(['tissue',KEY,'condition'],observed=True):
  ix=a.obs_names.get_indexer(o.index)
  for g in genes:
   if g not in a.raw.var_names:continue
   x=a.raw[ix,[g]].X
   dot.append({'tissue':str(t),'cluster':str(c),'condition':str(cond),'gene':g,'n_cells':len(ix),'mean_log1p':float(x.mean()),'pct_positive':float((x>0).mean())})
 csv(pd.DataFrame(dot),'marker_dotplot_values.csv')


def enrichment_stage(a):
 terms,libraries=go_library(set(a.var_names));all_results=[];status=[];manifest=[]
 annotated={ns:set().union(*lib.values()) for ns,lib in libraries.items()}
 def enrich(d,uid,tissue,cluster,analysis,definition):
  expressed=set(d.loc[(d.pct_expressing_tumor>0)|(d.pct_expressing_control>0),'gene'])
  for direction in (['cluster_up'] if analysis=='cluster_identity' else ['T2_up','Ctrl_up']):
   selected=set(d.loc[(d.pvals_adj<=.05)&((d.logfoldchanges>=.25) if direction!='Ctrl_up' else (d.logfoldchanges<=-.25)),'gene'])
   for ns,lib in libraries.items():
    universe=expressed&annotated[ns];r=ora(selected,universe,lib,terms)
    if len(r):
     for k,v in {'unit_id':uid,'tissue':tissue,'cluster':cluster,'analysis':analysis,'definition':definition,'direction':direction,'namespace':ns}.items():r[k]=v
     all_results.append(r)
    status.append({'unit_id':uid,'analysis':analysis,'definition':definition,'direction':direction,'namespace':ns,'n_expressed':len(expressed),'n_annotated_background':len(universe),'n_selected_before_annotation':len(selected),'n_query_annotated':len(selected&universe),'n_tested_terms':len(r)})
 # New within-tissue cluster-vs-rest DE for Borra Fig 4c, on the existing original cluster labels.
 for tissue in TISSUES:
  b=a[a.obs.tissue.astype(str)==tissue].copy();b.obs[KEY]=b.obs[KEY].cat.remove_unused_categories()
  groups=[c for c in ORDER if (b.obs[KEY]==c).sum()>=20 and (b.obs[KEY]!=c).sum()>=20]
  log('Within-tissue marker DE '+tissue)
  sc.tl.rank_genes_groups(b,KEY,groups=groups,reference='rest',method='wilcoxon',tie_correct=True,use_raw=True,pts=True,n_genes=b.raw.n_vars)
  for c in groups:
   d=sc.get.rank_genes_groups_df(b,group=c).rename(columns={'names':'gene','pct_nz_group':'pct_expressing_tumor','pct_nz_reference':'pct_expressing_control'})
   uid='identity__'+tissue+'__'+c;p=RES/'de'/(uid+'.csv.gz');d.to_csv(p,index=False)
   enrich(d,uid,tissue,c,'cluster_identity','original')
 # Reuse validated full DE rankings, copy into this self-contained results folder and hash.
 for definition,fn in [('original','de_status_original.csv'),('stable','de_status_stable.csv')]:
  ds=pd.read_csv(PREV/'tables'/fn)
  for row in ds.itertuples():
   if row.status!='computed' or not row.strict_min20:continue
   if definition=='original' and row.scope=='global':continue
   p=PREV/'de'/(row.unit_id+'.csv.gz');d=pd.read_csv(p)
   assert set(d.gene)==set(a.var_names),row.unit_id
   if definition=='stable':mask=(a.obs.stable_cluster.astype(str)==row.cluster)
   else:mask=(a.obs.tissue.astype(str)==row.tissue)&((a.obs[KEY].astype(str)==row.cluster) if row.scope=='cluster_tissue' else True)
   counts=a.obs.loc[mask,'condition'].value_counts()
   assert int(counts.get('Ctrl',0))==row.n_control and int(counts.get('T2',0))==row.n_tumor,row.unit_id
   shutil.copy2(p,RES/'de'/p.name);manifest.append({'file':str(p.relative_to(ROOT)),'sha256':sha(p),'copied_to':str((RES/'de'/p.name).relative_to(OUT))})
   enrich(d,row.unit_id,row.tissue,row.cluster,'condition','tissue' if getattr(row,'scope','')=='tissue' else definition)
 result=pd.concat(all_results,ignore_index=True)
 result['q_analysis_global']=result.groupby('analysis',sort=False).pvalue.transform(lambda x:bh(x.to_numpy()))
 csv(result,'GO_ORA_all.csv.gz');csv(result[result.q_family<=.05].sort_values(['analysis','q_family']),'GO_ORA_significant.csv')
 csv(pd.DataFrame(status),'GO_contrast_coverage.csv');js(SRC/'reused_DE_manifest.json',manifest)
 log(f'GO tests: {len(result)}, q_family <= .05: {(result.q_family<=.05).sum()}')
 # CC memberships for the actual cytotoxic genes, no invented PPI edges.
 cc=[]
 for g in sorted(set(sum(MODULES.values(),[]))):
  for t,gs in libraries['CC'].items():
   if g in gs and 5<=len(gs)<=500:cc.append({'gene':g,'go_id':t,'term':terms[t]['name'],'n_measured_genes':len(gs)})
 csv(pd.DataFrame(cc),'cytotoxic_gene_CC_membership.csv')


def trajectory_stage(a):
 audits=[];cells=[];bins=[];correlations=[]
 for tissue in TISSUES:
  log('Tissue-specific diffusion pseudotime '+tissue)
  b=ad.read_h5ad(PREV/f'objects/{tissue}_stable.h5ad')
  for module in MODULES:
   b.obs[module+'_mean']=a.obs.loc[b.obs_names,module+'_mean'].values
  # Do not reuse inherited global diffusion coordinates after within-tissue reclustering.
  ncomp,labels=connected_components(b.obsp['connectivities'],directed=False)
  largest=np.argmax(np.bincount(labels));selected=labels==largest;d=b[selected].copy()
  sc.tl.diffmap(d,n_comps=15,random_state=SEED)
  roots=['Cd24a','Egr2','Hivep3'];present=[g for g in roots if g in d.raw.var_names]
  sc.tl.score_genes(d,present,score_name='root_precursor_score',use_raw=True,random_state=SEED)
  cluster_means=d.obs.groupby('meeting_cluster',observed=True).root_precursor_score.mean()
  root_cluster=str(cluster_means.idxmax());ix=np.flatnonzero(d.obs.meeting_cluster.astype(str).to_numpy()==root_cluster)
  dc=d.obsm['X_diffmap'][:,1:10];center=np.median(dc[ix],axis=0);ordered=ix[np.argsort(np.linalg.norm(dc[ix]-center,axis=1))]
  root_ids=[ordered[0],ordered[min(len(ordered)-1,max(1,len(ordered)//4))],ordered[min(len(ordered)-1,max(2,len(ordered)//2))]]
  primary=None;root_corr=[]
  for j,r in enumerate(root_ids):
   d.uns['iroot']=int(r);sc.tl.dpt(d,n_dcs=10);pt=d.obs.dpt_pseudotime.to_numpy().copy()
   if j==0:primary=pt
   good=np.isfinite(primary)&np.isfinite(pt);root_corr.append(float(spearmanr(primary[good],pt[good]).statistic))
  d.obs['local_dpt']=primary
  audits.append({'tissue':tissue,'n_cells':b.n_obs,'n_connected_components':ncomp,'n_analyzed':d.n_obs,'n_excluded_disconnected':b.n_obs-d.n_obs,'root_cluster':root_cluster,'root_genes':';'.join(present),'root_cell':d.obs_names[root_ids[0]],'root_sensitivity_cells':';'.join(d.obs_names[root_ids]),'root_rho_1':root_corr[1],'root_rho_2':root_corr[2],'interpretation':'assumed precursor-rooted transcriptional ordering; not lineage proof or velocity'})
  for j,cell in enumerate(d.obs_names):cells.append({'cell':cell,'tissue':tissue,'condition':str(d.obs.condition.iloc[j]),'cluster':str(d.obs[KEY].iloc[j]),'stable_cluster':str(d.obs.meeting_cluster.iloc[j]),'local_dpt':float(primary[j]),'precursor_score':float(d.obs.root_precursor_score.iloc[j]),**{m:float(d.obs[m+'_mean'].iloc[j]) for m in MODULES}})
  for cond in ['Ctrl','T2']:
   sub=d.obs[(d.obs.condition==cond)&np.isfinite(d.obs.local_dpt)].copy();sub['bin']=pd.cut(sub.local_dpt,bins=np.linspace(0,1,21),include_lowest=True,labels=False)
   for module in MODULES:
    rho,p=spearmanr(sub.local_dpt,sub[module+'_mean']);correlations.append({'tissue':tissue,'condition':cond,'module':module,'n_cells':len(sub),'rho':rho,'pvalue_cell_level':p})
    for k,g in sub.groupby('bin',observed=True):bins.append({'tissue':tissue,'condition':cond,'module':module,'bin':int(k),'pseudotime':g.local_dpt.mean(),'n_cells':len(g),'mean_expression':g[module+'_mean'].mean(),'sd':g[module+'_mean'].std()})
  small=ad.AnnData(X=sparse.csr_matrix((d.n_obs,0)),obs=d.obs[[KEY,'meeting_cluster','condition','local_dpt','root_precursor_score']].copy());small.obsm['X_umap']=d.obsm['X_umap'];small.obsm['X_diffmap']=d.obsm['X_diffmap'];small.uns['root_audit']=audits[-1];small.write_h5ad(RES/'objects'/f'{tissue}_local_trajectory.h5ad')
 csv(pd.DataFrame(audits),'trajectory_root_audit.csv');csv(pd.DataFrame(cells),'trajectory_cell_values.csv.gz');csv(pd.DataFrame(bins),'trajectory_binned_curves.csv');csv(pd.DataFrame(correlations),'trajectory_correlations.csv')


def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['all','scores','go','trajectory']);args=parser.parse_args()
 start=time.time();a=ad.read_h5ad(PREV/'objects/scored_base.h5ad');assert a.shape==(15532,10670)
 a.obs['stable_cluster']=''
 for tissue in TISSUES:
  b=ad.read_h5ad(PREV/f'objects/{tissue}_stable.h5ad',backed='r');a.obs.loc[b.obs_names,'stable_cluster']=b.obs.meeting_cluster.astype(str).values;b.file.close()
 assert (a.obs.stable_cluster!='').all()
 js(SRC/'input_manifest.json',{'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PREV/'objects/scored_base.h5ad',*[PREV/f'objects/{t}_stable.h5ad' for t in TISSUES],ROOT/'docs/pager-scFGA.pdf',ROOT/'docs/Borra_et_al_2026_GAFA_CML_NK.pdf',PREV/'sources/MGI_MRK_ENSEMBL.rpt']},'packages':{m:importlib.metadata.version(m) for m in ['numpy','pandas','anndata','scanpy','scipy','matplotlib','python-pptx']},'shape':a.shape,'seed':SEED,'GPU_used':False,'source_dates':'20260905 analysis objects; 20260915 presentation only'})
 if args.stage in ['all','scores','trajectory']:log('Scoring modules');score_stage(a)
 if args.stage in ['all','go']:enrichment_stage(a)
 if args.stage in ['all','trajectory']:trajectory_stage(a)
 js(RES/'logs'/f'{args.stage}.completed.json',{'elapsed_seconds':time.time()-start,'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())});log('COMPLETE '+args.stage)
if __name__=='__main__':main()
