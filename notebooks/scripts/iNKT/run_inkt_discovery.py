"""Open discovery on the current iNKT pipeline, with cross-clustering and QC sensitivities.
No paper/PPT gene list or named pathway enters candidate selection.
"""
from __future__ import annotations
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import argparse,json,time,re,hashlib
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_inkt_meeting_followup import ROOT,INPUT,REFINED_CLUSTER_KEY,REFINED_ORDER,TISSUES,DEUnit,run_de_unit,ora_table,bh_adjust,sha256_file
from run_inkt_meeting_pathways import libraries,HEAT
PREV=ROOT/'output/iNKT_meeting_followup_20260905'
OUT=ROOT/'output/iNKT_discovery_20260905'
SEED=20260905
QC=['log_total_counts','log_n_genes','pct_counts_mt']

def log(s):print(time.strftime('%Y-%m-%d %H:%M:%S'),s,flush=True)
def js(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)))
def csv(d,name):d.to_csv(OUT/'tables'/name,index=False)
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def family(g):
    if re.match(r'^(Hsp|Dnaj)',g,re.I):return 'heat_shock'
    if re.match(r'^(Rpl|Rps|Mrpl|Mrps)',g):return 'ribosome'
    if re.match(r'^(mt-|Cox\d|Nduf|Uqcr|Atp5)',g,re.I):return 'respiratory_chain_or_mt'
    return 'other'
def setup():
    for x in ['tables','de_qc','enrichment','figures','logs']:(OUT/x).mkdir(parents=True,exist_ok=True)

def units_and_obs():
    a=ad.read_h5ad(PREV/'objects/scored_base.h5ad');a.obs['stable_cluster']=''
    for t in TISSUES:
        b=ad.read_h5ad(PREV/f'objects/{t}_stable.h5ad',backed='r');a.obs.loc[b.obs_names,'stable_cluster']=b.obs.meeting_cluster.astype(str);b.file.close()
    ori=pd.read_csv(PREV/'tables/de_status_original.csv');ori=ori[ori.strict_min20&ori.status.eq('computed')].copy();ori['definition']='original'
    sta=pd.read_csv(PREV/'tables/de_status_stable.csv');sta=sta[sta.strict_min20&sta.status.eq('computed')].copy();sta['definition']='stable';sta['scope']='stable_recluster'
    units=pd.concat([ori,sta],ignore_index=True);a.obs['log_total_counts']=np.log1p(a.obs.total_counts);a.obs['log_n_genes']=np.log1p(a.obs.n_genes_by_counts)
    return a,units

def mask_unit(obs,r):
    m=np.ones(len(obs),bool)
    if r.tissue!='all':m &= obs.tissue.astype(str).eq(r.tissue).to_numpy()
    if r.definition=='stable':m &= obs.stable_cluster.astype(str).eq(r.cluster).to_numpy()
    elif r.scope=='cluster_tissue':m &= obs[REFINED_CLUSTER_KEY].astype(str).eq(r.cluster).to_numpy()
    return m

def qc_strata(obs):
    """Condition-blind quartiles; equal sampling within each occupied QC stratum."""
    strata=[]
    for col in QC:
        x=obs[col].astype(float)
        if not np.isfinite(x).all():raise ValueError('Missing QC covariate')
        if x.nunique()<=1:code=np.zeros(len(x),int)
        else:code=pd.qcut(x,q=4,labels=False,duplicates='drop').fillna(0).to_numpy(int)
        strata.append(code)
    return np.ravel_multi_index(tuple(strata),(4,4,4))

def matched_indices(obs,seed=SEED):
    bins=qc_strata(obs);rng=np.random.default_rng(seed);cond=obs.condition.astype(str).to_numpy();ids=[];rows=[]
    for b in sorted(set(bins)):
        c=np.flatnonzero((bins==b)&(cond=='Ctrl'));t=np.flatnonzero((bins==b)&(cond=='T2'));n=min(len(c),len(t))
        if n:ids.extend(rng.choice(c,n,replace=False));ids.extend(rng.choice(t,n,replace=False))
        rows.append({'stratum':int(b),'n_Ctrl_available':len(c),'n_T2_available':len(t),'n_each_selected':n})
    return np.array(sorted(ids),int),pd.DataFrame(rows)

def balance(obs):
    rows={};cond=obs.condition.astype(str)
    for col in QC:
        c=obs.loc[cond.eq('Ctrl'),col].to_numpy(float);t=obs.loc[cond.eq('T2'),col].to_numpy(float)
        pooled=np.sqrt((np.var(c,ddof=1)+np.var(t,ddof=1))/2) if min(len(c),len(t))>1 else np.nan
        rows[col+'_SMD']=(np.mean(t)-np.mean(c))/pooled if pooled>0 else (0. if np.mean(t)==np.mean(c) else np.nan)
    return rows

def de_and_crosswalk():
    a,units=units_and_obs();csv(units,'analysis_units.csv');cw=[];freq=[];audit=[];members=[];strata=[]
    for t in TISSUES:
        o=a.obs[a.obs.tissue.astype(str).eq(t)];table=pd.crosstab(o[REFINED_CLUSTER_KEY],o.stable_cluster)
        for c in table.index:
            for s in table.columns:
                n=int(table.loc[c,s]);den=int(table.loc[c].sum());den_s=int(table[s].sum())
                if n:cw.append({'tissue':t,'original_cluster':c,'stable_cluster':s,'overlap':n,'fraction_original':n/den,'fraction_stable':n/den_s,'Jaccard':n/(den+den_s-n)})
    csv(pd.DataFrame(cw),'original_stable_crosswalk.csv')
    for r in units.itertuples():
        sub=a[mask_unit(a.obs,r)].copy();ns=sub.obs.condition.astype(str).value_counts();ids,st=matched_indices(sub.obs);st['unit_id']=r.unit_id;strata.append(st)
        match=sub[ids].copy();n=len(match)//2;row={'unit_id':r.unit_id,'definition':r.definition,'tissue':r.tissue,'cluster':r.cluster,'n_Ctrl':int(ns['Ctrl']),'n_T2':int(ns['T2']),'n_each_matched':n,'retained_minority_fraction':n/min(ns),'status':'computed' if n>=20 else 'insufficient_matched_cells'}
        row.update({'before_'+k:v for k,v in balance(sub.obs).items()});row.update({'after_'+k:v for k,v in balance(match.obs).items()})
        for cell in match.obs_names:members.append({'unit_id':r.unit_id,'cell_id':cell})
        if n>=20:
            log(f'QC-matched DE {r.unit_id}: {n} cells per condition');match.obs[REFINED_CLUSTER_KEY]=pd.Categorical([str(r.cluster)]*len(match));d=run_de_unit(match,DEUnit(r.unit_id,'QC_matched',str(r.cluster),None));d['tissue']=r.tissue;d.to_csv(OUT/f'de_qc/{r.unit_id}.csv.gz',index=False)
            row['n_robust_up']=int(d.robust_tumor_up.sum());row['n_robust_down']=int(d.robust_control_up.sum())
        audit.append(row);csv(pd.DataFrame(audit),'QC_matching_audit.csv')
        if r.scope in ['cluster_tissue','stable_recluster']:
            den=a.obs[a.obs.tissue.astype(str).eq(r.tissue)].condition.astype(str).value_counts();freq.append({'unit_id':r.unit_id,'tissue':r.tissue,'cluster':r.cluster,'definition':r.definition,'n_Ctrl':int(ns['Ctrl']),'n_T2':int(ns['T2']),'Ctrl_frequency':int(ns['Ctrl'])/den['Ctrl'],'T2_frequency':int(ns['T2'])/den['T2'],'delta_pp':100*(int(ns['T2'])/den['T2']-int(ns['Ctrl'])/den['Ctrl'])})
    csv(pd.concat(strata,ignore_index=True),'QC_stratum_counts.csv');csv(pd.DataFrame(members),'QC_selected_cell_ids.csv.gz');csv(pd.DataFrame(freq),'population_frequencies.csv')
    js(OUT/'logs/de.completed.json',{'units':len(units),'matching':'4 condition-blind quantiles each: logUMI, logGenes, pctMT; balance Ctrl/T2 within joint strata','seed':SEED,'inference':'descriptive QC sensitivity; post-condition covariates may include biological effects; not causal adjustment','input_SHA256':sha256_file(INPUT)})

def enrich():
    units=pd.read_csv(OUT/'tables/analysis_units.csv');lib=libraries();paths=[]
    for r in units.itertuples():
        for analysis in ['original_counts','QC_matched']:
            dp=PREV/f'de/{r.unit_id}.csv.gz' if analysis=='original_counts' else OUT/f'de_qc/{r.unit_id}.csv.gz'
            if not dp.exists():continue
            d=pd.read_csv(dp);universe=set(d.gene)
            for name in ['KEGG_Mouse_2019','Reactome_2022_ortholog']:
                for sens in ['full','without_heat_shock']:
                    u={g for g in universe if sens=='full' or not HEAT.match(g)}
                    p=OUT/f'enrichment/{r.unit_id}__{analysis}__{name}__{sens}.csv.gz'
                    previous=PREV/f'enrichment/{r.unit_id}__{sens}__{name}__ORA_robust.csv.gz'
                    if analysis=='original_counts' and previous.exists():q=pd.read_csv(previous)
                    else:q=ora_table(d,lib[name],u,'robust')
                    q['unit_id']=r.unit_id;q['analysis']=analysis;q['library']=name;q['sensitivity']=sens;q.to_csv(p,index=False);paths.append(p)
        log('Discovery ORA '+r.unit_id)
    frames=[pd.read_csv(p) for p in paths];allq=pd.concat(frames,ignore_index=True)
    # Additional family-wide correction, alongside original per-contrast correction.
    allq['FDR_across_units']=allq.groupby(['analysis','library','sensitivity','direction']).pvalue.transform(lambda p:bh_adjust(p.to_numpy()))
    allq.to_csv(OUT/'tables/functional_screen_all.csv.gz',index=False)
    js(OUT/'logs/enrichment.completed.json',{'n_tests':len(allq),'selection_scope':'all eligible KEGG/Reactome sets; no paper/PPT target restriction','ORA_threshold':'gene BH<=0.05, abs(log2FC)>=0.25; measured background; size3-500','additional_FDR':'BH across all units within library/direction/sensitivity/analysis','human_Reactome':'strict one-to-one mouse mapping; source version frozen from meeting run'})

def select_candidates():
    units=pd.read_csv(OUT/'tables/analysis_units.csv');cw=pd.read_csv(OUT/'tables/original_stable_crosswalk.csv');qa=pd.read_csv(OUT/'tables/QC_matching_audit.csv').set_index('unit_id');de={};rows=[]
    def get(uid,kind='original'):
        key=(uid,kind)
        if key not in de:
            p=(PREV/'de' if kind=='original' else OUT/'de_qc')/f'{uid}.csv.gz';de[key]=pd.read_csv(p).set_index('gene') if p.exists() else pd.DataFrame()
        return de[key]
    for r in units[units.scope.eq('cluster_tissue')].itertuples():
        mapping=cw[cw.tissue.eq(r.tissue)&cw.original_cluster.eq(r.cluster)].sort_values('overlap',ascending=False).iloc[0];sid='stable__'+mapping.stable_cluster;orig=get(r.unit_id);stable=get(sid);qc=get(r.unit_id,'qc');sqc=get(sid,'qc');tissue=get('tissue__'+r.tissue);global_de=get('global')
        if not len(qc) or not len(stable):continue
        cand=orig[(orig.pvals_adj<=.05)&(orig.logfoldchanges.abs()>=.25)&(orig[['pct_expressing_tumor','pct_expressing_control']].max(axis=1)>=.1)]
        for gene,o in cand.iterrows():
            s=stable.loc[gene];q=qc.loc[gene];b=tissue.loc[gene];g=global_de.loc[gene];z=sqc.loc[gene] if len(sqc) else None;sign=np.sign(o.logfoldchanges)
            effects=[o.logfoldchanges,s.logfoldchanges,q.logfoldchanges]+([z.logfoldchanges] if z is not None else [])
            all_sign=all(np.sign(e)==sign for e in effects);minfc=min(abs(e) for e in effects)
            smdcols=[c for c in qa.columns if c.startswith('after_')];maxsmd=float(np.abs(pd.to_numeric(qa.loc[r.unit_id,smdcols])).max());same_stable=bool(mapping.fraction_original>=.5 and np.sign(s.logfoldchanges)==sign and s.pvals_adj<=.05 and abs(s.logfoldchanges)>=.25)
            qc_support=bool(np.sign(q.logfoldchanges)==sign and q.pvals_adj<=.1 and abs(q.logfoldchanges)>=.25)
            local=bool(abs(o.logfoldchanges-b.logfoldchanges)>=.25 and (b.pvals_adj>.05 or sign!=np.sign(b.logfoldchanges) or abs(o.logfoldchanges)>abs(b.logfoldchanges)+.25))
            tier='supported_candidate' if same_stable and qc_support and all_sign and minfc>=.25 and min(r.n_tumor,r.n_control)>=50 and qa.loc[r.unit_id,'n_each_matched']>=50 and maxsmd<=.2 else 'sensitivity_or_small_group'
            row={'unit_id':r.unit_id,'tissue':r.tissue,'cluster':r.cluster,'gene':gene,'gene_family':family(gene),'original_log2FC':o.logfoldchanges,'original_FDR':o.pvals_adj,'detected_Ctrl':o.pct_expressing_control,'detected_T2':o.pct_expressing_tumor,'stable_unit':sid,'original_cells_in_stable_fraction':mapping.fraction_original,'stable_log2FC':s.logfoldchanges,'stable_FDR':s.pvals_adj,'QC_log2FC':q.logfoldchanges,'QC_FDR':q.pvals_adj,'stable_QC_log2FC':z.logfoldchanges if z is not None else np.nan,'stable_QC_FDR':z.pvals_adj if z is not None else np.nan,'tissue_log2FC':b.logfoldchanges,'tissue_FDR':b.pvals_adj,'global_log2FC':g.logfoldchanges,'global_FDR':g.pvals_adj,'same_direction_all':all_sign,'minimum_abs_log2FC':minfc,'stable_support':same_stable,'QC_support':qc_support,'local_effect_contrast':local,'matched_each_n':qa.loc[r.unit_id,'n_each_matched'],'max_postmatch_QC_SMD':maxsmd,'tier':tier}
            rows.append(row)
    out=pd.DataFrame(rows).sort_values(['tier','local_effect_contrast','minimum_abs_log2FC','original_FDR'],ascending=[False,False,False,True])
    # Annotation only after discovery; no effect on inclusion or rank.
    old=pd.read_csv(ROOT/'output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_de_tables_extracted.csv');out['present_in_legacy_DEG_tables']=out.gene.isin(set(old.names))
    csv(out,'gene_candidates_all_sensitivities.csv');core=out[out.tier.eq('supported_candidate')];csv(core,'gene_candidates_supported.csv');csv(core[core.local_effect_contrast&core.gene_family.eq('other')],'gene_candidates_local_non_dominant.csv')
    # Same dataset subsampling quantifies sensitivity to cell selection, not animal uncertainty.
    subsample(core)
    screen=pd.read_csv(OUT/'tables/functional_screen_all.csv.gz');sig=screen[(screen.analysis=='original_counts')&(screen.sensitivity=='full')&(screen.fdr<=.05)&screen.unit_id.str.startswith('cluster_tissue')];indexed=screen.set_index(['unit_id','analysis','library','sensitivity','direction','term']);prows=[]
    for r in sig.itertuples():
        meta=units.set_index('unit_id').loc[r.unit_id];maps=cw[cw.tissue.eq(meta.tissue)&cw.original_cluster.eq(meta.cluster)].sort_values('overlap',ascending=False);sid='stable__'+maps.iloc[0].stable_cluster
        def lookup(uid,analysis,sens):
            k=(uid,analysis,r.library,sens,r.direction,r.term)
            return indexed.loc[k] if k in indexed.index else None
        q=lookup(r.unit_id,'QC_matched','full');s=lookup(sid,'original_counts','full');h=lookup(r.unit_id,'original_counts','without_heat_shock');genes=set(str(r.overlap_genes).split(';'))-{''};non=[g for g in genes if family(g)=='other']
        prows.append({'unit_id':r.unit_id,'library':r.library,'term':r.term,'direction':r.direction,'FDR':r.fdr,'FDR_across_units':r.FDR_across_units,'QC_FDR':q.fdr if q is not None else np.nan,'stable_unit':sid,'stable_FDR':s.fdr if s is not None else np.nan,'no_heat_FDR':h.fdr if h is not None else np.nan,'n_drivers':len(genes),'n_other_drivers':len(non),'drivers':';'.join(sorted(genes)),'other_drivers':';'.join(sorted(non)),'stable_and_QC_supported':bool(q is not None and s is not None and q.fdr<=.05 and s.fdr<=.05),'source':'unrestricted library screen; pathway name alone does not identify mechanism'})
    pp=pd.DataFrame(prows);csv(pp,'pathway_candidates_all.csv');pp['driver_key']=pp.drivers
    grouped=[]
    for keys,d in pp.groupby(['unit_id','direction','driver_key']):
        best=d.sort_values('FDR').iloc[0];grouped.append({'unit_id':keys[0],'direction':keys[1],'driver_genes':keys[2],'representative_term':best.term,'n_repeated_terms':len(d),'terms':' | '.join(d.term),'min_FDR':d.FDR.min(),'n_other_drivers':int(best.n_other_drivers),'any_stable_QC_supported':bool(d.stable_and_QC_supported.any())})
    csv(pd.DataFrame(grouped),'pathway_exact_driver_groups.csv')
    js(OUT/'logs/candidates.completed.json',{'screened_DEG_records':len(out),'supported_gene_records':len(core),'local_non_dominant_records':int((core.local_effect_contrast&core.gene_family.eq('other')).sum()),'pathway_records':len(pp),'selection':'full expression universe; effect/detection/cluster/QC criteria; paper and PPT do not select candidates','inference':'all cell-level exploratory; sensitivity agreement is not independent biological validation'})

def subsample(core):
    a,units=units_and_obs();rows=[]
    for uid,g in core.groupby('unit_id'):
        r=units[units.unit_id.eq(uid)].iloc[0];sub=a[mask_unit(a.obs,r)];genes=g.gene.tolist();x=dense(sub.raw[:,genes].X)
        # Match Scanpy fold-change convention: back-transform mean log-expression.
        for seed in [11,23,37,53,71]:
            ids,_=matched_indices(sub.obs,seed);co=sub.obs.condition.astype(str).to_numpy()[ids];y=x[ids];c=np.mean(y[co=='Ctrl'],axis=0);t=np.mean(y[co=='T2'],axis=0);fc=np.log2((np.expm1(t)+1e-9)/(np.expm1(c)+1e-9))
            for gene,f in zip(genes,fc):rows.append({'unit_id':uid,'gene':gene,'seed':seed,'QC_matched_log2FC':f,'n_each':int((co=='Ctrl').sum())})
    z=pd.DataFrame(rows);csv(z,'candidate_cell_subsampling.csv');summ=z.groupby(['unit_id','gene']).QC_matched_log2FC.agg(['min','max','median']).reset_index();summ['same_sign_all_5']=summ['min']*summ['max']>0;csv(summ,'candidate_cell_subsampling_summary.csv')

def report():
    core=pd.read_csv(OUT/'tables/gene_candidates_supported.csv');local=pd.read_csv(OUT/'tables/gene_candidates_local_non_dominant.csv');qc=pd.read_csv(OUT/'tables/QC_matching_audit.csv');pp=pd.read_csv(OUT/'tables/pathway_candidates_all.csv');freq=pd.read_csv(OUT/'tables/population_frequencies.csv');sub=pd.read_csv(OUT/'tables/candidate_cell_subsampling_summary.csv')
    # Rank by weakest supported effect, keeping every passing result in the underlying table.
    local=local.merge(sub,on=['unit_id','gene'],how='left').sort_values(['same_sign_all_5','minimum_abs_log2FC','original_FDR'],ascending=[False,False,True]);csv(local,'discovery_shortlist.csv')
    top=local.head(15).copy();labels=[f'{r.gene} | {r.cluster} / {r.tissue}' for r in top.itertuples()];values=top[['original_log2FC','stable_log2FC','QC_log2FC','stable_QC_log2FC']].to_numpy()
    fig,ax=plt.subplots(figsize=(11,max(5,len(top)*.36)));im=ax.imshow(values,cmap='RdBu_r',vmin=-1.5,vmax=1.5,aspect='auto');ax.set_yticks(range(len(labels)),labels,fontsize=9);ax.set_xticks(range(4),['Original cluster','Tissue recluster','QC matched','Recluster + QC'],rotation=20,ha='right');fig.colorbar(im,ax=ax,label='log2FC T2 vs Ctrl');ax.set_title('Data-selected local candidates: agreement across sensitivities\nRanking uses all genes; no paper/PPT target restriction');fig.tight_layout();fig.savefig(OUT/'figures/discovery_shortlist.png',dpi=180);fig.savefig(OUT/'figures/discovery_shortlist.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,5));x=np.arange(len(qc));before=qc[[c for c in qc if c.startswith('before_')]].abs().max(axis=1);after=qc[[c for c in qc if c.startswith('after_')]].abs().max(axis=1);ax.plot(x,before,'o-',label='Before');ax.plot(x,after,'o-',label='After');ax.axhline(.2,color='gray',ls='--');ax.set_xticks(x,qc.unit_id.str.replace('cluster_tissue__','').str.replace('stable__',''),rotation=90,fontsize=6);ax.set_ylabel('Maximum absolute QC standardized mean difference');ax.legend();fig.tight_layout();fig.savefig(OUT/'figures/QC_balance.png',dpi=180);plt.close(fig)
    lines=['# iNKT：沿用最新pipeline寻找本数据中的候选发现','','本轮目标是从完整基因与功能库中筛选本数据里的变化。论文和旧PPT仅在筛选后提供已有结果标记；没有用其基因或通路清单限制候选。','','## 实际新增分析','','1. 复用最新原分群和三个组织独立重聚类的全部DE结果。','2. 在每个对照内按log UMI、log检测基因数、线粒体比例各四分位联合分层，在每层平衡Ctrl/T2细胞数，重新计算全基因DE。','3. 对所有可用对照做完整KEGG/Reactome ORA；保留去heat-shock敏感性和跨对照BH结果。','4. 用细胞交集连接旧分群与组织重聚类，检查效应方向与大小；对入选候选另外做5次细胞匹配抽样。','5. 以效应和检测率、重聚类一致性、QC敏感性筛选；用驱动基因合并重复的通路命名。','','## 优先候选（完整筛选规则见manifest）','','| 基因 | 组织/原群 | 原始log2FC / FDR | 重聚类log2FC | QC匹配log2FC / FDR | 组织整体log2FC | 5次抽样同号 |','|---|---|---|---:|---|---:|---|']
    for r in top.itertuples():lines.append(f'| {r.gene} | {r.tissue}/{r.cluster} | {r.original_log2FC:+.3f} / {r.original_FDR:.3g} | {r.stable_log2FC:+.3f} | {r.QC_log2FC:+.3f} / {r.QC_FDR:.3g} | {r.tissue_log2FC:+.3f} | {r.same_sign_all_5} |')
    lines+=['','这些是数据集内值得追踪的候选，不能称为已证明的新机制。筛选后的“局部效应”是原群与组织整体的效应对照，不是正式交互作用检验。','','## 输出','','- [候选图](figures/discovery_shortlist.png)','- [完整候选证据](tables/gene_candidates_all_sensitivities.csv)','- [优先候选](tables/discovery_shortlist.csv)','- [完整功能检验](tables/functional_screen_all.csv.gz)','- [功能候选的重聚类/QC对照](tables/pathway_candidates_all.csv)','- [重复驱动基因分组](tables/pathway_exact_driver_groups.csv)','- [细胞群组成](tables/population_frequencies.csv)','- [QC平衡情况](tables/QC_matching_audit.csv)','- [5次抽样结果](tables/candidate_cell_subsampling_summary.csv)','','## 解释边界','','当前每组织×条件只有一个sample标签；没有新增动物重复。匹配和抽样的一致性仅衡量同一批细胞的敏感性，不能提供动物总体置信区间。QC协变量也可能含真实生物变化，所以原始结果与匹配结果都保留。已有基因过滤保持不变，不能发现被过滤掉的低丰度基因。轨迹方向/velocity未在本轮声称完成。','']
    (OUT/'README.md').write_text('\n'.join(lines));js(OUT/'manifest.json',{'source_H5AD':str(INPUT),'source_SHA256':sha256_file(INPUT),'script_SHA256':sha256_file(Path(__file__)),'workflow':'latest iNKT analysis pipeline; open discovery plus clearly labelled QC/subsampling sensitivities','libraries':{k:len(v) for k,v in libraries().items() if k in ['KEGG_Mouse_2019','Reactome_2022_ortholog']},'selection':{'initial_gene_FDR':.05,'initial_abs_log2FC':.25,'max_detection_fraction_min':.1,'stable_FDR':.05,'QC_FDR':.1,'min_cells_per_condition':50,'minimum_original_overlap_in_stable':.5,'maximum_QC_SMD':.2,'minimum_abs_effect_each_sensitivity':.25,'all_effect_signs_consistent':True,'local_effect_difference_vs_whole_tissue':.25,'excluded_from_local_display':'Hsp/Dnaj, ribosomal, respiratory-chain/mt families; all retained in full tables','rank':'all-sign5 then minimum absolute effect across sensitivities then original FDR','paper_lists_used_for_selection':False},'counts':{'gene_records_supported':len(core),'local_other_gene_records':len(local),'pathway_records':len(pp)},'GPU_used':False,'inference':'exploratory, no independent animal replicates'})
    js(OUT/'logs/report.completed.json',{'status':'complete','shortlist_records':len(local)});log('DONE '+str(OUT))

if __name__=='__main__':
    setup();p=argparse.ArgumentParser();p.add_argument('stage',choices=['de','enrich','select','report','all']);args=p.parse_args()
    for stage,fun in [('de',de_and_crosswalk),('enrich',enrich),('select',select_candidates),('report',report)]:
        if args.stage in [stage,'all']:log('START '+stage);fun();log('COMPLETE '+stage)
