"""Fixed technical sensitivities, without replacing the main analysis or FDR."""
from pathlib import Path
import json,re
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import hypergeom,false_discovery_control
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/bio3';OUT=ROOT/'runs/meeting_followup';SRC=DATA/'iNKT_by_date/2026-09-19'
catalog=pd.read_csv(SRC/'results/tables/GO_term_catalog.csv').set_index('go_id');sets={}
for line in (SRC/'sources/GO_mouse_measured.gmt').read_text().splitlines():
    gid,name,*genes=line.split('\t');sets[gid]=set(genes)
de=pd.read_csv(OUT/'tables/all_20_condition_DE.csv.gz');go=pd.read_csv(OUT/'tables/all_original_condition_GO.csv.gz',keep_default_na=False)
ap1=go[go.go_id=='GO:0035976'];ap1.to_csv(OUT/'tables/AP1_complex_GO_exact_ID.csv',index=False)
# Export non-significant, tested annotations as well as significant coverage.
reverse={}
for gid,genes in sets.items():
    for gene in genes:reverse.setdefault(gene,[]).append(gid)
annotations=[]
for uid,d in de.groupby('unit_id'):
    for row in d[(d.pvals_adj<=.05)&(d.logfoldchanges.abs()>=.25)].itertuples():
        entries=reverse.get(row.gene,[])
        annotations.append(dict(unit_id=uid,gene=row.gene,direction='T2_up' if row.logfoldchanges>0 else 'Ctrl_up',
                                n_GO_annotations=len(entries),all_GO_ids=';'.join(sorted(entries)),
                                annotation_names=' | '.join(catalog.loc[g,'term'] for g in sorted(entries)),
                                interpretation='annotation only, not evidence of enrichment or activation'))
pd.DataFrame(annotations).to_csv(OUT/'tables/all_candidate_GO_annotations.csv.gz',index=False)
removed=sorted(g for g in de.gene.unique() if re.match(r'^(Rpl|Rps|Hsp|Dnaj)',g))
(OUT/'gene_lists/sensitivity_excluded_ribo_chaperone.txt').write_text('\n'.join(removed)+'\n')
results=[]
for uid in ['cluster_tissue__C0__bone_marrow','cluster_tissue__C3__spleen','cluster_tissue__C6__thymus']:
    d=de[de.unit_id==uid];expressed=set(d.loc[(d.pct_expressing_tumor>0)|(d.pct_expressing_control>0),'gene'])-set(removed)
    selected=set(d.loc[(d.pvals_adj<=.05)&(d.logfoldchanges>=.25),'gene'])-set(removed)
    for ns in ['BP','MF','CC']:
        lib={gid:sets[gid] for gid in catalog.index[catalog.namespace==ns]}
        universe=expressed&set().union(*lib.values());query=selected&universe;rows=[]
        for gid,members in lib.items():
            measured=members&universe
            if 5<=len(measured)<=500:
                hits=measured&query
                rows.append(dict(go_id=gid,term=catalog.loc[gid,'term'],M=len(universe),K=len(measured),n=len(query),k=len(hits),genes=';'.join(sorted(hits))))
        r=pd.DataFrame(rows);r['pvalue']=hypergeom.sf(r.k-1,len(universe),r.K,len(query)) if query else 1.
        r['q_sensitivity_family']=false_discovery_control(r.pvalue.to_numpy(),method='bh');r['unit_id']=uid;r['namespace']=ns
        results.append(r)
result=pd.concat(results,ignore_index=True);result.to_csv(OUT/'tables/GO_exclude_ribo_chaperone_sensitivity.csv.gz',index=False)
# Matched cell draws are technical sensitivity, not biological replicates or confidence intervals.
a=ad.read_h5ad(DATA/'output/iNKT_meeting_followup_20260905/objects/scored_base.h5ad')
panel='Fos Jun Junb Jund Fosb Fosl2 Dusp1 Nr4a1 Rorc Il23r'.split();raw=a.raw[:,panel].X.toarray()
qc=['total_counts','n_genes_by_counts','pct_counts_mt'];clusterkey='cluster_c5_split_20260830';runs=[]
for uid,d in de.groupby('unit_id'):
    tissue=d.tissue.iloc[0];cluster=d.cluster.iloc[0]
    mask=a.obs.tissue.astype(str).eq(tissue).to_numpy()
    if cluster!='ALL':mask &= a.obs[clusterkey].astype(str).eq(cluster).to_numpy()
    ix=np.flatnonzero(mask);obs=a.obs.iloc[ix];condition=obs.condition.astype(str).to_numpy()
    bins=np.column_stack([pd.qcut(obs[col].astype(float),4,labels=False,duplicates='drop').fillna(0).to_numpy(int) for col in qc])
    strata=np.ravel_multi_index(bins.T,(4,4,4));full_delta=raw[ix[condition=='T2']].mean(axis=0)-raw[ix[condition=='Ctrl']].mean(axis=0)
    for seed in range(20):
        rng=np.random.default_rng(seed);ctrl=[];t2=[]
        for group in np.unique(strata):
            c=np.flatnonzero((strata==group)&(condition=='Ctrl'));t=np.flatnonzero((strata==group)&(condition=='T2'));n=min(len(c),len(t))
            if n:ctrl.extend(rng.choice(c,n,replace=False));t2.extend(rng.choice(t,n,replace=False))
        ctrl=np.asarray(ctrl,int);t2=np.asarray(t2,int)
        if len(ctrl)==0:continue
        effect=raw[ix[t2]].mean(axis=0)-raw[ix[ctrl]].mean(axis=0);smd=[]
        for col in qc:
            c=obs.iloc[ctrl][col].to_numpy(float);t=obs.iloc[t2][col].to_numpy(float)
            pooled=np.sqrt((c.var()+t.var())/2);smd.append(abs(t.mean()-c.mean())/pooled if pooled>0 else 0.)
        for j,gene in enumerate(panel):
            runs.append(dict(unit_id=uid,gene=gene,seed=seed,n_each=len(ctrl),delta_mean_log1p=float(effect[j]),
                             primary_delta_mean_log1p=float(full_delta[j]),max_abs_QC_SMD=max(smd),same_direction=bool(np.sign(effect[j])==np.sign(full_delta[j]))))
pd.DataFrame(runs).to_csv(OUT/'tables/Rob_panel_QC_matched_runs.csv',index=False)
summary=pd.DataFrame(runs).groupby(['unit_id','gene']).agg(n_draws=('seed','size'),n_each=('n_each','min'),
    primary_delta=('primary_delta_mean_log1p','first'),min_delta=('delta_mean_log1p','min'),max_delta=('delta_mean_log1p','max'),
    same_direction_draws=('same_direction','sum'),max_QC_SMD=('max_abs_QC_SMD','max')).reset_index()
summary.to_csv(OUT/'tables/Rob_panel_QC_matched_summary.csv',index=False)
manifest=dict(passed=True,GO_comparisons=3,GO_tests=len(result),removed_genes=len(removed),
    GO_rule='T2-up robust list; remove ^(Rpl|Rps|Hsp|Dnaj) from query AND expressed background; rebuild term size5..500 and BH by unit/namespace',
    QC_rule='20 fixed seeds0..19, equal cell counts within combined quartile bins of total counts, detected genes, mitochondrial percentage',
    QC_effect='difference in mean primary log1p expression, not log2FC; draw range not biological CI; no new inferential P values',
    note='Sensitivity cannot separate true environmental stress from dissociation or batch; main analysis unchanged')
(OUT/'sensitivity_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
