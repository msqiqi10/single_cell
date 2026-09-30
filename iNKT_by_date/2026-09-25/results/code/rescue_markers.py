"""Read full-feature raw matrices for retained cells; preserve primary analysis."""
from pathlib import Path
import gzip,json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/bio3'
OUT=ROOT/'runs/meeting_followup';(OUT/'objects').mkdir(exist_ok=True)
base=ad.read_h5ad(DATA/'output/iNKT_meeting_followup_20260905/objects/scored_base.h5ad',backed='r')
obs=base.obs.copy();primary=set(base.var_names);base.file.close()
cluster='cluster_c5_split_20260830'
genes='Dusp1 Fos Jun Junb Jund Fosb Fosl1 Fosl2 Nr4a1 Il1r1 Il6ra Cd8a Ciita Rorc Il23r Ccr6 Il17a Il17f Tbx21 Gata3 Zbtb16 Nkg7 Prf1 Gzma Gzmb Ccl5'.split()
blocks=[];raw_audit=[]
for path in sorted((DATA/'input/iNKT/data').glob('*/sample_feature_bc_matrix')):
    sample=path.parent.name
    f=pd.read_csv(path/'features.tsv.gz',sep='\t',header=None)
    bc=pd.read_csv(path/'barcodes.tsv.gz',header=None)[0].astype(str)
    names=pd.Index(sample+'_'+bc)
    target=obs.index[obs['sample'].astype(str)==sample]
    positions=names.get_indexer(target);assert (positions>=0).all()
    with gzip.open(path/'matrix.mtx.gz','rb') as stream: matrix=mmread(stream).tocsr()
    feature_idx=np.flatnonzero(f[2].eq('Gene Expression').to_numpy())
    var=pd.DataFrame({'gene_ids':f.iloc[feature_idx,0].to_numpy()},index=f.iloc[feature_idx,1].astype(str).to_numpy())
    for gene in genes:
        ii=np.flatnonzero(f[1].eq(gene).to_numpy() & f[2].eq('Gene Expression').to_numpy())
        raw_audit.append(dict(sample=sample,gene=gene,n_matching_features=len(ii),
                              n_input_barcodes=len(bc),n_retained_cells=len(target),
                              input_cells_detected=int((matrix[ii].sum(axis=0)>0).sum()) if len(ii) else 0,
                              retained_cells_detected=int((matrix[ii][:,positions].sum(axis=0)>0).sum()) if len(ii) else 0))
    x=ad.AnnData(matrix[feature_idx][:,positions].T.tocsr(),obs=obs.loc[target].copy(),var=var)
    x.var_names_make_unique();blocks.append(x);print(sample,x.shape,flush=True)
a=ad.concat(blocks,join='inner',merge='same');a=a[obs.index].copy()
assert a.n_obs==15532 and a.obs_names.equals(obs.index)
assert (a.X.data>=0).all() and np.equal(a.X.data,np.floor(a.X.data)).all()
totals=np.asarray(a.X.sum(axis=1)).ravel();assert (totals>0).all()
a.obs['full_feature_total_counts']=totals
a.write_h5ad(OUT/'objects/retained_full_feature_counts.h5ad',compression='gzip')
check=ad.read_h5ad(DATA/'output/iNKT_meeting_followup_20260905/objects/scored_base.h5ad')
difference=a[check.obs_names,check.var_names].X-check.layers['counts']
difference.eliminate_zeros();assert difference.nnz==0
rows=[];global_rows=[]
for gene in genes:
    if gene not in a.var_names:
        global_rows.append(dict(gene=gene,status='not_in_input_feature_list',in_primary=gene in primary,n_detected=0,total_UMI=0));continue
    counts=a[:,[gene]].X.toarray().ravel()
    log=np.log1p(counts/totals*10000)
    global_rows.append(dict(gene=gene,status='detected' if (counts>0).any() else 'listed_but_zero_in_retained_cells',
                            in_primary=gene in primary,n_detected=int((counts>0).sum()),total_UMI=int(counts.sum())))
    frame=obs[['sample','condition','tissue',cluster]].copy();frame['count']=counts;frame['log1p_full']=log;frame['positive']=counts>0
    for key,part in frame.groupby(['sample','condition','tissue',cluster],observed=True):
        sample,condition,tissue,c=key
        rows.append(dict(gene=gene,sample=sample,condition=condition,tissue=tissue,cluster=c,n_cells=len(part),
                         n_detected=int(part.positive.sum()),pct_detected=100*part.positive.mean(),
                         total_UMI=int(part['count'].sum()),mean_UMI=part['count'].mean(),mean_log1p=part.log1p_full.mean(),in_primary=gene in primary))
pd.DataFrame(rows).to_csv(OUT/'tables/full_feature_marker_by_group.csv',index=False)
pd.DataFrame(global_rows).to_csv(OUT/'tables/full_feature_marker_status.csv',index=False)
pd.DataFrame(raw_audit).to_csv(OUT/'tables/raw_and_retained_marker_detection.csv',index=False)
report=dict(passed=True,shape=list(a.shape),cell_alignment='exact retained cell IDs and original order',
            full_feature_counts_match_primary_count_layer_for_all_primary_genes=True,
            integer_counts=True,normalization_for_marker_summary='full-feature count total per retained cell, target 10000, log1p',
            status='descriptive coverage/detection audit, not new differential-expression tests or secretion assay',
            markers=global_rows)
(OUT/'marker_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
