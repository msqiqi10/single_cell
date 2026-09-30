"""Coverage rescue and stability-selected tissue clustering, preserving primary outputs."""
from run_inkt_meeting_followup import *

def raw_coverage():
    base=ad.read_h5ad(OUT/'objects/scored_base.h5ad');blocks=[];audit=[]
    for path in sorted((ROOT/'input/iNKT/data').glob('*/sample_feature_bc_matrix')):
        sample=path.parent.name;f=pd.read_csv(path/'features.tsv.gz',sep='\t',header=None);bc=pd.read_csv(path/'barcodes.tsv.gz',header=None)[0].astype(str)
        names=pd.Index(sample+'_'+bc);positions=names.get_indexer(base.obs_names[base.obs['sample'].astype(str).eq(sample)])
        if (positions<0).any():raise ValueError('Unmatched retained cell barcode')
        with gzip.open(path/'matrix.mtx.gz','rb') as stream:m=mmread(stream).tocsr()
        idx=np.flatnonzero(f[2].eq('Gene Expression').to_numpy());m=m[idx][:,positions].T.tocsr()
        var=pd.DataFrame({'gene_ids':f.iloc[idx,0].to_numpy()},index=f.iloc[idx,1].astype(str).to_numpy());x=ad.AnnData(m,obs=base.obs.loc[names[positions]].copy(),var=var);x.var_names_make_unique();blocks.append(x)
    a=ad.concat(blocks,join='inner',merge='same');a=a[base.obs_names].copy();a.obsm['X_umap']=base.obsm['X_umap'].copy()
    detected=np.asarray((a.X>0).sum(axis=0)).ravel();total=np.asarray(a.X.sum(axis=0)).ravel();h2m,m2raw,_=mapping(a)
    sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a);gene_pool=a.var_names[detected>=10].tolist();summ=[]
    for state,(species,text) in STATE_GENES.items():
        used=[]
        for gene in text.split():
            candidates=([h2m[gene]] if gene in h2m else []) if species=='human' else ([gene] if gene in a.var_names else sorted(m2raw.get(gene,[])))
            g=candidates[0] if len(candidates)==1 else None;j=a.var_names.get_loc(g) if g else None;n=int(detected[j]) if j is not None else 0
            if g and n>=10:used.append(g)
            audit.append({'signature':state,'source_gene':gene,'measured_mouse_gene':g,'retained_cells_detected':n,'retained_total_UMI':int(total[j]) if j is not None else np.nan,'in_primary_filtered_genes':bool(g and g in base.var_names),'status':'used' if g and n>=10 else ('too_sparse_or_zero' if g else 'ambiguous_or_unmapped')})
        used=sorted(set(used));coverage=len(used)/len(text.split());ok=len(used)>=3 and coverage>=.4
        a.obs[state]=np.nan
        if ok:sc.tl.score_genes(a,used,score_name=state,ctrl_size=50,n_bins=25,random_state=SEED,use_raw=False,gene_pool=gene_pool)
        orig=base.obs[state].to_numpy(float);new=a.obs[state].to_numpy(float);rho=float(spearmanr(orig,new).statistic) if np.isfinite(orig).all() and np.isfinite(new).all() else np.nan
        summ.append({'signature':state,'n_used':len(used),'coverage':coverage,'genes':';'.join(used),'status':'scored_full_feature_sensitivity' if ok else 'insufficient_marker_coverage','spearman_to_primary':rho})
    csv(pd.DataFrame(audit),'full_feature_marker_coverage_audit.csv');csv(pd.DataFrame(summ),'full_feature_score_coverage.csv');a.obs[['sample','condition','tissue',REFINED_CLUSTER_KEY]+list(STATE_GENES)].to_csv(OUT/'tables/full_feature_state_scores_per_cell.csv.gz',index_label='cell_id')
    plot_score_maps(a,'09_full_feature_NK_state_sensitivity',FIVE_STATES)
    js(OUT/'logs/full_feature_scores.completed.json',{'n_cells':a.n_obs,'n_features':a.n_vars,'eligible_gene_pool':len(gene_pool),'min_detected_cells':10,'raw_filter':'Gene Expression only; same retained cells; normalize full-feature totals to10000 thenlog1p','scope':'marker coverage sensitivity, does not replace primary expression or DE matrices'})

def stable_clusters():
    stab=pd.read_csv(OUT/'tables/within_tissue_cluster_stability.csv');selected=[];status=[]
    for tissue in TISSUES:
        # Explicit data-dependent sensitivity, based on seed agreement, not condition results.
        s=stab[stab.tissue.eq(tissue)].sort_values(['minimum_ARI','resolution'],ascending=[False,True]).iloc[0];res=float(s.resolution)
        a=ad.read_h5ad(OUT/f'objects/{tissue}_reclustered.h5ad');a.obs['meeting_cluster_res05']=a.obs.meeting_cluster.copy()
        sc.tl.leiden(a,resolution=res,random_state=0,key_added='stable_cluster',flavor='igraph',n_iterations=2,directed=False)
        a.obs['meeting_cluster']=pd.Categorical([f'{tissue}_S{x}' for x in a.obs.stable_cluster.astype(str)])
        a.uns['meeting_stability_selection']={'resolution':res,'seed':0,'rule':'highest minimum pairwise seed ARI among0.3,0.5,0.8; post-audit sensitivity','minimum_ARI':float(s.minimum_ARI)}
        sc.tl.rank_genes_groups(a,'meeting_cluster',method='wilcoxon',use_raw=True,n_genes=30,tie_correct=True);sc.get.rank_genes_groups_df(a,group=None).to_csv(OUT/f'tables/{tissue}_stable_top30_markers.csv',index=False)
        plot_embedding(a,'10_stable_'+tissue,f'{tissue}: stability-selected resolution {res}; min seed ARI {s.minimum_ARI:.2f}')
        a.write_h5ad(OUT/f'objects/{tissue}_stable.h5ad',compression='gzip')
        pd.crosstab(a.obs[REFINED_CLUSTER_KEY],a.obs.meeting_cluster).to_csv(OUT/f'tables/{tissue}_stable_old_new_crosswalk.csv')
        for c in a.obs.meeting_cluster.cat.categories:
            b=a[a.obs.meeting_cluster.eq(c)].copy();ns=b.obs.condition.astype(str).value_counts();nt=int(ns.get('T2',0));nc=int(ns.get('Ctrl',0));uid='stable__'+c
            row={'unit_id':uid,'tissue':tissue,'cluster':c,'resolution':res,'n_tumor':nt,'n_control':nc,'strict_min20':min(nt,nc)>=20}
            if min(nt,nc)>=2:
                b.obs[REFINED_CLUSTER_KEY]=pd.Categorical([c]*len(b));d=run_de_unit(b,DEUnit(uid,'stable_recluster',c,None));d['tissue']=tissue;d.to_csv(OUT/f'de/{uid}.csv.gz',index=False);row.update(n_up=int(d.robust_tumor_up.sum()),n_down=int(d.robust_control_up.sum()),status='computed')
            else:row['status']='insufficient_cells'
            status.append(row)
        selected.append(s.to_dict());csv(pd.DataFrame(status),'de_status_stable.csv')
    csv(pd.DataFrame(selected),'selected_stability_resolution.csv');js(OUT/'logs/refinement.completed.json',{'stable_tissues':3,'clusters':len(status)})

if __name__=='__main__':
    raw_coverage();stable_clusters()
