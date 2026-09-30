#!/usr/bin/env python3
"""Execute the 2026-09-05 meeting follow-up, with resumable CPU-only stages."""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
import argparse, gzip, hashlib, json, re, sys, time
from pathlib import Path
from itertools import combinations
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from scipy.stats import hypergeom, spearmanr
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from inkt_palette import SAMPLE_PALETTE, TISSUE_PALETTE, CONDITION_ORDER
from run_inkt_c5_paper_followup import (
    REFINED_CLUSTER_KEY, REFINED_ORDER, REFINED_COLORS, DEUnit, run_de_unit,
    bh_adjust, read_gmt, write_gmt, sha256_file, build_cluster_proportion_tables,
    render_denominator_reconciliation, build_marker_complementarity,
    marker_complementarity_row, expression_and_counts,
)
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/iNKT_meeting_followup_20260905'
BASE=ROOT/'output/iNKT_reproduction_deck/20260830_C5_paper_Fig3DEF_followup'
INPUT=BASE/'20260830_inkt_selected_umap_C5_split.h5ad'
OLD=ROOT/'output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt'
SEED=20260905
TISSUES=('bone_marrow','spleen','thymus')
STATE_GENES={
 'NK_resting':('human','KLRC1 NCAM1 GZMK GZMA KLRB1'),
 'NK_adaptive':('human','KLRC2 GZMH LAG3 HLA-DRA'),
 'NK_activated':('human','TNFRSF18 TNFRSF9 TNFRSF4 CRTAM ENTPD1 HAVCR2 TIGIT TNFSF10 BCL2L11'),
 'NK_type_I_IFN':('human','ISG15 OAS3 MX1 IRF7 MX2 IRF9 OAS1 EIF2AK2 OAS2 LAG3'),
 'NK_cytokine':('human','CCL3 CCL4 TNF IFNG CD69'),
 'NK_early_activated':('human','GZMB XCL1 XCL2 TNFRSF9 CRTAM'),
 'iNKT1':('mouse','Il2rb Tbx21 Ifng'),
 'iNKT2':('mouse','Zbtb16 Gata3 Cd4 Il4 Il17rb'),
 'iNKT17':('mouse','Rorc Ccr6 Zbtb16 Il23r'),
 'iNKT_precursor':('mouse','Cd24a Egr2 Hivep3'),
 'Cycling':('mouse','Mki67 Top2a Cks1b'),
 'NK_maturation':('mouse','Itgam S1pr5 Gzma Zeb2 Klf2'),
 'NK_immaturity':('mouse','Cd27 Kit Il7r'),
 'TNF_SOCS_hypothesis':('mouse','Tnf Tnfrsf1b Cish Socs1 Socs2 Il2ra Il2rb Il7r'),
 'Cytotoxic_effector':('mouse','Prf1 Gzmb Gzma Nkg7'),
 'Heat_shock':('mouse','Hspa8 Hspa1a Hspa1b Hsp90aa1 Hsp90ab1 Dnaja1 Hsph1'),
}
FIVE_STATES=list(STATE_GENES)[:5]
PAPER_TERMS=['Hematopoietic cell lineage','Antigen processing and presentation','Type I diabetes mellitus','Cytokine-cytokine receptor interaction','Toxoplasmosis','Rheumatoid arthritis','Graft-versus-host disease','Asthma','Th17 cell differentiation','IL-17 signaling pathway','Autoimmune thyroid disease','NF-kappa B signaling pathway','Osteoclast differentiation','Cell adhesion molecules','Measles','Systemic lupus erythematosus','JAK-STAT signaling pathway']
PAPER_GENES='Il3ra Il2ra Il7r Tnf Ltb Il1b Ahr Il6st Tnfsf11 Bcl2a1d Socs1 Socs2 Cish Xcl1 Csf2 Tnfrsf9 Tnfsf8'.split()
WARN='Exploratory cell-level comparisons; one sample per tissue/condition. Tissue/sample Harmony are integration sensitivities, not identified technical batch effects.'

def log(msg):print(time.strftime('%Y-%m-%d %H:%M:%S'),msg,flush=True)
def slug(x):return re.sub(r'[^a-zA-Z0-9_-]+','_',str(x)).strip('_')
def js(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)))
def csv(df,name):
    p=OUT/'tables'/name;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False)
def save(fig,name):
    p=OUT/'figures'/name;p.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(p.with_suffix('.png'),dpi=170,bbox_inches='tight')
    fig.savefig(p.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig)
def vec(x):return np.asarray(x.toarray() if sparse.issparse(x) else x).ravel()
def set_source():
    OUT.mkdir(exist_ok=True)
    for d in ['tables','figures','objects','logs','sources','de','enrichment']: (OUT/d).mkdir(exist_ok=True)

def mapping(adata):
    hom=pd.read_csv(OUT/'sources/MGI_HOM_MouseHumanSequence.rpt',sep='\t',dtype=str).fillna('')
    ens=pd.read_csv(OUT/'sources/MGI_MRK_ENSEMBL.rpt',sep='\t',header=None,dtype=str).fillna('')
    ens_to_symbol=dict(zip(ens[5],ens[1]));current_to_raw={}
    for raw,eid in zip(adata.var_names,adata.var.gene_ids.astype(str)):
        current=ens_to_symbol.get(eid,str(raw));current_to_raw.setdefault(current,set()).add(str(raw))
    h2m={};mapping_rows=[]
    for _,g in hom.groupby('DB Class Key',sort=False):
        hs=set(g.loc[g['NCBI Taxon ID']=='9606','Symbol']);ms=set(g.loc[g['NCBI Taxon ID']=='10090','Symbol'])
        for h in hs:
            candidates=sorted({r for m in ms for r in current_to_raw.get(m,set())})
            status='one_to_one' if len(hs)==len(ms)==1 else 'ambiguous_orthology'
            if status=='one_to_one' and len(candidates)==1:h2m[h]=candidates[0]
            mapping_rows.append({'human':h,'mouse_symbols':';'.join(sorted(ms)),'measured_mouse_genes':';'.join(candidates),'orthology':status,'used':h in h2m})
    return h2m,current_to_raw,pd.DataFrame(mapping_rows)

def project_library(sets,h2m):
    return {name:set(h2m[g] for g in members if g in h2m) for name,members in sets.items()}

def audit_raw(targets):
    rows=[]
    for path in sorted((ROOT/'input/iNKT/data').glob('*/sample_feature_bc_matrix')):
        features=pd.read_csv(path/'features.tsv.gz',sep='\t',header=None)
        ids=np.flatnonzero(features[1].isin(targets).to_numpy())
        log('Raw marker count audit '+path.parent.name)
        with gzip.open(path/'matrix.mtx.gz','rb') as f: m=mmread(f).tocoo()
        totals=np.bincount(m.row,weights=m.data,minlength=m.shape[0])
        detected=np.bincount(m.row[m.data>0],minlength=m.shape[0])
        for i in ids:rows.append({'sample':path.parent.name,'gene':features.iloc[i,1],'total_counts':int(totals[i]),'detected_cells':int(detected[i]),'n_raw_cells':m.shape[1]})
    return pd.DataFrame(rows)

def score_and_audit():
    a=ad.read_h5ad(INPUT);h2m,m2raw,orth=mapping(a);csv(orth,'human_mouse_orthology_audit.csv')
    js(OUT/'sources/human_to_measured_mouse.json',h2m)
    library_audit=[]
    for name in ['Reactome_2022_Human','KEGG_2021_Human']:
        human=read_gmt(OUT/'sources'/f'{name}.gmt');mapped=project_library(human,h2m)
        write_gmt(OUT/'sources'/f'{name}_mapped_mouse.gmt',mapped)
        for term,genes in human.items():library_audit.append({'library':name,'term':term,'n_source_genes':len(genes),'n_mapped_measured':len(mapped[term]),'coverage':len(mapped[term])/max(1,len(genes))})
    csv(pd.DataFrame(library_audit),'library_gene_coverage.csv')
    gene_rows=[];summary=[];used_sets={}
    for state,(species,text) in STATE_GENES.items():
        mapped=[]
        for gene in text.split():
            candidates=([h2m[gene]] if gene in h2m else []) if species=='human' else ([gene] if gene in a.var_names else sorted(m2raw.get(gene,[])))
            measured=candidates[0] if len(candidates)==1 else None
            gene_rows.append({'signature':state,'species':species,'source_gene':gene,'measured_gene':measured,'status':'used' if measured else 'unmapped_ambiguous_or_filtered','source':'Dufva2023_Fig1B_G' if state.startswith('NK_') and state not in ['NK_maturation','NK_immaturity'] else ('Krovi2020_Fig1' if state.startswith('iNKT') else 'Blood2025_or_documented_control_panel')})
            if measured:mapped.append(measured)
        mapped=sorted(set(mapped));coverage=len(mapped)/len(text.split());can_score=len(mapped)>=3 and coverage>=.4
        used_sets[state]=mapped
        if can_score:
            sc.tl.score_genes(a,mapped,score_name=state,ctrl_size=50,n_bins=25,random_state=SEED,use_raw=True)
        else:a.obs[state]=np.nan
        summary.append({'signature':state,'n_original':len(text.split()),'n_used':len(mapped),'coverage':coverage,'score_status':'scored_continuous_program' if can_score else 'insufficient_marker_coverage','genes':';'.join(mapped)})
    csv(pd.DataFrame(gene_rows),'state_marker_gene_audit.csv');csv(pd.DataFrame(summary),'state_signature_coverage.csv');js(OUT/'sources/state_gene_sets_measured.json',used_sets)
    score_cols=list(STATE_GENES)
    csv(a.obs.reset_index().rename(columns={'index':'cell_id'})[['cell_id','sample','tissue','condition',REFINED_CLUSTER_KEY]+score_cols],'state_scores_per_cell.csv.gz')
    summary_frames=[]
    for keys in [['tissue','condition',REFINED_CLUSTER_KEY],['tissue','condition'],[REFINED_CLUSTER_KEY]]:
        g=a.obs.groupby(keys,observed=True)[score_cols].agg(['mean','median','count'])
        g.columns=['__'.join(x) for x in g.columns];g=g.reset_index();g['scope']='_'.join(keys);summary_frames.append(g)
    csv(pd.concat(summary_frames,ignore_index=True),'state_score_group_summary.csv')
    targets=set(PAPER_GENES+['Il4','Klrd1','Tmsb10','Cd52','Hspa8','Dnaja1','Dynll1','Mxd4','Jun','Fos','Ncam1','Itgam','Il17a','Il17f','Mx1'])
    targets.update(g for v in used_sets.values() for g in v)
    raw=audit_raw(targets);csv(raw,'raw_reference_gene_audit.csv')
    expr_rows=[]
    for gene in sorted(targets & set(a.var_names)):
        expr,cnt=expression_and_counts(a,gene)
        df=a.obs[['sample','tissue','condition',REFINED_CLUSTER_KEY]].copy();df['expression']=expr;df['detected']=cnt>0
        g=df.groupby(['sample','tissue','condition',REFINED_CLUSTER_KEY],observed=True).agg(n=('detected','size'),detection_fraction=('detected','mean'),mean_log1p=('expression','mean')).reset_index();g['gene']=gene;expr_rows.append(g)
    csv(pd.concat(expr_rows,ignore_index=True),'reference_gene_expression_by_group.csv')
    a.write_h5ad(OUT/'objects/scored_base.h5ad',compression='gzip')
    plot_score_maps(a,'01_reference_state_maps',FIVE_STATES)
    plot_score_maps(a,'02_iNKT_and_function_maps',['iNKT1','iNKT2','iNKT17','iNKT_precursor','Cycling','Cytotoxic_effector'])
    plot_score_group(a,'03_state_scores_by_tissue_condition')
    sources={str(p.relative_to(ROOT)):sha256_file(p) for p in (OUT/'sources').glob('*') if p.is_file()}
    js(OUT/'sources/provenance.json',{'input':str(INPUT),'input_sha256':sha256_file(INPUT),'sources_sha256':sources,'seed':SEED,'scoring':'Scanpy score_genes; expression-matched control genes; >=3 genes and >=40% coverage','urls':{'orthology':'https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt','mouse_symbols':'https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt','Reactome':'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022','Dufva':'https://www.sciencedirect.com/science/article/pii/S1074761323004909','Krovi':'https://www.nature.com/articles/s41467-020-20073-8','Blood':'https://doi.org/10.1182/bloodadvances.2024014592'}})
    js(OUT/'tables/input_missing_tasks.json',{'velocity':{'status':'blocked_missing_spliced_unspliced','layers':list(a.layers),'searched':str(ROOT/'input/iNKT'),'required':'matching spliced/unspliced counts or BAM/FASTQ plus reference annotations','existing_dpt_is_velocity':False},'toxicity':{'status':'blocked_missing_PMID_inventory','required':'emailed liver/kidney toxicity PMID list or local path'},'sample_metadata':{'condition_samples':pd.crosstab(a.obs['sample'],a.obs.condition).to_dict(),'technical_batch_column_present':False,'animal_replicate_column_present':False}})

def plot_score_maps(a,name,states):
    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained');xy=a.obsm['X_umap']
    for ax,state in zip(axes.flat,states):
        v=a.obs[state].to_numpy(float);ax.scatter(xy[:,0],xy[:,1],s=1,c='#dddddd',rasterized=True)
        if np.isfinite(v).any():
            lo,hi=np.nanquantile(v,[.02,.98]);im=ax.scatter(xy[:,0],xy[:,1],c=v,s=2,cmap='viridis',vmin=lo,vmax=hi,rasterized=True);fig.colorbar(im,ax=ax,shrink=.65)
        else:ax.text(.5,.5,'Insufficient marker coverage',transform=ax.transAxes,ha='center')
        for c in a.obs[REFINED_CLUSTER_KEY].astype(str).unique():
            z=xy[a.obs[REFINED_CLUSTER_KEY].astype(str).eq(c)];mid=np.median(z,axis=0);ax.text(*mid,c,fontsize=7,weight='bold')
        ax.set_title(state.replace('_',' '));ax.set_xticks([]);ax.set_yticks([])
    for ax in list(axes.flat)[len(states):]:ax.axis('off')
    fig.suptitle('Reference programs on current clusters | continuous scores; separate color scales',fontsize=14);save(fig,name)

def plot_score_group(a,name):
    keys=['tissue','condition'];g=a.obs.groupby(keys,observed=True)[list(STATE_GENES)].mean();v=g.T
    fig,ax=plt.subplots(figsize=(11,8));im=ax.imshow(v,cmap='coolwarm',aspect='auto');fig.colorbar(im,ax=ax,label='Mean expression-control score (not a cell fraction)');ax.set_xticks(range(len(g)),[' / '.join(x) for x in g.index],rotation=30,ha='right');ax.set_yticks(range(len(v)),v.index);ax.set_title('Reference score means by tissue and condition');save(fig,name)

def graph_metrics(a,rep,labels,method,scope):
    x=np.asarray(a.obsm[rep]);nbr=NearestNeighbors(n_neighbors=21).fit(x);idx=nbr.kneighbors(x,return_distance=False)[:,1:]
    row={'scope':scope,'method':method,'n_cells':len(x),'n_clusters':len(set(labels))}
    for key in ['tissue','condition',REFINED_CLUSTER_KEY]:
        y=a.obs[key].astype(str).to_numpy();row[key+'_neighbor_purity']=float(np.mean(y[idx]==y[:,None]))
        row[key+'_silhouette']=float(silhouette_score(x,y,sample_size=min(1800,len(x)),random_state=SEED)) if 1<len(set(y))<len(y) else np.nan
    for state in FIVE_STATES+['iNKT1','iNKT2','iNKT17','Cycling']:
        v=a.obs[state].to_numpy(float)
        row[state+'_neighbor_spearman']=float(spearmanr(v,v[idx].mean(axis=1)).statistic) if np.isfinite(v).all() else np.nan
    row['ARI_to_original_refined']=float(adjusted_rand_score(a.obs[REFINED_CLUSTER_KEY].astype(str),labels))
    return row

def embedding_stage():
    import harmonypy as hm
    a=ad.read_h5ad(OUT/'objects/scored_base.h5ad');metrics=[];stability=[]
    a.obsm['X_pca_baseline']=np.asarray(a.obsm['X_pca'])[:,:50].copy()
    for method in ['original','harmony_tissue','harmony_sample']:
        b=a.copy()
        if method=='original':rep='X_pca_baseline'
        else:
            key='tissue' if method=='harmony_tissue' else 'sample';log('Harmony sensitivity '+key)
            z=hm.run_harmony(b.obsm['X_pca_baseline'],b.obs,key,theta=2,max_iter_harmony=20,random_state=SEED,verbose=True).Z_corr
            b.obsm['X_harmony']=z.T if z.shape[0]==50 else z;rep='X_harmony'
        log('Building graph '+method);sc.pp.neighbors(b,n_neighbors=20,use_rep=rep,random_state=SEED)
        sc.tl.umap(b,min_dist=.5,random_state=SEED)
        sc.tl.leiden(b,resolution=.5,random_state=SEED,key_added='meeting_cluster',flavor='igraph',n_iterations=2,directed=False)
        metrics.append(graph_metrics(b,rep,b.obs.meeting_cluster.astype(str).to_numpy(),method,'all_tissues'))
        b.uns['meeting_analysis']={'method':method,'batch_key':method.replace('harmony_','') if method!='original' else 'none','warning':WARN}
        b.write_h5ad(OUT/f'objects/{method}.h5ad',compression='gzip');plot_embedding(b,'04_'+method,method)
        pd.DataFrame(b.obsm['X_umap'],index=b.obs_names,columns=['UMAP1','UMAP2']).join(b.obs[['tissue','condition',REFINED_CLUSTER_KEY,'meeting_cluster']]).to_csv(OUT/f'tables/{method}_coordinates.csv.gz',index_label='cell_id')
    for tissue in TISSUES:
        b=a[a.obs.tissue.astype(str).eq(tissue)].copy();b.X=b.raw.X.copy();b.uns.pop('iroot',None)
        log('Within-tissue HVG/PCA/graph '+tissue)
        sc.pp.highly_variable_genes(b,n_top_genes=3000,flavor='seurat',batch_key=None)
        sc.tl.pca(b,n_comps=50,mask_var='highly_variable',svd_solver='arpack',random_state=SEED)
        sc.pp.neighbors(b,n_neighbors=20,n_pcs=50,random_state=SEED);sc.tl.umap(b,min_dist=.5,random_state=SEED)
        for res in [.3,.5,.8]:
            labels=[]
            for seed in [0,17,42]:
                sc.tl.leiden(b,resolution=res,random_state=seed,key_added='temp_cluster',flavor='igraph',n_iterations=2,directed=False);labels.append(b.obs.temp_cluster.astype(str).to_numpy().copy())
            aris=[adjusted_rand_score(x,y) for x,y in combinations(labels,2)]
            stability.append({'tissue':tissue,'resolution':res,'n_clusters_seed0':len(set(labels[0])),'mean_ARI':np.mean(aris),'minimum_ARI':min(aris)})
            if res==.5:b.obs['meeting_cluster']=pd.Categorical([f'{tissue}_L{x}' for x in labels[0]])
        del b.obs['temp_cluster'];metrics.append(graph_metrics(b,'X_pca',b.obs.meeting_cluster.astype(str).to_numpy(),'within_tissue_recluster',tissue))
        # Marker genes describe each new cluster; these are identity contrasts, distinct from condition DE.
        sc.tl.rank_genes_groups(b,'meeting_cluster',method='wilcoxon',use_raw=True,n_genes=30,tie_correct=True)
        mk=sc.get.rank_genes_groups_df(b,group=None);mk.to_csv(OUT/f'tables/{tissue}_recluster_top30_markers.csv',index=False)
        b.write_h5ad(OUT/f'objects/{tissue}_reclustered.h5ad',compression='gzip');plot_embedding(b,'05_'+tissue,tissue+' independently reclustered')
        g=b.obs.groupby('meeting_cluster',observed=True)[list(STATE_GENES)].mean();g.to_csv(OUT/f'tables/{tissue}_recluster_state_scores.csv')
        pd.crosstab(b.obs[REFINED_CLUSTER_KEY],b.obs.meeting_cluster).to_csv(OUT/f'tables/{tissue}_old_new_crosswalk_counts.csv')
    csv(pd.DataFrame(metrics),'embedding_metrics.csv');csv(pd.DataFrame(stability),'within_tissue_cluster_stability.csv')

def plot_embedding(a,name,title):
    fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained');xy=a.obsm['X_umap']
    for ax,key in zip(axes.flat,['tissue','condition','meeting_cluster','NK_activated','NK_type_I_IFN','NK_cytokine']):
        if key.startswith('NK_'):
            v=a.obs[key].to_numpy(float);lo,hi=np.nanquantile(v,[.02,.98]);im=ax.scatter(xy[:,0],xy[:,1],s=2,c=v,cmap='viridis',vmin=lo,vmax=hi,rasterized=True);fig.colorbar(im,ax=ax,shrink=.7)
        else:
            cats=sorted(a.obs[key].astype(str).unique());colors=plt.get_cmap('tab20').colors
            for i,c in enumerate(cats):
                m=a.obs[key].astype(str).eq(c);color=TISSUE_PALETTE.get(c,colors[i%20]) if key=='tissue' else ({'Ctrl':'#2878B5','T2':'#D9534F'}.get(c,colors[i%20]));ax.scatter(xy[m,0],xy[m,1],s=2,c=[color],label=c,rasterized=True)
            ax.legend(markerscale=4,fontsize=6,loc='best',frameon=False)
        ax.set_title(key);ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(title+' | independent embeddings; axes not directly comparable',fontsize=14);save(fig,name)

def centroid_summary(xy,condition,positive):
    row={}
    for cond in ['Ctrl','T2']:
        m=(condition==cond)&positive;row['n_'+cond]=int(m.sum())
        center=np.median(xy[m],axis=0) if m.any() else np.array([np.nan,np.nan])
        row[cond+'_x'],row[cond+'_y']=center
    row['median_displacement']=float(np.linalg.norm(np.array([row['T2_x']-row['Ctrl_x'],row['T2_y']-row['Ctrl_y']])))
    return row

def marker_stage():
    a=ad.read_h5ad(OUT/'objects/scored_base.h5ad');base,arr=build_marker_complementarity(a);rows=[];shifts=[];xy=a.obsm['X_umap'];cond=a.obs.condition.astype(str).to_numpy();cl=a.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy();tis=a.obs.tissue.astype(str).to_numpy()
    for tissue in ['all']+list(TISSUES):
        for c in REFINED_ORDER:
            m=(cl==c)&((tis==tissue) if tissue!='all' else True)
            for co in CONDITION_ORDER:
                mm=m&(cond==co)
                rows.append(marker_complementarity_row(arr['Il4_expression'],arr['Il4_counts'],arr['Klrd1_expression'],arr['Klrd1_counts'],mm,tissue=tissue,cluster=c,condition=co))
            for marker in ['all_cells','Il4','Klrd1']:
                positive=np.ones(m.sum(),bool) if marker=='all_cells' else arr[marker+'_counts'][m]>0
                shifts.append({'tissue':tissue,'cluster':c,'marker':marker,**centroid_summary(xy[m],cond[m],positive)})
    counts=pd.DataFrame(rows);counts['Il4_positive']=counts.both_positive+counts.il4_only;counts['Klrd1_positive']=counts.both_positive+counts.klrd1_only
    csv(counts,'IL4_CD94_counts_by_cluster_tissue_condition.csv');csv(pd.DataFrame(shifts),'IL4_CD94_fixed_UMAP_median_shifts.csv')
    for tissue in ['all']+list(TISSUES):
        fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
        for ax,(marker,co) in zip(axes.flat,[('Il4','Ctrl'),('Il4','T2'),('Klrd1','Ctrl'),('Klrd1','T2')]):
            active=(cond==co)&((tis==tissue) if tissue!='all' else True);v=arr[marker+'_expression'];hi=np.quantile(v[v>0],.99)
            ax.scatter(xy[:,0],xy[:,1],s=1,c='#e3e3e3',rasterized=True);im=ax.scatter(xy[active,0],xy[active,1],s=3,c=v[active],cmap='viridis',vmin=0,vmax=hi,rasterized=True)
            for c in REFINED_ORDER:
                m=active&(cl==c)
                if m.any():
                    center=np.median(xy[m],axis=0);npos=int(np.sum(arr[marker+'_counts'][m]>0));ax.annotate(f'{c}\n{npos}/{m.sum()}',center,fontsize=7,weight='bold')
            ax.set_title(f'{tissue} / {co} / {marker}: positive / total');ax.set_xticks([]);ax.set_yticks([]);fig.colorbar(im,ax=ax,shrink=.7)
        save(fig,'06_counts_'+tissue)
    for tissue in ['all']+list(TISSUES):
        s=pd.DataFrame(shifts);s=s[(s.tissue==tissue)&(s.marker!='all_cells')]
        fig,ax=plt.subplots(figsize=(10,7));ax.scatter(xy[:,0],xy[:,1],s=1,c='#e3e3e3',rasterized=True)
        for _,r in s.iterrows():
            if min(r.n_Ctrl,r.n_T2)<5:continue
            color='#b2182b' if r.marker=='Il4' else '#2166ac';ax.annotate('',xy=(r.T2_x,r.T2_y),xytext=(r.Ctrl_x,r.Ctrl_y),arrowprops={'arrowstyle':'->','color':color,'lw':1.2});ax.text(r.T2_x,r.T2_y,f'{r.cluster} {r.marker}',fontsize=7,color=color)
        ax.set_title(f'{tissue}: positive-cell median shifts Ctrl → T2\nDescriptive UMAP positions; these arrows are NOT RNA velocity');save(fig,'07_median_shifts_'+tissue)
    props=build_cluster_proportion_tables(a.obs,a.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy(),REFINED_ORDER)
    for i,t in enumerate(props):
        if isinstance(t,pd.DataFrame):csv(t,f'frequency_reconciliation_{i}.csv')
    # Reuse verified denominator rendering with the function's documented return tuple.
    render_denominator_reconciliation(props[1],REFINED_ORDER,OUT/'figures/08_frequency_denominators.png',title_suffix='Meeting follow-up: explicit condition denominators')

def de_stage(local=False):
    a=ad.read_h5ad(OUT/'objects/scored_base.h5ad');jobs=[];statuses=[]
    if not local:
        jobs.append(('global','global','ALL','all',a))
        for tissue in TISSUES:
            sub=a[a.obs.tissue.astype(str).eq(tissue)].copy();jobs.append(('tissue__'+tissue,'tissue','ALL',tissue,sub))
            for c in REFINED_ORDER:
                sub=a[a.obs.tissue.astype(str).eq(tissue)&a.obs[REFINED_CLUSTER_KEY].astype(str).eq(c)].copy();jobs.append((f'cluster_tissue__{c}__{tissue}','cluster_tissue',c,tissue,sub))
    else:
        for tissue in TISSUES:
            a=ad.read_h5ad(OUT/f'objects/{tissue}_reclustered.h5ad')
            for c in sorted(a.obs.meeting_cluster.astype(str).unique()):
                sub=a[a.obs.meeting_cluster.astype(str).eq(c)].copy();jobs.append((f'local__{c}','local_recluster',c,tissue,sub))
    for unit_id,scope,c,tissue,sub in jobs:
        ns=sub.obs.condition.astype(str).value_counts();nt=int(ns.get('T2',0));nc=int(ns.get('Ctrl',0));row={'unit_id':unit_id,'scope':scope,'cluster':c,'tissue':tissue,'n_tumor':nt,'n_control':nc,'strict_min20':min(nt,nc)>=20,'status':'insufficient_cells' if min(nt,nc)<2 else 'computed'}
        if min(nt,nc)>=2:
            p=OUT/'de'/f'{unit_id}.csv.gz'
            if not p.exists():
                log('DE '+unit_id+f' T={nt} C={nc}');sub.obs[REFINED_CLUSTER_KEY]=pd.Categorical([c]*sub.n_obs)
                d=run_de_unit(sub,DEUnit(unit_id,scope,c,None));d['tissue']=tissue;d.to_csv(p,index=False)
                js(p.with_suffix('.receipt.json'),{'cell_sha256':hashlib.sha256('\n'.join(sub.obs_names).encode()).hexdigest(),'source_sha256':sha256_file(INPUT),'file_sha256':sha256_file(p),'n_genes':sub.n_vars,'method':'Wilcoxon tie-corrected BH; full log1p raw; T2 vs Ctrl'})
            else:
                receipt=json.loads(p.with_suffix('.receipt.json').read_text())
                if receipt['file_sha256']!=sha256_file(p) or receipt['source_sha256']!=sha256_file(INPUT):raise ValueError('DE cache mismatch '+str(p))
            d=pd.read_csv(p);row['n_up_fdr05_fc025']=int(d.robust_tumor_up.sum());row['n_down_fdr05_fc025']=int(d.robust_control_up.sum())
        statuses.append(row);csv(pd.DataFrame(statuses),'de_status_local.csv' if local else 'de_status_original.csv')

def ora_table(d,sets,universe,rule='robust'):
    rows=[];M=len(universe)
    for direction in ['Tumor','Control']:
        if rule=='robust':selected=set(d.loc[(d.pvals_adj<=.05)&((d.logfoldchanges>=.25) if direction=='Tumor' else (d.logfoldchanges<=-.25)),'gene'])&universe
        elif rule=='legacy':selected=set(d.loc[(d.pvals_adj<=.05)&(d.logfoldchanges.abs()>=1.2),'gene'])&universe
        elif rule=='paper':selected=set(d.loc[(d.pvals<=.05)&((d.logfoldchanges>=np.log2(1.5)) if direction=='Tumor' else (d.logfoldchanges<=-np.log2(1.5))),'gene'])&universe
        for term,genes in sets.items():
            present=genes&universe;K=len(present)
            if K<3 or K>500:continue
            hit=selected&present;k=len(hit);n=len(selected)
            p=float(hypergeom.sf(k-1,M,K,n)) if n and k else 1.
            rows.append({'term':term,'direction':direction,'n_query':n,'pathway_size':K,'overlap':k,'universe_size':M,'enrichment_ratio':k/max(n*K/M,1e-12),'pvalue':p,'overlap_genes':';'.join(sorted(hit)),'rule':rule})
        if rule=='legacy':break # legacy lists mix directions; do not label enrichment as activation
    r=pd.DataFrame(rows)
    if len(r):
        r['fdr']=r.groupby('direction')['pvalue'].transform(lambda x:bh_adjust(x.to_numpy()));r.loc[r.rule=='legacy','direction']='Mixed_legacy_DEG_union'
    return r

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['scores','embeddings','markers','de','de_local']);args=parser.parse_args();set_source()
    funcs={'scores':score_and_audit,'embeddings':embedding_stage,'markers':marker_stage,'de':de_stage,'de_local':lambda:de_stage(True)}
    start=time.time();log('START '+args.stage);funcs[args.stage]();js(OUT/f'logs/{args.stage}.completed.json',{'stage':args.stage,'elapsed_seconds':time.time()-start,'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())});log('COMPLETE '+args.stage)
if __name__=='__main__':main()
