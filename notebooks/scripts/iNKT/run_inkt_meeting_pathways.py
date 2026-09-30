"""Meeting-scoped KEGG/Reactome, signed DEG, legacy PPT and paper validation."""
from run_inkt_meeting_followup import *
import gseapy as gp
from run_inkt_c5_paper_followup import match_gene_sets_to_ranking

LIBRARIES={
 'KEGG_Mouse_2019':BASE/'tables/20260830_KEGG_2019_Mouse.gmt',
 'KEGG_Human_2021_ortholog':OUT/'sources/KEGG_2021_Human_mapped_mouse.gmt',
 'Reactome_2022_ortholog':OUT/'sources/Reactome_2022_Human_mapped_mouse.gmt',
 'PPT_Reactome_current':OUT/'sources/Legacy_named_Reactome_current_mapped_mouse.gmt',
}
HEAT=re.compile(r'^(Hsp|Dnaj)',re.I)

def libraries():
    sets={k:read_gmt(p) for k,p in LIBRARIES.items()}
    genes=pd.read_csv(OUT/'de/global.csv.gz',usecols=['gene']).gene.tolist()
    sets['KEGG_Mouse_2019']=match_gene_sets_to_ranking(sets['KEGG_Mouse_2019'],genes)
    audit=pd.read_csv(OUT/'tables/legacy_reactome_full_membership_audit.csv').fillna('')
    # Hypotheses defined by stable reaction ID, not repeated/truncated display names.
    sets['PPT_Reactome_current']={r.stable_id:set(r.mouse_genes.split(';'))-{''} for r in audit.itertuples() if r.status=='full_current_membership_retrieved'}
    return sets

def clean_rank(d,exclude_heat=False):
    d=d.copy();d['scores']=pd.to_numeric(d.scores,errors='coerce')
    if d.gene.duplicated().any():raise ValueError('Duplicated ranking genes')
    removed=d.loc[d.gene.map(lambda x:bool(HEAT.match(x))) if exclude_heat else np.zeros(len(d),dtype=bool),'gene']
    d=d[~d.gene.isin(removed)];finite=np.isfinite(d.scores)
    # Undefined scores from constant features carry no directional evidence.
    d.loc[~finite,'scores']=0.
    d=d.sort_values(['scores','gene'],ascending=[False,True],kind='mergesort')
    return d,{'n_ranked':len(d),'n_nonfinite_set_to_zero':int((~finite).sum()),'tie_fraction':float(d.scores.duplicated().mean()),'excluded_heat_genes':';'.join(removed)}

def enrichment_stage():
    libs=libraries();status=pd.read_csv(OUT/'tables/de_status_original.csv')
    # Original tissue x refined cluster contrasts form the fixed comparison grid.
    jobs=status[(status.status=='computed')&status.strict_min20]
    worker=int(os.environ.get('MEETING_WORKER','0'));workers=int(os.environ.get('MEETING_WORKERS','1'));jobs=jobs.iloc[worker::workers]
    audits=[];started=time.time()
    for row in jobs.itertuples():
        d=pd.read_csv(OUT/f'de/{row.unit_id}.csv.gz')
        for sensitivity in ['full','without_heat_shock']:
            clean,au=clean_rank(d,sensitivity!='full');universe=set(clean.gene);audits.append({'unit_id':row.unit_id,'sensitivity':sensitivity,**au})
            for lib,sets in libs.items():
                rules=['robust','paper'] if lib=='KEGG_Mouse_2019' else ['robust','legacy']
                if lib=='Reactome_2022_ortholog':rules=['robust']
                for rule in rules:
                    path=OUT/f'enrichment/{row.unit_id}__{sensitivity}__{lib}__ORA_{rule}.csv.gz'
                    if not path.exists():
                        r=ora_table(clean,sets,universe,rule);r['unit_id']=row.unit_id;r['library']=lib;r['sensitivity']=sensitivity
                        r.to_csv(path,index=False)
            # Same fixed family and signed ranking for full and heat-shock exclusion.
            path=OUT/f'enrichment/{row.unit_id}__{sensitivity}__KEGG_Mouse_2019__GSEA.csv.gz'
            if not path.exists():
                log(f'GSEA {row.unit_id} {sensitivity}; elapsed {time.time()-started:.0f}s')
                rnk=clean[['gene','scores']];sets={t:sorted(g&universe) for t,g in libs['KEGG_Mouse_2019'].items() if 5<=len(g&universe)<=500}
                res=gp.prerank(rnk=rnk,gene_sets=sets,outdir=None,min_size=5,max_size=500,permutation_num=5000,weight=1,ascending=None,threads=2,no_plot=True,seed=SEED,verbose=False,method='permutation').res2d
                res['unit_id']=row.unit_id;res['sensitivity']=sensitivity;res['library']='KEGG_Mouse_2019';res['permutation_num']=5000;res.to_csv(path,index=False)
        csv(pd.DataFrame(audits),f'enrichment_ranking_audit_worker{worker}.csv')
    js(OUT/f'logs/pathways_worker{worker}.completed.json',{'elapsed_seconds':time.time()-started,'n_original_min20_contrasts':len(jobs),'gsea_permutations':5000,'gsea_null':'gene-set permutations, not biological replicates','paper_ORA_threshold':'nominal p<=0.05, linear fold change>=1.5; supplementary methods p6; implemented log2FC>=log2(1.5)','background':'all 10670 DE-tested genes; Hsp*/Dnaj* removed from ranking, background and sets in sensitivity','FDR':'BH within complete eligible library/direction/contrast for ORA; GSEA empirical q within contrast, no cross-contrast control','limitation':'Mouse iNKT/T2 compared with paper NK/CML. No validated K1+K2 correspondence and no identified biological replicates.'})

def candidates(unit,status):
    if unit=='global':return ['global']
    if unit.startswith('tissue_'):return ['tissue__'+unit[7:]]
    tissue=next((t for t in TISSUES if unit.endswith(t)),None)
    return status.loc[(status.tissue==tissue)&(status.scope=='cluster_tissue')&(status.status=='computed'),'unit_id'].tolist()

def legacy_gene_stage():
    old=pd.read_csv(OLD/'legacy_ppt_de_tables_extracted.csv');status=pd.read_csv(OUT/'tables/de_status_original.csv');cache={}
    rows=[];summ=[]
    for lu,grp in old.groupby('legacy_unit'):
        for uid in candidates(lu,status):
            if uid not in cache:cache[uid]=pd.read_csv(OUT/f'de/{uid}.csv.gz').set_index('gene')
            d=cache[uid];s=status.set_index('unit_id').loc[uid];matched=[]
            for r in grp.itertuples():
                gene=str(r.names);q=d.loc[gene] if gene in d.index else None
                row={'legacy_unit':lu,'slide':r.slide,'current_unit':uid,'gene':gene,'legacy_log2FC':r.logFC,'legacy_fdr':r.pvals_adj,'n_tumor':s.n_tumor,'n_control':s.n_control,'strict_min20':s.strict_min20,'measured':q is not None}
                if q is not None:row.update(current_log2FC=q.logfoldchanges,current_fdr=q.pvals_adj,current_p=q.pvals,direction_concordant=bool(np.sign(r.logFC)==np.sign(q.logfoldchanges)),current_robust_significant=bool(q.pvals_adj<=.05 and abs(q.logfoldchanges)>=.25),current_legacy_significant=bool(q.pvals_adj<=.05 and abs(q.logfoldchanges)>=1.2));matched.append(row)
                rows.append(row)
            m=pd.DataFrame(matched);rho=float(spearmanr(m.legacy_log2FC,m.current_log2FC).statistic) if len(m)>2 else np.nan
            summ.append({'legacy_unit':lu,'current_unit':uid,'n_legacy':len(grp),'n_measured':len(m),'direction_concordance':m.direction_concordant.mean() if len(m) else np.nan,'n_concordant_robust':int((m.direction_concordant&m.current_robust_significant).sum()) if len(m) else 0,'n_concordant_legacy_threshold':int((m.direction_concordant&m.current_legacy_significant).sum()) if len(m) else 0,'spearman_logFC':rho,'strict_min20':s.strict_min20,'interpretation':'response-pattern comparison, NOT cell-identity crosswalk'})
    csv(pd.DataFrame(rows),'legacy_DEG_signed_validation_all_candidates.csv');csv(pd.DataFrame(summ),'legacy_DEG_response_similarity.csv')
    # Genes detected at cluster level but diluted at tissue/global level. This is descriptive.
    novel=[]
    global_d=cache.get('global',pd.read_csv(OUT/'de/global.csv.gz').set_index('gene'))
    for r in status[(status.scope=='cluster_tissue')&status.strict_min20&(status.status=='computed')].itertuples():
        d=pd.read_csv(OUT/f'de/{r.unit_id}.csv.gz').set_index('gene');t=pd.read_csv(OUT/f'de/tissue__{r.tissue}.csv.gz').set_index('gene')
        for gene,q in d[(d.pvals_adj<=.05)&(d.logfoldchanges.abs()>=.25)].iterrows():
            b=t.loc[gene];g=global_d.loc[gene]
            if b.pvals_adj>.05 or abs(b.logfoldchanges)<.25 or np.sign(b.logfoldchanges)!=np.sign(q.logfoldchanges):
                novel.append({'unit_id':r.unit_id,'tissue':r.tissue,'cluster':r.cluster,'gene':gene,'cluster_log2FC':q.logfoldchanges,'cluster_fdr':q.pvals_adj,'tissue_log2FC':b.logfoldchanges,'tissue_fdr':b.pvals_adj,'global_log2FC':g.logfoldchanges,'global_fdr':g.pvals_adj,'interpretation':'different detection/effect pattern; not a tested cluster-by-condition interaction'})
    csv(pd.DataFrame(novel),'cluster_signals_diluted_in_tissue.csv')

# Explicit KEGG IDs from slide text; do not fuzzy match truncated terms.
KEGG_ID_NAMES={'05417':'Lipid and atherosclerosis','05020':'Prion disease','04010':'MAPK signaling pathway','04915':'Estrogen signaling pathway','05162':'Measles','04141':'Protein processing in endoplasmic reticulum','05134':'Legionellosis','04213':'Longevity regulating pathway - multiple species','04612':'Antigen processing and presentation','05145':'Toxoplasmosis','03040':'Spliceosome'}
def normterm(s):return re.sub(r'[^a-z0-9]','',str(s).lower().replace('(cams)','').replace('(ibd)',''))
def matchterm(term,terms):
    match=[t for t in terms if normterm(t)==normterm(term)]
    return match[0] if len(match)==1 else None

def comparison_stage():
    legacy_gene_stage();libs=libraries();status=pd.read_csv(OUT/'tables/de_status_original.csv');h2m=json.loads((OUT/'sources/human_to_measured_mouse.json').read_text());cache={}
    def read(uid,sens,lib,method):
        k=(uid,sens,lib,method)
        if k not in cache:
            p=OUT/'enrichment'/('__'.join(k)+'.csv.gz');cache[k]=pd.read_csv(p) if p.exists() else pd.DataFrame()
        return cache[k]
    # Every original PPT row, every valid tissue-constrained candidate, and both sensitivities.
    old=pd.read_csv(OLD/'legacy_ppt_pathway_tables_extracted.csv').fillna('');rx=pd.read_csv(OUT/'tables/legacy_reactome_full_membership_audit.csv').fillna('').set_index('legacy_name');rows=[];drivers=[]
    de_cache={}
    def read_de(uid):
        if uid not in de_cache:de_cache[uid]=pd.read_csv(OUT/f'de/{uid}.csv.gz').set_index('gene')
        return de_cache[uid]
    id_names={gid:sorted(set(g.loc[~g.NAME.str.endswith('...'),'NAME'])) for gid,g in old.groupby('GS_ID')}
    for i,r in old.iterrows():
        lib=None;term=None;resolved_name=r.NAME
        fullnames=id_names.get(r.GS_ID,[])
        if r.NAME.endswith('...') and len(fullnames)==1:resolved_name=fullnames[0]
        if str(r.SOURCE).startswith('Reactome') and resolved_name in rx.index:
            lib='PPT_Reactome_current';term=rx.loc[resolved_name,'stable_id'] or None
        elif str(r.SOURCE).startswith('KEGG'):
            lib='KEGG_Human_2021_ortholog';m=re.search(r'hsa(\d{5})',r.NAME)
            term=matchterm(KEGG_ID_NAMES.get(m.group(1),'') if m else '',libs[lib])
        original_genes=[g.strip() for g in str(r.GENE_SYM).split(',') if g.strip()]
        for uid in candidates(r.legacy_unit,status):
            d=read_de(uid);valid=bool(status.set_index('unit_id').loc[uid,'strict_min20'])
            for g in original_genes:
                mg=h2m.get(g);q=d.loc[mg] if mg in d.index else None
                drivers.append({'ppt_row':i,'slide':r.slide,'legacy_unit':r.legacy_unit,'current_unit':uid,'term':r.NAME,'source_human_gene':g,'measured_mouse_gene':mg,'heat_shock':bool(mg and HEAT.match(mg)),'log2FC':q.logfoldchanges if q is not None else np.nan,'fdr':q.pvals_adj if q is not None else np.nan,'strict_min20':valid})
            for sens in ['full','without_heat_shock']:
                row={'ppt_row':i,'slide':r.slide,'legacy_unit':r.legacy_unit,'current_unit':uid,'legacy_name':r.NAME,'legacy_GS_ID':r.GS_ID,'resolved_name':resolved_name,'legacy_source':r.SOURCE,'legacy_pFDR':r.pFDR,'mapped_library':lib,'mapped_term':term,'sensitivity':sens,'strict_min20':valid,'status':'source_library_unavailable' if not lib else ('term_unresolved' if not term else ('below_min20' if not valid else 'below_minimum_set_size'))}
                if lib and term and valid:
                    for rule in ['legacy','robust']:
                        t=read(uid,sens,lib,'ORA_'+rule)
                        if len(t):
                            q=t[t.term.eq(term)]
                            for qrow in q.itertuples():
                                prefix=rule+'_'+qrow.direction;row.update({prefix+'_fdr':qrow.fdr,prefix+'_overlap':qrow.overlap,prefix+'_genes':qrow.overlap_genes});row['status']='tested_current_membership'
                rows.append(row)
    csv(pd.DataFrame(rows),'PPT_pathway_validation_every_row.csv');csv(pd.DataFrame(drivers),'PPT_pathway_driver_genes_signed.csv')
    paper=pd.read_csv(OUT/'sources/Blood_TableS3.csv');prows=[]
    for r in status[(status.status=='computed')&status.strict_min20].itertuples():
        for sens in ['full','without_heat_shock']:
            ora=read(r.unit_id,sens,'KEGG_Mouse_2019','ORA_paper');gsea=read(r.unit_id,sens,'KEGG_Mouse_2019','GSEA')
            for p in paper.itertuples(index=False,name=None):
                name=p[1].strip();term=matchterm(name,libs['KEGG_Mouse_2019']);row={'unit_id':r.unit_id,'sensitivity':sens,'paper_KEGG_ID':p[0],'paper_term':name,'paper_p':p[3],'paper_fdr':p[4],'in_Fig3E':name in PAPER_TERMS,'matched_term':term,'status':'term_not_in_frozen_2019_library' if term is None else 'not_testable_set_size'}
                if term:
                    for q in ora[ora.term.eq(term)].itertuples():row.update({q.direction+'_ORA_p':q.pvalue,q.direction+'_ORA_fdr':q.fdr,q.direction+'_ORA_drivers':q.overlap_genes});row['status']='tested'
                    for q in gsea[gsea.Term.eq(term)].to_dict('records'):row.update(GSEA_NES=q['NES'],GSEA_q=q['FDR q-val'],GSEA_p=q['NOM p-val'],GSEA_drivers=q['Lead_genes'])
                prows.append(row)
    csv(pd.DataFrame(prows),'Blood_TableS3_all48_pathway_validation.csv')
    # Exact 23-gene iNK-CML set, plus source-figure genes, signed in every contrast.
    pg=pd.read_csv(OUT/'sources/Blood_TableS5.csv');genes=sorted(set(pg['Gene name'])|set(PAPER_GENES)|{'Il4','Klrd1'})
    rows=[]
    for r in status[status.status=='computed'].itertuples():
        d=pd.read_csv(OUT/f'de/{r.unit_id}.csv.gz').set_index('gene')
        for gene in genes:
            q=d.loc[gene] if gene in d.index else None;paperrow=pg.loc[pg['Gene name']==gene]
            rows.append({'unit_id':r.unit_id,'tissue':r.tissue,'cluster':r.cluster,'gene':gene,'strict_min20':r.strict_min20,'in_TableS5':len(paperrow)>0,'paper_linear_FC':float(paperrow['Fold change (CML vs WT)'].iloc[0]) if len(paperrow) else np.nan,'measured':q is not None,'log2FC':q.logfoldchanges if q is not None else np.nan,'fdr':q.pvals_adj if q is not None else np.nan,'pvalue':q.pvals if q is not None else np.nan,'pct_T2':q.pct_expressing_tumor if q is not None else np.nan,'pct_Ctrl':q.pct_expressing_control if q is not None else np.nan})
    csv(pd.DataFrame(rows),'Blood_reference_driver_gene_validation.csv')
    # Display selection and redundancy are separate from statistical testing.
    top=[];redundancy=[]
    for r in status[(status.status=='computed')&status.strict_min20].itertuples():
        t=read(r.unit_id,'full','KEGG_Mouse_2019','GSEA').copy();t=t[pd.to_numeric(t['FDR q-val'])<=.05].sort_values(['FDR q-val','NES'],ascending=[True,False]);t['display_rank']=np.arange(1,len(t)+1);t['selected_top10']=t.display_rank<=10;t['selection_rule']='All q<=.05 terms ordered by q then NES; first10 displayed; no category grouping';top.append(t)
        for a,b in combinations(t.head(10).to_dict('records'),2):
            ga=set(str(a['Lead_genes']).split(';'));gb=set(str(b['Lead_genes']).split(';'));shared=ga&gb
            redundancy.append({'unit_id':r.unit_id,'term_a':a['Term'],'term_b':b['Term'],'leading_edge_jaccard':len(shared)/len(ga|gb),'shared_drivers':';'.join(sorted(shared)),'shared_heat_shock':';'.join(sorted(g for g in shared if HEAT.match(g)))})
    csv(pd.concat(top,ignore_index=True),'KEGG_display_selection_all_significant.csv');csv(pd.DataFrame(redundancy),'KEGG_top10_driver_redundancy.csv')
    js(OUT/'logs/comparisons.completed.json',{'PPT_pathway_rows':len(old),'PPT_gene_rows':int(pd.read_csv(OLD/'legacy_ppt_de_tables_extracted.csv').shape[0]),'paper_pathways':len(paper),'paper_TableS5_genes':len(pg)})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['enrichment','comparisons']);args=p.parse_args();set_source()
    {'enrichment':enrichment_stage,'comparisons':comparison_stage}[args.stage]()
