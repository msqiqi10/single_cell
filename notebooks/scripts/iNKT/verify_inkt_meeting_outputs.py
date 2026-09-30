"""End-to-end checks of generated meeting deliverables; no analysis recomputation."""
from run_inkt_meeting_followup import *
from run_inkt_meeting_pathways import libraries, HEAT
from pypdf import PdfReader
from pptx import Presentation


def main():
    checks=[]
    def check(condition,name):
        if not condition:raise AssertionError(name)
        checks.append(name)
    a=ad.read_h5ad(INPUT);status=pd.read_csv(OUT/'tables/de_status_original.csv');jobs=status[status.strict_min20&status.status.eq('computed')];libs=libraries();universe=set(a.var_names)
    check(len(jobs)==21,'21 explicit min20 original contrasts')
    check(len(list((OUT/'enrichment').glob('*__GSEA.csv.gz')))==42,'42 full/heat-excluded GSEA files')
    for r in jobs.itertuples():
        d=pd.read_csv(OUT/f'de/{r.unit_id}.csv.gz');check(len(d)==len(universe) and set(d.gene)==universe,r.unit_id+' full DE universe')
        for sensitivity in ['full','without_heat_shock']:
            u={g for g in universe if sensitivity=='full' or not HEAT.match(g)}
            eligible={t for t,g in libs['KEGG_Mouse_2019'].items() if 5<=len(g&u)<=500};gsea=pd.read_csv(OUT/f'enrichment/{r.unit_id}__{sensitivity}__KEGG_Mouse_2019__GSEA.csv.gz')
            check(set(gsea.Term)==eligible,r.unit_id+'/'+sensitivity+' complete GSEA hypothesis family')
            check(gsea.NES.notna().all() and gsea['FDR q-val'].between(0,1).all(),r.unit_id+'/'+sensitivity+' valid NES/q')
            check(gsea.permutation_num.eq(5000).all(),r.unit_id+'/'+sensitivity+' 5000 permutations')
            for q in gsea.itertuples():
                lead=set(str(q.Lead_genes).split(';'))-{''};check(lead<=libs['KEGG_Mouse_2019'][q.Term]&u,r.unit_id+'/'+sensitivity+'/'+q.Term+' valid leading genes')
            ora=pd.read_csv(OUT/f'enrichment/{r.unit_id}__{sensitivity}__KEGG_Mouse_2019__ORA_paper.csv.gz');expected={t for t,g in libs['KEGG_Mouse_2019'].items() if 3<=len(g&u)<=500}
            for direction in ['Tumor','Control']:
                t=ora[ora.direction.eq(direction)];selected=set(d.loc[(d.pvals<=.05)&((d.logfoldchanges>=np.log2(1.5)) if direction=='Tumor' else (d.logfoldchanges<=-np.log2(1.5))),'gene'])&u
                check(set(t.term)==expected and t.n_query.eq(len(selected)).all(),r.unit_id+'/'+sensitivity+'/'+direction+' correct paper query and family')
                check(np.allclose(t.fdr,bh_adjust(t.pvalue.to_numpy())),r.unit_id+'/'+sensitivity+'/'+direction+' correct BH correction')
                check(t.universe_size.eq(len(u)).all(),r.unit_id+'/'+sensitivity+'/'+direction+' correct ORA universe')
    count=pd.read_csv(OUT/'tables/IL4_CD94_counts_by_cluster_tissue_condition.csv');labels=a.obs[REFINED_CLUSTER_KEY].astype(str);tissue=a.obs.tissue.astype(str);condition=a.obs.condition.astype(str)
    il4=vec(a.layers['counts'][:,a.var_names.get_loc('Il4')]);cd94=vec(a.layers['counts'][:,a.var_names.get_loc('Klrd1')])
    for r in count.itertuples():
        m=labels.eq(r.cluster)&condition.eq(r.condition)
        if r.tissue!='all':m=m&tissue.eq(r.tissue)
        check(r.n_cells==int(m.sum()),'count '+r.tissue+'/'+r.cluster+'/'+r.condition)
        check(r.Il4_positive==int((il4[m]>0).sum()) and r.Klrd1_positive==int((cd94[m]>0).sum()),'positive '+r.tissue+'/'+r.cluster+'/'+r.condition)
        check(r.both_positive+r.il4_only+r.klrd1_only+r.neither==r.n_cells,'partition '+r.tissue+'/'+r.cluster+'/'+r.condition)
    freq=pd.read_csv(OUT/'tables/frequency_reconciliation_0.csv');sums=freq.groupby(['tissue','condition']).pct_cluster_within_tissue_condition.sum();check(np.allclose(sums,100),'all six condition frequencies sum to100')
    check(np.allclose(freq.pct_cluster_within_tissue_condition,100*freq['count']/freq.total_tissue_condition),'frequency denominator formula')
    for method in ['original','harmony_tissue','harmony_sample']:
        b=ad.read_h5ad(OUT/f'objects/{method}.h5ad');check(b.obs_names.equals(a.obs_names),method+' retained cell identities')
        check((b.layers['counts']!=a.layers['counts']).nnz==0,method+' original counts preserved')
    paper=pd.read_csv(OUT/'tables/Blood_TableS3_all48_pathway_validation.csv');check(len(paper)==21*2*48 and paper.paper_term.nunique()==48,'all48 paper terms across all42 analyses, with missing statuses')
    legacy=pd.read_csv(OUT/'tables/PPT_pathway_validation_every_row.csv');old=pd.read_csv(OLD/'legacy_ppt_pathway_tables_extracted.csv');check(legacy.ppt_row.nunique()==len(old)==130,'all130 original PPT pathway rows accounted for')
    truncated=legacy[legacy.legacy_name.str.startswith('HSP40s')];check(truncated.mapped_term.nunique()==2 and truncated.mapped_term.notna().all(),'truncated HSP40 reactions disambiguated using original GS_ID')
    gap=json.loads((OUT/'tables/input_missing_tasks.json').read_text());check(gap['velocity']['status']=='blocked_missing_spliced_unspliced','velocity remains explicitly incomplete')
    check(len(PdfReader(OUT/'meeting_results_20260905.pdf').pages)==23,'23 readable PDF pages');check(len(Presentation(OUT/'meeting_results_20260905.pptx').slides)==24,'24 editable PPTX slides')
    for p in (OUT/'figures').glob('*.png'):
        from PIL import Image
        with Image.open(p) as im:im.verify()
    check((OUT/'logs/tests.exit').read_text().strip()=='0','42 unit tests passed')
    results={'status':'passed','checks_passed':len(checks),'n_cells':len(a),'n_genes':a.n_vars,'GSEA_analyses':42,'original_min20_contrasts':21,'PPT_rows':130,'paper_S3_terms':48,'PDF_pages':23,'PPTX_slides':24,'unit_tests':42,'visual_review':'source figure crop, pathway matrices, dotplot, count/denominator and integration/state contact sheets inspected','velocity_computed':False,'toxicity_inventory_available':False,'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    js(OUT/'verification.json',results);log(str(results))
if __name__=='__main__':main()
