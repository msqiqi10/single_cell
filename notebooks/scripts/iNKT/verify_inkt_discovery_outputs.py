"""Audit real discovery outputs against source cells, complete test families and report claims."""
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from pypdf import PdfReader
from PIL import Image
from run_inkt_discovery import OUT,PREV,INPUT,QC,units_and_obs,mask_unit,qc_strata,sha256_file
from run_inkt_meeting_pathways import libraries,HEAT

def main():
    checks=[]
    def check(ok,label):
        if not ok:raise AssertionError(label)
        checks.append(label)
    print('Loading source and matching audit',flush=True)
    a,units=units_and_obs();obs=a.obs;qa=pd.read_csv(OUT/'tables/QC_matching_audit.csv').set_index('unit_id');members=pd.read_csv(OUT/'tables/QC_selected_cell_ids.csv.gz');strata=pd.read_csv(OUT/'tables/QC_stratum_counts.csv')
    old=json.loads((PREV/'manifest.json').read_text());manifest=json.loads((OUT/'manifest.json').read_text())
    check(sha256_file(INPUT)==old['source_sha256']==manifest['source_SHA256'],'source H5AD unchanged from meeting analysis')
    check(len(units)==35 and len(a)==15532 and a.n_vars==10670,'expected source and contrast inventory')
    check(not members.duplicated(['unit_id','cell_id']).any(),'selected cell IDs unique within each contrast')
    check(set(members.unit_id)==set(units.unit_id)==set(qa.index),'all contrast IDs accounted for')
    genes=set(a.raw.var_names);de_cache={}
    for r in units.itertuples():
        uid=r.unit_id;o=obs.loc[mask_unit(obs,r)].copy();ids=members.loc[members.unit_id.eq(uid),'cell_id'];sel=obs.loc[ids];q=qa.loc[uid]
        check(set(ids)<=set(o.index),uid+' selected cells belong to contrast')
        available=pd.crosstab(qc_strata(o),o.condition.astype(str)).reindex(columns=['Ctrl','T2'],fill_value=0)
        bins=pd.Series(qc_strata(o),index=o.index)
        chosen=pd.crosstab(bins.loc[ids].to_numpy(),sel.condition.astype(str).to_numpy()).reindex(index=available.index,columns=['Ctrl','T2'],fill_value=0)
        expected=available.min(axis=1)
        check((chosen.Ctrl==expected).all() and (chosen.T2==expected).all(),uid+' both conditions exactly balanced within every available QC stratum')
        check(int(expected.sum())==q.n_each_matched,uid+' matched counts agree with audit')
        recorded=strata[strata.unit_id.eq(uid)].set_index('stratum').sort_index()
        check(np.array_equal(recorded.n_each_selected,expected.sort_index()),uid+' saved stratum counts correct')
        for col in QC:
            c=sel.loc[sel.condition.astype(str).eq('Ctrl'),col].to_numpy(float);t=sel.loc[sel.condition.astype(str).eq('T2'),col].to_numpy(float)
            smd=(np.mean(t)-np.mean(c))/np.sqrt((np.var(t,ddof=1)+np.var(c,ddof=1))/2)
            check(np.isclose(smd,q['after_'+col+'_SMD'],atol=1e-10),uid+' independently checked SMD '+col)
        dp=OUT/f'de_qc/{uid}.csv.gz'
        check(dp.exists()==(q.n_each_matched>=20),uid+' DE availability follows explicit minimum cells')
        if not dp.exists():continue
        d=pd.read_csv(dp).set_index('gene');de_cache[(uid,'QC_matched')]=d
        check(len(d)==10670 and set(d.index)==genes,uid+' complete measured DE universe')
        check(d.n_control.eq(q.n_each_matched).all() and d.n_tumor.eq(q.n_each_matched).all(),uid+' DE sample counts agree with selected cells')
        check(d.pvals.between(0,1).all() and d.pvals_adj.between(0,1).all(),uid+' valid p and FDR')
        check(np.allclose(d.pvals_adj,multipletests(d.pvals,method='fdr_bh')[1]),uid+' independent BH recomputation')
        raw=a.raw[ids,d.index].X.astype(np.float64);co=sel.condition.astype(str).to_numpy();c=np.asarray(raw[co=='Ctrl'].mean(axis=0)).ravel();t=np.asarray(raw[co=='T2'].mean(axis=0)).ravel();fc=np.log2((np.expm1(t)+1e-9)/(np.expm1(c)+1e-9))
        check(np.allclose(fc,d.logfoldchanges,atol=3e-6),uid+' effects recomputed from selected source cells')
        check(np.allclose(c,d.mean_log1p_expression_control,atol=3e-5) and np.allclose(t,d.mean_log1p_expression_tumor,atol=3e-5),uid+' group expression means checked')
    print('Checking full ORA universe, direction, driver sets and corrections',flush=True)
    screen=pd.read_csv(OUT/'tables/functional_screen_all.csv.gz',low_memory=False);libs=libraries()
    grouping=['unit_id','analysis','library','sensitivity','direction']
    check(not screen.duplicated(grouping+['term']).any(),'no duplicated functional hypotheses')
    group_count=0
    for keys,d in screen.groupby(grouping,sort=False):
        uid,analysis,lib,sens,direction=keys;group_count+=1
        if (uid,analysis) not in de_cache:de_cache[(uid,analysis)]=pd.read_csv(PREV/f'de/{uid}.csv.gz').set_index('gene')
        deg=de_cache[(uid,analysis)];u={g for g in genes if sens=='full' or not HEAT.match(g)}
        selected=set(deg.index[(deg.pvals_adj<=.05)&((deg.logfoldchanges>=.25) if direction=='Tumor' else (deg.logfoldchanges<=-.25))])&u
        eligible={term:g&u for term,g in libs[lib].items() if 3<=len(g&u)<=500}
        label='/'.join(keys)
        check(set(d.term)==set(eligible),label+' complete eligible gene-set family')
        check(d.n_query.eq(len(selected)).all() and d.universe_size.eq(len(u)).all(),label+' correct query and background')
        exp_size=np.array([len(eligible[t]) for t in d.term]);hits=[';'.join(sorted(eligible[t]&selected)) for t in d.term];exp_overlap=np.array([len(eligible[t]&selected) for t in d.term])
        check(np.array_equal(d.pathway_size,exp_size) and np.array_equal(d.overlap,exp_overlap),label+' pathway membership and overlap counts')
        check(d.overlap_genes.fillna('').tolist()==hits,label+' exact driver identities')
        expected_p=hypergeom.sf(exp_overlap-1,len(u),exp_size,len(selected))
        check(np.allclose(expected_p,d.pvalue,rtol=1e-9,atol=1e-100),label+' independent hypergeometric p values')
        check(np.allclose(d.fdr,multipletests(d.pvalue,method='fdr_bh')[1],rtol=1e-9,atol=1e-100),label+' per-family BH')
    for keys,d in screen.groupby(['analysis','library','sensitivity','direction']):
        check(np.allclose(d.FDR_across_units,multipletests(d.pvalue,method='fdr_bh')[1],rtol=1e-9,atol=1e-100),str(keys)+' cross-contrast BH')
    check(group_count==(35+29)*2*2*2,'all expected functional analysis families present')
    pri=pd.read_csv(OUT/'tables/discovery_priority_all4_FDR05.csv');qs=['original_FDR','stable_FDR','QC_FDR','stable_QC_FDR'];fcs=['original_log2FC','stable_log2FC','QC_log2FC','stable_QC_log2FC']
    check(len(pri)==17 and (pri[qs]<=.05).all().all(),'17 priority records pass all four FDR thresholds')
    check((pri[fcs].abs()>=.25).all().all() and (np.sign(pri[fcs]).nunique(axis=1)==1).all(),'priority effects all sufficiently large and same direction')
    check((pri.original_cells_in_stable_fraction>=.5).all() and (pri.max_postmatch_QC_SMD<=.2).all(),'priority cluster overlap and QC balance thresholds')
    draws=pd.read_csv(OUT/'tables/candidate_cell_subsampling.csv')
    for r in pri.itertuples():
        s=draws[draws.unit_id.eq(r.unit_id)&draws.gene.eq(r.gene)]
        check(len(s)==5 and set(s.seed)=={11,23,37,53,71} and np.all(np.sign(s.QC_matched_log2FC)==np.sign(r.original_log2FC)),r.unit_id+'/'+r.gene+' five actual selections agree with original sign')
        for k,(uid,kind) in enumerate([(r.unit_id,'original_counts'),(r.stable_unit,'original_counts'),(r.unit_id,'QC_matched'),(r.stable_unit,'QC_matched')]):
            src=de_cache[(uid,kind)].loc[r.gene]
            check(np.isclose(getattr(r,fcs[k]),src.logfoldchanges) and np.isclose(getattr(r,qs[k]),src.pvals_adj),r.unit_id+'/'+r.gene+'/'+str(k)+' priority evidence matches DE')
    focal=pd.read_csv(OUT/'tables/discovery_BM_C4_focal_genes.csv')
    check((~focal.present_in_legacy_DEG_tables).all(),'six focal genes absent from extracted legacy DEG lists')
    funcs=pd.read_csv(OUT/'tables/discovery_functional_priority.csv');sp=funcs[funcs.term.eq('Spliceosome')].iloc[0]
    check(sp.original_GSEA_FDR>.05 and sp.no_heat_GSEA_FDR>.05 and sp.no_heat_FDR<=.05,'splicing ORA versus GSEA distinction preserved')
    fig_pdf=OUT/'discovery_brief_20260905.pdf';pdf=PdfReader(fig_pdf)
    check(len(pdf.pages)==6 and all(len(p.extract_text())>100 for p in pdf.pages),'six readable evidence report pages')
    for p in (OUT/'figures').glob('*.png'):
        with Image.open(p) as im:im.verify()
        check(True,p.name+' readable image')
    check((OUT/'logs/run.exit').read_text().strip()=='0' and (OUT/'logs/brief.exit').read_text().strip()=='0','analysis and brief jobs finished successfully')
    check((OUT/'logs/verification_tests.exit').read_text().strip()=='0','related unit tests passed')
    result={'status':'passed','checks_passed':len(checks),'source_cells':len(a),'source_genes':a.n_vars,'audited_contrasts':len(units),'QC_matched_DE_contrasts':int(qa.status.eq('computed').sum()),'functional_hypotheses_checked':len(screen),'functional_families_checked':group_count,'strict_priority_records':len(pri),'unit_test_log':'logs/verification_tests.log','PDF_pages':len(pdf.pages),'visual_review':'pending','completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'script_SHA256':sha256_file(Path(__file__))}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
