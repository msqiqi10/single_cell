"""Independently verify data integrity, statistics, figures and presentation."""
from pathlib import Path
import json, hashlib, zipfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import hypergeom,mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pypdf import PdfReader
from PIL import Image
from analyze import OUT,ROOT,PREV,RES,SRC,KEY,MODULES,sha
checks=[]
def check(ok,label):
 if not bool(ok):raise AssertionError(label)
 checks.append(label);print('PASS',label,flush=True)
def main():
 a=ad.read_h5ad(PREV/'objects/scored_base.h5ad');o=pd.read_csv(RES/'tables/cell_metadata_scores.csv.gz').set_index('cell');check(o.index.equals(a.obs_names),'Cell identities and order unchanged')
 manifest=json.loads((SRC/'input_manifest.json').read_text())
 for p,h in manifest['inputs'].items():check(sha(ROOT/p)==h,'Input unchanged: '+p)
 x=np.asarray(a.raw[:,MODULES['PAGER_core']].X.mean(axis=1)).ravel();check(np.allclose(x,o.PAGER_core_mean),'Module mean independently recomputed from raw normalized expression')
 c=pd.read_csv(RES/'tables/cytotoxicity_contrasts.csv');ok=c.status.eq('tested');check(c.loc[ok,['n_Ctrl','n_T2']].min().min()>=20,'Minimum condition counts enforced');check(c.loc[~ok,'pvalue'].isna().all(),'Insufficient groups have no statistical claims');check(np.allclose(c.loc[ok,'q_all'],multipletests(c.loc[ok,'pvalue'],method='fdr_bh')[1]),'Module BH independently reproduced')
 r=c[(c.definition=='original')&(c.tissue=='bone_marrow')&(c.cluster=='C4')&(c.module=='PAGER_core')].iloc[0];m=(o.tissue=='bone_marrow')&(o[KEY]=='C4');t=o.loc[m&(o.condition=='T2'),'PAGER_core_mean'];ctrl=o.loc[m&(o.condition=='Ctrl'),'PAGER_core_mean'];check(np.isclose(t.mean()-ctrl.mean(),r.delta_mean),'T2-minus-Ctrl effect direction verified');check(np.isclose(mannwhitneyu(t,ctrl).pvalue,r.pvalue,rtol=1e-8),'C4 cytotoxic test independently reproduced')
 g=pd.read_csv(RES/'tables/GO_ORA_all.csv.gz',low_memory=False);check(g[['pvalue','q_family','q_analysis_global']].notna().all().all(),'GO statistical values finite');check(g[['pvalue','q_family','q_analysis_global']].ge(0).all().all() and g[['pvalue','q_family','q_analysis_global']].le(1).all().all(),'GO probabilities bounded');check(((g.K>=5)&(g.K<=500)&(g.k<=g.K)&(g.k<=g.n)&(g.n<=g.M)).all(),'GO hypergeometric counts valid')
 keys=['analysis','definition','unit_id','direction','namespace']
 for key,d in g.groupby(keys,sort=False):
  check(np.allclose(d.q_family,multipletests(d.pvalue,method='fdr_bh')[1]),'GO BH: '+str(key))
 sample=g.sample(200,random_state=42);expected=hypergeom.sf(sample.k.to_numpy()-1,sample.M.to_numpy(),sample.K.to_numpy(),sample.n.to_numpy());check(np.allclose(expected,sample.pvalue,rtol=1e-8,atol=1e-300),'200 independent hypergeometric recalculations')
 for analysis,d in g.groupby('analysis'):check(np.allclose(d.q_analysis_global,multipletests(d.pvalue,method='fdr_bh')[1]),'Cross-contrast GO BH '+analysis)
 direct=pd.read_csv(RES/'tables/GO_direct_annotation_audit.csv.gz');check(~direct.qualifier.fillna('').str.split('|').apply(lambda x:'NOT' in x).any(),'NOT annotations excluded');check(~direct.evidence.eq('ND').any(),'ND annotations excluded')
 r=pd.read_csv(RES/'tables/trajectory_root_audit.csv');tr=pd.read_csv(RES/'tables/trajectory_cell_values.csv.gz');check(len(tr)==r.n_analyzed.sum(),'Trajectory counts reconcile');check(np.isfinite(tr.local_dpt).all() and tr.local_dpt.between(0,1).all(),'DPT values finite and normalized');check(set(tr.cell).issubset(set(o.index)),'Trajectory cells match source')
 index=pd.read_csv(OUT/'presentation/slide_index.csv');pdf=PdfReader(OUT/'presentation/iNKT_cytotoxicity_GO_20260919.pdf');check(len(pdf.pages)==len(index),'PDF page count')
 with zipfile.ZipFile(OUT/'presentation/iNKT_cytotoxicity_GO_20260919.pptx') as z:
  slides=[s for s in z.namelist() if s.startswith('ppt/slides/slide') and s.endswith('.xml')];check(len(slides)==len(index),'PPTX page count')
 for p in index.figure:
  with Image.open(OUT/p) as im:check(im.width>=1500 and im.height>=900,'Figure resolution '+p);im.verify()
 focus=pd.read_csv(RES/'tables/focus_stable_cluster_crosswalk.csv');check((focus.original_fraction_in_stable>.98).all(),'Focused cluster crosswalk coverage independently checked')
 kill=pd.read_csv(RES/'tables/GO_cytotoxicity_targeted_audit.csv');check(not (kill.q_family<=.05).any(),'No significant killing-related GO condition result overclaimed')
 e1=c[(c.definition=='original')&(c.module=='PAGER_E1')&(((c.tissue=='bone_marrow')&(c.cluster=='C4'))|((c.tissue=='spleen')&(c.cluster=='C3')))];check((e1.q_all>.05).all(),'Broader E1 signature limitation confirmed')
 check(index.title.str.contains('Core signature changes do not imply').any(),'Signature and GO limitations have a dedicated slide')
 for p in (OUT/'code').glob('*.py'):compile(p.read_text(),str(p),'exec')
 check(True,'All dated Python files compile')
 js={'status':'passed','checks':checks,'n_checks':len(checks),'n_slides':len(index),'GO_test_rows':len(g),'GPU_used':False}
 (RES/'verification.json').write_text(json.dumps(js,indent=2))
 files={str(p.relative_to(OUT)):sha(p) for directory in ['code','presentation'] for p in (OUT/directory).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
 (RES/'artifact_manifest.json').write_text(json.dumps(files,indent=2));print('VERIFIED',len(checks),'checks')
if __name__=='__main__':main()
