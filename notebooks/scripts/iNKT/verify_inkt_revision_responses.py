"""Verify document completeness, copied evidence and key numerical claims."""
from pathlib import Path
import json,hashlib,re,zipfile,time
import pandas as pd,numpy as np
from lxml import etree
from pypdf import PdfReader
from PIL import Image
from fontTools.ttLib import TTFont
ROOT=Path('/home/zzz0054/bio3');OUT=ROOT/'output/iNKT_revision_responses_20260906';M=ROOT/'output/iNKT_meeting_followup_20260905';D=ROOT/'output/iNKT_discovery_20260905'

def main():
 checks=[]
 def check(ok,label):
  if not ok:raise AssertionError(label)
  checks.append(label)
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 s=json.loads((OUT/'response_sections.json').read_text());m=json.loads((OUT/'source_manifest.json').read_text())
 check([x['id'] for x in s]==[f'R{i:02d}' for i in range(1,14)],'R01-R13 present exactly once and in order')
 for x in s:
  for k in ['status','done','method','result','pending','talk','meeting']:check(bool(x[k]),x['id']+' nonempty '+k)
  for f in x['files']:check((OUT/f['path']).exists() and sha(OUT/f['path'])==f['sha256'],x['id']+' evidence copy '+f['name'])
  for f in x['figures']:
   with Image.open(OUT/f['path']) as im:im.verify()
   check(True,x['id']+' readable figure '+f['path'])
  t=x['tables'][0];saved=pd.read_csv(OUT/f'data/{x["id"]}_report_table_1.csv',dtype=str,keep_default_na=False)
  check(saved.columns.tolist()==t['headers'] and saved.values.tolist()==t['rows'],x['id']+' report table equals exported CSV')
 for f,v in m['source_SHA256'].items():check(sha(ROOT/f)==v,'unchanged source '+f)
 z=zipfile.ZipFile(OUT/'iNKT_revision_responses_20260906.docx');check(z.testzip() is None,'Word zip integrity')
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};xml=etree.fromstring(z.read('word/document.xml'));text=''.join(xml.xpath('//w:t/text()',namespaces=ns))
 for x in s:
  check(x['id']+'　'+x['title'] in text,x['id']+' Word section heading')
  for k in ['done','method','result','pending','talk']:check(x[k] in text,x['id']+' complete Word '+k)
 check(len(z.namelist())>10 and len([n for n in z.namelist() if n.startswith('word/media/')])==12,'12 embedded Word figures')
 html=(OUT/'iNKT_revision_responses_20260906.html').read_text();md=(OUT/'iNKT_revision_responses_20260906.md').read_text();font=TTFont(OUT/'assets/report_zh_subset.otf');cm=font.getBestCmap()
 plain=re.sub('<[^>]*>','',html);needed={ord(ch) for ch in plain+md if '\u4e00'<=ch<='\u9fff'}
 check(needed<=set(cm),'all Chinese text covered by embedded subset font')
 for x in s:check(f'id="{x["id"]}"' in html and '## '+x['id']+' ' in md,x['id']+' HTML and Markdown headings')
 pdf=PdfReader(OUT/'iNKT_revision_responses_20260906.pdf');pt=[p.extract_text() for p in pdf.pages];combined=''.join(pt)
 for x in s:check(x['id'] in combined and x['title'].replace(' ','') in combined.replace(' ','').replace('\n',''),x['id']+' PDF readable heading')
 check(all(len(t)>200 for t in pt),'no empty or link-only PDF pages')
 def mt(name):return pd.read_csv(M/'tables'/name)
 emb=mt('embedding_metrics.csv').set_index('method');check(round(emb.loc['harmony_tissue','tissue_neighbor_purity'],3)==.504,'R01 integration metric')
 cov=mt('full_feature_score_coverage.csv').set_index('signature');check(cov.loc['NK_adaptive','n_used']==1 and cov.loc['NK_resting','n_used']==3,'R02 incomplete signature coverage')
 st=mt('selected_stability_resolution.csv');check(st.n_clusters_seed0.tolist()==[5,4,6] and st.resolution.tolist()==[.3,.3,.5],'R03 selected resolutions')
 gap=json.loads((M/'tables/input_missing_tasks.json').read_text());check(gap['velocity']['layers']==['counts'] and s[3]['status'].startswith('未运行'),'R04 velocity explicitly incomplete')
 legacy=mt('legacy_best_response_candidates_NOT_identity.csv').set_index('legacy_unit');check(legacy.loc['c5_bone_marrow','n_measured']==61 and legacy.loc['c5_bone_marrow','direction_concordance']==1,'R05 DEG response example')
 de=mt('de_status_original.csv');check(len(de)==31 and de.status.eq('computed').sum()==26 and de.strict_min20.sum()==21,'R06 computed and eligible comparisons')
 ppt=mt('PPT_pathway_validation_every_row.csv');f=ppt[ppt.sensitivity.eq('full')&ppt.strict_min20];check(ppt.ppt_row.nunique()==130 and f[f.status.eq('tested_current_membership')].ppt_row.nunique()==120,'R07 130 registered and 120 tested PPT records')
 b=mt('Blood_TableS3_all48_pathway_validation.csv');check(len(b)==2016 and b[b.in_Fig3E].paper_term.nunique()==17 and b[b.status.eq('tested')].paper_term.nunique()==44,'R08 main and supplementary coverage')
 disp=mt('KEGG_display_selection_all_significant.csv');check(len(disp)==145,'R09 significant contrast-term record count')
 check(len(list((M/'enrichment').glob('*__KEGG_Mouse_2019__GSEA.csv.gz')))==42,'R10 full and heat-excluded GSEA files')
 c=mt('IL4_CD94_counts_by_cluster_tissue_condition.csv');cc=c[c.tissue.eq('bone_marrow')&c.cluster.eq('C0')].set_index('condition');check(len(c)==72 and cc.loc['Ctrl','Il4_positive']==354 and cc.loc['T2','Il4_positive']==286,'R11 per-cluster counts')
 fr=mt('frequency_reconciliation_1.csv');check(len(fr)==27 and fr.signs_reconcile.all(),'R12 all denominator signs reconciled')
 pri=pd.read_csv(D/'tables/discovery_priority_all4_FDR05.csv');check(len(pri)==17 and pri.gene.nunique()==16 and (pri[['original_FDR','stable_FDR','QC_FDR','stable_QC_FDR']]<=.05).all().all(),'R13 strict priority count and four FDR checks')
 check((OUT/'logs/build.exit').read_text().strip()=='0' and (OUT/'logs/render.exit').read_text().strip()=='0','document build and rendering succeeded')
 result={'status':'passed','checks_passed':len(checks),'sections':13,'main_figures':12,'report_tables':13,'PDF_pages':len(pdf.pages),'visual_review':'pending','analysis_recomputed':False,'verification_script_SHA256':sha(Path(__file__)),'finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
 (OUT/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
