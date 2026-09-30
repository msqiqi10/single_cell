"""Fetch current full Reactome memberships for reactions named in the legacy PPT."""
from pathlib import Path
import concurrent.futures,json,re,time,xml.etree.ElementTree as ET
import pandas as pd
import requests
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/iNKT_meeting_followup_20260905';SRC=OUT/'sources';CACHE=SRC/'reactome_api';CACHE.mkdir(exist_ok=True)
legacy=pd.read_csv(ROOT/'output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_pathway_tables_extracted.csv').fillna('')
# Source DOI/XML was independently verified with Europe PMC.
root=ET.parse(SRC/'Blood2025_fulltext.xml').getroot()
for el in root.iter():
 if el.tag in ['supplementary-material','media']:
  print('SUPPLEMENT',ET.tostring(el,encoding='unicode')[:1300],flush=True)

def norm(s):return re.sub(r'[^a-z0-9]','',re.sub('<[^>]+>','',str(s)).lower())
def one(name):
 fname=CACHE/(re.sub(r'[^A-Za-z0-9]+','_',name)[:160]+'.json')
 if fname.exists():return json.loads(fname.read_text())
 result={'legacy_name':name,'status':'unresolved'}
 try:
  q=name.replace('...','').strip();r=requests.get('https://reactome.org/ContentService/search/query',params={'query':q,'species':'Homo sapiens','rows':100},timeout=30);r.raise_for_status();data=r.json()
  entries=[e for g in data.get('results',[]) for e in g.get('entries',[])]
  entries=[e for e in entries if str(e.get('stId','')).startswith('R-HSA-') and e.get('type') in ['Reaction','Pathway']]
  matches=[e for e in entries if norm(e.get('name',''))==norm(q) or (name.endswith('...') and norm(e.get('name','')).startswith(norm(q)))]
  matches={e['stId']:e for e in matches}
  if len(matches)!=1:result.update(status='no_unique_name_match',candidates=[{'id':e.get('stId'),'name':re.sub('<[^>]+>','',e.get('name',''))} for e in entries[:8]])
  else:
   rid,e=next(iter(matches.items()));r=requests.get(f'https://reactome.org/ContentService/data/participants/{rid}/referenceEntities',timeout=30);r.raise_for_status();refs=r.json()
   genes=sorted({g for ref in refs for g in ref.get('geneName',[])})
   result.update(status='full_current_membership_retrieved' if genes else 'no_gene_names',stable_id=rid,current_name=re.sub('<[^>]+>','',e['name']),human_genes=genes,n_human=len(genes),reference_entities=refs,version_note='Current Reactome retrieval 2026-09-05, not original Reactome_2021')
 except Exception as e:result.update(status='request_failed',error=str(e))
 fname.write_text(json.dumps(result,indent=2));print(name,result['status'],result.get('n_human'),flush=True);return result
names=sorted(set(legacy.loc[legacy.SOURCE.str.startswith('Reactome'),'NAME']))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,names))
h2m=json.loads((SRC/'human_to_measured_mouse.json').read_text())
sets={};rows=[]
for r in results:
 genes=sorted({h2m[g] for g in r.get('human_genes',[]) if g in h2m})
 if r['status']=='full_current_membership_retrieved':sets[r['legacy_name']]=genes
 rows.append({k:v for k,v in r.items() if k not in ['reference_entities','human_genes','candidates']}|{'human_genes':';'.join(r.get('human_genes',[])),'mouse_genes':';'.join(genes),'n_measured_mouse':len(genes)})
with (SRC/'Legacy_named_Reactome_current_mapped_mouse.gmt').open('w') as f:
 for name,genes in sets.items():f.write(name+'\tcurrent_full_membership\t'+'\t'.join(genes)+'\n')
pd.DataFrame(rows).to_csv(OUT/'tables/legacy_reactome_full_membership_audit.csv',index=False)
print('COMPLETE',len(results),len(sets),flush=True)
