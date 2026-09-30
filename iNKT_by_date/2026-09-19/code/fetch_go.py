"""Freeze matched GO ontology and mouse annotations, recording download provenance."""
from pathlib import Path
import requests, hashlib, json, datetime
out=Path(__file__).resolve().parents[1]/'sources'
urls={'go-basic.obo':'https://current.geneontology.org/ontology/go-basic.obo','mgi.gaf.gz':'https://current.geneontology.org/annotations/mgi.gaf.gz'}
manifest={}
for name,url in urls.items():
 p=out/name
 if not p.exists():
  print('Downloading',url,flush=True)
  with requests.get(url,stream=True,timeout=(30,120)) as r:
   r.raise_for_status()
   with p.with_suffix(p.suffix+'.part').open('wb') as f:
    for chunk in r.iter_content(1024*1024):f.write(chunk)
   p.with_suffix(p.suffix+'.part').rename(p)
 manifest[name]={'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 print(name,manifest[name],flush=True)
(out/'download_manifest.json').write_text(json.dumps(manifest,indent=2))
