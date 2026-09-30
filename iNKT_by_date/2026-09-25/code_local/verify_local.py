from pathlib import Path
import json,hashlib,importlib.metadata,time
import anndata as ad
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'data/bio3'; manifest=json.loads((ROOT/'provenance/remote_manifest.json').read_text())
checked=[];bad=[]
for rel,meta in manifest['files'].items():
 p=DATA/rel
 if not p.is_file():bad.append({'file':rel,'issue':'missing'});continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 ok=p.stat().st_size==meta['bytes'] and h.hexdigest()==meta['sha256'];checked.append({'file':rel,'matches':ok})
 if not ok:bad.append({'file':rel,'issue':'hash or size mismatch'})
print('File checks',len(checked),'mismatches',len(bad),flush=True)
if bad:
 (ROOT/'provenance/local_verification.json').write_text(json.dumps({'pass':False,'issues':bad},indent=2));raise RuntimeError(bad[:5])
p=DATA/'output/iNKT_meeting_followup_20260905/objects/scored_base.h5ad';a=ad.read_h5ad(p,backed='r')
assert a.shape==(15532,10670)
assert a.obs_names.is_unique and a.var_names.is_unique
sample={str(k):int(v) for k,v in a.obs['sample'].value_counts().items()}
assert sum(sample.values())==15532
shape=list(a.shape);rawshape=list(a.raw.shape);a.file.close()
inputs=[]
for p in sorted((DATA/'input/iNKT/data').glob('*/sample_feature_bc_matrix')):
 assert all((p/n).is_file() for n in ['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']);inputs.append(p.parent.name)
assert len(inputs)==6
packages={k:importlib.metadata.version(k) for k in ['numpy','pandas','anndata','scanpy','scipy','matplotlib','networkx','scikit-learn','igraph','leidenalg','gseapy']}
report={'pass':True,'files_verified':len(checked),'mismatches':0,'snapshot_bytes':manifest['total_bytes'],'shape':shape,'raw_shape':rawshape,'sample_counts':sample,'raw_matrix_samples':inputs,'packages':packages,'verification':'offline local SHA256 + AnnData backed read + input matrix triplet completeness','no_remote_access_in_verification':True}
(ROOT/'provenance/local_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
