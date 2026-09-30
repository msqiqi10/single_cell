"""Mouse GO adapter for GOLDEN's description/adjacency fusion, with explicit IDs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sentence_transformers import SentenceTransformer
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, adjusted_rand_score

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup'
nodes=pd.read_csv(OUT/'network/mouse_GO_nodes.csv');bundle=np.load(OUT/'network/overlap_matrices.npz')
ids=nodes.GOID.to_numpy();assert np.array_equal(ids,bundle['ID'])
model=SentenceTransformer(str(ROOT/'models/all-MiniLM-L6-v2'),device='cpu',local_files_only=True)
semantic=model.encode(nodes.DESCRIPTION.tolist(),normalize_embeddings=True,show_progress_bar=False)
assert semantic.shape==(len(ids),384) and np.isfinite(semantic).all()
np.savez_compressed(OUT/'network/semantic_embeddings.npz',ID=ids.astype(str),embeddings=semantic)
adj=bundle['adjacency'];scores=[];solutions={}
for alpha in [0.,.25,.5,.75,1.]:
    features=np.concatenate([alpha*semantic,(1-alpha)*adj],axis=1)
    for k in range(2,13):
        labels=AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(features)
        score=silhouette_score(features,labels,metric='euclidean')
        scores.append(dict(alpha=alpha,k=k,silhouette=score,min_size=int(np.bincount(labels).min())))
        solutions[alpha,k]=labels
scores=pd.DataFrame(scores);scores.to_csv(OUT/'network/fusion_k_sensitivity.csv',index=False)
chosen=scores[scores.alpha==.5].sort_values(['silhouette','k'],ascending=[False,True]).iloc[0]
k=int(chosen.k);reference=solutions[.5,k]
# Canonicalize module labels by smallest GO ID, making labels independent of sklearn integers.
order=sorted(set(reference),key=lambda label:min(ids[reference==label]));mapping={label:f'F{i+1:02d}' for i,label in enumerate(order)}
labels=[mapping[x] for x in reference]
nodes['fusion_module']=labels;nodes.to_csv(OUT/'network/fusion_nodes.csv',index=False)
pd.DataFrame([dict(alpha=alpha,k=k,ARI_to_alpha05=float(adjusted_rand_score(reference,solutions[alpha,k]))) for alpha in [0.,.25,.5,.75,1.]]).to_csv(OUT/'network/fusion_alpha_stability.csv',index=False)
manifest=dict(method='Local mouse GO adaptation of GOLDEN description plus adjacency fusion',
    upstream_commit='94f9d367483ccb4872630ebaa3b611b11220d9a6',
    formula='concat(alpha * L2(description embedding), (1-alpha) * raw weighted adjacency)',
    alpha=.5,clustering='Ward agglomerative, Euclidean',k=k,k_selection='maximum silhouette for k=2..12 at fixed alpha=0.5; descriptive exploratory selection, no biological label tuning',
    semantic_model=json.loads((ROOT/'provenance/semantic_model.json').read_text()),nodes=len(ids),embedding_dimensions=384,
    adaptations=['mouse frozen GO definitions and measured gene members replace unavailable PAGER metadata',
                 'Jaccard>=0.25 and shared>=3 replaces unavailable m-type network; not mathematically identical to original m-type',
                 'explicit common GO order, retained isolates, correct label-to-node mapping',
                 'k sweep and alpha sensitivity recorded, rather than externally assumed number of modules'],
    not_claimed='Unmodified official GOLDEN reproduction, gold-standard biological modules, pathway activation, or causal edges')
(OUT/'fusion_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
