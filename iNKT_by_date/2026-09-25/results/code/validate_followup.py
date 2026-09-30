"""Independent invariants for source reuse, direction, coverage, IDs and holdouts."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/meeting_followup';DATA=ROOT/'data/bio3';checks={}
de=pd.read_csv(OUT/'tables/all_20_condition_DE.csv.gz');assert len(de)==213400 and de.unit_id.nunique()==20
assert not de.duplicated(['unit_id','gene']).any();checks['full_DE_20x10670']=True
all_de=pd.read_csv(OUT/'tables/all_26_existing_DE_with_eligibility.csv.gz');assert len(all_de)==277420
checks['all_26_existing_DE_including_low_cell_and_global']=True
go=pd.read_csv(OUT/'tables/all_original_condition_GO.csv.gz',keep_default_na=False)
for uid,d in de.groupby('unit_id'):
    source=pd.read_csv(DATA/'iNKT_by_date/2026-09-19/results/de'/f'{uid}.csv.gz')
    pd.testing.assert_frame_equal(d[source.columns].reset_index(drop=True),source,check_exact=False,rtol=1e-12,atol=1e-300)
checks['DE_matches_source']=True
cov=pd.read_csv(OUT/'tables/full_candidate_GO_coverage.csv')
assert len(cov)==120
for row in cov.itertuples():
    d=de[de.unit_id==row.unit_id];mask=d.logfoldchanges>=.25 if row.direction=='T2_up' else d.logfoldchanges<=-.25
    query=set(d.loc[(d.pvals_adj<=.05)&mask,'gene']);assert len(query)==row.n_selected
    f=go[(go.unit_id==row.unit_id)&(go.direction==row.direction)&(go.namespace==row.namespace)]
    covered={g for text in f.loc[f.q_family<=.05,'genes'] for g in text.split(';') if g}
    assert covered<=query and len(covered)==row.n_in_significant_terms
    assert row.n_without_branch_annotation+row.n_annotated==row.n_selected
    assert row.n_not_in_significant_terms+row.n_in_significant_terms==row.n_selected
    np.testing.assert_allclose(false_discovery_control(f.pvalue.to_numpy(),method='bh'),f.q_family,rtol=1e-10,atol=1e-300)
checks['coverage_direction_conservation_and_family_BH']=True
nodes=pd.read_csv(OUT/'network/annotated_nodes.csv');ids=nodes.GOID.to_numpy();assert len(ids)==97 and len(set(ids))==97
expected=set(go.loc[(go.namespace=='BP')&(go.q_family<=.05),'go_id']);assert set(ids)==expected
bundle=np.load(OUT/'network/overlap_matrices.npz');assert np.array_equal(ids,bundle['ID'])
edges=pd.read_csv(OUT/'network/mouse_GO_edges.tsv',sep='\t');assert len(edges)==577
index={gid:i for i,gid in enumerate(ids)};expected_edges=set()
for row in edges.itertuples():
    assert row.GS_A_ID!=row.GS_B_ID and row.JACCARD>=.25 and row.n_shared>=3
    expected_edges.add(tuple(sorted((index[row.GS_A_ID],index[row.GS_B_ID]))))
assert len(expected_edges)==577
assert int((bundle['adjacency']>0).sum())==2*len(edges)
for seed in [42,43,44]:
    split=np.load(OUT/'models'/f'gnn_seed{seed}.csv.splits.npz');assert np.array_equal(split['ID'],ids)
    groups={key:{tuple(x) for x in split[key]} for key in ['train_positive','val_positive','test_positive','val_negative','test_negative']}
    pos=groups['train_positive']|groups['val_positive']|groups['test_positive'];assert pos==expected_edges
    assert sum(len(groups[k]) for k in ['train_positive','val_positive','test_positive'])==len(pos)
    assert not ((groups['val_negative']|groups['test_negative']) & pos)
    assert not (groups['val_negative']&groups['test_negative'])
    for key in groups:assert len(groups[key])==len(split[key])
checks['network_node_edge_and_holdout_alignment']=True
link=pd.read_csv(OUT/'tables/module_gene_evidence.csv.gz');assert np.isfinite(link.log2FC).all()
lookup=de.set_index(['unit_id','gene'])
for row in link.itertuples():
    v=lookup.loc[row.unit_id,row.gene];assert np.isclose(v.logfoldchanges,row.log2FC)
    assert v.pvals_adj<=.05
    assert (v.logfoldchanges>=.25) if row.direction=='T2_up' else (v.logfoldchanges<=-.25)
checks['module_to_gene_effect_and_direction']=True
marker=json.loads((OUT/'marker_validation.json').read_text());assert marker['passed'] and marker['shape']==[15532,32285]
raw=pd.read_csv(OUT/'tables/full_feature_marker_by_group.csv');glob=pd.read_csv(OUT/'tables/full_feature_marker_status.csv').set_index('gene')
for gene,group in raw.groupby('gene'):
    assert group.n_cells.sum()==15532
    assert group.n_detected.sum()==glob.loc[gene,'n_detected']
    assert group.total_UMI.sum()==glob.loc[gene,'total_UMI']
    assert ((group.n_detected>=0)&(group.n_detected<=group.n_cells)).all()
checks['full_feature_marker_counts_and_denominators']=True
files={str(p.relative_to(OUT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='validation.json' and 'objects' not in p.parts}
report=dict(passed=all(checks.values()),checks=checks,source_hash_manifest='../../provenance/remote_manifest.json',files=files,
    external_replication=False,verified='Computational invariants and reproducibility, not biological validation')
(OUT/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps({'passed':report['passed'],'checks':checks,'files_hashed':len(files)},indent=2))
