"""Offline reproduction of BM C0 GO ORA from frozen DE and GO membership.

Reuses upstream DE and ontology annotations; independently recomputes ORA,
including zero-hit terms in each multiple-testing family. Does not rerun DE.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/bio3/iNKT_by_date/2026-09-19'
OUT = ROOT / 'runs/BM_C0_GO_local_validation'
OUT.mkdir(parents=True, exist_ok=True)
catalog = pd.read_csv(SOURCE / 'results/tables/GO_term_catalog.csv').set_index('go_id')
membership = {}
for line in (SOURCE / 'sources/GO_mouse_measured.gmt').read_text().splitlines():
    go_id, name, *genes = line.split('\t')
    membership[go_id] = set(genes)
de = pd.read_csv(SOURCE / 'results/de/cluster_tissue__C0__bone_marrow.csv.gz')
expressed = set(de.loc[(de.pct_expressing_tumor > 0) | (de.pct_expressing_control > 0), 'gene'])
results, coverage = [], []
for direction in ['T2_up', 'Ctrl_up']:
    effect = de.logfoldchanges >= .25 if direction == 'T2_up' else de.logfoldchanges <= -.25
    selected = set(de.loc[(de.pvals_adj <= .05) & effect, 'gene'])
    for namespace in ['BP', 'MF', 'CC']:
        ids = catalog.index[catalog.namespace == namespace]
        universe = expressed & set().union(*(membership[g] for g in ids))
        query = selected & universe
        rows = []
        for go_id in ids:
            members = membership[go_id] & universe
            if not 5 <= len(members) <= 500:
                continue
            hits = members & query
            rows.append(dict(go_id=go_id, term=catalog.loc[go_id, 'term'],
                             M=len(universe), K=len(members), n=len(query), k=len(hits),
                             genes=';'.join(sorted(hits))))
        frame = pd.DataFrame(rows)
        frame['pvalue'] = hypergeom.sf(frame.k - 1, len(universe), frame.K, len(query)) if query else 1.
        frame['q_family'] = false_discovery_control(frame.pvalue.to_numpy(), method='bh')
        frame['direction'], frame['namespace'] = direction, namespace
        results.append(frame)
        coverage.append(dict(direction=direction, namespace=namespace,
                             selected_genes=len(selected), annotated_query_genes=len(query),
                             tested_terms=len(frame), significant_terms=int((frame.q_family <= .05).sum())))
result = pd.concat(results, ignore_index=True)
old = pd.read_csv(SOURCE / 'results/tables/GO_ORA_all.csv.gz', keep_default_na=False)
old = old[(old.tissue == 'bone_marrow') & (old.cluster == 'C0') &
          (old.analysis == 'condition') & (old.definition == 'original')]
keys = ['go_id', 'direction', 'namespace']
joined = result.merge(old, on=keys, suffixes=('_local', '_reference'), validate='one_to_one', how='outer', indicator=True)
assert joined['_merge'].eq('both').all()
for column in ['M', 'K', 'n', 'k', 'genes']:
    assert joined[column + '_local'].eq(joined[column + '_reference']).all(), column
errors = {}
for column in ['pvalue', 'q_family']:
    local, reference = joined[column + '_local'], joined[column + '_reference']
    np.testing.assert_allclose(local, reference, rtol=1e-10, atol=1e-300)
    errors[column] = float(np.max(np.abs(local-reference)))
result.to_csv(OUT / 'GO_ORA_recomputed.csv.gz', index=False)
pd.DataFrame(coverage).to_csv(OUT / 'coverage.csv', index=False)
report = dict(passed=True, comparison='bone marrow C0, T2 versus Ctrl, original clusters',
              tests_recomputed=len(result), exact_match_columns=['M', 'K', 'n', 'k', 'genes'],
              max_absolute_errors=errors, coverage=coverage,
              scope='Recomputed ORA from existing full DE and frozen mouse GO membership; did not rerun DE or ontology parsing. q_analysis_global requires other contrasts and is not recomputed here.',
              offline=True)
(OUT / 'validation.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
