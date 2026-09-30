"""Meaningful numerical tests for enrichment and QC matching."""
import tempfile, unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from analyze import bh, ora, read_obo, qc_match
class AnalysisTests(unittest.TestCase):
 def test_bh_matches_known_values(self):
  np.testing.assert_allclose(bh([.01,.04,.03,.002]),[.02,.04,.04,.008])
 def test_ora_uses_measured_universe_and_tests_zero_hits(self):
  genes=set('abcdefghij');sets={'T1':set('abcdeZ'),'T2':set('fghijZ')};terms={t:{'name':t} for t in sets}
  d=ora({'a','b','c','Z'},genes,sets,terms).set_index('go_id')
  self.assertEqual(d.loc['T1','n'],3);self.assertEqual(d.loc['T1','K'],5)
  self.assertAlmostEqual(d.loc['T1','pvalue'],hypergeom.sf(2,10,5,3))
  self.assertEqual(d.loc['T2','pvalue'],1);self.assertEqual(len(d),2)
  self.assertAlmostEqual(d.loc['T1','q_family'],min(1,2*d.loc['T1','pvalue']))
 def test_empty_query_has_no_false_positive(self):
  d=ora(set(),set('abcdefghij'),{'T':set('abcde')},{'T':{'name':'test'}})
  self.assertEqual(d.pvalue.iloc[0],1);self.assertEqual(d.q_family.iloc[0],1)
 def test_ontology_avoids_regulates_and_obsolete(self):
  txt='data-version: releases/test\n[Term]\nid: GO:1\nname: one\nnamespace: biological_process\nis_a: GO:2 ! two\nrelationship: part_of GO:3 ! three\nrelationship: regulates GO:4 ! four\nalt_id: GO:5\n[Term]\nid: GO:6\nname: old\nis_obsolete: true\n'
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'go.obo';p.write_text(txt);terms,alt,version=read_obo(p)
  self.assertEqual(terms['GO:1']['parents'],['GO:2','GO:3']);self.assertNotIn('GO:6',terms);self.assertEqual(alt['GO:5'],'GO:1')
 def test_qc_match_balanced_unique_reproducible(self):
  obs=pd.DataFrame({'condition':['Ctrl']*50+['T2']*70,'total_counts':np.tile(np.arange(20),6),'n_genes_by_counts':np.tile(np.arange(20),6),'pct_counts_mt':np.ones(120)})
  ix=qc_match(obs,42);self.assertEqual(len(set(ix)),len(ix));c=obs.iloc[ix].condition.value_counts();self.assertEqual(c['Ctrl'],c['T2']);np.testing.assert_array_equal(ix,qc_match(obs,42))
if __name__=='__main__':unittest.main(verbosity=2)
