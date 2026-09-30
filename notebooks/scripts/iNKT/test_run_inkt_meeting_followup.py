"""Regression checks for inference direction, gene matching and sensitivity controls."""
import unittest
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from run_inkt_meeting_followup import ora_table, centroid_summary, bh_adjust
from run_inkt_meeting_pathways import clean_rank, libraries, normterm

class MeetingScientificChecks(unittest.TestCase):
    def setUp(self):
        self.d=pd.DataFrame({'gene':['A','B','C','D','E','F','Hspa8','Dnaja1'],'scores':[2,1,-2,-1,0,np.nan,3,-3],'logfoldchanges':[.60,.55,-.60,-.55,0,0,2,-2],'pvals':[.01]*8,'pvals_adj':[.02]*8})
    def test_paper_linear_fold_is_not_log2_one(self):
        r=ora_table(self.d,{'set':set(self.d.gene)},set(self.d.gene),'paper').set_index('direction')
        self.assertEqual(set(r.loc['Tumor','overlap_genes'].split(';')),{'A','Hspa8'})
        self.assertEqual(set(r.loc['Control','overlap_genes'].split(';')),{'C','Dnaja1'})
    def test_heat_removal_changes_universe_and_membership(self):
        d,a=clean_rank(self.d,True);self.assertEqual(len(d),6);self.assertTrue(np.isfinite(d.scores).all())
        r=ora_table(d,{'set':{'A','B','E','Hspa8','Dnaja1'}},set(d.gene),'robust').set_index('direction')
        self.assertEqual(r.loc['Tumor','universe_size'],6);self.assertEqual(r.loc['Tumor','pathway_size'],3)
        self.assertAlmostEqual(r.loc['Tumor','pvalue'],hypergeom.sf(1,6,3,2))
        self.assertNotIn('Hspa8',r.loc['Tumor','overlap_genes'])
    def test_legacy_union_has_no_activation_direction(self):
        r=ora_table(self.d,{'set':set(self.d.gene)},set(self.d.gene),'legacy')
        self.assertEqual(r.direction.tolist(),['Mixed_legacy_DEG_union']);self.assertEqual(r.iloc[0].overlap,2)
    def test_all_hypotheses_including_zero_hits_enter_BH(self):
        r=ora_table(self.d,{'hit':{'A','B','Hspa8'},'none':{'E','F','D'}},set(self.d.gene),'paper')
        t=r[r.direction.eq('Tumor')];self.assertEqual(len(t),2);self.assertEqual(t.loc[t.term.eq('none'),'pvalue'].iloc[0],1)
        self.assertTrue(np.allclose(t.fdr,bh_adjust(t.pvalue.to_numpy())))
    def test_empty_centroid_is_missing_not_zero(self):
        r=centroid_summary(np.array([[2,4],[4,8]]),np.array(['Ctrl','T2']),np.array([True,False]))
        self.assertEqual(r['n_T2'],0);self.assertTrue(np.isnan(r['T2_x']));self.assertTrue(np.isnan(r['median_displacement']))
    def test_duplicate_ranking_is_rejected(self):
        with self.assertRaises(ValueError):clean_rank(pd.concat([self.d,self.d.iloc[:1]]))
    def test_frozen_mouse_gmt_case_is_restored(self):
        sets=libraries()['KEGG_Mouse_2019']
        self.assertIn('Tnf',sets['TNF signaling pathway']);self.assertGreater(len(sets['TNF signaling pathway']),50)
        self.assertIn('Hspa8',sets['Antigen processing and presentation']);self.assertGreater(len(sets['Ribosome']),100)
    def test_term_matching_CAMS_alias(self):
        self.assertEqual(normterm('Cell adhesion molecules (CAMs)'),normterm('Cell adhesion molecules'))

if __name__=='__main__':unittest.main()
