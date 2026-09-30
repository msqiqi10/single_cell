"""Meaningful invariants for condition-blind QC matching and discovery inputs."""
import unittest
import numpy as np,pandas as pd
from run_inkt_discovery import qc_strata,matched_indices,balance,family

class DiscoveryMatchingTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(18);n=160
        self.obs=pd.DataFrame({'condition':['Ctrl']*100+['T2']*60,'log_total_counts':rng.normal(8.5,.2,n),'log_n_genes':rng.normal(7.5,.1,n),'pct_counts_mt':rng.uniform(.5,4,n)})
    def test_quantiles_do_not_use_condition(self):
        flipped=self.obs.copy();flipped.condition=flipped.condition.map({'Ctrl':'T2','T2':'Ctrl'})
        np.testing.assert_array_equal(qc_strata(self.obs),qc_strata(flipped))
    def test_exact_condition_balance_within_selected_strata(self):
        ids,t=matched_indices(self.obs);self.assertEqual(len(ids),len(set(ids)));bins=qc_strata(self.obs)
        for b in set(bins[ids]):
            cond=self.obs.iloc[ids[bins[ids]==b]].condition.value_counts();self.assertEqual(cond['Ctrl'],cond['T2'])
        self.assertEqual(len(ids),2*t.n_each_selected.sum())
    def test_constant_covariates_are_supported(self):
        obs=pd.DataFrame({'condition':['Ctrl']*7+['T2']*3,'log_total_counts':[1]*10,'log_n_genes':[1]*10,'pct_counts_mt':[1]*10})
        ids,_=matched_indices(obs);self.assertEqual(len(ids),6);self.assertTrue(all(x==0 for x in balance(obs.iloc[ids]).values()))
    def test_matching_seed_is_reproducible(self):
        a,_=matched_indices(self.obs,41);b,_=matched_indices(self.obs,41);np.testing.assert_array_equal(a,b)
    def test_missing_qc_is_not_silently_assigned(self):
        obs=self.obs.copy();obs.loc[0,'pct_counts_mt']=np.nan
        with self.assertRaises(ValueError):matched_indices(obs)
    def test_dominant_driver_classes_are_explicit(self):
        self.assertEqual(family('Hspa8'),'heat_shock');self.assertEqual(family('Rpl34'),'ribosome');self.assertEqual(family('Ndufa13'),'respiratory_chain_or_mt');self.assertEqual(family('Il7r'),'other')

if __name__=='__main__':unittest.main()
