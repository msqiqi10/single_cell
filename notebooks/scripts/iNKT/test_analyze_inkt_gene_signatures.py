from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).with_name("analyze_inkt_gene_signatures.py")
SPEC = importlib.util.spec_from_file_location("analyze_inkt_gene_signatures", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class GeneSignatureTests(unittest.TestCase):
    def test_clean_gene_token_and_alias(self):
        self.assertEqual(MODULE.clean_gene_token(" Phgdh \u00a0"), "Phgdh")
        self.assertEqual(MODULE.GENE_ALIASES["SEPP1"], "Selenop")
        self.assertEqual(MODULE.GENE_ALIASES["ITGAB4"], "Itgb4")

    def test_describe_all_signature_families(self):
        cases = {
            "iNKT1": ("inkt1", "literature"),
            "My list iNKT2": ("inkt2", "inhouse"),
            "iNKT17 (Wang et al., 2022)": ("inkt17", "wang_2022"),
            "Circulatory markers (Wang et al": ("circulatory", "wang_2022"),
            "Direct TCR activation ": ("direct_tcr_activation", "curated"),
            "Residency markers (Wang et al.,": ("residency", "wang_2022"),
        }
        for sheet, expected in cases.items():
            with self.subTest(sheet=sheet):
                descriptor = MODULE.describe_signature_sheet(sheet)
                self.assertIsNotNone(descriptor)
                self.assertEqual((descriptor.subtype, descriptor.source), expected)

    def test_condition_effects(self):
        summary = pd.DataFrame(
            [
                {
                    "stratum_type": "global",
                    "tissue": None,
                    "cluster": None,
                    "set_id": "test",
                    "condition": "Ctrl",
                    "n_cells": 10,
                    "mean_score": 1.0,
                },
                {
                    "stratum_type": "global",
                    "tissue": None,
                    "cluster": None,
                    "set_id": "test",
                    "condition": "T2",
                    "n_cells": 12,
                    "mean_score": 1.75,
                },
            ]
        )
        effects = MODULE.condition_effects(summary, ["mean_score"])
        self.assertEqual(len(effects), 1)
        self.assertAlmostEqual(effects.iloc[0]["delta_mean_score_t2_minus_ctrl"], 0.75)
        self.assertEqual(effects.iloc[0]["n_ctrl"], 10)
        self.assertEqual(effects.iloc[0]["n_t2"], 12)

    def test_marker_summary_uses_position_indices(self):
        obs = pd.DataFrame(
            {
                "condition": ["Ctrl", "Ctrl", "T2", "T2"],
                "tissue": ["spleen"] * 4,
                "cluster": ["0", "0", "0", "1"],
            },
            index=["a", "b", "c", "d"],
        )
        log_expression = np.array([[0.0], [1.0], [2.0], [0.0]])
        normalized = np.expm1(log_expression)
        result = MODULE.marker_expression_summary(
            obs, log_expression, normalized, ["Xcl1"], "cluster"
        )
        global_ctrl = result.loc[
            result["stratum_type"].eq("global") & result["condition"].eq("Ctrl")
        ].iloc[0]
        self.assertEqual(global_ctrl["n_cells"], 2)
        self.assertAlmostEqual(global_ctrl["mean_log1p_expression"], 0.5)
        self.assertAlmostEqual(global_ctrl["fraction_expressing"], 0.5)

    def test_standardized_mean_z_scores_puts_unequal_sets_on_shared_scale(self):
        expression = np.array(
            [
                [0.0, 0.0, 1.0, 2.0],
                [1.0, 1.0, 2.0, 3.0],
                [2.0, 2.0, 3.0, 4.0],
                [3.0, 3.0, 4.0, 5.0],
            ]
        )
        scores, gene_stats, signature_stats = MODULE.standardized_mean_z_scores(
            expression,
            ["g1", "g2", "g3", "g4"],
            {
                "inkt1": ["g1"],
                "inkt2": ["g2", "g3"],
                "inkt17": ["g3", "g4"],
            },
        )
        self.assertEqual(list(scores.columns), ["inkt1", "inkt2", "inkt17"])
        for subtype in scores:
            self.assertAlmostEqual(float(scores[subtype].mean()), 0.0)
            self.assertAlmostEqual(float(scores[subtype].std(ddof=0)), 1.0)
        self.assertEqual(len(gene_stats), 4)
        self.assertEqual(signature_stats.set_index("subtype").loc["inkt2", "n_variable_genes"], 2)

    def test_bootstrap_cluster_calls_can_return_unclassified(self):
        metadata = pd.DataFrame(
            {
                "condition": ["Ctrl"] * 12,
                "sample": ["s1", "s1", "s2", "s2"] * 3,
                "cluster": ["0"] * 4 + ["1"] * 4 + ["2"] * 4,
            }
        )
        scores = pd.DataFrame(
            {
                "inkt1": [2.0] * 4 + [-0.2] * 4 + [1.0] * 4,
                "inkt2": [0.5] * 4 + [-0.3] * 4 + [1.0] * 4,
                "inkt17": [0.0] * 4 + [-0.4] * 4 + [0.0] * 4,
            }
        )
        result = MODULE.bootstrap_standardized_cluster_calls(
            metadata,
            scores,
            "cluster",
            n_bootstrap=20,
            seed=0,
            batch_size=5,
            strata_key="sample",
        ).set_index("cluster")
        self.assertEqual(result.loc["0", "supported_call"], "legacy-list inkt1-enriched")
        self.assertEqual(result.loc["1", "supported_call"], "unclassified")
        self.assertEqual(result.loc["2", "supported_call"], "mixed")
        self.assertEqual(result.loc["0", "bootstrap_strata"], "sample")


if __name__ == "__main__":
    unittest.main()
