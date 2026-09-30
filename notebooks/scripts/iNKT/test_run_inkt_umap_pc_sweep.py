from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from run_inkt_umap_pc_sweep import (
    density_mass_threshold,
    label_knn_purity,
    normalize_metric_ranks,
    procrustes_concordance,
    select_recommended_pc,
    validate_pcs,
)


class UmapPcSweepTests(unittest.TestCase):
    def test_validate_pcs_requires_sorted_unique_values(self) -> None:
        self.assertEqual(validate_pcs([50, 100, 200]), (50, 100, 200))
        with self.assertRaisesRegex(ValueError, "sorted"):
            validate_pcs([100, 50])
        with self.assertRaisesRegex(ValueError, "unique"):
            validate_pcs([50, 50])

    def test_procrustes_concordance_is_rotation_and_scale_invariant(self) -> None:
        reference = np.array([[0.0, 0.0], [2.0, 0.0], [0.5, 1.5], [-1.0, 0.5]])
        rotation = np.array([[0.0, -1.0], [1.0, 0.0]])
        candidate = reference @ rotation * 3.5 + np.array([10.0, -4.0])
        self.assertAlmostEqual(procrustes_concordance(reference, candidate), 1.0, places=12)

    def test_label_purity_uses_neighbor_labels(self) -> None:
        labels = np.array(["a", "a", "b", "b"])
        neighbors = np.array([[1], [0], [3], [2]])
        self.assertEqual(label_knn_purity(labels, neighbors), 1.0)

    def test_density_threshold_retains_requested_high_density_mass(self) -> None:
        density = np.array([[9.0, 4.0], [2.0, 1.0]])
        threshold = density_mass_threshold(density, mass=0.50)
        self.assertEqual(threshold, 9.0)

    def test_selection_score_is_weighted_and_selects_best_candidate(self) -> None:
        metrics = pd.DataFrame(
            {
                "n_pcs": [50, 100, 200],
                "tissue_silhouette": [0.1, 0.4, 0.2],
                "leiden_silhouette": [0.1, 0.5, 0.2],
                "trustworthiness": [0.8, 0.9, 0.85],
                "mean_stability": [0.7, 0.8, 0.75],
                "legacy_centroid_concordance": [0.6, 0.9, 0.7],
                "program_knn_smoothness": [0.5, 0.8, 0.6],
            }
        )
        selected, weights = select_recommended_pc(metrics)
        self.assertEqual(selected, 100)
        self.assertAlmostEqual(sum(weights.values()), 1.0)
        self.assertEqual(int(metrics.loc[metrics["selection_score"].idxmax(), "n_pcs"]), 100)

    def test_equal_metric_ranks_are_neutral(self) -> None:
        metrics = pd.DataFrame({"x": [1.0, 1.0, 1.0]})
        ranks = normalize_metric_ranks(metrics, {"x": 1.0})
        self.assertTrue(np.allclose(ranks, 0.5))


if __name__ == "__main__":
    unittest.main()
