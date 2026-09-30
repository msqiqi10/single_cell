from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from run_inkt_c5_paper_followup import (
    C5_SUBCLUSTER_KEY,
    FIGURE_F_MAIN_SET_ORDER,
    GSEA_CLUSTERS,
    GSEA_HEATMAP_MODULES,
    ORIGINAL_CLUSTER_KEY,
    PAPER_VENN_THRESHOLD,
    REFINED_CLUSTER_KEY,
    ROBUST_THRESHOLD,
    RUN_DATE,
    SAMPLE_DISPLAY_ORDER,
    audit_raw_ncam1,
    bh_adjust,
    build_cluster_proportion_tables,
    build_gsea_bundle_receipt,
    build_figure_f_main_membership,
    build_gsea_heatmap_matrices,
    build_gsea_ranking,
    exact_intersection_regions,
    gsea_computation_identity,
    load_validated_gsea_bundle,
    install_refined_labels,
    marker_complementarity_row,
    odds_ratio_haldane,
    run_requested_de,
    select_gsea_heatmap_pathways,
    select_gsea_heatmap_modules,
    select_pathways_for_plot,
    select_two_cluster_resolution,
    validate_full_gene_input,
    validate_h5ad_roundtrip,
    venn_region_counts,
)


def make_synthetic_full_adata(
    var_names: tuple[str, ...] = ("Il4", "Klrd1", "Other"),
) -> ad.AnnData:
    obs_names = ["cell_c0", "cell_c5_a", "cell_c6", "cell_c5_b"]
    obs = pd.DataFrame(
        {
            "sample": ["Ctrl_BM", "Ctrl_BM", "T2_Spleen", "T2_Spleen"],
            "condition": ["Ctrl", "Ctrl", "T2", "T2"],
            "tissue": ["bone_marrow", "bone_marrow", "spleen", "spleen"],
            ORIGINAL_CLUSTER_KEY: ["0", "5", "6", "5"],
        },
        index=obs_names,
    )
    values = np.arange(
        1,
        len(obs_names) * len(var_names) + 1,
        dtype=np.float32,
    ).reshape(len(obs_names), len(var_names))
    adata = ad.AnnData(
        X=values,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(var_names, name="gene")),
    )
    adata.layers["counts"] = values.astype(np.int32)
    adata.raw = adata.copy()
    adata.obsm["X_pca"] = np.arange(12, dtype=np.float32).reshape(4, 3)
    adata.obsm["X_umap"] = np.array(
        [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]],
        dtype=np.float32,
    )
    return adata


def make_scrambled_c5(adata: ad.AnnData) -> ad.AnnData:
    c5 = adata[["cell_c5_b", "cell_c5_a"]].copy()
    c5.obs[C5_SUBCLUSTER_KEY] = pd.Categorical(
        ["C5-2", "C5-1"],
        categories=["C5-1", "C5-2"],
        ordered=True,
    )
    c5.obsm["X_umap_original_full"] = np.asarray(c5.obsm["X_umap"]).copy()
    c5.obsm["X_umap"] = np.array(
        [[20.0, 21.0], [10.0, 11.0]],
        dtype=np.float32,
    )
    return c5


def write_raw_sample(root: Path, sample: str, audit_status: str) -> None:
    matrix_dir = root / sample / "sample_feature_bc_matrix"
    matrix_dir.mkdir(parents=True)
    if audit_status == "feature_absent":
        feature_rows = (
            "ENSMUSG_OTHER1\tOther1\tGene Expression\n"
            "ENSMUSG_OTHER2\tOther2\tGene Expression\n"
            "ENSMUSG_OTHER3\tOther3\tGene Expression\n"
        )
    else:
        feature_rows = (
            "ENSMUSG_OTHER1\tOther1\tGene Expression\n"
            "ENSMUSG_NCAM1\tNcam1\tGene Expression\n"
            "ENSMUSG_OTHER3\tOther3\tGene Expression\n"
        )
    with gzip.open(matrix_dir / "features.tsv.gz", "wt") as handle:
        handle.write(feature_rows)

    entries = "2 1 3\n" if audit_status == "present_with_counts" else ""
    n_entries = 1 if entries else 0
    matrix_text = (
        "%%MatrixMarket matrix coordinate integer general\n"
        "% synthetic regression fixture\n"
        f"3 2 {n_entries}\n"
        f"{entries}"
    )
    with gzip.open(matrix_dir / "matrix.mtx.gz", "wt") as handle:
        handle.write(matrix_text)


class InktC5PaperFollowupTests(unittest.TestCase):
    def test_bh_adjust_preserves_order_and_caps_at_one(self) -> None:
        adjusted = bh_adjust([0.01, 0.04, 0.03, np.nan])
        self.assertTrue(np.allclose(adjusted[:3], [0.03, 0.04, 0.04]))
        self.assertTrue(np.isnan(adjusted[3]))

    def test_resolution_selection_requires_two_clusters_for_every_seed(self) -> None:
        metrics = pd.DataFrame(
            {
                "resolution": [0.08, 0.10, 0.15],
                "all_seeds_two_clusters": [True, True, False],
                "mean_pairwise_ari": [0.988, 0.992, 1.0],
                "min_pairwise_ari": [0.967, 0.967, 1.0],
            }
        )
        self.assertEqual(select_two_cluster_resolution(metrics), 0.10)

    def test_denominator_reconciliation_uses_tissue_specific_baseline(self) -> None:
        rows: list[dict[str, str]] = []
        for _ in range(7):
            rows.append({"tissue": "bone_marrow", "condition": "Ctrl", "cluster": "C0"})
        for _ in range(3):
            rows.append({"tissue": "bone_marrow", "condition": "Ctrl", "cluster": "C1"})
        for _ in range(4):
            rows.append({"tissue": "bone_marrow", "condition": "T2", "cluster": "C0"})
        rows.append({"tissue": "bone_marrow", "condition": "T2", "cluster": "C1"})
        for tissue in ("thymus", "spleen"):
            for condition in ("Ctrl", "T2"):
                rows.append({"tissue": tissue, "condition": condition, "cluster": "C0"})
                rows.append({"tissue": tissue, "condition": condition, "cluster": "C1"})
        frame = pd.DataFrame(rows)
        table, delta = build_cluster_proportion_tables(
            frame,
            frame["cluster"],
            ["C0", "C1"],
        )
        bm_c0 = delta[
            delta["tissue"].eq("bone_marrow") & delta["cluster"].eq("C0")
        ].iloc[0]
        self.assertAlmostEqual(bm_c0["tumor_share_zero_delta_baseline_pct"], 100 / 3)
        self.assertAlmostEqual(bm_c0["tumor_share_within_tissue_cluster_pct"], 400 / 11)
        self.assertAlmostEqual(bm_c0["tumor_minus_control_percentage_points"], 10.0)
        self.assertLess(bm_c0["tumor_share_within_tissue_cluster_pct"], 50.0)
        self.assertTrue(bool(bm_c0["signs_reconcile"]))
        self.assertEqual(int(table["count"].sum()), len(frame))

    def test_venn_regions_are_mutually_exclusive(self) -> None:
        sets = {
            "bone_marrow": {"a", "ab", "ac", "abc"},
            "spleen": {"b", "ab", "bc", "abc"},
            "thymus": {"c", "ac", "bc", "abc"},
        }
        self.assertEqual(
            venn_region_counts(sets),
            {
                "a_only": 1,
                "b_only": 1,
                "c_only": 1,
                "ab_only": 1,
                "ac_only": 1,
                "bc_only": 1,
                "abc": 1,
            },
        )

    def test_four_set_figure_f_membership_is_exact_and_union_conserving(self) -> None:
        expected_sets = {
            "C5-1__bone_marrow": {"a", "ab", "all"},
            "C5-1__spleen": {"b", "ab", "all"},
            "C5-2__bone_marrow": {"c", "all"},
            "C5-2__spleen": {"d", "all"},
        }
        universe = sorted(set().union(*expected_sets.values()))
        rows: list[dict[str, object]] = []
        for set_id in FIGURE_F_MAIN_SET_ORDER:
            cluster, tissue = set_id.split("__", maxsplit=1)
            for gene in universe:
                selected = gene in expected_sets[set_id]
                rows.append(
                    {
                        "unit_id": set_id,
                        "cluster": cluster,
                        "tissue": tissue,
                        "gene": gene,
                        "paper_venn_tumor_up": selected,
                        "robust_tumor_up": selected,
                    }
                )
        de = pd.DataFrame(rows)

        membership, summary, observed_sets = build_figure_f_main_membership(
            de,
            PAPER_VENN_THRESHOLD,
        )

        self.assertEqual(observed_sets, expected_sets)
        self.assertEqual(set(membership["gene"]), set(universe))
        indexed = summary.set_index("metric")
        self.assertEqual(int(indexed.loc["all_four", "n_genes"]), 1)
        self.assertEqual(str(indexed.loc["all_four", "genes"]), "all")
        ab_metric = (
            "intersection_exact__C5-1__bone_marrow&C5-1__spleen"
        )
        self.assertEqual(int(indexed.loc[ab_metric, "n_genes"]), 1)
        self.assertEqual(str(indexed.loc[ab_metric, "genes"]), "ab")

        regions = exact_intersection_regions(observed_sets, FIGURE_F_MAIN_SET_ORDER)
        region_union = set().union(*regions.values())
        self.assertEqual(region_union, set(universe))
        self.assertEqual(sum(map(len, regions.values())), len(universe))

    def test_four_set_figure_f_rejects_missing_de_unit(self) -> None:
        incomplete = pd.DataFrame(
            {
                "unit_id": ["C5-1__bone_marrow"],
                "cluster": ["C5-1"],
                "tissue": ["bone_marrow"],
                "gene": ["A"],
                "paper_venn_tumor_up": [True],
                "robust_tumor_up": [True],
            }
        )
        with self.assertRaisesRegex(RuntimeError, "Missing DE rows"):
            build_figure_f_main_membership(incomplete, PAPER_VENN_THRESHOLD)

    def test_complementarity_phi_is_negative_for_mutually_exclusive_detection(self) -> None:
        il4_counts = np.array([1, 1, 0, 0], dtype=float)
        klrd1_counts = np.array([0, 0, 1, 1], dtype=float)
        row = marker_complementarity_row(
            il4_counts,
            il4_counts,
            klrd1_counts,
            klrd1_counts,
            np.ones(4, dtype=bool),
            scope="test",
        )
        self.assertEqual(row["both_positive"], 0)
        self.assertEqual(row["il4_only"], 2)
        self.assertEqual(row["klrd1_only"], 2)
        self.assertAlmostEqual(row["detection_phi"], -1.0)

    def test_haldane_odds_ratio_is_finite_with_zero_cells(self) -> None:
        value = odds_ratio_haldane(
            overlap=0,
            selected_size=5,
            pathway_size=10,
            universe_size=100,
        )
        self.assertTrue(np.isfinite(value))
        self.assertGreater(value, 0)
    def test_build_gsea_ranking_is_full_signed_and_deterministic(self) -> None:
        frame = pd.DataFrame(
            {
                "gene": ["z", "b", "A", "a", "c", "d"],
                "scores": [2.0, 2.0, -3.0, 1.0, 0.0, np.nan],
                "logfoldchanges": [0.2, 0.5, -1.0, 0.1, 0.0, 0.0],
            }
        )

        ranking, audit = build_gsea_ranking(frame)

        self.assertEqual(ranking["gene"].tolist(), ["b", "z", "c", "A"])
        self.assertEqual(ranking["rank_position"].tolist(), [1, 2, 3, 4])
        self.assertEqual(ranking["rank_score"].tolist(), [2.0, 2.0, 0.0, -3.0])
        self.assertEqual(audit["n_dropped_blank_or_nonfinite"], 1)
        self.assertEqual(audit["n_casefold_collision_keys"], 1)
        self.assertEqual(audit["n_casefold_collision_rows"], 2)
        self.assertEqual(audit["n_genes_in_tied_score_groups"], 2)
        self.assertTrue(audit["has_positive_scores"])
        self.assertTrue(audit["has_negative_scores"])

    def test_gsea_heatmap_selection_keeps_shared_specific_and_na(self) -> None:
        gsea = pd.DataFrame(
            [
                {
                    "cluster": "C0",
                    "term": "Shared tumor",
                    "nes": 2.0,
                    "fdr": 0.001,
                    "direction": "Tumor",
                },
                {
                    "cluster": "C1",
                    "term": "Shared tumor",
                    "nes": 1.5,
                    "fdr": 0.01,
                    "direction": "Tumor",
                },
                {
                    "cluster": "C2",
                    "term": "Shared mixed",
                    "nes": 1.4,
                    "fdr": 0.02,
                    "direction": "Tumor",
                },
                {
                    "cluster": "C3",
                    "term": "Shared mixed",
                    "nes": -1.6,
                    "fdr": 0.03,
                    "direction": "Control",
                },
                {
                    "cluster": "C5-1",
                    "term": "Specific C5-1",
                    "nes": -2.1,
                    "fdr": 0.005,
                    "direction": "Control",
                },
                {
                    "cluster": "C7",
                    "term": "Specific C7",
                    "nes": 1.8,
                    "fdr": 0.02,
                    "direction": "Tumor",
                },
                {
                    "cluster": "C0",
                    "term": "Not significant",
                    "nes": 0.4,
                    "fdr": 0.8,
                    "direction": "Tumor",
                },
            ]
        )

        first = select_gsea_heatmap_pathways(gsea, max_pathways=4)
        second = select_gsea_heatmap_pathways(gsea, max_pathways=4)
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(len(first), 4)
        self.assertEqual(int(first["classification"].eq("shared").sum()), 2)
        self.assertEqual(
            int(
                first["classification"]
                .astype(str)
                .str.startswith("cluster-specific:")
                .sum()
            ),
            2,
        )
        self.assertNotIn("Not significant", set(first["term"]))

        nes, fdr = build_gsea_heatmap_matrices(gsea, first)
        self.assertEqual(list(nes.columns), list(GSEA_CLUSTERS))
        self.assertEqual(list(fdr.columns), list(GSEA_CLUSTERS))
        self.assertAlmostEqual(float(nes.loc["Shared tumor", "C0"]), 2.0)
        self.assertTrue(np.isnan(nes.loc["Shared tumor", "C7"]))
        self.assertTrue(np.isnan(fdr.loc["Shared tumor", "C7"]))


    def test_gsea_module_map_accounts_for_all_significant_terms_once(self) -> None:
        rows: list[dict[str, object]] = []
        representative_nes: dict[str, float] = {}
        for module_index, module in enumerate(GSEA_HEATMAP_MODULES):
            for term_index, term in enumerate(module["source_terms"]):
                nes = 1.0 + module_index / 10 + term_index / 100
                rows.append(
                    {
                        "cluster": "C0",
                        "term": term,
                        "nes": nes,
                        "fdr": 0.01,
                        "direction": "Tumor",
                    }
                )
                if term == module["representative_term"]:
                    representative_nes[str(module["display_label"])] = nes
        gsea = pd.DataFrame(rows)

        selection = select_gsea_heatmap_modules(gsea, clusters=("C0", "C1"))

        self.assertEqual(len(selection), 10)
        self.assertEqual(int(selection["n_source_terms"].sum()), 20)
        source_terms = [
            term
            for value in selection["source_terms"].astype(str)
            for term in value.split(";")
        ]
        self.assertEqual(len(source_terms), len(set(source_terms)))
        self.assertEqual(set(source_terms), set(gsea["term"]))
        self.assertTrue(
            selection["statistic_source"]
            .eq("representative_KEGG_term_NES_and_FDR; no_NES_averaging")
            .all()
        )
        nes, fdr = build_gsea_heatmap_matrices(
            gsea,
            selection,
            clusters=("C0", "C1"),
        )
        self.assertEqual(nes.index.name, "module")
        self.assertEqual(list(nes.index), selection["display_label"].tolist())
        for label, expected in representative_nes.items():
            self.assertAlmostEqual(float(nes.loc[label, "C0"]), expected)
            self.assertAlmostEqual(float(fdr.loc[label, "C0"]), 0.01)
            self.assertTrue(np.isnan(nes.loc[label, "C1"]))

        unmapped = pd.concat(
            [
                gsea,
                pd.DataFrame(
                    [
                        {
                            "cluster": "C0",
                            "term": "Unexpected significant term",
                            "nes": 1.5,
                            "fdr": 0.01,
                            "direction": "Tumor",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
        with self.assertRaisesRegex(RuntimeError, "unmapped"):
            select_gsea_heatmap_modules(unmapped, clusters=("C0", "C1"))

    def test_reusable_gsea_bundle_validates_current_ranking_and_receipt(self) -> None:
        de = pd.DataFrame(
            {
                "scope": ["pooled_cluster"] * 3,
                "cluster": ["C0"] * 3,
                "gene": ["A", "B", "C"],
                "scores": np.array([10.0, 0.0, -5.0], dtype=np.float32),
                "logfoldchanges": np.array([0.5, 0.0, -0.25], dtype=np.float32),
                "n_tumor": [12] * 3,
                "n_control": [14] * 3,
            }
        )
        ranking, _ = build_gsea_ranking(de)
        ranking.insert(0, "cluster", "C0")
        ranking["score_direction"] = np.select(
            [ranking["rank_score"].gt(0), ranking["rank_score"].lt(0)],
            ["Tumor", "Control"],
            default="Zero",
        )
        gene_sets = {"Term A": {"A", "B"}, "Term B": {"C"}}
        current_de = de.copy()
        current_de["scores"] = current_de["scores"].astype(np.float64)
        current_de["logfoldchanges"] = current_de["logfoldchanges"].astype(np.float64)
        current_de.loc[current_de["gene"].eq("A"), "scores"] += 4e-7
        current_de.loc[current_de["gene"].eq("C"), "scores"] -= 2e-7
        current_de.loc[current_de["gene"].eq("A"), "logfoldchanges"] += 2e-8

        gsea = pd.DataFrame(
            [
                {
                    "cluster": "C0",
                    "term": "Term A",
                    "es": 0.5,
                    "nes": 1.2,
                    "nominal_pvalue": 0.02,
                    "fdr": 0.03,
                    "fwer_pvalue": 0.04,
                    "leading_edge_genes": "A",
                    "direction": "Tumor",
                    "significant_fdr_0_05": True,
                    "pathway_size_in_ranking": 2,
                    "n_ranked_genes": 3,
                    "n_tumor": 12,
                    "n_control": 14,
                    "rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
                    "permutation_type": "gene_set",
                    "permutation_num": 50,
                    "weight": 1.0,
                    "seed": 7,
                },
                {
                    "cluster": "C0",
                    "term": "Term B",
                    "es": -0.4,
                    "nes": -1.1,
                    "nominal_pvalue": 0.03,
                    "fdr": 0.04,
                    "fwer_pvalue": 0.05,
                    "leading_edge_genes": "C",
                    "direction": "Control",
                    "significant_fdr_0_05": True,
                    "pathway_size_in_ranking": 1,
                    "n_ranked_genes": 3,
                    "n_tumor": 12,
                    "n_control": 14,
                    "rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
                    "permutation_type": "gene_set",
                    "permutation_num": 50,
                    "weight": 1.0,
                    "seed": 7,
                },
            ]
        )
        identity = gsea_computation_identity(
            permutation_num=50,
            seed=7,
            min_size=1,
            max_size=3,
            clusters=("C0",),
        )
        audit = {
            "parameters": {**identity, "threads": 1, "fdr_threshold": 0.05},
            "gene_set_library": {"frozen_copy_sha256": "synthetic-gmt-sha"},
            "n_results": 2,
            "n_pathways_by_cluster": {"C0": 2},
        }
        with tempfile.TemporaryDirectory() as temporary_directory:
            tables_dir = Path(temporary_directory)
            result_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz"
            ranking_path = (
                tables_dir
                / f"{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz"
            )
            audit_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_audit.json"
            gsea.to_csv(result_path, index=False, compression="gzip")
            ranking.to_csv(ranking_path, index=False, compression="gzip")
            audit_path.write_text(json.dumps(audit), encoding="utf-8")
            receipt_path = (
                tables_dir / f"{RUN_DATE}_KEGG_GSEA_bundle_receipt.json"
            )
            with self.assertRaisesRegex(RuntimeError, "bundle is incomplete"):
                load_validated_gsea_bundle(
                    current_de,
                    gene_sets,
                    {"frozen_copy_sha256": "synthetic-gmt-sha"},
                    tables_dir,
                    permutation_num=50,
                    seed=7,
                    threads=2,
                    min_size=1,
                    max_size=3,
                    clusters=("C0",),
                )
            receipt = build_gsea_bundle_receipt(
                result_path,
                ranking_path,
                gsea=gsea,
                rankings=ranking,
                kegg_audit={"frozen_copy_sha256": "synthetic-gmt-sha"},
                identity=identity,
                provenance="synthetic_trusted_fixture",
            )
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")


            loaded, loaded_rankings, loaded_audit = load_validated_gsea_bundle(
                current_de,
                gene_sets,
                {"frozen_copy_sha256": "synthetic-gmt-sha"},
                tables_dir,
                permutation_num=50,
                seed=7,
                threads=2,
                min_size=1,
                max_size=3,
                clusters=("C0",),
            )
            self.assertEqual(len(loaded), 2)
            self.assertEqual(len(loaded_rankings), 3)
            self.assertEqual(loaded_audit["reuse_validation"]["status"], "passed")
            self.assertTrue(
                receipt_path.is_file()
            )
            score_validation = loaded_audit["reuse_validation"][
                "ranking_validation"
            ][0]["numeric_validation"]["rank_score"]
            self.assertGreater(
                score_validation["max_abs_delta_before_float32_canonicalization"], 0
            )

            tampered = ranking.copy()
            tampered.loc[0, "rank_score"] += 0.1
            tampered.to_csv(ranking_path, index=False, compression="gzip")
            with self.assertRaisesRegex(RuntimeError, "ranking rank_score mismatch"):
                load_validated_gsea_bundle(
                    current_de,
                    gene_sets,
                    {"frozen_copy_sha256": "synthetic-gmt-sha"},
                    tables_dir,
                    permutation_num=50,
                    seed=7,
                    threads=1,
                    min_size=1,
                    max_size=3,
                    clusters=("C0",),
                )

    def test_install_refined_labels_maps_scrambled_c5_by_obs_name(self) -> None:
        adata = make_synthetic_full_adata()
        c5 = make_scrambled_c5(adata)

        mapping = install_refined_labels(adata, c5).set_index("cell_id")

        expected_labels = {
            "cell_c0": "C0",
            "cell_c5_a": "C5-1",
            "cell_c6": "C6",
            "cell_c5_b": "C5-2",
        }
        self.assertEqual(
            adata.obs[REFINED_CLUSTER_KEY].astype(str).to_dict(),
            expected_labels,
        )
        self.assertEqual(
            mapping[REFINED_CLUSTER_KEY].astype(str).to_dict(),
            expected_labels,
        )

        split_key = f"X_umap_c5_only_{RUN_DATE}"
        split_coordinates = np.asarray(adata.obsm[split_key])
        self.assertTrue(np.isnan(split_coordinates[[0, 2]]).all())
        np.testing.assert_allclose(split_coordinates[1], [10.0, 11.0])
        np.testing.assert_allclose(split_coordinates[3], [20.0, 21.0])

    def test_validate_full_gene_input_accepts_valid_and_rejects_invalid_inputs(
        self,
    ) -> None:
        valid = make_synthetic_full_adata()
        summary = validate_full_gene_input(valid)
        self.assertTrue(summary["var_raw_names_equal"])
        self.assertEqual(
            summary["required_markers"],
            {"Il4": True, "Klrd1": True},
        )

        mismatched = make_synthetic_full_adata()
        mismatched.raw = mismatched[:, ["Klrd1", "Il4", "Other"]].copy()
        with self.assertRaisesRegex(RuntimeError, "same full-gene universe"):
            validate_full_gene_input(mismatched)

        for missing_marker in ("Il4", "Klrd1"):
            with self.subTest(missing_marker=missing_marker):
                genes = tuple(
                    gene
                    for gene in ("Il4", "Klrd1", "Other")
                    if gene != missing_marker
                )
                missing = make_synthetic_full_adata(genes)
                with self.assertRaisesRegex(RuntimeError, missing_marker):
                    validate_full_gene_input(missing)

    def test_select_pathways_for_plot_excludes_zero_overlap_terms(self) -> None:
        ora = pd.DataFrame(
            [
                {
                    "threshold": ROBUST_THRESHOLD,
                    "cluster": "C0",
                    "term": "zero-overlap",
                    "overlap_size": 0,
                    "significant_fdr_0_05": True,
                    "fdr": 0.001,
                    "pvalue": 0.001,
                },
                {
                    "threshold": ROBUST_THRESHOLD,
                    "cluster": "C0",
                    "term": "positive-overlap",
                    "overlap_size": 3,
                    "significant_fdr_0_05": False,
                    "fdr": 0.2,
                    "pvalue": 0.02,
                },
                {
                    "threshold": "paper-only",
                    "cluster": "C0",
                    "term": "wrong-threshold",
                    "overlap_size": 5,
                    "significant_fdr_0_05": True,
                    "fdr": 0.001,
                    "pvalue": 0.001,
                },
            ]
        )

        selected = select_pathways_for_plot(ora)

        self.assertEqual(selected, ["positive-overlap"])
        self.assertNotIn("zero-overlap", selected)

    def test_validate_h5ad_roundtrip_preserves_mapping_and_obs_order(self) -> None:
        full = make_synthetic_full_adata()
        c5 = make_scrambled_c5(full)
        install_refined_labels(full, c5)
        split_key = f"X_umap_c5_only_{RUN_DATE}"
        self.assertIn(split_key, full.obsm)

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            full_path = temporary_path / "full.h5ad"
            c5_path = temporary_path / "c5.h5ad"
            full.write_h5ad(full_path)
            c5.write_h5ad(c5_path)

            result = validate_h5ad_roundtrip(
                full_path,
                c5_path,
                full,
                c5,
            )

        self.assertEqual(result["status"], "passed")
        self.assertTrue(result["files"]["full"]["checks"]["obs_names_order"])
        self.assertTrue(result["files"]["c5"]["checks"]["obs_names_order"])
        self.assertTrue(
            result["files"]["full"]["checks"]["required_obsm"]
        )
        self.assertEqual(
            list(full.obs_names),
            ["cell_c0", "cell_c5_a", "cell_c6", "cell_c5_b"],
        )
        self.assertEqual(list(c5.obs_names), ["cell_c5_b", "cell_c5_a"])
        self.assertEqual(
            full.obs[REFINED_CLUSTER_KEY].astype(str).to_dict(),
            {
                "cell_c0": "C0",
                "cell_c5_a": "C5-1",
                "cell_c6": "C6",
                "cell_c5_b": "C5-2",
            },
        )

    def test_run_requested_de_requires_all_figure_f_units_and_propagates_size_flags(
        self,
    ) -> None:
        def minimal_de_frame(unit: object) -> pd.DataFrame:
            tissue = getattr(unit, "tissue")
            passes_min20 = tissue != "thymus"
            n_cells = 25 if passes_min20 else 11
            return pd.DataFrame(
                {
                    "unit_id": [getattr(unit, "unit_id")],
                    "n_tumor": [n_cells],
                    "n_control": [n_cells],
                    "paper_venn_tumor_up": [False],
                    "robust_tumor_up": [False],
                    "robust_control_up": [False],
                    "passes_recommended_min20": [passes_min20],
                    "strict_estimable_min20": [passes_min20],
                    "sample_size_note": [
                        "passes_recommended_min20"
                        if passes_min20
                        else "exploratory_below_recommended_min20"
                    ],
                }
            )

        def fail_one_tissue(_adata: ad.AnnData, unit: object) -> pd.DataFrame:
            if getattr(unit, "unit_id") == "C5-2__spleen":
                raise RuntimeError("synthetic tissue failure")
            return minimal_de_frame(unit)

        def fail_one_c5_1(_adata: ad.AnnData, unit: object) -> pd.DataFrame:
            if getattr(unit, "unit_id") == "C5-1__spleen":
                raise RuntimeError("synthetic C5-1 failure")
            return minimal_de_frame(unit)
        def fail_one_pooled(_adata: ad.AnnData, unit: object) -> pd.DataFrame:
            if getattr(unit, "unit_id") == "pooled_cluster__C3":
                raise RuntimeError("synthetic pooled failure")
            return minimal_de_frame(unit)

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            with patch(
                "run_inkt_c5_paper_followup.run_de_unit",
                side_effect=fail_one_pooled,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "Pooled DE units incomplete",
                ):
                    run_requested_de(
                        ad.AnnData(), root / "failed_pooled"
                    )

            with patch(
                "run_inkt_c5_paper_followup.run_de_unit",
                side_effect=fail_one_tissue,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "C5-2 tissue DE units incomplete",
                ):
                    run_requested_de(ad.AnnData(), root / "failed")

            with patch(
                "run_inkt_c5_paper_followup.run_de_unit",
                side_effect=fail_one_c5_1,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "C5-1 BM/Spleen DE units incomplete",
                ):
                    run_requested_de(
                        ad.AnnData(), root / "failed_c5_1"
                    )

            with patch(
                "run_inkt_c5_paper_followup.run_de_unit",
                side_effect=lambda _adata, unit: minimal_de_frame(unit),
            ):
                combined, status = run_requested_de(
                    ad.AnnData(),
                    root / "completed",
                )

        quality_columns = {
            "passes_recommended_min20",
            "strict_estimable_min20",
            "sample_size_note",
        }
        self.assertTrue(quality_columns.issubset(combined.columns))
        self.assertTrue(quality_columns.issubset(status.columns))
        tissue_status = status[status["scope"].eq("C5-2_tissue")].set_index("tissue")
        pooled_status = status[status["scope"].eq("pooled_cluster")]
        self.assertEqual(
            set(pooled_status["cluster"]),
            set(GSEA_CLUSTERS),
        )
        self.assertTrue(pooled_status["status"].eq("completed").all())
        self.assertFalse(bool(tissue_status.loc["thymus", "passes_recommended_min20"]))
        self.assertFalse(bool(tissue_status.loc["thymus", "strict_estimable_min20"]))
        self.assertEqual(
            tissue_status.loc["thymus", "sample_size_note"],
            "exploratory_below_recommended_min20",
        )
        self.assertTrue(
            tissue_status.loc[
                ["bone_marrow", "spleen"],
                "passes_recommended_min20",
            ].all()
        )

        main_status = status[
            status["unit_id"].isin(FIGURE_F_MAIN_SET_ORDER)
        ].set_index("unit_id")
        self.assertEqual(
            set(main_status.index),
            set(FIGURE_F_MAIN_SET_ORDER),
        )
        self.assertTrue(main_status["passes_recommended_min20"].all())
        self.assertTrue(main_status["strict_estimable_min20"].all())
    def test_audit_raw_ncam1_reports_all_three_observed_states(self) -> None:
        statuses = ("feature_absent", "present_no_counts", "present_with_counts")
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for expected_status in statuses:
                with self.subTest(expected_status=expected_status):
                    status_root = root / expected_status
                    for sample in SAMPLE_DISPLAY_ORDER:
                        write_raw_sample(status_root, sample, expected_status)

                    result = audit_raw_ncam1(status_root)

                    self.assertEqual(
                        result["sample"].tolist(),
                        list(SAMPLE_DISPLAY_ORDER),
                    )
                    self.assertEqual(
                        set(result["audit_status"]),
                        {expected_status},
                    )
                    if expected_status == "feature_absent":
                        self.assertFalse(result["feature_present"].any())
                    else:
                        self.assertTrue(result["feature_present"].all())
                    expected_entries = (
                        1 if expected_status == "present_with_counts" else 0
                    )
                    expected_counts = 3.0 if expected_entries else 0.0
                    self.assertTrue(
                        result["raw_matrix_nonzero_entries"].eq(expected_entries).all()
                    )
                    self.assertTrue(
                        result["raw_matrix_total_counts"].eq(expected_counts).all()
                    )


if __name__ == "__main__":
    unittest.main()
