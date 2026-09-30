#!/usr/bin/env python3
"""Run the dated C5 split and paper-adapted Figure 3D--F follow-up.

The workflow is intentionally non-destructive: the audited 50-PC H5AD is read
as input, the original ``leiden_res_0_5`` labels are retained, and a new
``cluster_c5_split_20260830`` column is written to a new H5AD.  Differential
expression p-values are exploratory cell-level statistics because each
tissue-condition combination contains one sample.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

# Keep the run on CPU and prevent BLAS oversubscription.  In particular, never
# allow a library to select the known-faulty GPU 6 on this host.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"
os.environ["CUDA_VISIBLE_DEVICES"] = ""

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from PIL import Image
from scipy import sparse
from scipy.stats import fisher_exact, hypergeom, spearmanr
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.cluster import KMeans

from inkt_palette import (
    CONDITION_ORDER,
    SAMPLE_DISPLAY_ORDER,
    SAMPLE_PALETTE,
    TISSUE_ORDER,
    TISSUE_PALETTE,
)


ROOT = Path(__file__).resolve().parents[3]
RUN_DATE = "20260830"
DEFAULT_H5AD = (
    ROOT
    / "output/iNKT_reproduction_deck/20260825_UMAP_PC_sweep/umap_pc_sweep/inkt_selected_umap.h5ad"
)
DEFAULT_RAW_ROOT = ROOT / "input/iNKT/data"
DEFAULT_PAPER = ROOT / "docs/blooda_adv-2024-014592-main.pdf"
DEFAULT_OUT_DIR = (
    ROOT / f"output/iNKT_reproduction_deck/{RUN_DATE}_C5_paper_Fig3DEF_followup"
)
ORIGINAL_CLUSTER_KEY = "leiden_res_0_5"
REFINED_CLUSTER_KEY = f"cluster_c5_split_{RUN_DATE}"
C5_SUBCLUSTER_KEY = f"c5_subcluster_{RUN_DATE}"
TISSUE_LABELS = {
    "thymus": "Thymus",
    "bone_marrow": "Bone marrow",
    "spleen": "Spleen",
}
CONDITION_LABELS = {"Ctrl": "Control", "T2": "Tumor"}
REFINED_ORDER = ("C0", "C1", "C2", "C3", "C4", "C5-1", "C5-2", "C6", "C7")
REFINED_COLORS = {
    "C0": "#1F77B4",
    "C1": "#FF7F0E",
    "C2": "#2CA02C",
    "C3": "#D62728",
    "C4": "#9467BD",
    "C5-1": "#8C564B",
    "C5-2": "#E377C2",
    "C6": "#17BECF",
    "C7": "#7F7F7F",
}
TARGET_CLUSTERS = ("C5-1", "C5-2", "C0", "C6")
GSEA_CLUSTERS = REFINED_ORDER
GSEA_MIN_SIZE = 15
GSEA_MAX_SIZE = 500
GSEA_DEFAULT_PERMUTATIONS = 5000
GSEA_DEFAULT_SEED = 20260830
GSEA_WEIGHT = 1.0
GSEA_FDR_THRESHOLD = 0.05
GSEA_RANK_SCORE_ATOL = 1e-6
GSEA_RANK_LOGFC_ATOL = 1e-6
GSEA_HEATMAP_MAX_PATHWAYS = 24
GSEA_HEATMAP_MODULES = (
    {
        "display_label": "Ribosomal gene program",
        "representative_term": "Ribosome",
        "source_terms": ("Ribosome",),
        "interpretation": "Rpl/Rps-led ribosomal gene program",
    },
    {
        "display_label": "Mitochondrial ETC / OXPHOS",
        "representative_term": "Oxidative phosphorylation",
        "source_terms": (
            "Oxidative phosphorylation",
            "Parkinson disease",
            "Cardiac muscle contraction",
            "Alzheimer disease",
            "Thermogenesis",
            "Huntington disease",
            "Non-alcoholic fatty liver disease (NAFLD)",
            "Retrograde endocannabinoid signaling",
        ),
        "interpretation": "Cox/Nduf/Atp5/Uqcr-led mitochondrial ETC program",
    },
    {
        "display_label": "Proteasome",
        "representative_term": "Proteasome",
        "source_terms": ("Proteasome",),
        "interpretation": "Proteasome gene program",
    },
    {
        "display_label": "RNA transport",
        "representative_term": "RNA transport",
        "source_terms": ("RNA transport",),
        "interpretation": "Broad KEGG RNA transport program",
    },
    {
        "display_label": "RNA polymerase",
        "representative_term": "RNA polymerase",
        "source_terms": ("RNA polymerase",),
        "interpretation": "Polr1/2/3-subunit program",
    },
    {
        "display_label": "Antigen processing/presentation (MHC-I-led)",
        "representative_term": "Antigen processing and presentation",
        "source_terms": (
            "Antigen processing and presentation",
            "Type I diabetes mellitus",
            "Graft-versus-host disease",
        ),
        "interpretation": "H2/B2m/Tap-led antigen-presentation program",
    },
    {
        "display_label": "PI3K-Akt signaling",
        "representative_term": "PI3K-Akt signaling pathway",
        "source_terms": ("PI3K-Akt signaling pathway",),
        "interpretation": "KEGG PI3K-Akt transcriptional program; not kinase activity",
    },
    {
        "display_label": "DNA repair / homologous recombination",
        "representative_term": "Homologous recombination",
        "source_terms": ("Homologous recombination",),
        "interpretation": "Brca/Rad/Rpa-led homologous-recombination program",
    },
    {
        "display_label": "Lysine degradation (KMT/SETD-led edge)",
        "representative_term": "Lysine degradation",
        "source_terms": ("Lysine degradation",),
        "interpretation": "Original KEGG term; leading edge is Kmt/Setd-enriched",
    },
    {
        "display_label": "AP-1 immediate-early proxy‡",
        "representative_term": "Amphetamine addiction",
        "source_terms": ("Amphetamine addiction", "Cocaine addiction"),
        "interpretation": "Post-hoc Jun/Fos leading-edge proxy; not a formal AP-1 gene set",
        "representative_rationale": "broader Jun/Fos/Fosb leading edge among the two significant proxy pathways",
    },
)
FIGURE_F_MAIN_CLUSTERS = ("C5-1", "C5-2")
FIGURE_F_MAIN_TISSUES = ("bone_marrow", "spleen")
FIGURE_F_MAIN_SET_ORDER = tuple(
    f"{cluster}__{tissue}"
    for cluster in FIGURE_F_MAIN_CLUSTERS
    for tissue in FIGURE_F_MAIN_TISSUES
)
FIGURE_F_MAIN_SET_LABELS = {
    f"{cluster}__{tissue}": f"{cluster} · {TISSUE_LABELS[tissue]}"
    for cluster in FIGURE_F_MAIN_CLUSTERS
    for tissue in FIGURE_F_MAIN_TISSUES
}
PAPER_VENN_THRESHOLD = "paper_nominal_p_lt_0.05_log2fc_ge_1"
ROBUST_THRESHOLD = "robust_bh_fdr_le_0.05_log2fc_ge_0.25"
STATISTICAL_WARNING = (
    "exploratory_cell_level_comparison; one_sample_per_tissue_condition; "
    "p_values_do_not_establish_biological_replicate_significance"
)
EXPLORATORY_MIN_CELLS_PER_CONDITION = 2
RECOMMENDED_MIN_CELLS_PER_CONDITION = 20


@dataclass(frozen=True)
class DEUnit:
    unit_id: str
    scope: str
    cluster: str
    tissue: str | None = None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--paper", type=Path, default=DEFAULT_PAPER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-hvg", type=int, default=3000)
    parser.add_argument("--min-dist", type=float, default=0.5)
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=float,
        default=[0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00],
    )
    parser.add_argument("--stability-seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    parser.add_argument("--kegg-library", default="KEGG_2019_Mouse")
    parser.add_argument("--kegg-gmt", type=Path)
    parser.add_argument(
        "--gsea-permutations",
        type=int,
        default=GSEA_DEFAULT_PERMUTATIONS,
    )
    parser.add_argument("--gsea-seed", type=int, default=GSEA_DEFAULT_SEED)
    parser.add_argument("--gsea-threads", type=int, default=1)
    parser.add_argument(
        "--gsea-max-pathways",
        type=int,
        default=GSEA_HEATMAP_MAX_PATHWAYS,
        help="Legacy raw-pathway selector cap; ignored by the curated main heatmap.",
    )
    parser.add_argument(
        "--reuse-gsea-tables",
        action="store_true",
        help=(
            "Reuse the dated GSEA result/ranking tables only after fail-closed "
            "validation against the current DE rankings, frozen GMT, and audit."
        ),
    )
    parser.add_argument("--skip-h5ad", action="store_true")
    args = parser.parse_args(argv)
    if args.n_neighbors < 2:
        parser.error("--n-neighbors must be >= 2")
    if args.n_pcs < 2 or args.n_hvg < 2:
        parser.error("--n-pcs and --n-hvg must be >= 2")
    if not args.resolutions or any(value <= 0 for value in args.resolutions):
        parser.error("--resolutions must contain positive values")
    if not args.stability_seeds:
        parser.error("--stability-seeds cannot be empty")
    if args.gsea_permutations < 10:
        parser.error("--gsea-permutations must be >= 10")
    if args.gsea_threads < 1:
        parser.error("--gsea-threads must be >= 1")
    if args.gsea_max_pathways < 2:
        parser.error("--gsea-max-pathways must be >= 2")
    return args


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    adjusted = np.full(array.shape, np.nan, dtype=float)
    finite = np.isfinite(array)
    if not finite.any():
        return adjusted
    raw = array[finite]
    order = np.argsort(raw)
    ranked = raw[order]
    ranked_adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    ranked_adjusted = np.minimum.accumulate(ranked_adjusted[::-1])[::-1]
    restored = np.empty_like(ranked_adjusted)
    restored[order] = np.clip(ranked_adjusted, 0.0, 1.0)
    adjusted[finite] = restored
    return adjusted


def save_figure(figure: Any, png_path: Path, *, dpi: int = 220) -> tuple[Path, Path]:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = png_path.with_suffix(".pdf")
    figure.savefig(png_path, dpi=dpi, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def axis_extent(coordinates: np.ndarray) -> tuple[float, float, float, float]:
    minimum = np.asarray(coordinates).min(axis=0)
    maximum = np.asarray(coordinates).max(axis=0)
    padding = np.maximum((maximum - minimum) * 0.035, 0.1)
    return (
        float(minimum[0] - padding[0]),
        float(maximum[0] + padding[0]),
        float(minimum[1] - padding[1]),
        float(maximum[1] + padding[1]),
    )


def style_umap_axis(
    axis: Any,
    extent: tuple[float, float, float, float],
    *,
    xlabel: str = "UMAP1",
    ylabel: str = "UMAP2",
) -> None:
    axis.set_xlim(extent[0], extent[1])
    axis.set_ylim(extent[2], extent[3])
    axis.set_xlabel(xlabel, fontsize=8)
    axis.set_ylabel(ylabel, fontsize=8)
    axis.tick_params(labelsize=6, length=2, colors="#596674")
    axis.grid(color="#D8DEE5", linewidth=0.35, alpha=0.55, zorder=0)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color("#94A3B8")


def mean_pairwise_ari(label_sets: Sequence[np.ndarray]) -> tuple[float, float]:
    values = [adjusted_rand_score(label_sets[left], label_sets[right]) for left, right in combinations(range(len(label_sets)), 2)]
    if not values:
        return 1.0, 1.0
    return float(np.mean(values)), float(np.min(values))


def select_two_cluster_resolution(metrics: pd.DataFrame) -> float:
    required = {"resolution", "all_seeds_two_clusters", "mean_pairwise_ari", "min_pairwise_ari"}
    missing = required - set(metrics.columns)
    if missing:
        raise KeyError(f"Missing resolution metrics: {sorted(missing)}")
    candidates = metrics[metrics["all_seeds_two_clusters"].astype(bool)].copy()
    if candidates.empty:
        raise RuntimeError("No candidate resolution produced exactly two clusters for every seed")
    candidates["distance_to_0_1"] = (candidates["resolution"] - 0.1).abs()
    candidates = candidates.sort_values(
        ["mean_pairwise_ari", "min_pairwise_ari", "distance_to_0_1", "resolution"],
        ascending=[False, False, True, True],
        kind="stable",
    )
    return float(candidates.iloc[0]["resolution"])


def validate_full_gene_input(adata: ad.AnnData) -> dict[str, Any]:
    if adata.raw is None:
        raise RuntimeError("Full-gene adata.raw is required")
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] is required")
    if not adata.obs_names.is_unique:
        raise RuntimeError("Cell index must be unique before C5 labels are mapped back")
    var_names = pd.Index(adata.var_names.astype(str))
    raw_var_names = pd.Index(adata.raw.var_names.astype(str))
    if adata.n_vars != adata.raw.n_vars or not var_names.equals(raw_var_names):
        raise RuntimeError(
            "C5 HVG selection requires X and raw to contain the same full-gene universe "
            f"in the same order; X={adata.n_vars}, raw={adata.raw.n_vars}"
        )
    marker_names = ("Il4", "Klrd1")
    missing_markers = [gene for gene in marker_names if gene not in var_names]
    if missing_markers:
        raise RuntimeError(f"Required marker genes are absent: {missing_markers}")
    if adata.layers["counts"].shape != adata.shape:
        raise RuntimeError("layers['counts'] shape does not match adata.X")
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "raw_n_vars": int(adata.raw.n_vars),
        "var_raw_names_equal": True,
        "counts_shape": list(map(int, adata.layers["counts"].shape)),
        "obs_names_unique": True,
        "required_markers": {gene: True for gene in marker_names},
    }


def run_c5_subclustering(
    adata: ad.AnnData,
    *,
    n_hvg: int,
    n_pcs: int,
    n_neighbors: int,
    min_dist: float,
    resolutions: Sequence[float],
    seeds: Sequence[int],
    random_state: int,
) -> tuple[ad.AnnData, pd.DataFrame, float]:
    original_labels = adata.obs[ORIGINAL_CLUSTER_KEY].astype(str)
    if len({int(seed) for seed in seeds}) < 2:
        raise ValueError("C5 stability/medoid selection requires at least two unique seeds")
    mask = original_labels.eq("5").to_numpy()
    if not mask.any():
        raise RuntimeError("No C5 cells were found")
    c5 = adata[mask].copy()
    c5.obsm["X_umap_original_full"] = np.asarray(adata.obsm["X_umap"])[mask].copy()
    c5.uns.pop("iroot", None)
    for column in [name for name in c5.var.columns if name.startswith("highly_variable")]:
        del c5.var[column]
    for column in ("means", "dispersions", "dispersions_norm"):
        if column in c5.var:
            del c5.var[column]
    sc.pp.highly_variable_genes(
        c5,
        n_top_genes=min(int(n_hvg), c5.n_vars),
        batch_key="sample",
        flavor="seurat",
    )
    available_hvg = int(c5.var["highly_variable"].sum())
    effective_pcs = min(int(n_pcs), c5.n_obs - 1, available_hvg - 1)
    if effective_pcs < 2:
        raise RuntimeError(f"Insufficient C5 dimensions for PCA: {effective_pcs}")
    sc.tl.pca(
        c5,
        n_comps=effective_pcs,
        svd_solver="arpack",
        mask_var="highly_variable",
        random_state=random_state,
    )
    sc.pp.neighbors(
        c5,
        n_neighbors=min(int(n_neighbors), c5.n_obs - 1),
        n_pcs=effective_pcs,
        use_rep="X_pca",
        metric="euclidean",
        random_state=random_state,
    )
    sc.tl.umap(c5, min_dist=min_dist, spread=1.0, random_state=random_state)

    aggregate_rows: list[dict[str, Any]] = []
    labels_by_resolution: dict[float, list[np.ndarray]] = {}
    original_c5_umap = np.asarray(c5.obsm["X_umap_original_full"])
    geometry_labels = KMeans(n_clusters=2, n_init=50, random_state=random_state).fit_predict(
        original_c5_umap
    )
    for resolution in sorted({float(value) for value in resolutions}):
        seed_labels: list[np.ndarray] = []
        cluster_counts_by_seed: list[int] = []
        for seed in seeds:
            key = f"_c5_leiden_res_{resolution:g}_seed_{seed}"
            sc.tl.leiden(
                c5,
                resolution=resolution,
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=int(seed),
            )
            labels = c5.obs[key].astype(str).to_numpy()
            seed_labels.append(labels)
            cluster_counts_by_seed.append(int(len(np.unique(labels))))
        labels_by_resolution[resolution] = seed_labels
        mean_ari, minimum_ari = mean_pairwise_ari(seed_labels)
        reference = seed_labels[0]
        n_reference = len(np.unique(reference))
        pca_silhouette = (
            float(silhouette_score(np.asarray(c5.obsm["X_pca"]), reference))
            if n_reference > 1
            else float("nan")
        )
        umap_silhouette = (
            float(silhouette_score(np.asarray(c5.obsm["X_umap"]), reference))
            if n_reference > 1
            else float("nan")
        )
        original_umap_silhouette = (
            float(silhouette_score(original_c5_umap, reference))
            if n_reference > 1
            else float("nan")
        )
        geometry_ari = (
            float(adjusted_rand_score(geometry_labels, reference))
            if n_reference == 2
            else float("nan")
        )
        aggregate_rows.append(
            {
                "resolution": resolution,
                "cluster_counts_by_seed": ";".join(map(str, cluster_counts_by_seed)),
                "n_clusters_seed0": n_reference,
                "all_seeds_two_clusters": all(value == 2 for value in cluster_counts_by_seed),
                "mean_pairwise_ari": mean_ari,
                "min_pairwise_ari": minimum_ari,
                "pca_silhouette_seed0": pca_silhouette,
                "umap_silhouette_seed0": umap_silhouette,
                "original_umap_silhouette_seed0": original_umap_silhouette,
                "original_umap_kmeans_ari_seed0": geometry_ari,
                "cluster_sizes_seed0": ";".join(
                    f"{label}:{count}"
                    for label, count in pd.Series(reference).value_counts().sort_index().items()
                ),
            }
        )
    metrics = pd.DataFrame(aggregate_rows).sort_values("resolution").reset_index(drop=True)
    selected_resolution = select_two_cluster_resolution(metrics)
    selected_seed_labels = labels_by_resolution[selected_resolution]
    mean_ari_by_seed = [
        float(
            np.mean(
                [
                    adjusted_rand_score(labels, other)
                    for other_index, other in enumerate(selected_seed_labels)
                    if other_index != label_index
                ]
            )
        )
        for label_index, labels in enumerate(selected_seed_labels)
    ]
    medoid_index = int(np.argmax(mean_ari_by_seed))
    selected_labels = selected_seed_labels[medoid_index]
    selected_partition_seed = int(seeds[medoid_index])

    # Make C5-2 the right-hand island in the audited full-data UMAP, as requested.
    original_c5_umap = np.asarray(c5.obsm["X_umap_original_full"])
    label_median_x = {
        label: float(np.median(original_c5_umap[selected_labels == label, 0]))
        for label in np.unique(selected_labels)
    }
    ordered_labels = sorted(label_median_x, key=label_median_x.get)
    if len(ordered_labels) != 2:
        raise RuntimeError(f"Selected C5 solution has {len(ordered_labels)} clusters, expected 2")
    mapping = {ordered_labels[0]: "C5-1", ordered_labels[1]: "C5-2"}
    c5.obs[C5_SUBCLUSTER_KEY] = pd.Categorical(
        [mapping[label] for label in selected_labels],
        categories=["C5-1", "C5-2"],
        ordered=True,
    )
    metrics["selected"] = np.isclose(metrics["resolution"], selected_resolution)
    metrics["n_c5_cells"] = c5.n_obs
    metrics["n_hvg"] = available_hvg
    metrics["n_pcs"] = effective_pcs
    selected_mask = metrics["selected"]
    left_x = original_c5_umap[selected_labels == ordered_labels[0], 0]
    right_x = original_c5_umap[selected_labels == ordered_labels[1], 0]
    metrics["selected_partition_seed"] = np.nan
    metrics["selected_seed_mean_ari"] = np.nan
    metrics["selected_original_umap_kmeans_ari"] = np.nan
    metrics["selected_original_umap_silhouette"] = np.nan
    metrics["selected_original_umap_median_x_gap"] = np.nan
    metrics["selected_original_umap_q05_q95_gap"] = np.nan
    metrics.loc[selected_mask, "selected_partition_seed"] = selected_partition_seed
    metrics.loc[selected_mask, "selected_seed_mean_ari"] = mean_ari_by_seed[medoid_index]
    metrics.loc[selected_mask, "selected_original_umap_kmeans_ari"] = adjusted_rand_score(
        geometry_labels, selected_labels
    )
    metrics.loc[selected_mask, "selected_original_umap_silhouette"] = silhouette_score(
        original_c5_umap, selected_labels
    )
    metrics.loc[selected_mask, "selected_original_umap_median_x_gap"] = (
        float(np.median(right_x)) - float(np.median(left_x))
    )
    metrics.loc[selected_mask, "selected_original_umap_q05_q95_gap"] = (
        float(np.quantile(right_x, 0.05)) - float(np.quantile(left_x, 0.95))
    )

    metrics["n_neighbors"] = min(int(n_neighbors), c5.n_obs - 1)
    return c5, metrics, selected_resolution


def install_refined_labels(adata: ad.AnnData, c5: ad.AnnData) -> pd.DataFrame:
    original = adata.obs[ORIGINAL_CLUSTER_KEY].astype(str)
    expected_c5_names = adata.obs_names[original.eq("5").to_numpy()]
    if not c5.obs_names.is_unique:
        raise RuntimeError("C5 cell index must be unique")
    missing_c5 = expected_c5_names.difference(c5.obs_names)
    extra_c5 = c5.obs_names.difference(expected_c5_names)
    if len(missing_c5) or len(extra_c5):
        raise RuntimeError(
            f"C5 index mismatch: missing={len(missing_c5)}, extra={len(extra_c5)}"
        )
    refined = pd.Series("C" + original, index=adata.obs_names, dtype="string")
    c5_labels = c5.obs[C5_SUBCLUSTER_KEY].astype(str)
    refined.loc[c5.obs_names] = c5_labels.to_numpy()
    unknown = sorted(set(refined.dropna()) - set(REFINED_ORDER))
    if unknown:
        raise RuntimeError(f"Unexpected refined cluster labels: {unknown}")
    adata.obs[REFINED_CLUSTER_KEY] = pd.Categorical(refined, categories=REFINED_ORDER, ordered=True)
    c5_column = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
    c5_column.loc[c5.obs_names] = c5_labels.to_numpy()
    adata.obs[C5_SUBCLUSTER_KEY] = pd.Categorical(
        c5_column,
        categories=["C5-1", "C5-2"],
        ordered=True,
    )
    coordinates = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    positions = adata.obs_names.get_indexer(c5.obs_names)
    if np.any(positions < 0):
        raise RuntimeError("C5 cell names could not be mapped back to the full object")
    coordinates[positions] = np.asarray(c5.obsm["X_umap"], dtype=np.float32)
    adata.obsm[f"X_umap_c5_only_{RUN_DATE}"] = coordinates
    adata.uns[f"{REFINED_CLUSTER_KEY}_colors"] = [REFINED_COLORS[label] for label in REFINED_ORDER]

    mapping = adata.obs[["sample", "condition", "tissue", ORIGINAL_CLUSTER_KEY, REFINED_CLUSTER_KEY, C5_SUBCLUSTER_KEY]].copy()
    mapping.insert(0, "cell_id", mapping.index.astype(str))
    mapping["original_umap1"] = np.asarray(adata.obsm["X_umap"])[:, 0]
    mapping["original_umap2"] = np.asarray(adata.obsm["X_umap"])[:, 1]
    mapping["c5_only_umap1"] = coordinates[:, 0]
    mapping["c5_only_umap2"] = coordinates[:, 1]
    mapping["is_original_c5"] = original.eq("5").to_numpy()
    return mapping.reset_index(drop=True)


def build_cluster_proportion_tables(
    obs: pd.DataFrame,
    labels: Sequence[str],
    cluster_order: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.DataFrame(
        {
            "tissue": obs["tissue"].astype(str).to_numpy(),
            "condition": obs["condition"].astype(str).to_numpy(),
            "cluster": np.asarray(labels, dtype=str),
        },
        index=obs.index,
    )
    grouped = frame.groupby(["tissue", "condition", "cluster"], observed=True).size()
    complete = pd.MultiIndex.from_product(
        [TISSUE_ORDER, CONDITION_ORDER, list(cluster_order)],
        names=["tissue", "condition", "cluster"],
    )
    table = grouped.reindex(complete, fill_value=0).rename("count").reset_index()
    table["count"] = table["count"].astype(int)
    table["total_tissue_condition"] = table.groupby(["tissue", "condition"], observed=True)["count"].transform("sum")
    table["pct_cluster_within_tissue_condition"] = table["count"] / table["total_tissue_condition"] * 100.0
    table["total_tissue_cluster"] = table.groupby(["tissue", "cluster"], observed=True)["count"].transform("sum")
    table["pct_condition_within_tissue_cluster"] = np.where(
        table["total_tissue_cluster"] > 0,
        table["count"] / table["total_tissue_cluster"] * 100.0,
        np.nan,
    )
    table["total_tissue_all_conditions"] = table.groupby("tissue", observed=True)["count"].transform("sum")
    tissue_condition_total = (
        table[["tissue", "condition", "total_tissue_condition"]]
        .drop_duplicates()
        .pivot(index="tissue", columns="condition", values="total_tissue_condition")
    )
    tumor_baseline = (
        tissue_condition_total["T2"] / (tissue_condition_total["Ctrl"] + tissue_condition_total["T2"]) * 100.0
    )
    table["tumor_share_zero_delta_baseline_pct"] = table["tissue"].map(tumor_baseline)

    frequency = table.pivot(index=["tissue", "cluster"], columns="condition", values="pct_cluster_within_tissue_condition").reset_index()
    counts = table.pivot(index=["tissue", "cluster"], columns="condition", values="count").reset_index()
    shares = table.pivot(index=["tissue", "cluster"], columns="condition", values="pct_condition_within_tissue_cluster").reset_index()
    delta = frequency.merge(counts, on=["tissue", "cluster"], suffixes=("_pct", "_count"))
    delta = delta.merge(shares, on=["tissue", "cluster"], suffixes=("", "_share"))
    delta = delta.rename(
        columns={
            "Ctrl_pct": "control_cluster_frequency_pct",
            "T2_pct": "tumor_cluster_frequency_pct",
            "Ctrl_count": "control_count",
            "T2_count": "tumor_count",
            "Ctrl": "control_share_within_tissue_cluster_pct",
            "T2": "tumor_share_within_tissue_cluster_pct",
        }
    )
    delta["tumor_minus_control_percentage_points"] = delta["tumor_cluster_frequency_pct"] - delta["control_cluster_frequency_pct"]
    delta["tumor_share_zero_delta_baseline_pct"] = delta["tissue"].map(tumor_baseline)
    delta["tumor_share_minus_zero_delta_baseline_pp"] = (
        delta["tumor_share_within_tissue_cluster_pct"] - delta["tumor_share_zero_delta_baseline_pct"]
    )
    nonempty = delta["control_count"].add(delta["tumor_count"]).gt(0)
    raw_a = delta.loc[nonempty, "tumor_minus_control_percentage_points"].to_numpy(dtype=float)
    raw_b = delta.loc[nonempty, "tumor_share_minus_zero_delta_baseline_pp"].to_numpy(dtype=float)
    sign_a = np.where(np.isclose(raw_a, 0.0, atol=1e-12), 0.0, np.sign(raw_a))
    sign_b = np.where(np.isclose(raw_b, 0.0, atol=1e-12), 0.0, np.sign(raw_b))
    delta["signs_reconcile"] = True
    delta.loc[nonempty, "signs_reconcile"] = sign_a == sign_b

    sums = table.groupby(["tissue", "condition"], observed=True)["pct_cluster_within_tissue_condition"].sum()
    if not np.allclose(sums.to_numpy(), 100.0, atol=1e-9):
        raise RuntimeError(f"Cluster frequencies do not sum to 100: {sums.to_dict()}")
    if int(table["count"].sum()) != len(frame):
        raise RuntimeError("Proportion table count does not match the source cells")
    if not bool(delta["signs_reconcile"].all()):
        raise RuntimeError("Conditional-share baseline and cluster-frequency delta signs disagree")
    return table, delta
def add_cluster_centroid_labels(
    axis: Any,
    coordinates: np.ndarray,
    labels: Sequence[str],
    order: Sequence[str],
    *,
    fontsize: float = 6.0,
) -> None:
    values = np.asarray(labels, dtype=str)
    for label in order:
        points = coordinates[values == label]
        if not len(points):
            continue
        center = np.median(points, axis=0)
        axis.text(
            center[0],
            center[1],
            label,
            ha="center",
            va="center",
            fontsize=fontsize,
            fontweight="bold",
            color="white",
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": REFINED_COLORS.get(label, "#475569"),
                "edgecolor": "white",
                "linewidth": 0.5,
                "alpha": 0.90,
            },
            zorder=6,
        )


def render_c5_validation(
    adata: ad.AnnData,
    c5: ad.AnnData,
    metrics: pd.DataFrame,
    selected_resolution: float,
    path: Path,
) -> tuple[Path, Path]:
    full_coordinates = np.asarray(adata.obsm["X_umap"])
    full_extent = axis_extent(full_coordinates)
    refined = adata.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy()
    c5_mask = np.isin(refined, ["C5-1", "C5-2"])
    c5_labels = c5.obs[C5_SUBCLUSTER_KEY].astype(str).to_numpy()
    c5_coordinates = np.asarray(c5.obsm["X_umap"])
    c5_extent = axis_extent(c5_coordinates)

    figure, axes = plt.subplots(2, 2, figsize=(13.5, 10.2), layout="constrained")
    axis = axes[0, 0]
    axis.scatter(
        full_coordinates[~c5_mask, 0],
        full_coordinates[~c5_mask, 1],
        s=2.0,
        c="#D7DEE7",
        alpha=0.28,
        edgecolors="none",
        rasterized=True,
    )
    for label in ("C5-1", "C5-2"):
        mask = refined == label
        axis.scatter(
            full_coordinates[mask, 0],
            full_coordinates[mask, 1],
            s=5.0,
            c=REFINED_COLORS[label],
            alpha=0.78,
            edgecolors="none",
            rasterized=True,
            label=f"{label} (n={int(mask.sum()):,})",
        )
    add_cluster_centroid_labels(axis, full_coordinates, refined, ["C5-1", "C5-2"], fontsize=7)
    style_umap_axis(axis, full_extent)
    axis.set_title("A. C5 split mapped back to the audited full UMAP", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=8, loc="upper right")

    axis = axes[0, 1]
    for label in ("C5-1", "C5-2"):
        mask = c5_labels == label
        axis.scatter(
            c5_coordinates[mask, 0],
            c5_coordinates[mask, 1],
            s=8.0,
            c=REFINED_COLORS[label],
            alpha=0.82,
            edgecolors="none",
            rasterized=True,
            label=f"{label} (n={int(mask.sum()):,})",
        )
    add_cluster_centroid_labels(axis, c5_coordinates, c5_labels, ["C5-1", "C5-2"], fontsize=8)
    style_umap_axis(axis, c5_extent, xlabel="C5 UMAP1", ylabel="C5 UMAP2")
    axis.set_title("B. C5-only UMAP after within-C5 HVG/PCA/20-NN", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[1, 0]
    composition = pd.crosstab(
        c5.obs[C5_SUBCLUSTER_KEY].astype(str),
        c5.obs["sample"].astype(str),
    ).reindex(index=["C5-1", "C5-2"], columns=SAMPLE_DISPLAY_ORDER, fill_value=0)
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100.0
    bottom = np.zeros(len(composition_pct))
    x_positions = np.arange(len(composition_pct))
    for sample in SAMPLE_DISPLAY_ORDER:
        values = composition_pct[sample].to_numpy(dtype=float)
        axis.bar(
            x_positions,
            values,
            bottom=bottom,
            color=SAMPLE_PALETTE[sample],
            width=0.68,
            edgecolor="white",
            linewidth=0.5,
            label=sample,
        )
        bottom += values
    axis.set_xticks(x_positions, composition_pct.index)
    axis.set_ylim(0, 100)
    axis.set_ylabel("Cells within refined C5 cluster (%)")
    axis.set_title("C. Sample composition (descriptive only)", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#E2E8F0", linewidth=0.5)

    axis = axes[1, 1]
    resolution = metrics["resolution"].to_numpy(dtype=float)
    axis.plot(
        resolution,
        metrics["n_clusters_seed0"].to_numpy(dtype=float),
        marker="o",
        color="#334155",
        label="clusters (seed 0)",
    )
    axis.set_xlabel("Leiden resolution")
    axis.set_ylabel("Number of clusters", color="#334155")
    axis.tick_params(axis="y", labelcolor="#334155")
    axis.axvline(selected_resolution, color="#C2410C", linestyle="--", linewidth=1.2)
    axis_ari = axis.twinx()
    axis_ari.plot(
        resolution,
        metrics["mean_pairwise_ari"].to_numpy(dtype=float),
        marker="s",
        color="#0F766E",
        label="mean pairwise ARI",
    )
    axis_ari.set_ylim(0, 1.05)
    axis_ari.set_ylabel("Cross-seed ARI", color="#0F766E")
    axis_ari.tick_params(axis="y", labelcolor="#0F766E")
    selected_row = metrics.loc[np.isclose(metrics["resolution"], selected_resolution)].iloc[0]
    axis.text(
        selected_resolution,
        1.05,
        f" selected {selected_resolution:g}\n"
        f"seed ARI={selected_row['mean_pairwise_ari']:.3f}\n"
        f"geometry ARI={selected_row['selected_original_umap_kmeans_ari']:.3f}",
        color="#9A3412",
        fontsize=8,
        va="bottom",
    )
    axis.set_title("D. Resolution and random-seed stability", loc="left", fontweight="bold")
    axis.spines["top"].set_visible(False)
    axis.grid(axis="x", color="#E2E8F0", linewidth=0.5)
    handles = [
        Line2D([0], [0], color="#334155", marker="o", label="clusters (seed 0)"),
        Line2D([0], [0], color="#0F766E", marker="s", label="mean pairwise ARI"),
    ]
    axis.legend(handles=handles, frameon=False, fontsize=8, loc="upper right")

    figure.suptitle(
        f"{RUN_DATE} C5 split validation — selected resolution {selected_resolution:g}\n"
        "C5-2 is defined as the right-hand island in the original 50-PC UMAP",
        fontsize=15,
        fontweight="bold",
    )
    return save_figure(figure, path)


def render_refined_full_umap(adata: ad.AnnData, path: Path) -> tuple[Path, Path]:
    coordinates = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy()
    figure, axis = plt.subplots(figsize=(10.8, 8.4), layout="constrained")
    for label in REFINED_ORDER:
        mask = labels == label
        if not mask.any():
            continue
        axis.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=2.4 if not label.startswith("C5-") else 3.6,
            c=REFINED_COLORS[label],
            alpha=0.68 if not label.startswith("C5-") else 0.84,
            edgecolors="none",
            rasterized=True,
        )
    add_cluster_centroid_labels(axis, coordinates, labels, REFINED_ORDER, fontsize=6.5)
    style_umap_axis(axis, axis_extent(coordinates))
    axis.set_title(
        f"{RUN_DATE} refined full UMAP: original C5 replaced by C5-1 / C5-2",
        fontsize=14,
        fontweight="bold",
    )
    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=REFINED_COLORS[label],
            markeredgecolor="none",
            label=label,
            markersize=6,
        )
        for label in REFINED_ORDER
    ]
    axis.legend(handles=legend, frameon=False, ncol=3, loc="upper right", fontsize=8)
    axis.text(
        0.01,
        0.01,
        "Coordinates unchanged from selected 50-PC UMAP; only C5 labels are refined.",
        transform=axis.transAxes,
        fontsize=7.5,
        color="#475569",
    )
    return save_figure(figure, path)


def render_denominator_reconciliation(
    delta: pd.DataFrame,
    cluster_order: Sequence[str],
    path: Path,
    *,
    title_suffix: str,
) -> tuple[Path, Path]:
    figure, axes = plt.subplots(2, 3, figsize=(17.2, 8.7), sharex="col", layout="constrained")
    for column, tissue in enumerate(TISSUE_ORDER):
        subset = (
            delta[delta["tissue"].eq(tissue)]
            .set_index("cluster")
            .reindex(cluster_order)
            .reset_index()
        )
        x_positions = np.arange(len(cluster_order))
        baseline = float(subset["tumor_share_zero_delta_baseline_pct"].dropna().iloc[0])
        top = axes[0, column]
        tumor_share = subset["tumor_share_within_tissue_cluster_pct"].to_numpy(dtype=float)
        control_share = 100.0 - tumor_share
        top.bar(
            x_positions,
            control_share,
            color="#457B9D",
            width=0.72,
            label="Control share",
        )
        top.bar(
            x_positions,
            tumor_share,
            bottom=control_share,
            color="#B23A48",
            width=0.72,
            label="Tumor share",
        )
        top.axhline(
            100.0 - baseline,
            color="#111827",
            linestyle="--",
            linewidth=1.2,
            label=f"zero-delta boundary (Tumor {baseline:.2f}%)",
        )
        top.axhline(50.0, color="#94A3B8", linestyle=":", linewidth=1.0, label="50%")
        top.set_ylim(0, 100)
        top.set_ylabel("Condition share inside tissue×cluster (%)" if column == 0 else "")
        top.set_title(TISSUE_LABELS[tissue], fontweight="bold", color=TISSUE_PALETTE[tissue])
        top.legend(frameon=False, fontsize=7, loc="upper right")
        top.grid(axis="y", color="#E2E8F0", linewidth=0.45)
        top.spines[["top", "right"]].set_visible(False)

        bottom = axes[1, column]
        changes = subset["tumor_minus_control_percentage_points"].to_numpy(dtype=float)
        bottom.bar(
            x_positions,
            changes,
            color=np.where(changes >= 0, "#B23A48", "#457B9D"),
            width=0.72,
        )
        bottom.axhline(0.0, color="#111827", linewidth=0.8)
        bottom.set_ylabel("Tumor − Control cluster frequency (pp)" if column == 0 else "")
        bottom.set_xticks(x_positions, cluster_order, rotation=45, ha="right")
        bottom.grid(axis="y", color="#E2E8F0", linewidth=0.45)
        bottom.spines[["top", "right"]].set_visible(False)
    figure.suptitle(
        f"{RUN_DATE} denominator reconciliation — {title_suffix}\n"
        "The stacked-share zero-change reference is tissue-specific, not necessarily 50%",
        fontsize=14,
        fontweight="bold",
    )
    return save_figure(figure, path)


def vector(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def expression_and_counts(adata: ad.AnnData, symbol: str) -> tuple[np.ndarray, np.ndarray]:
    if adata.raw is None or symbol not in adata.raw.var_names:
        raise KeyError(f"{symbol} is unavailable in adata.raw")
    if symbol not in adata.var_names or "counts" not in adata.layers:
        raise KeyError(f"{symbol} is unavailable in the counts layer")
    expression = vector(adata.raw[:, symbol].X).astype(float, copy=False)
    counts = vector(adata[:, symbol].layers["counts"]).astype(float, copy=False)
    return expression, counts


def marker_complementarity_row(
    il4_expression: np.ndarray,
    il4_counts: np.ndarray,
    klrd1_expression: np.ndarray,
    klrd1_counts: np.ndarray,
    mask: np.ndarray,
    **metadata: Any,
) -> dict[str, Any]:
    il4_positive = il4_counts[mask] > 0
    klrd1_positive = klrd1_counts[mask] > 0
    both = int(np.sum(il4_positive & klrd1_positive))
    il4_only = int(np.sum(il4_positive & ~klrd1_positive))
    klrd1_only = int(np.sum(~il4_positive & klrd1_positive))
    neither = int(np.sum(~il4_positive & ~klrd1_positive))
    total = both + il4_only + klrd1_only + neither
    contingency = np.array([[both, il4_only], [klrd1_only, neither]], dtype=int)
    odds_ratio, fisher_p = fisher_exact(contingency)
    denominator = math.sqrt(
        max(
            (both + il4_only)
            * (klrd1_only + neither)
            * (both + klrd1_only)
            * (il4_only + neither),
            0,
        )
    )
    phi = (
        (both * neither - il4_only * klrd1_only) / denominator
        if denominator
        else float("nan")
    )
    correlation = (
        spearmanr(il4_expression[mask], klrd1_expression[mask]).statistic
        if total > 1
        else float("nan")
    )
    return {
        **metadata,
        "n_cells": total,
        "both_positive": both,
        "il4_only": il4_only,
        "klrd1_only": klrd1_only,
        "neither": neither,
        "both_positive_pct": both / total * 100.0 if total else np.nan,
        "il4_only_pct": il4_only / total * 100.0 if total else np.nan,
        "klrd1_only_pct": klrd1_only / total * 100.0 if total else np.nan,
        "neither_pct": neither / total * 100.0 if total else np.nan,
        "positive_set_jaccard": (
            both / (both + il4_only + klrd1_only)
            if both + il4_only + klrd1_only
            else np.nan
        ),
        "detection_odds_ratio": float(odds_ratio),
        "fisher_exact_pvalue": float(fisher_p),
        "detection_phi": float(phi),
        "spearman_log1p_expression": float(correlation),
    }


def build_marker_complementarity(
    adata: ad.AnnData,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    il4_expression, il4_counts = expression_and_counts(adata, "Il4")
    klrd1_expression, klrd1_counts = expression_and_counts(adata, "Klrd1")
    condition = adata.obs["condition"].astype(str).to_numpy()
    tissue = adata.obs["tissue"].astype(str).to_numpy()
    refined = adata.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy()
    rows = [
        marker_complementarity_row(
            il4_expression,
            il4_counts,
            klrd1_expression,
            klrd1_counts,
            np.ones(adata.n_obs, dtype=bool),
            scope="all_cells",
            condition="all",
            tissue="all",
            cluster="all",
        )
    ]
    for condition_value in CONDITION_ORDER:
        rows.append(
            marker_complementarity_row(
                il4_expression,
                il4_counts,
                klrd1_expression,
                klrd1_counts,
                condition == condition_value,
                scope="condition",
                condition=condition_value,
                tissue="all",
                cluster="all",
            )
        )
    for condition_value in CONDITION_ORDER:
        for tissue_value in TISSUE_ORDER:
            rows.append(
                marker_complementarity_row(
                    il4_expression,
                    il4_counts,
                    klrd1_expression,
                    klrd1_counts,
                    (condition == condition_value) & (tissue == tissue_value),
                    scope="condition_tissue",
                    condition=condition_value,
                    tissue=tissue_value,
                    cluster="all",
                )
            )
    for condition_value in CONDITION_ORDER:
        for cluster_value in REFINED_ORDER:
            rows.append(
                marker_complementarity_row(
                    il4_expression,
                    il4_counts,
                    klrd1_expression,
                    klrd1_counts,
                    (condition == condition_value) & (refined == cluster_value),
                    scope="condition_cluster",
                    condition=condition_value,
                    tissue="all",
                    cluster=cluster_value,
                )
            )
    arrays = {
        "Il4_expression": il4_expression,
        "Il4_counts": il4_counts,
        "Klrd1_expression": klrd1_expression,
        "Klrd1_counts": klrd1_counts,
    }
    return pd.DataFrame(rows), arrays
def render_condition_masked_markers(
    adata: ad.AnnData,
    arrays: Mapping[str, np.ndarray],
    path: Path,
) -> tuple[Path, Path]:
    coordinates = np.asarray(adata.obsm["X_umap"])
    extent = axis_extent(coordinates)
    conditions = adata.obs["condition"].astype(str).to_numpy()
    refined = adata.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy()
    markers = [
        ("Il4", "IL-4 / Il4", "magma"),
        ("Klrd1", "CD94 / Klrd1", "viridis"),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(13.6, 10.0), sharex=True, sharey=True, layout="constrained")
    for row, (symbol, display, cmap) in enumerate(markers):
        values = np.asarray(arrays[f"{symbol}_expression"], dtype=float)
        positive_values = values[values > 0]
        vmax = max(float(np.quantile(positive_values, 0.99)), 1e-6) if len(positive_values) else 1.0
        norm = Normalize(vmin=0.0, vmax=vmax)
        for column, condition_value in enumerate(CONDITION_ORDER):
            axis = axes[row, column]
            active = conditions == condition_value
            faded = ~active
            positive = active & (values > 0)
            zero = active & ~positive
            axis.scatter(
                coordinates[faded, 0],
                coordinates[faded, 1],
                s=1.8,
                c="#E2E8F0",
                alpha=0.10,
                edgecolors="none",
                rasterized=True,
                zorder=1,
            )
            axis.scatter(
                coordinates[zero, 0],
                coordinates[zero, 1],
                s=2.1,
                c="#CBD5E1",
                alpha=0.48,
                edgecolors="none",
                rasterized=True,
                zorder=2,
            )
            order = np.argsort(values[positive])
            positive_indices = np.flatnonzero(positive)[order]
            scatter = axis.scatter(
                coordinates[positive_indices, 0],
                coordinates[positive_indices, 1],
                s=4.0,
                c=values[positive_indices],
                cmap=cmap,
                norm=norm,
                alpha=0.86,
                edgecolors="none",
                rasterized=True,
                zorder=3,
            )
            add_cluster_centroid_labels(
                axis,
                coordinates,
                refined,
                ["C5-1", "C5-2"],
                fontsize=5.5,
            )
            style_umap_axis(
                axis,
                extent,
                xlabel="UMAP1" if row == 1 else "",
                ylabel="UMAP2" if column == 0 else "",
            )
            detected = int(positive.sum())
            axis.set_title(
                f"{display} — {CONDITION_LABELS[condition_value]} active\n"
                f"{detected:,}/{int(active.sum()):,} cells detected; other condition faded",
                fontsize=10,
                fontweight="bold",
            )
            colorbar = figure.colorbar(scatter, ax=axis, fraction=0.040, pad=0.02)
            colorbar.set_label("log1p normalized expression", fontsize=7)
            colorbar.ax.tick_params(labelsize=6)
    figure.suptitle(
        f"{RUN_DATE} condition-masked UMAP expression landscapes\n"
        "These are transcriptomic UMAP coordinates, not physical tissue spatial coordinates",
        fontsize=15,
        fontweight="bold",
    )
    return save_figure(figure, path)


def render_marker_overlap_bars(
    complementarity: pd.DataFrame,
    path: Path,
) -> tuple[Path, Path]:
    subset = complementarity[complementarity["scope"].eq("condition_tissue")].copy()
    order = [
        (condition, tissue)
        for tissue in TISSUE_ORDER
        for condition in CONDITION_ORDER
    ]
    subset["_order"] = [
        order.index((str(condition), str(tissue)))
        for condition, tissue in zip(subset["condition"], subset["tissue"], strict=True)
    ]
    subset = subset.sort_values("_order")
    labels = [
        f"{CONDITION_LABELS[row.condition]}\n{TISSUE_LABELS[row.tissue]}"
        for row in subset.itertuples(index=False)
    ]
    categories = [
        ("il4_only_pct", "IL-4 only", "#D95F02"),
        ("both_positive_pct", "Both", "#7B2CBF"),
        ("klrd1_only_pct", "CD94 only", "#1B9E77"),
        ("neither_pct", "Neither", "#CBD5E1"),
    ]
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(12.8, 8.7),
        gridspec_kw={"height_ratios": [2.0, 1.0]},
        layout="constrained",
    )
    axis = axes[0]
    x_positions = np.arange(len(subset))
    bottom = np.zeros(len(subset))
    for column, label, color in categories:
        values = subset[column].to_numpy(dtype=float)
        axis.bar(
            x_positions,
            values,
            bottom=bottom,
            color=color,
            label=label,
            width=0.72,
            edgecolor="white",
            linewidth=0.45,
        )
        bottom += values
    axis.set_xticks(x_positions, labels)
    axis.set_ylim(0, 100)
    axis.set_ylabel("Cells (%)")
    axis.set_title(
        "A. IL-4/CD94 detection categories within each condition×tissue",
        loc="left",
        fontweight="bold",
    )
    axis.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    axis.grid(axis="y", color="#E2E8F0", linewidth=0.5)
    axis.spines[["top", "right"]].set_visible(False)

    axis = axes[1]
    condition_rows = (
        complementarity[complementarity["scope"].eq("condition")]
        .set_index("condition")
        .reindex(CONDITION_ORDER)
    )
    x_positions = np.arange(len(condition_rows))
    phi = condition_rows["detection_phi"].to_numpy(dtype=float)
    spearman = condition_rows["spearman_log1p_expression"].to_numpy(dtype=float)
    width = 0.34
    axis.bar(x_positions - width / 2, phi, width, color="#475569", label="Detection phi")
    axis.bar(
        x_positions + width / 2,
        spearman,
        width,
        color="#0F766E",
        label="Expression Spearman",
    )
    axis.axhline(0.0, color="#111827", linewidth=0.8)
    axis.set_xticks(x_positions, [CONDITION_LABELS[value] for value in CONDITION_ORDER])
    axis.set_ylabel("Association")
    axis.set_title(
        "B. Quantitative association (negative supports complementarity)",
        loc="left",
        fontweight="bold",
    )
    axis.legend(frameon=False, ncol=2)
    axis.grid(axis="y", color="#E2E8F0", linewidth=0.5)
    axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(
        f"{RUN_DATE} IL-4 / CD94 complementarity audit",
        fontsize=15,
        fontweight="bold",
    )
    return save_figure(figure, path)


def _matrix_summary(matrix: Any, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    selected = matrix[mask]
    means = np.asarray(selected.mean(axis=0)).ravel()
    if sparse.issparse(selected):
        fractions = np.asarray((selected > 0).mean(axis=0)).ravel()
    else:
        fractions = np.mean(np.asarray(selected) > 0, axis=0)
    return means, fractions


def run_de_unit(adata: ad.AnnData, unit: DEUnit) -> pd.DataFrame:
    cluster_values = adata.obs[REFINED_CLUSTER_KEY].astype(str).to_numpy()
    mask = cluster_values == unit.cluster
    if unit.tissue is not None:
        mask &= adata.obs["tissue"].astype(str).to_numpy() == unit.tissue
    sub = adata[mask].copy()
    condition = sub.obs["condition"].astype(str).to_numpy()
    n_tumor = int(np.sum(condition == "T2"))
    n_control = int(np.sum(condition == "Ctrl"))
    if min(n_tumor, n_control) < EXPLORATORY_MIN_CELLS_PER_CONDITION:
        raise RuntimeError(
            f"{unit.unit_id} has insufficient cells: tumor={n_tumor}, control={n_control}"
        )
    if sub.raw is None:
        raise RuntimeError("Full-gene adata.raw is required for DE")
    sub.obs["_de_condition"] = pd.Categorical(condition, categories=["Ctrl", "T2"])
    key = f"_de_{re.sub('[^A-Za-z0-9_]+', '_', unit.unit_id)}"
    sc.tl.rank_genes_groups(
        sub,
        groupby="_de_condition",
        groups=["T2"],
        reference="Ctrl",
        method="wilcoxon",
        use_raw=True,
        n_genes=sub.raw.n_vars,
        corr_method="benjamini-hochberg",
        tie_correct=True,
        pts=True,
        key_added=key,
    )
    frame = sc.get.rank_genes_groups_df(sub, group="T2", key=key).rename(
        columns={"names": "gene"}
    )
    raw_names = pd.Index(sub.raw.var_names.astype(str))
    mean_tumor, pct_tumor = _matrix_summary(sub.raw.X, condition == "T2")
    mean_control, pct_control = _matrix_summary(sub.raw.X, condition == "Ctrl")
    summaries = pd.DataFrame(
        {
            "gene": raw_names,
            "mean_log1p_expression_tumor": mean_tumor,
            "mean_log1p_expression_control": mean_control,
            "pct_expressing_tumor": pct_tumor,
            "pct_expressing_control": pct_control,
        }
    )
    frame = frame.merge(summaries, on="gene", how="left", validate="one_to_one")
    frame.insert(0, "unit_id", unit.unit_id)
    frame.insert(1, "scope", unit.scope)
    frame.insert(2, "cluster", unit.cluster)
    frame.insert(3, "tissue", unit.tissue if unit.tissue is not None else "pooled_all_tissues")
    frame.insert(4, "contrast", "Tumor_T2_vs_Control_Ctrl")
    frame.insert(5, "n_tumor", n_tumor)
    frame.insert(6, "n_control", n_control)
    frame["paper_venn_tumor_up"] = (
        pd.to_numeric(frame["logfoldchanges"], errors="coerce").ge(1.0)
        & pd.to_numeric(frame["pvals"], errors="coerce").lt(0.05)
    )
    frame["robust_tumor_up"] = (
        pd.to_numeric(frame["logfoldchanges"], errors="coerce").ge(0.25)
        & pd.to_numeric(frame["pvals_adj"], errors="coerce").le(0.05)
    )
    frame["robust_control_up"] = (
        pd.to_numeric(frame["logfoldchanges"], errors="coerce").le(-0.25)
        & pd.to_numeric(frame["pvals_adj"], errors="coerce").le(0.05)
    )
    passes_recommended_min20 = (
        min(n_tumor, n_control) >= RECOMMENDED_MIN_CELLS_PER_CONDITION
    )
    frame["passes_recommended_min20"] = passes_recommended_min20
    frame["strict_estimable_min20"] = passes_recommended_min20
    frame["sample_size_note"] = (
        "passes_recommended_min20"
        if passes_recommended_min20
        else "exploratory_below_recommended_min20"
    )
    frame["statistical_unit_warning"] = STATISTICAL_WARNING
    return frame.sort_values(["pvals_adj", "pvals"], kind="stable", na_position="last")


def run_requested_de(adata: ad.AnnData, unit_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    units = [
        DEUnit(
            unit_id=f"pooled_cluster__{cluster}",
            scope="pooled_cluster",
            cluster=cluster,
        )
        for cluster in GSEA_CLUSTERS
    ]
    units.extend(
        DEUnit(
            unit_id=f"C5-1__{tissue}",
            scope="C5-1_tissue",
            cluster="C5-1",
            tissue=tissue,
        )
        for tissue in FIGURE_F_MAIN_TISSUES
    )
    units.extend(
        DEUnit(
            unit_id=f"C5-2__{tissue}",
            scope="C5-2_tissue",
            cluster="C5-2",
            tissue=tissue,
        )
        for tissue in TISSUE_ORDER
    )
    unit_dir.mkdir(parents=True, exist_ok=True)
    frames: list[pd.DataFrame] = []
    status_rows: list[dict[str, Any]] = []
    for unit in units:
        print(f"Running DE: {unit.unit_id}", flush=True)
        try:
            frame = run_de_unit(adata, unit)
        except Exception as exc:
            status_rows.append(
                {
                    **asdict(unit),
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                    "passes_recommended_min20": False,
                    "strict_estimable_min20": False,
                    "sample_size_note": "failed_not_estimable",
                }
            )
            continue
        output_path = unit_dir / f"{RUN_DATE}_{unit.unit_id.replace('__', '_')}_DE.csv.gz"
        frame.to_csv(output_path, index=False, compression="gzip")
        frames.append(frame)
        status_rows.append(
            {
                **asdict(unit),
                "status": "completed",
                "n_tumor": int(frame["n_tumor"].iloc[0]),
                "n_control": int(frame["n_control"].iloc[0]),
                "n_genes": int(len(frame)),
                "n_paper_tumor_up": int(frame["paper_venn_tumor_up"].sum()),
                "n_robust_tumor_up": int(frame["robust_tumor_up"].sum()),
                "passes_recommended_min20": bool(frame["passes_recommended_min20"].iloc[0]),
                "strict_estimable_min20": bool(frame["strict_estimable_min20"].iloc[0]),
                "sample_size_note": str(frame["sample_size_note"].iloc[0]),
                "n_robust_control_up": int(frame["robust_control_up"].sum()),
                "output": str(output_path),
            }
        )
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    status = pd.DataFrame(status_rows)
    requested_pooled = set(GSEA_CLUSTERS)
    completed_pooled = set(
        status.loc[
            status["scope"].eq("pooled_cluster") & status["status"].eq("completed"),
            "cluster",
        ]
    )
    if completed_pooled != requested_pooled:
        raise RuntimeError(
            f"Pooled DE units incomplete: requested={sorted(requested_pooled)}, "
            f"completed={sorted(completed_pooled)}"
        )
    requested_c5_2_tissues = set(TISSUE_ORDER)
    completed_c5_2_tissues = set(
        status.loc[
            status["scope"].eq("C5-2_tissue") & status["status"].eq("completed"),
            "tissue",
        ]
    )
    if completed_c5_2_tissues != requested_c5_2_tissues:
        raise RuntimeError(
            "C5-2 tissue DE units incomplete; Venn must not convert failures "
            f"to empty sets: requested={sorted(requested_c5_2_tissues)}, "
            f"completed={sorted(completed_c5_2_tissues)}"
        )
    requested_c5_1_tissues = set(FIGURE_F_MAIN_TISSUES)
    completed_c5_1_tissues = set(
        status.loc[
            status["scope"].eq("C5-1_tissue") & status["status"].eq("completed"),
            "tissue",
        ]
    )
    if completed_c5_1_tissues != requested_c5_1_tissues:
        raise RuntimeError(
            "C5-1 BM/Spleen DE units incomplete; main Figure 3F must not "
            f"convert failures to empty sets: requested={sorted(requested_c5_1_tissues)}, "
            f"completed={sorted(completed_c5_1_tissues)}"
        )
    return combined, status


def render_volcano_figure_d(de: pd.DataFrame, path: Path) -> tuple[Path, Path]:
    pooled = de[de["scope"].eq("pooled_cluster")].copy()
    figure, axes = plt.subplots(2, 2, figsize=(14.2, 11.0), layout="constrained")
    for axis, cluster in zip(axes.flat, TARGET_CLUSTERS, strict=True):
        frame = pooled[pooled["cluster"].eq(cluster)].copy()
        x_values = pd.to_numeric(frame["logfoldchanges"], errors="coerce").to_numpy(dtype=float)
        raw_p = pd.to_numeric(frame["pvals"], errors="coerce").to_numpy(dtype=float)
        raw_p = np.clip(raw_p, np.finfo(float).tiny, 1.0)
        y_values = np.minimum(-np.log10(raw_p), 50.0)
        finite = np.isfinite(x_values) & np.isfinite(y_values)
        control_up = frame["robust_control_up"].to_numpy(dtype=bool) & finite
        tumor_up = frame["robust_tumor_up"].to_numpy(dtype=bool) & finite
        other = finite & ~control_up & ~tumor_up
        axis.scatter(
            x_values[other],
            y_values[other],
            s=5,
            c="#CBD5E1",
            alpha=0.55,
            edgecolors="none",
            rasterized=True,
        )
        axis.scatter(
            x_values[control_up],
            y_values[control_up],
            s=8,
            c="#457B9D",
            alpha=0.72,
            edgecolors="none",
            rasterized=True,
            label="Control enriched",
        )
        axis.scatter(
            x_values[tumor_up],
            y_values[tumor_up],
            s=8,
            c="#B23A48",
            alpha=0.72,
            edgecolors="none",
            rasterized=True,
            label="Tumor enriched",
        )
        axis.axvline(0.0, color="#64748B", linewidth=0.6)
        axis.axvline(-0.25, color="#94A3B8", linestyle=":", linewidth=0.8)
        axis.axvline(0.25, color="#94A3B8", linestyle=":", linewidth=0.8)
        label_frames: list[pd.DataFrame] = []
        for direction_mask in (control_up, tumor_up):
            candidate_direction = frame.loc[direction_mask].copy()
            candidate_direction["_label_score"] = (
                -np.log10(
                    np.clip(
                        pd.to_numeric(
                            candidate_direction["pvals_adj"], errors="coerce"
                        ).to_numpy(dtype=float),
                        np.finfo(float).tiny,
                        1.0,
                    )
                )
                + pd.to_numeric(
                    candidate_direction["logfoldchanges"], errors="coerce"
                ).abs().to_numpy(dtype=float)
            )
            label_frames.append(candidate_direction.nlargest(4, "_label_score"))
        candidates = pd.concat(label_frames, ignore_index=True).drop_duplicates("gene")
        texts = []
        for row in candidates.itertuples(index=False):
            x_value = float(row.logfoldchanges)
            y_value = min(-math.log10(max(float(row.pvals), np.finfo(float).tiny)), 50.0)
            texts.append(
                axis.text(
                    np.clip(x_value, -7.8, 7.8),
                    y_value,
                    str(row.gene),
                    fontsize=6.5,
                    color="#1E293B",
                )
            )
        axis.set_xlim(-8, 8)
        axis.set_ylim(0, 52)
        if texts:
            adjust_text(
                texts,
                ax=axis,
                ensure_inside_axes=True,
                arrowprops={"arrowstyle": "-", "color": "#64748B", "lw": 0.4},
            )
        axis.set_xlabel("log2 fold change (Tumor / Control)")
        axis.set_ylabel("−log10 nominal P")
        n_tumor = int(frame["n_tumor"].iloc[0])
        n_control = int(frame["n_control"].iloc[0])
        axis.set_title(
            f"{cluster}: Tumor n={n_tumor:,}, Control n={n_control:,}\n"
            f"robust up: Tumor {int(tumor_up.sum()):,}, Control {int(control_up.sum()):,}",
            fontweight="bold",
        )
        axis.grid(color="#E2E8F0", linewidth=0.4)
        axis.spines[["top", "right"]].set_visible(False)
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#457B9D", markeredgecolor="none", label="Control enriched"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#B23A48", markeredgecolor="none", label="Tumor enriched"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#CBD5E1", markeredgecolor="none", label="Other"),
    ]
    figure.legend(handles=handles, frameon=False, ncol=3, loc="lower center")
    figure.suptitle(
        f"{RUN_DATE} paper-adapted Figure 3D — pooled tissue, cell-level Wilcoxon\n"
        "Color: BH-FDR≤0.05 and |log2FC|≥0.25; dotted guides: log2FC=±0.25; y-axis: nominal P",
        fontsize=15,
        fontweight="bold",
    )
    return save_figure(figure, path)


def read_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = {gene for gene in fields[2:] if gene}
    return gene_sets


def write_gmt(path: Path, gene_sets: Mapping[str, Iterable[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for term in sorted(gene_sets):
            members = sorted({str(gene) for gene in gene_sets[term] if str(gene)})
            handle.write("\t".join([str(term), "Enrichr", *members]) + "\n")


def load_kegg_gene_sets(
    library: str,
    output_gmt: Path,
    supplied_gmt: Path | None,
) -> tuple[dict[str, set[str]], dict[str, Any]]:
    if supplied_gmt is not None:
        if not supplied_gmt.exists():
            raise FileNotFoundError(supplied_gmt)
        gene_sets = read_gmt(supplied_gmt)
        source = supplied_gmt
        status = "supplied_offline_gmt"
    else:
        import gseapy as gp

        downloaded = gp.get_library(name=library, organism="Mouse")
        gene_sets = {str(term): set(map(str, genes)) for term, genes in downloaded.items()}
        write_gmt(output_gmt, gene_sets)
        source = output_gmt
        status = "downloaded_from_enrichr_and_frozen"
    if not gene_sets:
        raise RuntimeError("No KEGG gene sets were loaded")
    audit = {
        "library": library,
        "status": status,
        "source": str(source),
        "sha256": sha256_file(source),
        "n_terms": len(gene_sets),
        "n_unique_genes": len({gene.casefold() for genes in gene_sets.values() for gene in genes}),
    }
    return gene_sets, audit


def odds_ratio_haldane(
    overlap: int,
    selected_size: int,
    pathway_size: int,
    universe_size: int,
) -> float:
    a = overlap + 0.5
    b = selected_size - overlap + 0.5
    c = pathway_size - overlap + 0.5
    d = universe_size - selected_size - pathway_size + overlap + 0.5
    return float((a * d) / (b * c))


def run_kegg_ora(
    de: pd.DataFrame,
    tested_genes: Sequence[str],
    gene_sets: Mapping[str, set[str]],
) -> pd.DataFrame:
    universe_map: dict[str, str] = {}
    for gene in tested_genes:
        universe_map.setdefault(str(gene).casefold(), str(gene))
    universe = set(universe_map)
    pathway_case = {
        term: {str(gene).casefold() for gene in members} & universe
        for term, members in gene_sets.items()
    }
    thresholds = {
        ROBUST_THRESHOLD: lambda frame: frame["robust_tumor_up"].astype(bool),
        PAPER_VENN_THRESHOLD: lambda frame: frame["paper_venn_tumor_up"].astype(bool),
    }
    rows: list[dict[str, Any]] = []
    pooled = de[de["scope"].eq("pooled_cluster")]
    for cluster in TARGET_CLUSTERS:
        frame = pooled[pooled["cluster"].eq(cluster)]
        if frame.empty:
            continue
        for threshold_name, selector in thresholds.items():
            selected_original = set(frame.loc[selector(frame), "gene"].dropna().astype(str))
            selected = {gene.casefold() for gene in selected_original if gene.casefold() in universe}
            for term, members in pathway_case.items():
                if not 5 <= len(members) <= 500:
                    continue
                overlap = selected & members
                pvalue = (
                    float(
                        hypergeom.sf(
                            len(overlap) - 1,
                            len(universe),
                            len(members),
                            len(selected),
                        )
                    )
                    if selected
                    else 1.0
                )
                rows.append(
                    {
                        "cluster": cluster,
                        "threshold": threshold_name,
                        "term": term,
                        "universe_size": len(universe),
                        "selected_tumor_up_size": len(selected),
                        "pathway_size_in_universe": len(members),
                        "overlap_size": len(overlap),
                        "overlap_genes": ";".join(
                            sorted(universe_map[gene] for gene in overlap)
                        ),
                        "fold_enrichment": (
                            (len(overlap) / len(selected)) / (len(members) / len(universe))
                            if selected and members
                            else np.nan
                        ),
                        "odds_ratio_haldane": odds_ratio_haldane(
                            len(overlap),
                            len(selected),
                            len(members),
                            len(universe),
                        ),
                        "pvalue": pvalue,
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr"] = result.groupby(
        ["cluster", "threshold"], observed=True
    )["pvalue"].transform(lambda values: bh_adjust(values.to_numpy(dtype=float)))
    result["significant_fdr_0_05"] = result["fdr"].le(0.05)
    return result.sort_values(
        ["threshold", "cluster", "fdr", "pvalue"],
        kind="stable",
    ).reset_index(drop=True)


def select_pathways_for_plot(ora: pd.DataFrame, *, max_per_cluster: int = 5) -> list[str]:
    robust = ora[
        ora["threshold"].eq(ROBUST_THRESHOLD) & ora["overlap_size"].gt(0)
    ]
    chosen: list[str] = []
    for cluster in TARGET_CLUSTERS:
        frame = robust[robust["cluster"].eq(cluster)].sort_values(
            ["significant_fdr_0_05", "fdr", "pvalue", "overlap_size"],
            ascending=[False, True, True, False],
            kind="stable",
        )
        chosen.extend(frame.head(max_per_cluster)["term"].astype(str))
    deduplicated = list(dict.fromkeys(chosen))
    return deduplicated[:20]


def render_pathway_figure_e(ora: pd.DataFrame, path: Path) -> tuple[Path, Path]:
    terms = select_pathways_for_plot(ora)
    robust = ora[
        ora["threshold"].eq(ROBUST_THRESHOLD)
        & ora["term"].isin(terms)
        & ora["overlap_size"].gt(0)
    ].copy()
    figure, axis = plt.subplots(
        figsize=(11.5, max(6.2, 0.42 * max(len(terms), 1) + 2.2)),
        layout="constrained",
    )
    if not terms or robust.empty:
        axis.text(
            0.5,
            0.5,
            "No KEGG pathways available for plotting",
            ha="center",
            va="center",
            transform=axis.transAxes,
        )
        axis.set_axis_off()
    else:
        y_lookup = {term: index for index, term in enumerate(reversed(terms))}
        x_lookup = {cluster: index for index, cluster in enumerate(TARGET_CLUSTERS)}
        x_values = robust["cluster"].map(x_lookup).to_numpy(dtype=float)
        y_values = robust["term"].map(y_lookup).to_numpy(dtype=float)
        color_values = np.minimum(
            -np.log10(
                np.clip(
                    robust["fdr"].to_numpy(dtype=float),
                    np.finfo(float).tiny,
                    1.0,
                )
            ),
            12.0,
        )
        sizes = 24 + 16 * np.sqrt(robust["overlap_size"].to_numpy(dtype=float))
        scatter = axis.scatter(
            x_values,
            y_values,
            s=sizes,
            c=color_values,
            cmap="magma_r",
            vmin=0,
            vmax=max(3.0, float(np.nanmax(color_values))),
            edgecolors=np.where(robust["significant_fdr_0_05"], "#111827", "#94A3B8"),
            linewidths=np.where(robust["significant_fdr_0_05"], 0.9, 0.35),
        )
        axis.set_xticks(range(len(TARGET_CLUSTERS)), TARGET_CLUSTERS)
        axis.set_yticks(range(len(terms)), list(reversed(terms)), fontsize=8)
        axis.set_xlabel("Refined cluster")
        axis.set_ylabel("KEGG pathway")
        axis.grid(color="#E2E8F0", linewidth=0.5)
        axis.spines[["top", "right"]].set_visible(False)
        colorbar = figure.colorbar(scatter, ax=axis, fraction=0.035, pad=0.02)
        colorbar.set_label("−log10 pathway BH-FDR")
        size_handles = [
            axis.scatter([], [], s=24 + 16 * math.sqrt(value), c="#CBD5E1", edgecolors="#64748B", label=f"{value} genes")
            for value in (2, 5, 10)
        ]
        axis.legend(handles=size_handles, title="Overlap", frameon=False, loc="lower right")
    axis.set_title(
        f"{RUN_DATE} paper-adapted Figure 3E — Tumor-up mouse KEGG ORA\n"
        "Input genes: BH-FDR≤0.05 and log2FC≥0.25; black edge = pathway BH-FDR≤0.05",
        fontsize=13,
        fontweight="bold",
    )
    return save_figure(figure, path)


def build_gsea_ranking(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a deterministic, full-gene prerank table for one cluster."""
    required = {"gene", "scores", "logfoldchanges"}
    missing = required - set(frame.columns)
    if missing:
        raise KeyError(f"GSEA ranking is missing columns: {sorted(missing)}")
    ranking = frame.loc[:, ["gene", "scores", "logfoldchanges"]].copy()
    ranking["gene"] = ranking["gene"].fillna("").astype(str).str.strip()
    ranking["rank_score"] = pd.to_numeric(ranking["scores"], errors="coerce")
    ranking["logfoldchanges"] = pd.to_numeric(
        ranking["logfoldchanges"], errors="coerce"
    )
    ranking["gene_key"] = ranking["gene"].str.casefold()
    valid = ranking["gene"].ne("") & np.isfinite(
        ranking["rank_score"].to_numpy(dtype=float)
    )
    n_input_rows = int(len(ranking))
    ranking = ranking.loc[valid].copy()
    if ranking.empty:
        raise RuntimeError("No finite genes remain for GSEA ranking")

    collision_sizes = ranking.groupby("gene_key", observed=True).size()
    collision_keys = set(collision_sizes[collision_sizes.gt(1)].index.astype(str))
    ranking["_abs_score"] = ranking["rank_score"].abs()
    ranking = ranking.sort_values(
        ["gene_key", "_abs_score", "rank_score", "logfoldchanges", "gene"],
        ascending=[True, False, False, False, True],
        kind="stable",
        na_position="last",
    ).drop_duplicates("gene_key", keep="first")
    ranking = ranking.sort_values(
        ["rank_score", "logfoldchanges", "gene"],
        ascending=[False, False, True],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)
    if ranking["gene_key"].duplicated().any():
        raise RuntimeError("Case-insensitive duplicate genes remain in GSEA ranking")

    tie_sizes = ranking.groupby("rank_score", observed=True).size()
    ranking.insert(0, "rank_position", np.arange(1, len(ranking) + 1))
    audit = {
        "n_input_rows": n_input_rows,
        "n_ranked_genes": int(len(ranking)),
        "n_dropped_blank_or_nonfinite": int(n_input_rows - int(valid.sum())),
        "n_casefold_collision_keys": int(len(collision_keys)),
        "n_casefold_collision_rows": int(
            collision_sizes[collision_sizes.gt(1)].sum()
        ),
        "n_unique_scores": int(ranking["rank_score"].nunique()),
        "n_tied_score_groups": int(tie_sizes.gt(1).sum()),
        "n_genes_in_tied_score_groups": int(tie_sizes[tie_sizes.gt(1)].sum()),
        "min_rank_score": float(ranking["rank_score"].min()),
        "max_rank_score": float(ranking["rank_score"].max()),
        "has_positive_scores": bool(ranking["rank_score"].gt(0).any()),
        "has_negative_scores": bool(ranking["rank_score"].lt(0).any()),
        "sort_rule": "rank_score_desc; log2fc_desc; gene_asc; stable",
        "rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
    }
    return (
        ranking.loc[
            :,
            [
                "rank_position",
                "gene",
                "gene_key",
                "rank_score",
                "logfoldchanges",
            ],
        ],
        audit,
    )


def match_gene_sets_to_ranking(
    gene_sets: Mapping[str, set[str]],
    ranked_genes: Sequence[str],
) -> dict[str, set[str]]:
    """Case-insensitively map GMT members to the measured gene spelling."""
    measured: dict[str, str] = {}
    for gene in ranked_genes:
        measured.setdefault(str(gene).casefold(), str(gene))
    matched: dict[str, set[str]] = {}
    for term in sorted(gene_sets):
        members = {
            measured[key]
            for key in (str(gene).casefold() for gene in gene_sets[term])
            if key in measured
        }
        matched[str(term)] = members
    return matched


def run_preranked_gsea(
    de: pd.DataFrame,
    gene_sets: Mapping[str, set[str]],
    *,
    permutation_num: int,
    seed: int,
    threads: int,
    min_size: int = GSEA_MIN_SIZE,
    max_size: int = GSEA_MAX_SIZE,
    weight: float = GSEA_WEIGHT,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Run full-gene KEGG preranked GSEA for all requested refined clusters."""
    import gseapy as gp

    pooled = de[de["scope"].eq("pooled_cluster")].copy()
    observed = set(pooled["cluster"].dropna().astype(str))
    missing = set(map(str, clusters)) - observed
    if missing:
        raise RuntimeError(f"Missing pooled DE rankings for GSEA: {sorted(missing)}")

    result_frames: list[pd.DataFrame] = []
    ranking_frames: list[pd.DataFrame] = []
    ranking_audits: list[dict[str, Any]] = []
    for cluster in clusters:
        print(f"Running preranked KEGG GSEA: {cluster}", flush=True)
        cluster_frame = pooled[pooled["cluster"].eq(cluster)].copy()
        ranking, ranking_audit = build_gsea_ranking(cluster_frame)
        ranking.insert(0, "cluster", str(cluster))
        ranking["score_direction"] = np.select(
            [ranking["rank_score"].gt(0), ranking["rank_score"].lt(0)],
            ["Tumor", "Control"],
            default="Zero",
        )
        ranking_frames.append(ranking)
        ranking_audits.append({"cluster": str(cluster), **ranking_audit})

        matched_gene_sets = match_gene_sets_to_ranking(
            gene_sets,
            ranking["gene"].astype(str).tolist(),
        )
        prerank = gp.prerank(
            rnk=ranking.loc[:, ["gene", "rank_score"]],
            gene_sets=matched_gene_sets,
            outdir=None,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutation_num,
            weight=weight,
            ascending=None,
            threads=threads,
            no_plot=True,
            seed=seed,
            verbose=False,
            method="permutation",
        )
        cluster_result = prerank.res2d.copy()
        required_columns = {
            "Term",
            "ES",
            "NES",
            "NOM p-val",
            "FDR q-val",
            "FWER p-val",
            "Tag %",
            "Gene %",
            "Lead_genes",
        }
        missing_columns = required_columns - set(cluster_result.columns)
        if missing_columns:
            raise RuntimeError(
                "Unexpected GSEApy output; missing columns: "
                f"{sorted(missing_columns)}"
            )
        cluster_result = cluster_result.rename(
            columns={
                "Term": "term",
                "ES": "es",
                "NES": "nes",
                "NOM p-val": "nominal_pvalue",
                "FDR q-val": "fdr",
                "FWER p-val": "fwer_pvalue",
                "Tag %": "tag_fraction",
                "Gene %": "gene_fraction",
                "Lead_genes": "leading_edge_genes",
            }
        )
        for column in ("es", "nes", "nominal_pvalue", "fdr", "fwer_pvalue"):
            cluster_result[column] = pd.to_numeric(
                cluster_result[column], errors="coerce"
            )
        cluster_result["term"] = cluster_result["term"].astype(str)
        cluster_result["leading_edge_genes"] = (
            cluster_result["leading_edge_genes"].fillna("").astype(str)
        )
        pathway_sizes = {
            term: len(members) for term, members in matched_gene_sets.items()
        }
        cluster_result.insert(0, "cluster", str(cluster))
        cluster_result["direction"] = np.select(
            [cluster_result["nes"].gt(0), cluster_result["nes"].lt(0)],
            ["Tumor", "Control"],
            default="Zero",
        )
        cluster_result["significant_fdr_0_05"] = cluster_result["fdr"].le(
            GSEA_FDR_THRESHOLD
        )
        cluster_result["pathway_size_in_ranking"] = (
            cluster_result["term"].map(pathway_sizes).astype("Int64")
        )
        cluster_result["leading_edge_size"] = cluster_result[
            "leading_edge_genes"
        ].map(lambda value: len([gene for gene in value.split(";") if gene]))
        cluster_result["n_ranked_genes"] = int(len(ranking))
        cluster_result["n_tumor"] = int(cluster_frame["n_tumor"].iloc[0])
        cluster_result["n_control"] = int(cluster_frame["n_control"].iloc[0])
        cluster_result["rank_metric"] = "scanpy_wilcoxon_score_T2_vs_Ctrl"
        cluster_result["permutation_type"] = "gene_set"
        cluster_result["permutation_num"] = int(permutation_num)
        cluster_result["weight"] = float(weight)
        cluster_result["seed"] = int(seed)
        cluster_result["statistical_unit_warning"] = STATISTICAL_WARNING

        finite = np.isfinite(
            cluster_result[["nes", "fdr"]].to_numpy(dtype=float)
        ).all(axis=1)
        if not bool(finite.all()):
            raise RuntimeError(f"{cluster} GSEA returned non-finite NES/FDR values")
        if not cluster_result["fdr"].between(0.0, 1.0).all():
            raise RuntimeError(f"{cluster} GSEA returned FDR outside [0, 1]")
        result_frames.append(cluster_result)

    result = pd.concat(result_frames, ignore_index=True)
    rankings = pd.concat(ranking_frames, ignore_index=True)
    completed = set(result["cluster"].astype(str))
    if completed != set(map(str, clusters)):
        raise RuntimeError(
            f"GSEA incomplete: requested={list(clusters)}, completed={sorted(completed)}"
        )
    result["_cluster_order"] = pd.Categorical(
        result["cluster"], categories=list(clusters), ordered=True
    )
    result = result.sort_values(
        ["_cluster_order", "fdr", "nes", "term"],
        ascending=[True, True, False, True],
        kind="stable",
    ).drop(columns="_cluster_order").reset_index(drop=True)
    audit = {
        "software": {"gseapy": str(gp.__version__)},
        "parameters": {
            "clusters": list(map(str, clusters)),
            "rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
            "positive_direction": "Tumor_T2",
            "negative_direction": "Control_Ctrl",
            "gene_sets": "KEGG_2019_Mouse_frozen_GMT",
            "min_size": int(min_size),
            "max_size": int(max_size),
            "permutation_num": int(permutation_num),
            "permutation_type": "gene_set",
            "weight": float(weight),
            "seed": int(seed),
            "threads": int(threads),
            "ascending": None,
            "fdr_threshold": float(GSEA_FDR_THRESHOLD),
        },
        "rankings": ranking_audits,
        "n_results": int(len(result)),
        "n_pathways_by_cluster": {
            str(cluster): int(result["cluster"].eq(cluster).sum())
            for cluster in clusters
        },
        "statistical_warning": STATISTICAL_WARNING,
    }
    return result, rankings, audit


def gsea_computation_identity(
    *,
    permutation_num: int,
    seed: int,
    min_size: int = GSEA_MIN_SIZE,
    max_size: int = GSEA_MAX_SIZE,
    weight: float = GSEA_WEIGHT,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> dict[str, Any]:
    """Return the result-defining GSEA parameters used for cache identity."""
    return {
        "clusters": list(map(str, clusters)),
        "rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
        "positive_direction": "Tumor_T2",
        "negative_direction": "Control_Ctrl",
        "gene_sets": "KEGG_2019_Mouse_frozen_GMT",
        "min_size": int(min_size),
        "max_size": int(max_size),
        "permutation_num": int(permutation_num),
        "permutation_type": "gene_set",
        "weight": float(weight),
        "seed": int(seed),
        "ascending": None,
    }


def build_gsea_bundle_receipt(
    result_path: Path,
    ranking_path: Path,
    *,
    gsea: pd.DataFrame,
    rankings: pd.DataFrame,
    kegg_audit: Mapping[str, Any],
    identity: Mapping[str, Any],
    provenance: str,
) -> dict[str, Any]:
    """Bind the exact result and ranking bytes to their GMT and parameters."""
    gmt_sha256 = str(
        kegg_audit.get("frozen_copy_sha256", kegg_audit.get("sha256", ""))
    )
    if not gmt_sha256:
        raise RuntimeError("Cannot create GSEA receipt without frozen GMT SHA256")
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": str(provenance),
        "computation_identity": dict(identity),
        "frozen_gmt_sha256": gmt_sha256,
        "files": {
            "gsea_results": {
                "path": str(result_path.resolve()),
                "size_bytes": int(result_path.stat().st_size),
                "sha256": sha256_file(result_path),
                "n_rows": int(len(gsea)),
            },
            "ranked_genes": {
                "path": str(ranking_path.resolve()),
                "size_bytes": int(ranking_path.stat().st_size),
                "sha256": sha256_file(ranking_path),
                "n_rows": int(len(rankings)),
            },
        },
    }


def load_validated_gsea_bundle(
    de: pd.DataFrame,
    gene_sets: Mapping[str, set[str]],
    kegg_audit: Mapping[str, Any],
    tables_dir: Path,
    *,
    permutation_num: int,
    seed: int,
    threads: int,
    min_size: int = GSEA_MIN_SIZE,
    max_size: int = GSEA_MAX_SIZE,
    weight: float = GSEA_WEIGHT,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Load a dated GSEA bundle only after current-input, GMT, and schema checks."""
    result_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz"
    ranking_path = (
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz"
    )
    audit_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_audit.json"
    receipt_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_bundle_receipt.json"
    required_paths = (result_path, ranking_path, audit_path, receipt_path)
    missing_paths = [str(path) for path in required_paths if not path.is_file()]
    empty_paths = [
        str(path) for path in required_paths if path.is_file() and path.stat().st_size == 0
    ]
    if missing_paths or empty_paths:
        raise RuntimeError(
            f"Reusable GSEA bundle is incomplete; missing={missing_paths}, "
            f"empty={empty_paths}"
        )

    gsea = pd.read_csv(result_path)
    rankings = pd.read_csv(ranking_path)
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    expected_identity = gsea_computation_identity(
        permutation_num=permutation_num,
        seed=seed,
        min_size=min_size,
        max_size=max_size,
        weight=weight,
        clusters=clusters,
    )
    audit_parameters = audit.get("parameters")
    if not isinstance(audit_parameters, dict):
        raise RuntimeError("Reusable GSEA audit has no parameter dictionary")
    identity_mismatches = {
        key: {"expected": expected, "observed": audit_parameters.get(key)}
        for key, expected in expected_identity.items()
        if audit_parameters.get(key) != expected
    }
    if identity_mismatches:
        raise RuntimeError(
            "Reusable GSEA computation identity mismatch: "
            f"{identity_mismatches}"
        )

    current_gmt_sha256 = str(
        kegg_audit.get("frozen_copy_sha256", kegg_audit.get("sha256", ""))
    )
    audited_library = audit.get("gene_set_library", {})
    audited_gmt_sha256 = str(
        audited_library.get(
            "frozen_copy_sha256", audited_library.get("sha256", "")
        )
    )
    if not current_gmt_sha256 or current_gmt_sha256 != audited_gmt_sha256:
        raise RuntimeError(
            "Reusable GSEA frozen GMT SHA256 mismatch: "
            f"current={current_gmt_sha256}, audited={audited_gmt_sha256}"
        )

    ranking_required = {
        "cluster",
        "rank_position",
        "gene",
        "gene_key",
        "rank_score",
        "score_direction",
        "logfoldchanges",
    }
    missing_ranking_columns = ranking_required - set(rankings.columns)
    if missing_ranking_columns:
        raise RuntimeError(
            f"Reusable GSEA rankings missing columns: {sorted(missing_ranking_columns)}"
        )
    observed_ranking_order = rankings["cluster"].astype(str).drop_duplicates().tolist()
    if observed_ranking_order != list(map(str, clusters)):
        raise RuntimeError(
            "Reusable GSEA ranking cluster order mismatch: "
            f"{observed_ranking_order}"
        )

    pooled = de[de["scope"].eq("pooled_cluster")].copy()
    current_ranking_frames: list[pd.DataFrame] = []
    ranking_validation: list[dict[str, Any]] = []
    for cluster in clusters:
        cluster_de = pooled[pooled["cluster"].eq(cluster)].copy()
        if cluster_de.empty:
            raise RuntimeError(f"Current DE lacks pooled ranking for {cluster}")
        current, _ = build_gsea_ranking(cluster_de)
        current.insert(0, "cluster", str(cluster))
        current["score_direction"] = np.select(
            [current["rank_score"].gt(0), current["rank_score"].lt(0)],
            ["Tumor", "Control"],
            default="Zero",
        )
        cached = rankings[rankings["cluster"].eq(cluster)].reset_index(drop=True)
        current = current.reset_index(drop=True)
        if cached.empty:
            raise RuntimeError(f"Reusable GSEA ranking is empty for {cluster}")
        expected_positions = np.arange(1, len(cached) + 1)
        if not np.array_equal(
            pd.to_numeric(cached["rank_position"]).to_numpy(), expected_positions
        ):
            raise RuntimeError(f"Reusable GSEA rank positions are invalid for {cluster}")
        if cached["gene_key"].astype(str).duplicated().any():
            raise RuntimeError(f"Reusable GSEA gene keys are not unique for {cluster}")
        if len(cached) != len(current):
            raise RuntimeError(
                f"Reusable GSEA ranking length mismatch for {cluster}: "
                f"cached={len(cached)}, current={len(current)}"
            )
        for column in ("rank_position", "gene", "gene_key", "score_direction"):
            cached_values = cached[column].astype(str).to_numpy()
            current_values = current[column].astype(str).to_numpy()
            if not np.array_equal(cached_values, current_values):
                first = int(np.flatnonzero(cached_values != current_values)[0])
                raise RuntimeError(
                    f"Reusable GSEA ranking {column} mismatch for {cluster} "
                    f"at zero-based row {first}"
                )
        numeric_validation: dict[str, dict[str, Any]] = {}
        for column, tolerance in (
            ("rank_score", GSEA_RANK_SCORE_ATOL),
            ("logfoldchanges", GSEA_RANK_LOGFC_ATOL),
        ):
            cached_numeric = pd.to_numeric(
                cached[column], errors="coerce"
            ).to_numpy(dtype=float)
            current_numeric = pd.to_numeric(
                current[column], errors="coerce"
            ).to_numpy(dtype=float)
            cached_finite = np.isfinite(cached_numeric)
            current_finite = np.isfinite(current_numeric)
            finite_pattern_equal = np.array_equal(cached_finite, current_finite)
            absolute_delta = np.full(cached_numeric.shape, np.nan, dtype=float)
            jointly_finite = cached_finite & current_finite
            absolute_delta[jointly_finite] = np.abs(
                cached_numeric[jointly_finite] - current_numeric[jointly_finite]
            )
            max_abs_delta = (
                float(np.nanmax(absolute_delta)) if jointly_finite.any() else 0.0
            )
            canonical_equal = np.array_equal(
                cached_numeric.astype(np.float32),
                current_numeric.astype(np.float32),
                equal_nan=True,
            )
            bounded_delta = max_abs_delta <= tolerance
            if not finite_pattern_equal or not canonical_equal or not bounded_delta:
                mismatch_mask = np.zeros(cached_numeric.shape, dtype=bool)
                mismatch_mask[jointly_finite] = (
                    absolute_delta[jointly_finite] > tolerance
                    ) | (
                    cached_numeric[jointly_finite].astype(np.float32)
                    != current_numeric[jointly_finite].astype(np.float32)
                )
                mismatch_mask |= cached_finite != current_finite
                first = (
                    int(np.flatnonzero(mismatch_mask)[0])
                    if mismatch_mask.any()
                    else -1
                )
                raise RuntimeError(
                    f"Reusable GSEA ranking {column} mismatch for {cluster}; "
                    f"max_abs_delta={max_abs_delta:.12g}, "
                    f"tolerance={tolerance:.1g}, first_row={first}"
                )
            numeric_validation[column] = {
                "max_abs_delta_before_float32_canonicalization": max_abs_delta,
                "absolute_tolerance": float(tolerance),
                "float32_canonical_values_equal": True,
            }
        current_ranking_frames.append(current)
        ranking_validation.append(
            {
                "cluster": str(cluster),
                "n_ranked_genes": int(len(current)),
                "status": (
                    "exact_gene_order_sign_and_float32_canonical_numeric_match"
                ),
                "numeric_validation": numeric_validation,
            }
        )
    current_rankings = pd.concat(current_ranking_frames, ignore_index=True)

    result_required = {
        "cluster",
        "term",
        "es",
        "nes",
        "nominal_pvalue",
        "fdr",
        "fwer_pvalue",
        "leading_edge_genes",
        "direction",
        "significant_fdr_0_05",
        "pathway_size_in_ranking",
        "n_ranked_genes",
        "n_tumor",
        "n_control",
        "rank_metric",
        "permutation_type",
        "permutation_num",
        "weight",
        "seed",
    }
    missing_result_columns = result_required - set(gsea.columns)
    if missing_result_columns:
        raise RuntimeError(
            f"Reusable GSEA results missing columns: {sorted(missing_result_columns)}"
        )
    observed_result_order = gsea["cluster"].astype(str).drop_duplicates().tolist()
    if observed_result_order != list(map(str, clusters)):
        raise RuntimeError(
            f"Reusable GSEA result cluster order mismatch: {observed_result_order}"
        )
    if gsea.duplicated(["cluster", "term"]).any():
        raise RuntimeError("Reusable GSEA results contain duplicate cluster×term rows")
    for column in ("es", "nes", "nominal_pvalue", "fdr", "fwer_pvalue"):
        gsea[column] = pd.to_numeric(gsea[column], errors="coerce")
    if not np.isfinite(
        gsea[["es", "nes", "nominal_pvalue", "fdr", "fwer_pvalue"]].to_numpy(
            dtype=float
        )
    ).all():
        raise RuntimeError("Reusable GSEA results contain non-finite statistics")
    for column in ("nominal_pvalue", "fdr", "fwer_pvalue"):
        if not gsea[column].between(0.0, 1.0).all():
            raise RuntimeError(f"Reusable GSEA {column} is outside [0, 1]")
    expected_direction = np.select(
        [gsea["nes"].gt(0), gsea["nes"].lt(0)],
        ["Tumor", "Control"],
        default="Zero",
    )
    if not np.array_equal(gsea["direction"].astype(str).to_numpy(), expected_direction):
        raise RuntimeError("Reusable GSEA direction is inconsistent with signed NES")

    for column, expected in (
        ("rank_metric", expected_identity["rank_metric"]),
        ("permutation_type", expected_identity["permutation_type"]),
        ("permutation_num", expected_identity["permutation_num"]),
        ("weight", expected_identity["weight"]),
        ("seed", expected_identity["seed"]),
    ):
        values = set(gsea[column].dropna().tolist())
        if values != {expected}:
            raise RuntimeError(
                f"Reusable GSEA result column {column} mismatch: {values}"
            )

    for cluster in clusters:
        current = current_rankings[current_rankings["cluster"].eq(cluster)]
        cluster_result = gsea[gsea["cluster"].eq(cluster)]
        matched = match_gene_sets_to_ranking(
            gene_sets, current["gene"].astype(str).tolist()
        )
        eligible = {
            str(term): len(members)
            for term, members in matched.items()
            if min_size <= len(members) <= max_size
        }
        observed_terms = set(cluster_result["term"].astype(str))
        if observed_terms != set(eligible):
            raise RuntimeError(
                f"Reusable GSEA eligible-term mismatch for {cluster}; "
                f"missing={sorted(set(eligible) - observed_terms)}, "
                f"extra={sorted(observed_terms - set(eligible))}"
            )
        sizes = pd.to_numeric(
            cluster_result["pathway_size_in_ranking"], errors="coerce"
        )
        expected_sizes = cluster_result["term"].astype(str).map(eligible)
        if not np.array_equal(sizes.to_numpy(), expected_sizes.to_numpy()):
            raise RuntimeError(f"Reusable GSEA pathway sizes mismatch for {cluster}")
        if not cluster_result["n_ranked_genes"].eq(len(current)).all():
            raise RuntimeError(f"Reusable GSEA ranked-gene count mismatch for {cluster}")
        expected_n_tumor = int(
            pd.to_numeric(
                pooled.loc[pooled["cluster"].eq(cluster), "n_tumor"]
            ).iloc[0]
        )
        expected_n_control = int(
            pd.to_numeric(
                pooled.loc[pooled["cluster"].eq(cluster), "n_control"]
            ).iloc[0]
        )
        if (
            not cluster_result["n_tumor"].eq(expected_n_tumor).all()
            or not cluster_result["n_control"].eq(expected_n_control).all()
        ):
            raise RuntimeError(f"Reusable GSEA condition counts mismatch for {cluster}")
        for row in cluster_result.itertuples(index=False):
            leading_edge = {
                gene
                for gene in str(row.leading_edge_genes).split(";")
                if gene and gene.lower() != "nan"
            }
            if not leading_edge.issubset(matched[str(row.term)]):
                raise RuntimeError(
                    f"Reusable GSEA leading edge is invalid for {cluster}/{row.term}"
                )

    expected_results = sum(
        1
        for cluster in clusters
        for members in match_gene_sets_to_ranking(
            gene_sets,
            current_rankings.loc[
                current_rankings["cluster"].eq(cluster), "gene"
            ].astype(str),
        ).values()
        if min_size <= len(members) <= max_size
    )
    if len(gsea) != expected_results or int(audit.get("n_results", -1)) != len(gsea):
        raise RuntimeError(
            "Reusable GSEA total result count mismatch: "
            f"table={len(gsea)}, expected={expected_results}, "
            f"audit={audit.get('n_results')}"
        )
    audited_counts = audit.get("n_pathways_by_cluster", {})
    observed_counts = {
        str(cluster): int(gsea["cluster"].eq(cluster).sum()) for cluster in clusters
    }
    if audited_counts != observed_counts:
        raise RuntimeError(
            f"Reusable GSEA per-cluster result counts mismatch: {audited_counts}"
        )

    gsea["direction"] = expected_direction
    gsea["significant_fdr_0_05"] = gsea["fdr"].le(GSEA_FDR_THRESHOLD)
    result_sha256 = sha256_file(result_path)
    ranking_sha256 = sha256_file(ranking_path)
    receipt_status = "validated_existing_receipt"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt_mismatches: list[str] = []
    if receipt.get("schema_version") != 1:
        receipt_mismatches.append("schema_version")
    if receipt.get("computation_identity") != expected_identity:
        receipt_mismatches.append("computation_identity")
    if receipt.get("frozen_gmt_sha256") != current_gmt_sha256:
        receipt_mismatches.append("frozen_gmt_sha256")
    receipt_files = receipt.get("files", {})
    if receipt_files.get("gsea_results", {}).get("sha256") != result_sha256:
        receipt_mismatches.append("gsea_results.sha256")
    if receipt_files.get("ranked_genes", {}).get("sha256") != ranking_sha256:
        receipt_mismatches.append("ranked_genes.sha256")
    if receipt_mismatches:
        raise RuntimeError(
            f"Reusable GSEA bundle receipt mismatch: {receipt_mismatches}"
        )

    audit["reuse_validation"] = {
        "status": "passed",
        "reused": True,
        "receipt_status": receipt_status,
        "receipt_path": str(receipt_path.resolve()),
        "gsea_results_sha256": result_sha256,
        "ranked_genes_sha256": ranking_sha256,
        "current_frozen_gmt_sha256": current_gmt_sha256,
        "current_threads": int(threads),
        "current_fdr_threshold": float(GSEA_FDR_THRESHOLD),
        "ranking_validation": ranking_validation,
        "validated_fields": [
            "computation_identity",
            "frozen_GMT_sha256",
            "current_DE_gene_order_and_rank_score",
            "cluster_and_term_completeness",
            "eligible_term_universe",
            "numeric_ranges_and_signed_direction",
            "leading_edge_membership",
            "condition_cell_counts",
            "bundle_receipt_sha256_binding",
        ],
    }
    audit["gene_set_library"] = dict(kegg_audit)
    return gsea, rankings, audit


def select_gsea_heatmap_pathways(
    gsea: pd.DataFrame,
    *,
    max_pathways: int = GSEA_HEATMAP_MAX_PATHWAYS,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> pd.DataFrame:
    """Select shared and cluster-specific significant pathways deterministically."""
    required = {"cluster", "term", "nes", "fdr", "direction"}
    missing = required - set(gsea.columns)
    if missing:
        raise KeyError(f"GSEA selection is missing columns: {sorted(missing)}")
    work = gsea.copy()
    work["nes"] = pd.to_numeric(work["nes"], errors="coerce")
    work["fdr"] = pd.to_numeric(work["fdr"], errors="coerce")
    work = work[
        np.isfinite(work["nes"].to_numpy(dtype=float))
        & np.isfinite(work["fdr"].to_numpy(dtype=float))
    ].copy()
    significant = work[work["fdr"].le(GSEA_FDR_THRESHOLD)].copy()
    if significant.empty:
        return pd.DataFrame(
            columns=[
                "display_order",
                "term",
                "selection_reason",
                "classification",
                "n_significant_clusters",
                "n_tumor_significant_clusters",
                "n_control_significant_clusters",
                "tumor_significant_clusters",
                "control_significant_clusters",
                "min_fdr",
                "max_abs_nes",
            ]
        )

    cluster_order = {str(cluster): index for index, cluster in enumerate(clusters)}
    rows: list[dict[str, Any]] = []
    for term, frame in significant.groupby("term", sort=True, observed=True):
        tumor_clusters = tuple(
            cluster
            for cluster in clusters
            if (
                frame["cluster"].eq(cluster)
                & frame["direction"].eq("Tumor")
            ).any()
        )
        control_clusters = tuple(
            cluster
            for cluster in clusters
            if (
                frame["cluster"].eq(cluster)
                & frame["direction"].eq("Control")
            ).any()
        )
        significant_clusters = tuple(
            cluster
            for cluster in clusters
            if frame["cluster"].eq(cluster).any()
        )
        n_significant = len(significant_clusters)
        rows.append(
            {
                "term": str(term),
                "classification": (
                    "shared"
                    if n_significant >= 2
                    else f"cluster-specific:{significant_clusters[0]}"
                ),
                "n_significant_clusters": n_significant,
                "n_tumor_significant_clusters": len(tumor_clusters),
                "n_control_significant_clusters": len(control_clusters),
                "tumor_significant_clusters": ";".join(tumor_clusters),
                "control_significant_clusters": ";".join(control_clusters),
                "significant_clusters": ";".join(significant_clusters),
                "min_fdr": float(frame["fdr"].min()),
                "max_abs_nes": float(frame["nes"].abs().max()),
                "_specific_cluster_order": (
                    cluster_order[significant_clusters[0]]
                    if n_significant == 1
                    else len(cluster_order)
                ),
            }
        )
    metadata = pd.DataFrame(rows).sort_values(
        [
            "n_significant_clusters",
            "min_fdr",
            "max_abs_nes",
            "term",
        ],
        ascending=[False, True, False, True],
        kind="stable",
    )
    selected: list[str] = []
    reasons: dict[str, list[str]] = {}

    def add_term(term: str, reason: str) -> None:
        if term in reasons:
            if reason not in reasons[term]:
                reasons[term].append(reason)
            return
        if len(selected) >= max_pathways:
            return
        selected.append(term)
        reasons[term] = [reason]

    shared_limit = min(6, max_pathways)
    shared = metadata[metadata["n_significant_clusters"].ge(2)]
    for term in shared.head(shared_limit)["term"].astype(str):
        add_term(term, "shared_top")

    specific = metadata[metadata["n_significant_clusters"].eq(1)]
    for cluster in clusters:
        frame = specific[
            specific["significant_clusters"].eq(str(cluster))
        ].sort_values(
            ["min_fdr", "max_abs_nes", "term"],
            ascending=[True, False, True],
            kind="stable",
        )
        if not frame.empty:
            add_term(str(frame.iloc[0]["term"]), f"specific_{cluster}")

    for cluster in clusters:
        for direction in ("Tumor", "Control"):
            frame = significant[
                significant["cluster"].eq(cluster)
                & significant["direction"].eq(direction)
            ].copy()
            if frame.empty:
                continue
            frame["_abs_nes"] = frame["nes"].abs()
            frame = frame.sort_values(
                ["fdr", "_abs_nes", "term"],
                ascending=[True, False, True],
                kind="stable",
            )
            add_term(
                str(frame.iloc[0]["term"]),
                f"top_{cluster}_{direction.lower()}",
            )

    for term in metadata["term"].astype(str):
        add_term(term, "global_fill")
    output = metadata.set_index("term").loc[selected].reset_index()
    output.insert(0, "display_order", np.arange(1, len(output) + 1))
    output.insert(
        2,
        "selection_reason",
        [";".join(reasons[term]) for term in output["term"].astype(str)],
    )
    return output.drop(columns="_specific_cluster_order")


def select_gsea_heatmap_modules(
    gsea: pd.DataFrame,
    *,
    modules: Sequence[Mapping[str, Any]] = GSEA_HEATMAP_MODULES,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> pd.DataFrame:
    """Map every significant raw KEGG term once to a representative module row.

    The displayed NES and FDR always come from the representative term.  Terms
    with overlapping biological labels are documented but never averaged.
    """
    required = {"cluster", "term", "nes", "fdr", "direction"}
    missing = required - set(gsea.columns)
    if missing:
        raise KeyError(f"GSEA module selection is missing columns: {sorted(missing)}")

    work = gsea.copy()
    work["fdr"] = pd.to_numeric(work["fdr"], errors="coerce")
    significant_terms = set(
        work.loc[
            np.isfinite(work["fdr"].to_numpy(dtype=float))
            & work["fdr"].le(GSEA_FDR_THRESHOLD),
            "term",
        ].astype(str)
    )
    if not significant_terms:
        return pd.DataFrame(
            columns=[
                "display_order",
                "display_label",
                "term",
                "representative_term",
                "source_terms",
                "omitted_redundant_terms",
                "n_source_terms",
                "interpretation",
                "statistic_source",
            ]
        )

    mapped_terms: list[str] = []
    labels: list[str] = []
    representatives: list[str] = []
    for module in modules:
        label = str(module["display_label"])
        representative = str(module["representative_term"])
        source_terms = tuple(map(str, module["source_terms"]))
        if representative not in source_terms:
            raise RuntimeError(
                f"Representative term {representative!r} is absent from module {label!r}"
            )
        labels.append(label)
        representatives.append(representative)
        mapped_terms.extend(source_terms)
    duplicate_terms = sorted(
        term for term in set(mapped_terms) if mapped_terms.count(term) > 1
    )
    if duplicate_terms:
        raise RuntimeError(f"KEGG terms mapped to multiple modules: {duplicate_terms}")
    if len(labels) != len(set(labels)):
        raise RuntimeError("GSEA module display labels must be unique")
    if len(representatives) != len(set(representatives)):
        raise RuntimeError("GSEA module representative terms must be unique")

    mapped_set = set(mapped_terms)
    unmapped = sorted(significant_terms - mapped_set)
    mapped_but_not_significant = sorted(mapped_set - significant_terms)
    if unmapped or mapped_but_not_significant:
        raise RuntimeError(
            "Curated GSEA module map does not exactly account for significant "
            f"raw KEGG terms; unmapped={unmapped}, "
            f"mapped_but_not_significant={mapped_but_not_significant}"
        )

    raw_metadata = select_gsea_heatmap_pathways(
        gsea,
        max_pathways=len(significant_terms),
        clusters=clusters,
    )
    if set(raw_metadata["term"].astype(str)) != significant_terms:
        raise RuntimeError("Raw significant-pathway metadata is incomplete")
    metadata = raw_metadata.set_index("term")
    rows: list[dict[str, Any]] = []
    for display_order, module in enumerate(modules, start=1):
        representative = str(module["representative_term"])
        source_terms = tuple(map(str, module["source_terms"]))
        row = metadata.loc[representative].to_dict()
        row["representative_selection_reason"] = row.pop("selection_reason")
        row.update(
            {
                "display_order": display_order,
                "display_label": str(module["display_label"]),
                "term": representative,
                "representative_term": representative,
                "source_terms": ";".join(source_terms),
                "omitted_redundant_terms": ";".join(
                    term for term in source_terms if term != representative
                ),
                "n_source_terms": len(source_terms),
                "n_omitted_redundant_terms": len(source_terms) - 1,
                "interpretation": str(module["interpretation"]),
                "representative_rationale": str(
                    module.get(
                        "representative_rationale", "direct canonical KEGG term"
                    )
                ),
                "statistic_source": (
                    "representative_KEGG_term_NES_and_FDR; no_NES_averaging"
                ),
                "selection_reason": "curated_de_redundant_module",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def build_gsea_heatmap_matrices(
    gsea: pd.DataFrame,
    selection: pd.DataFrame,
    *,
    clusters: Sequence[str] = GSEA_CLUSTERS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ordered = selection.sort_values("display_order")
    terms = ordered["term"].astype(str).tolist()
    if not terms:
        index_name = "module" if "display_label" in selection.columns else "term"
        empty = pd.DataFrame(
            index=pd.Index([], name=index_name), columns=list(clusters)
        )
        return empty.copy(), empty.copy()
    if len(terms) != len(set(terms)):
        raise RuntimeError("GSEA heatmap representative terms must be unique")
    selected = gsea[gsea["term"].isin(terms)].copy()
    if selected.duplicated(["term", "cluster"]).any():
        raise RuntimeError("Duplicate term×cluster rows prevent GSEA heatmap pivot")
    nes = selected.pivot(index="term", columns="cluster", values="nes").reindex(
        index=terms,
        columns=list(clusters),
    )
    fdr = selected.pivot(index="term", columns="cluster", values="fdr").reindex(
        index=terms,
        columns=list(clusters),
    )
    if "display_label" in ordered.columns:
        labels = ordered["display_label"].astype(str).tolist()
        if len(labels) != len(set(labels)):
            raise RuntimeError("GSEA heatmap display labels must be unique")
        term_to_label = dict(zip(terms, labels, strict=True))
        nes = nes.rename(index=term_to_label)
        fdr = fdr.rename(index=term_to_label)
        nes.index.name = "module"
        fdr.index.name = "module"
    return nes, fdr


def render_gsea_heatmap_figure_e(
    gsea: pd.DataFrame,
    selection: pd.DataFrame,
    de_status: pd.DataFrame,
    path: Path,
) -> tuple[Path, Path]:
    nes, fdr = build_gsea_heatmap_matrices(gsea, selection)
    n_modules = int(len(nes))
    figure = plt.figure(
        figsize=(18.5, max(11.5, 0.68 * max(n_modules, 1) + 4.0)),
        layout="constrained",
    )
    grid = figure.add_gridspec(1, 2, width_ratios=[10.0, 4.5])
    axis = figure.add_subplot(grid[0, 0])
    annotation_axis = figure.add_subplot(grid[0, 1])
    if nes.empty:
        axis.text(
            0.5,
            0.5,
            "No KEGG pathways met GSEA FDR≤0.05",
            ha="center",
            va="center",
            transform=axis.transAxes,
        )
        axis.set_axis_off()
        annotation_axis.set_axis_off()
    else:
        values = nes.to_numpy(dtype=float)
        vmax = max(1.0, float(np.nanmax(np.abs(values))))
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#E5E7EB")
        image = axis.imshow(
            values,
            cmap=cmap,
            norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
            aspect="auto",
            interpolation="none",
        )
        pooled_status = (
            de_status[
                de_status["scope"].eq("pooled_cluster")
                & de_status["status"].eq("completed")
            ]
            .set_index("cluster")
            .reindex(GSEA_CLUSTERS)
        )
        labels: list[str] = []
        for cluster in GSEA_CLUSTERS:
            row = pooled_status.loc[cluster]
            dagger = "†" if cluster == "C7" else ""
            labels.append(
                f"{cluster}{dagger}\nT {int(row['n_tumor']):,} / C {int(row['n_control']):,}"
            )
        axis.set_xticks(np.arange(len(GSEA_CLUSTERS)), labels, fontsize=9)
        axis.set_yticks(
            np.arange(n_modules),
            nes.index.astype(str).tolist(),
            fontsize=9.2,
        )
        axis.set_xlabel(
            "Refined cluster (pooled tissues)\n"
            "Exploratory cell-level ranking (one sample per tissue×condition); "
            "NES is not a biological-replicate effect size. †C7 is small.\n"
            "Rows use one representative KEGG term (no NES averaging); "
            "AP-1‡ is a post-hoc Jun/Fos leading-edge proxy.",
            fontsize=8.4,
            color="#475569",
            labelpad=8,
        )
        axis.set_ylabel("De-redundant presentation group\n(representative KEGG NES/FDR)")
        axis.set_xticks(
            np.arange(-0.5, len(GSEA_CLUSTERS), 1),
            minor=True,
        )
        axis.set_yticks(np.arange(-0.5, n_modules, 1), minor=True)
        axis.grid(which="minor", color="white", linewidth=0.7)
        axis.tick_params(which="minor", bottom=False, left=False)
        for row_index in range(n_modules):
            for column_index in range(len(GSEA_CLUSTERS)):
                nes_value = values[row_index, column_index]
                if not np.isfinite(nes_value):
                    continue
                fdr_value = float(fdr.iloc[row_index, column_index])
                marker = "*" if fdr_value <= GSEA_FDR_THRESHOLD else ""
                text_color = "white" if abs(nes_value) >= 0.60 * vmax else "#111827"
                axis.text(
                    column_index,
                    row_index,
                    f"{nes_value:.2f}{marker}",
                    ha="center",
                    va="center",
                    fontsize=7.6,
                    color=text_color,
                    fontweight="bold" if marker else "normal",
                )
        colorbar = figure.colorbar(image, ax=axis, fraction=0.034, pad=0.02)
        colorbar.set_label("Signed normalized enrichment score (NES)")

        annotation_axis.set_xlim(0, 1)
        annotation_axis.set_ylim(axis.get_ylim())
        annotation_axis.set_xticks([])
        annotation_axis.set_yticks([])
        for spine in annotation_axis.spines.values():
            spine.set_visible(False)
        annotation_axis.set_title(
            "Representative KEGG term\nFDR≤0.05 pattern",
            fontsize=10,
            fontweight="bold",
        )
        metadata_key = (
            "display_label" if "display_label" in selection.columns else "term"
        )
        selected_meta = selection.set_index(metadata_key).reindex(nes.index)
        for row_index, row in enumerate(selected_meta.itertuples()):
            source_term = str(
                getattr(row, "representative_term", getattr(row, "Index"))
            )
            omitted_count = int(getattr(row, "n_omitted_redundant_terms", 0))
            if int(row.n_significant_clusters) >= 2:
                pattern = f"shared ({int(row.n_significant_clusters)} clusters)"
            else:
                cluster = str(row.significant_clusters)
                direction = (
                    "Tumor"
                    if str(row.tumor_significant_clusters)
                    else "Control"
                )
                pattern = f"cluster-specific · {cluster} {direction}"
            label_lines = [
                pattern,
                f"source KEGG: {source_term}",
            ]
            if omitted_count:
                label_lines.append(f"redundant KEGG labels omitted: {omitted_count}")
            label_lines.extend(
                [
                    f"T: {row.tumor_significant_clusters or '—'}",
                    f"C: {row.control_significant_clusters or '—'}",
                ]
            )
            annotation_axis.text(
                0.0,
                row_index,
                "\n".join(label_lines),
                ha="left",
                va="center",
                fontsize=8.0,
                color="#334155",
            )
    figure.suptitle(
        f"{RUN_DATE} method-adapted Figure 3E — all-cluster preranked KEGG GSEA\n"
        "full-gene Wilcoxon rank; signed NES: red/positive = Tumor, "
        "blue/negative = Control; * GSEA FDR≤0.05; "
        "10 de-redundant presentation groups",
        fontsize=14,
        fontweight="bold",
    )
    return save_figure(figure, path)


VENN_TISSUES = ("bone_marrow", "spleen", "thymus")


def tumor_up_mask(frame: pd.DataFrame, threshold: str) -> pd.Series:
    if threshold == PAPER_VENN_THRESHOLD:
        return frame["paper_venn_tumor_up"].astype(bool)
    if threshold == ROBUST_THRESHOLD:
        return frame["robust_tumor_up"].astype(bool)
    raise ValueError(f"Unknown overlap threshold: {threshold}")


def exact_intersection_regions(
    sets: Mapping[str, set[str]],
    set_order: Sequence[str],
) -> dict[tuple[str, ...], set[str]]:
    missing = [set_id for set_id in set_order if set_id not in sets]
    if missing:
        raise KeyError(f"Missing overlap sets: {missing}")
    regions: dict[tuple[str, ...], set[str]] = {}
    for bitmask in range(1, 1 << len(set_order)):
        active = tuple(
            set_id
            for index, set_id in enumerate(set_order)
            if bitmask & (1 << index)
        )
        inactive = tuple(set_id for set_id in set_order if set_id not in active)
        genes = set.intersection(*(sets[set_id] for set_id in active)).copy()
        for set_id in inactive:
            genes.difference_update(sets[set_id])
        regions[active] = genes
    return regions


def build_figure_f_main_membership(
    de: pd.DataFrame,
    threshold: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    sets: dict[str, set[str]] = {}
    for set_id in FIGURE_F_MAIN_SET_ORDER:
        cluster, tissue = set_id.split("__", maxsplit=1)
        frame = de[
            de["unit_id"].eq(set_id)
            & de["cluster"].eq(cluster)
            & de["tissue"].eq(tissue)
        ]
        if frame.empty:
            raise RuntimeError(f"Missing DE rows for main Figure 3F set {set_id}")
        selected = tumor_up_mask(frame, threshold)
        sets[set_id] = set(frame.loc[selected, "gene"].dropna().astype(str))

    universe = sorted(set().union(*(sets[set_id] for set_id in FIGURE_F_MAIN_SET_ORDER)))
    membership = pd.DataFrame(
        {
            "gene": universe,
            **{
                set_id: [gene in sets[set_id] for gene in universe]
                for set_id in FIGURE_F_MAIN_SET_ORDER
            },
        }
    )
    if membership.empty:
        membership["n_sets"] = pd.Series(dtype=int)
        membership["intersection"] = pd.Series(dtype=str)
    else:
        membership["n_sets"] = membership[list(FIGURE_F_MAIN_SET_ORDER)].sum(axis=1)
        membership["intersection"] = [
            "&".join(
                set_id
                for set_id in FIGURE_F_MAIN_SET_ORDER
                if bool(row[set_id])
            )
            for row in membership.to_dict(orient="records")
        ]

    regions = exact_intersection_regions(sets, FIGURE_F_MAIN_SET_ORDER)
    summary_rows = [
        {
            "threshold": threshold,
            "metric": f"set_size__{set_id}",
            "n_genes": len(sets[set_id]),
            "n_sets": 1,
            "genes": ";".join(sorted(sets[set_id])),
        }
        for set_id in FIGURE_F_MAIN_SET_ORDER
    ]
    summary_rows.extend(
        {
            "threshold": threshold,
            "metric": "intersection_exact__" + "&".join(active),
            "n_genes": len(genes),
            "n_sets": len(active),
            "genes": ";".join(sorted(genes)),
        }
        for active, genes in regions.items()
    )
    all_four = set.intersection(*(sets[set_id] for set_id in FIGURE_F_MAIN_SET_ORDER))
    summary_rows.extend(
        [
            {
                "threshold": threshold,
                "metric": "all_four",
                "n_genes": len(all_four),
                "n_sets": len(FIGURE_F_MAIN_SET_ORDER),
                "genes": ";".join(sorted(all_four)),
            },
            {
                "threshold": threshold,
                "metric": "union",
                "n_genes": len(universe),
                "n_sets": 0,
                "genes": ";".join(universe),
            },
        ]
    )
    return membership, pd.DataFrame(summary_rows), sets


def render_figure_f_main_upset(
    sets: Mapping[str, set[str]],
    threshold: str,
    path: Path,
    de_status: pd.DataFrame,
) -> tuple[Path, Path]:
    regions = exact_intersection_regions(sets, FIGURE_F_MAIN_SET_ORDER)
    all_four_pattern = tuple(FIGURE_F_MAIN_SET_ORDER)
    patterns = [all_four_pattern]
    patterns.extend(
        sorted(
            (
                pattern
                for pattern, genes in regions.items()
                if genes and pattern != all_four_pattern
            ),
            key=lambda pattern: (
                -len(regions[pattern]),
                -len(pattern),
                tuple(FIGURE_F_MAIN_SET_ORDER.index(item) for item in pattern),
            ),
        )
    )
    counts = np.asarray([len(regions[pattern]) for pattern in patterns], dtype=int)
    x_positions = np.arange(len(patterns), dtype=float)

    status = de_status.set_index("unit_id").reindex(FIGURE_F_MAIN_SET_ORDER)
    if status["status"].ne("completed").any():
        raise RuntimeError("Main Figure 3F contains incomplete DE units")
    strict_estimable = bool(
        status["passes_recommended_min20"].fillna(False).astype(bool).all()
    )
    sample_text = " | ".join(
        f"{FIGURE_F_MAIN_SET_LABELS[set_id]} T={int(status.loc[set_id, 'n_tumor'])}"
        f"/C={int(status.loc[set_id, 'n_control'])}"
        for set_id in FIGURE_F_MAIN_SET_ORDER
    )
    threshold_label = (
        "paper Figure 3F rule: log2FC≥1 and nominal P<0.05"
        if threshold == PAPER_VENN_THRESHOLD
        else "sensitivity rule: log2FC≥0.25 and BH-FDR≤0.05"
    )

    figure = plt.figure(figsize=(16.0, 9.6), layout="constrained")
    grid = figure.add_gridspec(2, 1, height_ratios=(2.6, 1.25), hspace=0.04)
    bar_axis = figure.add_subplot(grid[0])
    matrix_axis = figure.add_subplot(grid[1], sharex=bar_axis)
    colors = ["#7C3AED"] + ["#475569"] * (len(patterns) - 1)
    bars = bar_axis.bar(x_positions, counts, color=colors, width=0.72)
    bar_axis.bar_label(
        bars,
        labels=[str(value) for value in counts],
        padding=3,
        fontsize=9,
    )
    bar_axis.set_ylabel("Genes in exact intersection")
    bar_axis.grid(axis="y", color="#E2E8F0", linewidth=0.7)
    bar_axis.spines[["top", "right"]].set_visible(False)
    bar_axis.tick_params(axis="x", bottom=False, labelbottom=False)

    shared_genes = sorted(regions[all_four_pattern])
    shared_text = ", ".join(shared_genes) if shared_genes else "none"
    wrapped_shared = "\n".join(
        shared_text[index : index + 92]
        for index in range(0, len(shared_text), 92)
    )
    bar_axis.text(
        0.99,
        0.98,
        f"All-four shared genes (n={len(shared_genes)})\n{wrapped_shared}",
        transform=bar_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "#F5F3FF",
            "edgecolor": "#8B5CF6",
        },
    )

    y_positions = np.arange(len(FIGURE_F_MAIN_SET_ORDER))[::-1]
    for column, pattern in enumerate(patterns):
        matrix_axis.scatter(
            np.full(len(y_positions), column),
            y_positions,
            s=38,
            color="#CBD5E1",
            zorder=1,
        )
        active_y = [
            y_positions[index]
            for index, set_id in enumerate(FIGURE_F_MAIN_SET_ORDER)
            if set_id in pattern
        ]
        if len(active_y) > 1:
            matrix_axis.plot(
                [column, column],
                [min(active_y), max(active_y)],
                color="#111827",
                linewidth=1.7,
                zorder=2,
            )
        matrix_axis.scatter(
            np.full(len(active_y), column),
            active_y,
            s=58,
            color="#111827" if pattern != all_four_pattern else "#7C3AED",
            zorder=3,
        )
    matrix_axis.set_yticks(
        y_positions,
        [
            f"{FIGURE_F_MAIN_SET_LABELS[set_id]}  (set n={len(sets[set_id])})"
            for set_id in FIGURE_F_MAIN_SET_ORDER
        ],
    )
    matrix_axis.set_xticks(
        x_positions,
        [str(index + 1) for index in range(len(patterns))],
    )
    matrix_axis.set_xlabel("Exact intersection pattern (all-four shown first)")
    matrix_axis.set_ylim(-0.7, len(y_positions) - 0.3)
    matrix_axis.spines[["top", "right", "left"]].set_visible(False)
    matrix_axis.tick_params(axis="y", length=0)
    matrix_axis.grid(axis="x", color="#F1F5F9", linewidth=0.5)

    figure.suptitle(
        f"{RUN_DATE} paper-adapted Figure 3F — C5 state × tissue Tumor-up overlap\n"
        f"C5-1/C5-2 × bone marrow/spleen; {threshold_label}\n"
        f"{sample_text}\n"
        f"strict min20 all-four intersection: "
        f"{'estimable' if strict_estimable else 'NA (not estimable)'}",
        fontsize=13,
        fontweight="bold",
    )
    return save_figure(figure, path)


def build_venn_membership(
    de: pd.DataFrame,
    threshold: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    tissue_de = de[de["scope"].eq("C5-2_tissue")].copy()
    sets: dict[str, set[str]] = {}
    for tissue in VENN_TISSUES:
        frame = tissue_de[tissue_de["tissue"].eq(tissue)]
        selected = tumor_up_mask(frame, threshold)
        sets[tissue] = set(frame.loc[selected, "gene"].dropna().astype(str))
    universe = sorted(set().union(*sets.values())) if sets else []
    membership = pd.DataFrame(
        {
            "gene": universe,
            **{
                tissue: [gene in sets[tissue] for gene in universe]
                for tissue in VENN_TISSUES
            },
        }
    )
    if membership.empty:
        membership["region"] = pd.Series(dtype=str)
    else:
        membership["region"] = [
            "&".join(tissue for tissue in VENN_TISSUES if bool(row[tissue]))
            for row in membership.to_dict(orient="records")
        ]
    a, b, c = (sets[tissue] for tissue in VENN_TISSUES)
    regions = {
        "bone_marrow_only": a - b - c,
        "spleen_only": b - a - c,
        "thymus_only": c - a - b,
        "bone_marrow&spleen_only": (a & b) - c,
        "bone_marrow&thymus_only": (a & c) - b,
        "spleen&thymus_only": (b & c) - a,
        "all_three": a & b & c,
    }
    summary_rows = [
        {
            "threshold": threshold,
            "metric": f"set_size__{tissue}",
            "n_genes": len(sets[tissue]),
            "genes": ";".join(sorted(sets[tissue])),
        }
        for tissue in VENN_TISSUES
    ]
    summary_rows.extend(
        {
            "threshold": threshold,
            "metric": region,
            "n_genes": len(genes),
            "genes": ";".join(sorted(genes)),
        }
        for region, genes in regions.items()
    )
    return membership, pd.DataFrame(summary_rows), sets


def venn_region_counts(sets: Mapping[str, set[str]]) -> dict[str, int]:
    a, b, c = (sets[tissue] for tissue in VENN_TISSUES)
    return {
        "a_only": len(a - b - c),
        "b_only": len(b - a - c),
        "c_only": len(c - a - b),
        "ab_only": len((a & b) - c),
        "ac_only": len((a & c) - b),
        "bc_only": len((b & c) - a),
        "abc": len(a & b & c),
    }


def render_venn_figure_f(
    sets: Mapping[str, set[str]],
    threshold: str,
    path: Path,
    de_status: pd.DataFrame,
) -> tuple[Path, Path]:
    counts = venn_region_counts(sets)
    figure, axis = plt.subplots(figsize=(11.2, 8.7), layout="constrained")
    circle_specs = [
        ((-0.62, 0.24), "#F59E0B", TISSUE_LABELS["bone_marrow"]),
        ((0.62, 0.24), "#DC2626", TISSUE_LABELS["spleen"]),
        ((0.0, -0.56), "#2563EB", TISSUE_LABELS["thymus"]),
    ]
    for center, color, _ in circle_specs:
        axis.add_patch(
            Circle(
                center,
                radius=1.0,
                facecolor=color,
                edgecolor=color,
                linewidth=2.0,
                alpha=0.22,
            )
        )
    text_positions = {
        "a_only": (-1.08, 0.48),
        "b_only": (1.08, 0.48),
        "c_only": (0.0, -1.18),
        "ab_only": (0.0, 0.76),
        "ac_only": (-0.48, -0.40),
        "bc_only": (0.48, -0.40),
        "abc": (0.0, -0.02),
    }
    for region, position in text_positions.items():
        axis.text(
            position[0],
            position[1],
            str(counts[region]),
            ha="center",
            va="center",
            fontsize=15 if region == "abc" else 12,
            fontweight="bold" if region == "abc" else "normal",
            color="#111827",
        )
    axis.text(-1.18, 1.19, f"Bone marrow\nn={len(sets['bone_marrow']):,}", ha="center", fontweight="bold")
    axis.text(1.18, 1.19, f"Spleen\nn={len(sets['spleen']):,}", ha="center", fontweight="bold")
    axis.text(0.0, -1.73, f"Thymus\nn={len(sets['thymus']):,}", ha="center", fontweight="bold")
    shared = sorted(set.intersection(*(sets[tissue] for tissue in VENN_TISSUES)))
    shared_text = ", ".join(shared) if shared else "none"
    axis.text(
        2.0,
        0.15,
        "All-three shared genes\n" + "\n".join(
            [shared_text[index : index + 54] for index in range(0, len(shared_text), 54)]
        ),
        ha="left",
        va="center",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "#F8FAFC",
            "edgecolor": "#94A3B8",
        },
    )
    status = (
        de_status[de_status["scope"].eq("C5-2_tissue")]
        .set_index("tissue")
        .reindex(VENN_TISSUES)
    )
    strict_three_tissue_estimable = bool(
        status["passes_recommended_min20"].fillna(False).astype(bool).all()
    )
    estimability_text = (
        "strict min20 all-three intersection: estimable"
        if strict_three_tissue_estimable
        else "strict min20 all-three intersection: NA (not estimable)"
    )
    sample_text = " | ".join(
        f"{TISSUE_LABELS[tissue]} T={int(status.loc[tissue, 'n_tumor'])}/C={int(status.loc[tissue, 'n_control'])}"
        f"{' [<20; exploratory]' if not bool(status.loc[tissue, 'passes_recommended_min20']) else ''}"
        for tissue in VENN_TISSUES
        if tissue in status.index and status.loc[tissue, "status"] == "completed"
    )
    if threshold == PAPER_VENN_THRESHOLD:
        threshold_label = "paper Figure 3F rule: log2FC≥1 and nominal P<0.05"
    else:
        threshold_label = "sensitivity rule: log2FC≥0.25 and BH-FDR≤0.05"
    axis.set_title(
        f"{RUN_DATE} supplemental Figure 3F — C5-2 Tumor-up genes by tissue\n"
        f"{threshold_label}\n{sample_text}\n{estimability_text}",
        fontsize=13,
        fontweight="bold",
    )
    axis.set_xlim(-2.0, 3.25)
    axis.set_ylim(-1.95, 1.55)
    axis.set_aspect("equal")
    axis.axis("off")
    return save_figure(figure, path)


def feature_rows(features_path: Path, symbol: str) -> list[tuple[int, str, str]]:
    matches: list[tuple[int, str, str]] = []
    with gzip.open(features_path, "rt") as handle:
        for row_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[1] == symbol and fields[2] == "Gene Expression":
                matches.append((row_number, fields[0], fields[2]))
    return matches


def audit_raw_ncam1(raw_root: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for sample in SAMPLE_DISPLAY_ORDER:
        matrix_dir = raw_root / sample / "sample_feature_bc_matrix"
        features_path = matrix_dir / "features.tsv.gz"
        matrix_path = matrix_dir / "matrix.mtx.gz"
        for required in (features_path, matrix_path):
            if not required.exists():
                raise FileNotFoundError(required)
        matches = feature_rows(features_path, "Ncam1")
        target_rows = [row_number for row_number, _, _ in matches]
        entries: list[tuple[int, int, float]] = []
        if target_rows:
            row_pattern = "^(" + "|".join(map(str, target_rows)) + ")[[:space:]]"
            completed = subprocess.run(
                ["zgrep", "-E", row_pattern, str(matrix_path)],
                check=False,
                capture_output=True,
                text=True,
            )
            if completed.returncode not in (0, 1):
                raise RuntimeError(
                    f"zgrep failed for {sample}/Ncam1: {completed.stderr.strip()}"
                )
            for line in completed.stdout.splitlines():
                row_number, column_number, count = line.split()
                entries.append((int(row_number), int(column_number), float(count)))
        rows.append(
            {
                "sample": sample,
                "human_marker": "CD56",
                "mouse_symbol": "Ncam1",
                "feature_present": bool(matches),
                "feature_rows_1_based": ";".join(str(row) for row, _, _ in matches),
                "feature_ids": ";".join(feature_id for _, feature_id, _ in matches),
                "raw_matrix_nonzero_entries": len(entries),
                "raw_matrix_total_counts": float(sum(value for _, _, value in entries)),
                "features_path": str(features_path.resolve()),
                "matrix_path": str(matrix_path.resolve()),
                "features_size_bytes": features_path.stat().st_size,
                "matrix_size_bytes": matrix_path.stat().st_size,
            }
        )
    result = pd.DataFrame(rows)
    present_with_counts = result["feature_present"].astype(bool) & (
        result["raw_matrix_nonzero_entries"].gt(0)
        | result["raw_matrix_total_counts"].gt(0)
    )
    result["audit_status"] = np.select(
        [~result["feature_present"].astype(bool), present_with_counts],
        ["feature_absent", "present_with_counts"],
        default="present_no_counts",
    )
    return result


def validate_h5ad_roundtrip(
    full_path: Path,
    c5_path: Path,
    expected_full: ad.AnnData,
    expected_c5: ad.AnnData,
) -> dict[str, Any]:
    result: dict[str, Any] = {"status": "passed", "files": {}}
    specifications = (
        (
            "full",
            full_path,
            expected_full,
            ("X_pca", "X_umap", f"X_umap_c5_only_{RUN_DATE}"),
            (ORIGINAL_CLUSTER_KEY, REFINED_CLUSTER_KEY),
        ),
        (
            "c5",
            c5_path,
            expected_c5,
            ("X_pca", "X_umap", "X_umap_original_full"),
            (ORIGINAL_CLUSTER_KEY, C5_SUBCLUSTER_KEY),
        ),
    )
    for label, path, expected, required_obsm, required_obs in specifications:
        loaded = ad.read_h5ad(path, backed="r")
        try:
            checks = {
                "n_obs": loaded.n_obs == expected.n_obs,
                "n_vars": loaded.n_vars == expected.n_vars,
                "obs_names_order": loaded.obs_names.equals(expected.obs_names),
                "raw_present": loaded.raw is not None,
                "raw_n_vars": loaded.raw is not None
                and loaded.raw.n_vars == expected.raw.n_vars,
                "counts_layer": "counts" in loaded.layers
                and loaded.layers["counts"].shape == expected.shape,
                "required_obsm": set(required_obsm).issubset(loaded.obsm.keys()),
                "required_obs": set(required_obs).issubset(loaded.obs.columns),
            }
            if checks["required_obs"]:
                for key in required_obs:
                    checks[f"obs_values__{key}"] = np.array_equal(
                        loaded.obs[key].astype(str).to_numpy(),
                        expected.obs[key].astype(str).to_numpy(),
                    )
            failed = sorted(name for name, passed in checks.items() if not bool(passed))
            if failed:
                raise RuntimeError(f"{label} H5AD round-trip checks failed: {failed}")
            result["files"][label] = {
                "path": str(path),
                "size_bytes": int(path.stat().st_size),
                "checks": checks,
            }
        finally:
            loaded.file.close()
    return result


def build_overview_pdf(image_paths: Sequence[Path], output_path: Path) -> None:
    pages: list[Image.Image] = []
    for path in image_paths:
        if not path.exists():
            continue
        with Image.open(path) as source:
            image = source.convert("RGB")
            if image.width > 1800:
                height = int(round(image.height * 1800 / image.width))
                image = image.resize((1800, height), Image.Resampling.LANCZOS)
            pages.append(image.copy())
    if not pages:
        raise RuntimeError("No figure pages were available for the overview PDF")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    first, remaining = pages[0], pages[1:]
    first.save(output_path, "PDF", resolution=160.0, save_all=True, append_images=remaining)
    for page in pages:
        page.close()


def top_gene_summary(de: pd.DataFrame, *, n: int = 8) -> dict[str, dict[str, list[str]]]:
    output: dict[str, dict[str, list[str]]] = {}
    pooled = de[de["scope"].eq("pooled_cluster")]
    for cluster in TARGET_CLUSTERS:
        frame = pooled[pooled["cluster"].eq(cluster)]
        tumor = (
            frame[frame["robust_tumor_up"].astype(bool)]
            .sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False])
            .head(n)["gene"]
            .astype(str)
            .tolist()
        )
        control = (
            frame[frame["robust_control_up"].astype(bool)]
            .sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, True])
            .head(n)["gene"]
            .astype(str)
            .tolist()
        )
        output[cluster] = {"tumor_up": tumor, "control_up": control}
    return output


def write_readme(
    path: Path,
    *,
    c5: ad.AnnData,
    metrics: pd.DataFrame,
    selected_resolution: float,
    complementarity: pd.DataFrame,
    de_status: pd.DataFrame,
    de: pd.DataFrame,
    ora: pd.DataFrame,
    gsea: pd.DataFrame,
    gsea_selection: pd.DataFrame,
    gsea_audit: Mapping[str, Any],
    main_f_summaries: Mapping[str, pd.DataFrame],
    supplemental_venn_summaries: Mapping[str, pd.DataFrame],
    ncam1_audit: pd.DataFrame,
    original_delta: pd.DataFrame,
    kegg_audit: Mapping[str, Any],
    h5ad_written: bool,
) -> None:
    selected = metrics.loc[np.isclose(metrics["resolution"], selected_resolution)].iloc[0]
    c5_counts = c5.obs[C5_SUBCLUSTER_KEY].astype(str).value_counts()
    original_coordinates = np.asarray(c5.obsm["X_umap_original_full"])
    labels = c5.obs[C5_SUBCLUSTER_KEY].astype(str).to_numpy()
    medians = {
        label: float(np.median(original_coordinates[labels == label, 0]))
        for label in ("C5-1", "C5-2")
    }
    condition_rows = (
        complementarity[complementarity["scope"].eq("condition")]
        .set_index("condition")
        .reindex(CONDITION_ORDER)
    )
    bm_marker_rows = (
        complementarity[
            complementarity["scope"].eq("condition_tissue")
            & complementarity["tissue"].eq("bone_marrow")
        ]
        .set_index("condition")
        .reindex(CONDITION_ORDER)
    )
    main_f_status = (
        de_status.set_index("unit_id").reindex(FIGURE_F_MAIN_SET_ORDER)
    )
    strict_main_f_estimable = bool(
        main_f_status["passes_recommended_min20"]
        .fillna(False)
        .astype(bool)
        .all()
    )
    c5_tissue_status = (
        de_status[de_status["scope"].eq("C5-2_tissue")]
        .set_index("tissue")
        .reindex(VENN_TISSUES)
    )
    strict_three_tissue_estimable = bool(
        c5_tissue_status["passes_recommended_min20"].fillna(False).astype(bool).all()
    )
    main_f_sample_summary = " | ".join(
        f"{FIGURE_F_MAIN_SET_LABELS[set_id]} "
        f"T={int(main_f_status.loc[set_id, 'n_tumor'])}"
        f"/C={int(main_f_status.loc[set_id, 'n_control'])}"
        f"{' [<20; exploratory]' if not bool(main_f_status.loc[set_id, 'passes_recommended_min20']) else ''}"
        for set_id in FIGURE_F_MAIN_SET_ORDER
    )
    tissue_sample_summary = " | ".join(
        f"{TISSUE_LABELS[tissue]} T={int(c5_tissue_status.loc[tissue, 'n_tumor'])}"
        f"/C={int(c5_tissue_status.loc[tissue, 'n_control'])}"
        f"{' [<20; exploratory]' if not bool(c5_tissue_status.loc[tissue, 'passes_recommended_min20']) else ''}"
        for tissue in VENN_TISSUES
    )
    main_f_estimability_suffix = (
        ""
        if strict_main_f_estimable
        else " [strict min20 all-four = NA, not 0]"
    )
    supplemental_estimability_suffix = (
        ""
        if strict_three_tissue_estimable
        else " [exploratory; strict min20 all-three = NA, not 0]"
    )
    feature_present_samples = int(ncam1_audit["feature_present"].sum())
    raw_sample_count = int(len(ncam1_audit))
    ncam1_nonzero_entries = int(ncam1_audit["raw_matrix_nonzero_entries"].sum())
    ncam1_total_counts = float(ncam1_audit["raw_matrix_total_counts"].sum())
    if ncam1_nonzero_entries > 0 or ncam1_total_counts > 0:
        ncam1_conclusion = "Raw Ncam1 reads are present; use the audit table for sample-level counts."
    elif feature_present_samples > 0:
        ncam1_conclusion = "Ncam1 is listed as a feature but has zero raw reads; no CD56/Ncam1 distribution can be analyzed."
    else:
        ncam1_conclusion = "Ncam1 is absent from the raw feature tables; no CD56/Ncam1 distribution can be analyzed."
    frozen_gmt_sha256 = str(
        kegg_audit.get("frozen_copy_sha256", kegg_audit["sha256"])
    )

    top_genes = top_gene_summary(de)
    gsea_parameters = gsea_audit["parameters"]
    shared_selected = int(
        gsea_selection["classification"].eq("shared").sum()
    )
    specific_selected = int(
        gsea_selection["classification"].astype(str).str.startswith(
            "cluster-specific:"
        ).sum()
    )
    n_source_terms = int(
        pd.to_numeric(gsea_selection["n_source_terms"]).sum()
    )
    n_omitted_terms = int(
        pd.to_numeric(gsea_selection["n_omitted_redundant_terms"]).sum()
    )
    reuse_status = str(gsea_audit.get("reuse_validation", {}).get("status", "unknown"))
    reuse_receipt_status = str(gsea_audit.get("reuse_validation", {}).get("receipt_status", "fresh"))
    gsea_pathway_lines: list[str] = []
    for cluster in GSEA_CLUSTERS:
        frame = gsea[
            gsea["cluster"].eq(cluster)
            & gsea["significant_fdr_0_05"].astype(bool)
        ].copy()
        tumor = frame[frame["nes"].gt(0)].sort_values(
            ["fdr", "nes", "term"],
            ascending=[True, False, True],
            kind="stable",
        )
        control = frame[frame["nes"].lt(0)].sort_values(
            ["fdr", "nes", "term"],
            ascending=[True, True, True],
            kind="stable",
        )
        tumor_terms = tumor.head(3)["term"].astype(str).tolist()
        control_terms = control.head(3)["term"].astype(str).tolist()
        gsea_pathway_lines.append(
            f"- {cluster}: Tumor-enriched {len(tumor)}, "
            f"top {', '.join(tumor_terms) if tumor_terms else 'none'}; "
            f"Control-enriched {len(control)}, "
            f"top {', '.join(control_terms) if control_terms else 'none'}."
        )
    de_lines: list[str] = []
    for cluster in TARGET_CLUSTERS:
        status = de_status[
            de_status["scope"].eq("pooled_cluster")
            & de_status["cluster"].eq(cluster)
        ].iloc[0]
        genes = top_genes[cluster]
        de_lines.append(
            f"- {cluster}: Tumor/Control = {int(status.n_tumor)}/{int(status.n_control)}; "
            f"robust Tumor-up {int(status.n_robust_tumor_up)}, "
            f"Control-up {int(status.n_robust_control_up)}; "
            f"top Tumor-up: {', '.join(genes['tumor_up']) if genes['tumor_up'] else 'none'}."
        )
    main_f_lines: list[str] = []
    for threshold, frame in main_f_summaries.items():
        indexed = frame.set_index("metric")
        shared = str(indexed.loc["all_four", "genes"])
        shared = shared if shared else "none"
        set_sizes = "/".join(
            str(int(indexed.loc[f"set_size__{set_id}", "n_genes"]))
            for set_id in FIGURE_F_MAIN_SET_ORDER
        )
        main_f_lines.append(
            f"- {threshold}: C5-1_BM/C5-1_Spleen/C5-2_BM/C5-2_Spleen "
            f"set sizes {set_sizes}; all-four "
            f"{int(indexed.loc['all_four', 'n_genes'])}: {shared}."
            f"{main_f_estimability_suffix}"
        )
    supplemental_venn_lines: list[str] = []
    for threshold, frame in supplemental_venn_summaries.items():
        indexed = frame.set_index("metric")
        shared = str(indexed.loc["all_three", "genes"])
        shared = shared if shared else "none"
        supplemental_venn_lines.append(
            f"- {threshold}: BM/Spleen/Thymus set sizes "
            f"{int(indexed.loc['set_size__bone_marrow', 'n_genes'])}/"
            f"{int(indexed.loc['set_size__spleen', 'n_genes'])}/"
            f"{int(indexed.loc['set_size__thymus', 'n_genes'])}; "
            f"all-three {int(indexed.loc['all_three', 'n_genes'])}: {shared}.{supplemental_estimability_suffix}"
        )
    baseline = (
        original_delta[["tissue", "tumor_share_zero_delta_baseline_pct"]]
        .drop_duplicates()
        .set_index("tissue")
    )
    bm_c0 = original_delta[
        original_delta["tissue"].eq("bone_marrow")
        & original_delta["cluster"].eq("C0")
    ].iloc[0]
    h5ad_line = (
        f"- {RUN_DATE}_inkt_selected_umap_C5_split.h5ad retains original clusters "
        f"and adds {REFINED_CLUSTER_KEY}; dated round-trip validation passed."
        if h5ad_written
        else "- Full H5AD was skipped; the CSV cell mapping remains complete."
    )
    text = f"""# {RUN_DATE} C5 split + paper Figure 3D–F follow-up

## 1. C5 split

- Original C5: {c5.n_obs:,} cells; C5-1 = {int(c5_counts.get('C5-1', 0)):,}, C5-2 = {int(c5_counts.get('C5-2', 0)):,}.
- Within-C5: {int(selected.n_hvg):,} HVGs, {int(selected.n_pcs)} PCs, {int(selected.n_neighbors)}-NN, selected Leiden resolution {selected_resolution:g}.
- Five-seed mean pairwise ARI = {selected.mean_pairwise_ari:.4f}; minimum ARI = {selected.min_pairwise_ari:.4f}.
- Selected cross-seed medoid seed = {int(selected.selected_partition_seed)}; its mean ARI to other seeds = {selected.selected_seed_mean_ari:.4f}.
- Original-UMAP geometry: KMeans agreement ARI = {selected.selected_original_umap_kmeans_ari:.4f}; silhouette = {selected.selected_original_umap_silhouette:.4f}.
- Original-UMAP1 median: C5-1 {medians['C5-1']:.3f}, C5-2 {medians['C5-2']:.3f}. C5-2 is the right-hand island.
- Median UMAP1 gap = {selected.selected_original_umap_median_x_gap:.3f}; q05(C5-2) − q95(C5-1) = {selected.selected_original_umap_q05_q95_gap:.3f}. A positive quantile gap directly supports the visible separation.
- Input feature-universe validation passed: X/raw have the same {c5.n_vars:,} genes in the same order; Il4 and Klrd1 are present.
- Original {ORIGINAL_CLUSTER_KEY} is retained. Refined labels use {REFINED_CLUSTER_KEY}.

## 2. IL-4 / CD94 condition masks

- Control: detection phi {condition_rows.loc['Ctrl', 'detection_phi']:.4f}; expression Spearman {condition_rows.loc['Ctrl', 'spearman_log1p_expression']:.4f}; double-positive {condition_rows.loc['Ctrl', 'both_positive_pct']:.2f}%.
- Tumor: detection phi {condition_rows.loc['T2', 'detection_phi']:.4f}; expression Spearman {condition_rows.loc['T2', 'spearman_log1p_expression']:.4f}; double-positive {condition_rows.loc['T2', 'both_positive_pct']:.2f}%.
- Bone marrow detection phi: Control {bm_marker_rows.loc['Ctrl', 'detection_phi']:.4f}; Tumor {bm_marker_rows.loc['T2', 'detection_phi']:.4f}.
- Negative phi/Spearman values support a weak complementary pattern, strongest in bone marrow; the non-zero double-positive fraction shows it is not strict mutual exclusivity.
- Each mask panel highlights one condition and fades the other. The table also contains condition×tissue and condition×cluster Fisher/phi/Spearman statistics.
- “Spatial” here means the transcriptomic UMAP landscape, not physical spatial-transcriptomics coordinates.

## 3. Paper-adapted Figure 3D

{chr(10).join(de_lines)}

The x/y axes follow the paper: log2FC and −log10 nominal P. Color uses the robust project rule BH-FDR≤0.05 and |log2FC|≥0.25. Tissues are pooled, so tissue-composition effects can contribute.

## 4. Method-adapted Figure 3E — all-cluster preranked GSEA

- The paper supplies only the display/method framework; all biological results below come from this iNKT dataset.
- Every refined cluster ({', '.join(GSEA_CLUSTERS)}) uses the full {int(gsea['n_ranked_genes'].min()) if not gsea.empty else 0:,}-gene Wilcoxon-score ranking for Tumor T2 vs Control Ctrl.
- Signed NES convention: positive = Tumor enriched; negative = Control enriched. An asterisk in the heatmap marks the native GSEA FDR q≤{GSEA_FDR_THRESHOLD:g}.
- Frozen gene-set library: Enrichr {kegg_audit['library']}, GMT SHA256 {frozen_gmt_sha256[:16]}….
- Parameters: gene-set permutations={int(gsea_parameters['permutation_num']):,}, weight={float(gsea_parameters['weight']):g}, size {int(gsea_parameters['min_size'])}–{int(gsea_parameters['max_size'])}, seed={int(gsea_parameters['seed'])}.
- Heatmap display: {len(gsea_selection)} manually de-redundant presentation groups ({shared_selected} shared, {specific_selected} cluster-specific representative terms), accounting exactly once for all {n_source_terms} FDR-significant raw KEGG terms; {n_omitted_terms} redundant labels are omitted from the main y-axis.
- These are presentation groups, not independently tested module-level gene sets. Each plotted NES/FDR is from the explicitly listed representative raw KEGG term; no NES values are averaged. The module-selection table records every source term, representative rationale, and omitted label; the full unfiltered cluster×pathway table is retained.
- “Lysine degradation” keeps its original KEGG identity and is only annotated as KMT/SETD-enriched at the leading edge. “AP-1 immediate-early proxy” is a post-hoc Jun/Fos leading-edge interpretation of the addiction pathways, not a formal AP-1 gene-set test.
- GSEA materialization status: {reuse_status}; receipt: {reuse_receipt_status}. Reuse is accepted only after current DE ranking, GMT SHA256, parameter, eligible-term, direction/range, leading-edge, and bundle-hash validation.

{chr(10).join(gsea_pathway_lines)}

The former thresholded Tumor-up KEGG ORA is retained as a supplemental sensitivity analysis, not used for the main Figure E. Pooled-tissue rankings can reflect tissue composition, and gene-set permutation FDR does not provide biological-replicate inference.

## 5. Paper-adapted Figure 3F

### Main: C5 state × tissue overlap

{chr(10).join(main_f_lines)}

- Main sets: C5-1_BM, C5-1_Spleen, C5-2_BM, and C5-2_Spleen, matching the paper's two-state × two-tissue structure.
- Main-set cell counts: {main_f_sample_summary}.
- Strict ≥20 cells per condition in all four main sets: {'estimable' if strict_main_f_estimable else 'NA (not estimable)'}.
- A four-set UpSet is used instead of hand-drawn four-circle geometry so every exact intersection remains unambiguous.
- The paper rule and robust sensitivity rule are plotted separately; all DE P values remain exploratory cell-level results.

### Supplemental: C5-2 across three tissues

{chr(10).join(supplemental_venn_lines)}

- Supplemental tissue cell counts: {tissue_sample_summary}.
- Strict ≥20 cells per condition in all three tissues: {'estimable' if strict_three_tissue_estimable else 'NA (not estimable), because at least one tissue fails the threshold; this is not a zero intersection'}.
- Thymus remains exploratory and does not determine the main Figure 3F conclusion.

## 6. Why the two cluster-frequency figures differ

- Stacked share is P(condition | tissue, cluster).
- Delta is P(cluster | tissue, condition)Tumor minus P(cluster | tissue, condition)Control.
- The zero-delta Tumor-share baseline is tissue-specific, not fixed at 50%: Thymus {baseline.loc['thymus'].iloc[0]:.3f}%, BM {baseline.loc['bone_marrow'].iloc[0]:.3f}%, Spleen {baseline.loc['spleen'].iloc[0]:.3f}%.
- BM C0 illustrates this: Tumor stacked share {bm_c0.tumor_share_within_tissue_cluster_pct:.3f}% exceeds its {bm_c0.tumor_share_zero_delta_baseline_pct:.3f}% baseline, so the frequency delta is {bm_c0.tumor_minus_control_percentage_points:+.3f} pp. All row-level invariants pass.

## 7. Raw-matrix CD56 audit

- Mouse Ncam1 is present in {feature_present_samples}/{raw_sample_count} raw 10x feature tables.
- Across {raw_sample_count} matrix.mtx.gz files: {ncam1_nonzero_entries} non-zero entries and {ncam1_total_counts:.0f} total counts.
- {ncam1_conclusion} No surrogate gene is relabeled as CD56.

## Paper mapping

Source: docs/blooda_adv-2024-014592-main.pdf.

- PDF p6 / journal p764: Figure 3D–F.
- PDF p7 / journal p765: results text and continued legend.
- Figure 3F specifies fold change≥2 and P<.05, implemented as log2FC≥1 and nominal P<0.05.
- Full scRNA-seq methods are deferred to an unavailable supplement.

## Key outputs

- figures/{RUN_DATE}_C5_subclustering_validation.png
- figures/{RUN_DATE}_full_UMAP_C5_replaced.png
- figures/{RUN_DATE}_IL4_CD94_condition_masks.png
- figures/{RUN_DATE}_IL4_CD94_complementarity.png
- figures/{RUN_DATE}_paper_adapted_Fig3D_volcano.png
- figures/{RUN_DATE}_Fig3E_all_refined_clusters_KEGG_GSEA_signed_NES_heatmap.png
- figures/{RUN_DATE}_paper_adapted_Fig3E_KEGG_ORA.png (supplemental)
- figures/{RUN_DATE}_paper_adapted_Fig3F_C5split_BM_Spleen_UpSet_paper.png
- figures/{RUN_DATE}_paper_adapted_Fig3F_C5split_BM_Spleen_UpSet_robust.png
- figures/{RUN_DATE}_paper_adapted_Fig3F_C5-2_tissue_Venn_paper.png
- figures/{RUN_DATE}_paper_adapted_Fig3F_C5-2_tissue_Venn_robust.png
- {RUN_DATE}_iNKT_C5_Fig3DEF_results.pdf
- tables/{RUN_DATE}_cell_index_C5_split.csv.gz
- tables/{RUN_DATE}_DE_all_units.csv.gz and tables/de/
- tables/{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz
- tables/{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz
- tables/{RUN_DATE}_KEGG_GSEA_heatmap_module_selection.csv
- tables/{RUN_DATE}_KEGG_GSEA_heatmap_pathway_selection.csv (compatibility copy)
- tables/{RUN_DATE}_KEGG_GSEA_heatmap_signed_NES.csv and *_FDR.csv
- tables/{RUN_DATE}_KEGG_GSEA_bundle_receipt.json
- tables/{RUN_DATE}_KEGG_GSEA_audit.json
- tables/{RUN_DATE}_KEGG_ORA_all.csv.gz (supplemental)
- tables/{RUN_DATE}_Fig3F_C5split_BM_Spleen_*
- tables/{RUN_DATE}_C5-2_Venn_*
- tables/{RUN_DATE}_Ncam1_CD56_raw_matrix_audit.csv
- tables/{RUN_DATE}_input_feature_universe_validation.json
- {RUN_DATE}_H5AD_roundtrip_validation.json
{h5ad_line}

## Statistical limit

Each tissue×condition has one sample. All P values are exploratory cell-level/pseudoreplicate comparisons and do not establish biological-replicate significance. Composition describes QC-retained captured cells, not absolute tissue abundance.
"""
    path.write_text(text, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.h5ad = args.h5ad.resolve()
    args.raw_root = args.raw_root.resolve()
    args.paper = args.paper.resolve()
    args.out_dir = args.out_dir.resolve()
    for required in (args.h5ad, args.raw_root, args.paper):
        if not required.exists():
            raise FileNotFoundError(required)
    figures_dir = args.out_dir / "figures"
    tables_dir = args.out_dir / "tables"
    de_unit_dir = tables_dir / "de"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 1

    print(f"Reading {args.h5ad}", flush=True)
    adata = sc.read_h5ad(args.h5ad)
    required_obs = {"sample", "condition", "tissue", ORIGINAL_CLUSTER_KEY}
    missing_obs = required_obs - set(adata.obs.columns)
    if missing_obs:
        raise KeyError(f"Missing obs columns: {sorted(missing_obs)}")
    if adata.raw is None or "counts" not in adata.layers:
        raise RuntimeError("The selected H5AD must contain adata.raw and layers['counts']")
    if "X_pca" not in adata.obsm or "X_umap" not in adata.obsm:
        raise RuntimeError("The selected H5AD must contain X_pca and X_umap")
    input_validation = validate_full_gene_input(adata)
    write_json(
        tables_dir / f"{RUN_DATE}_input_feature_universe_validation.json",
        input_validation,
    )

    c5, c5_metrics, selected_resolution = run_c5_subclustering(
        adata,
        n_hvg=args.n_hvg,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        resolutions=args.resolutions,
        seeds=args.stability_seeds,
        random_state=args.random_state,
    )
    mapping = install_refined_labels(adata, c5)
    c5_metrics.to_csv(tables_dir / f"{RUN_DATE}_C5_resolution_stability.csv", index=False)
    mapping.to_csv(
        tables_dir / f"{RUN_DATE}_cell_index_C5_split.csv.gz",
        index=False,
        compression="gzip",
    )
    (
        c5.obs.groupby(
            [C5_SUBCLUSTER_KEY, "sample", "condition", "tissue"],
            observed=True,
        )
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(tables_dir / f"{RUN_DATE}_C5_subcluster_counts.csv", index=False)
    )
    render_c5_validation(
        adata,
        c5,
        c5_metrics,
        selected_resolution,
        figures_dir / f"{RUN_DATE}_C5_subclustering_validation.png",
    )
    render_refined_full_umap(
        adata,
        figures_dir / f"{RUN_DATE}_full_UMAP_C5_replaced.png",
    )

    original_order = tuple(f"C{index}" for index in range(8))
    original_labels = "C" + adata.obs[ORIGINAL_CLUSTER_KEY].astype(str)
    original_proportions, original_delta = build_cluster_proportion_tables(
        adata.obs,
        original_labels,
        original_order,
    )
    refined_proportions, refined_delta = build_cluster_proportion_tables(
        adata.obs,
        adata.obs[REFINED_CLUSTER_KEY].astype(str),
        REFINED_ORDER,
    )
    original_proportions.to_csv(
        tables_dir / f"{RUN_DATE}_original_cluster_counts_and_proportions.csv",
        index=False,
    )
    original_delta.to_csv(
        tables_dir / f"{RUN_DATE}_original_denominator_reconciliation.csv",
        index=False,
    )
    refined_proportions.to_csv(
        tables_dir / f"{RUN_DATE}_C5split_cluster_counts_and_proportions.csv",
        index=False,
    )
    refined_delta.to_csv(
        tables_dir / f"{RUN_DATE}_C5split_denominator_reconciliation.csv",
        index=False,
    )
    render_denominator_reconciliation(
        original_delta,
        original_order,
        figures_dir / f"{RUN_DATE}_denominator_reconciliation_original_clusters.png",
        title_suffix="original C0–C7",
    )
    render_denominator_reconciliation(
        refined_delta,
        REFINED_ORDER,
        figures_dir / f"{RUN_DATE}_denominator_reconciliation_C5split.png",
        title_suffix="C5 replaced by C5-1/C5-2",
    )

    complementarity, marker_arrays = build_marker_complementarity(adata)
    complementarity.to_csv(
        tables_dir / f"{RUN_DATE}_IL4_CD94_complementarity.csv",
        index=False,
    )
    render_condition_masked_markers(
        adata,
        marker_arrays,
        figures_dir / f"{RUN_DATE}_IL4_CD94_condition_masks.png",
    )
    render_marker_overlap_bars(
        complementarity,
        figures_dir / f"{RUN_DATE}_IL4_CD94_complementarity.png",
    )

    de, de_status = run_requested_de(adata, de_unit_dir)
    de.to_csv(
        tables_dir / f"{RUN_DATE}_DE_all_units.csv.gz",
        index=False,
        compression="gzip",
    )
    de_status.to_csv(tables_dir / f"{RUN_DATE}_DE_unit_status.csv", index=False)
    render_volcano_figure_d(
        de,
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3D_volcano.png",
    )

    frozen_gmt = tables_dir / f"{RUN_DATE}_{args.kegg_library}.gmt"
    gene_sets, kegg_audit = load_kegg_gene_sets(
        args.kegg_library,
        frozen_gmt,
        args.kegg_gmt.resolve() if args.kegg_gmt is not None else None,
    )
    if args.kegg_gmt is not None:
        write_gmt(frozen_gmt, gene_sets)
        kegg_audit = {
            **kegg_audit,
            "frozen_copy": str(frozen_gmt),
            "frozen_copy_sha256": sha256_file(frozen_gmt),
        }
    ora = run_kegg_ora(de, list(adata.raw.var_names.astype(str)), gene_sets)
    ora.to_csv(
        tables_dir / f"{RUN_DATE}_KEGG_ORA_all.csv.gz",
        index=False,
        compression="gzip",
    )
    write_json(
        tables_dir / f"{RUN_DATE}_KEGG_library_audit.json",
        json_ready(kegg_audit),
    )
    render_pathway_figure_e(
        ora,
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3E_KEGG_ORA.png",
    )
    gsea_result_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz"
    gsea_ranking_path = (
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz"
    )
    gsea_receipt_path = tables_dir / f"{RUN_DATE}_KEGG_GSEA_bundle_receipt.json"
    if args.reuse_gsea_tables:
        print("Validating reusable dated GSEA bundle", flush=True)
        gsea, gsea_rankings, gsea_audit = load_validated_gsea_bundle(
            de,
            gene_sets,
            kegg_audit,
            tables_dir,
            permutation_num=args.gsea_permutations,
            seed=args.gsea_seed,
            threads=args.gsea_threads,
        )
    else:
        gsea, gsea_rankings, gsea_audit = run_preranked_gsea(
            de,
            gene_sets,
            permutation_num=args.gsea_permutations,
            seed=args.gsea_seed,
            threads=args.gsea_threads,
        )
        gsea_audit["gene_set_library"] = kegg_audit
        gsea.to_csv(
            gsea_result_path,
            index=False,
            compression="gzip",
        )
        gsea_rankings.to_csv(
            gsea_ranking_path,
            index=False,
            compression="gzip",
        )
        fresh_identity = gsea_computation_identity(
            permutation_num=args.gsea_permutations,
            seed=args.gsea_seed,
        )
        fresh_receipt = build_gsea_bundle_receipt(
            gsea_result_path,
            gsea_ranking_path,
            gsea=gsea,
            rankings=gsea_rankings,
            kegg_audit=kegg_audit,
            identity=fresh_identity,
            provenance="fresh_preranked_gsea_run",
        )
        write_json(gsea_receipt_path, json_ready(fresh_receipt))
        gsea_audit["reuse_validation"] = {
            "status": "not_reused_fresh_run",
            "reused": False,
            "receipt_path": str(gsea_receipt_path.resolve()),
            "gsea_results_sha256": fresh_receipt["files"]["gsea_results"]["sha256"],
            "ranked_genes_sha256": fresh_receipt["files"]["ranked_genes"]["sha256"],
        }
    gsea_selection = select_gsea_heatmap_modules(gsea)
    gsea_selection.to_csv(
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_heatmap_module_selection.csv",
        index=False,
    )
    gsea_selection.to_csv(
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_heatmap_pathway_selection.csv",
        index=False,
    )
    gsea_nes_matrix, gsea_fdr_matrix = build_gsea_heatmap_matrices(
        gsea,
        gsea_selection,
    )
    gsea_nes_matrix.to_csv(
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_heatmap_signed_NES.csv"
    )
    gsea_fdr_matrix.to_csv(
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_heatmap_FDR.csv"
    )
    write_json(
        tables_dir / f"{RUN_DATE}_KEGG_GSEA_audit.json",
        json_ready(gsea_audit),
    )
    render_gsea_heatmap_figure_e(
        gsea,
        gsea_selection,
        de_status,
        figures_dir
        / f"{RUN_DATE}_Fig3E_all_refined_clusters_KEGG_GSEA_signed_NES_heatmap.png",
    )

    main_f_status = (
        de_status.set_index("unit_id").reindex(FIGURE_F_MAIN_SET_ORDER).reset_index()
    )
    strict_main_f_estimable = bool(
        main_f_status["passes_recommended_min20"].fillna(False).astype(bool).all()
    )
    c5_tissue_status = de_status[de_status["scope"].eq("C5-2_tissue")].copy()
    strict_three_tissue_estimable = bool(
        c5_tissue_status["passes_recommended_min20"].fillna(False).astype(bool).all()
    )
    main_f_summaries: dict[str, pd.DataFrame] = {}
    supplemental_venn_summaries: dict[str, pd.DataFrame] = {}
    for threshold, suffix in (
        (PAPER_VENN_THRESHOLD, "paper"),
        (ROBUST_THRESHOLD, "robust"),
    ):
        main_membership, main_summary, main_sets = build_figure_f_main_membership(
            de,
            threshold,
        )
        main_all_four = int(
            main_summary.set_index("metric").loc["all_four", "n_genes"]
        )
        main_summary["strict_min20_all_four_estimable"] = strict_main_f_estimable
        main_summary["strict_min20_all_four_n_genes"] = (
            float(main_all_four) if strict_main_f_estimable else np.nan
        )
        main_summary["strict_min20_interpretation"] = (
            "estimable" if strict_main_f_estimable else "NA_not_estimable_not_zero"
        )
        main_membership.to_csv(
            tables_dir
            / f"{RUN_DATE}_Fig3F_C5split_BM_Spleen_membership_{suffix}.csv",
            index=False,
        )
        main_summary.to_csv(
            tables_dir / f"{RUN_DATE}_Fig3F_C5split_BM_Spleen_summary_{suffix}.csv",
            index=False,
        )
        main_f_summaries[threshold] = main_summary
        render_figure_f_main_upset(
            main_sets,
            threshold,
            figures_dir
            / f"{RUN_DATE}_paper_adapted_Fig3F_C5split_BM_Spleen_UpSet_{suffix}.png",
            de_status,
        )

        membership, summary, sets = build_venn_membership(de, threshold)
        exploratory_all_three = int(
            summary.set_index("metric").loc["all_three", "n_genes"]
        )
        summary["strict_min20_all_three_estimable"] = strict_three_tissue_estimable
        summary["strict_min20_all_three_n_genes"] = (
            float(exploratory_all_three) if strict_three_tissue_estimable else np.nan
        )
        summary["strict_min20_interpretation"] = (
            "estimable" if strict_three_tissue_estimable else "NA_not_estimable_not_zero"
        )
        membership.to_csv(
            tables_dir / f"{RUN_DATE}_C5-2_Venn_membership_{suffix}.csv",
            index=False,
        )
        summary.to_csv(
            tables_dir / f"{RUN_DATE}_C5-2_Venn_summary_{suffix}.csv",
            index=False,
        )
        supplemental_venn_summaries[threshold] = summary
        render_venn_figure_f(
            sets,
            threshold,
            figures_dir
            / f"{RUN_DATE}_paper_adapted_Fig3F_C5-2_tissue_Venn_{suffix}.png",
            de_status,
        )

    ncam1_audit = audit_raw_ncam1(args.raw_root)
    ncam1_audit.to_csv(
        tables_dir / f"{RUN_DATE}_Ncam1_CD56_raw_matrix_audit.csv",
        index=False,
    )
    selected_metric = c5_metrics.loc[c5_metrics["selected"]].iloc[0]

    adata.uns[f"c5_subclustering_{RUN_DATE}"] = {
        "source_cluster_key": ORIGINAL_CLUSTER_KEY,
        "source_cluster": "5",
        "new_cluster_key": REFINED_CLUSTER_KEY,
        "c5_subcluster_key": C5_SUBCLUSTER_KEY,
        "n_c5_cells": int(c5.n_obs),
        "n_hvg": int(c5.var["highly_variable"].sum()),
        "n_pcs": int(c5.obsm["X_pca"].shape[1]),
        "n_neighbors": int(args.n_neighbors),
        "min_dist": float(args.min_dist),
        "selected_resolution": float(selected_resolution),
        "selected_partition_seed": int(selected_metric["selected_partition_seed"]),
        "original_umap_kmeans_ari": float(
            selected_metric["selected_original_umap_kmeans_ari"]
        ),
        "original_umap_silhouette": float(
            selected_metric["selected_original_umap_silhouette"]
        ),
        "original_umap_median_x_gap": float(selected_metric["selected_original_umap_median_x_gap"]),
        "original_umap_q05_q95_gap": float(selected_metric["selected_original_umap_q05_q95_gap"]),
        "selected_seed_mean_ari": float(selected_metric["selected_seed_mean_ari"]),
        "stability_seeds": list(map(int, args.stability_seeds)),
        "label_rule": "lower_original_UMAP1_median=C5-1; higher_original_UMAP1_median=C5-2",
    }
    full_h5ad = args.out_dir / f"{RUN_DATE}_inkt_selected_umap_C5_split.h5ad"
    c5_h5ad = args.out_dir / f"{RUN_DATE}_C5_only_reclustered.h5ad"
    if not args.skip_h5ad:
        print(f"Writing {full_h5ad}", flush=True)
        adata.write_h5ad(full_h5ad, compression="gzip")
        c5.write_h5ad(c5_h5ad, compression="gzip")
        h5ad_roundtrip = validate_h5ad_roundtrip(full_h5ad, c5_h5ad, adata, c5)
    else:
        h5ad_roundtrip = {
            "status": "skipped_by_flag",
            "files": {},
        }
    write_json(
        args.out_dir / f"{RUN_DATE}_H5AD_roundtrip_validation.json",
        json_ready(h5ad_roundtrip),
    )


    strict_main_f_intersections = (
        {
            threshold: int(frame.set_index("metric").loc["all_four", "n_genes"])
            for threshold, frame in main_f_summaries.items()
        }
        if strict_main_f_estimable
        else None
    )
    strict_supplemental_intersections = (
        {
            threshold: int(
                frame.set_index("metric").loc["all_three", "n_genes"]
            )
            for threshold, frame in supplemental_venn_summaries.items()
        }
        if strict_three_tissue_estimable
        else None
    )

    summary_payload = {
        "run_date": RUN_DATE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "input_validation": input_validation,
        "c5_split": {
            "n_original_c5": int(c5.n_obs),
            "counts": c5.obs[C5_SUBCLUSTER_KEY].astype(str).value_counts().to_dict(),
            "selected_resolution": float(selected_resolution),
            "mean_pairwise_ari": float(selected_metric["mean_pairwise_ari"]),
            "min_pairwise_ari": float(selected_metric["min_pairwise_ari"]),
            "selected_partition_seed": int(selected_metric["selected_partition_seed"]),
            "selected_seed_mean_ari": float(selected_metric["selected_seed_mean_ari"]),
            "original_umap_kmeans_ari": float(
                selected_metric["selected_original_umap_kmeans_ari"]
            ),
            "original_umap_silhouette": float(
                selected_metric["selected_original_umap_silhouette"]
            ),
            "original_umap_median_x_gap": float(selected_metric["selected_original_umap_median_x_gap"]),
            "original_umap_q05_q95_gap": float(selected_metric["selected_original_umap_q05_q95_gap"]),
            "label_rule": "C5-2 is the right-hand original-UMAP island",
        },
        "de_status": de_status.to_dict(orient="records"),
        "top_genes": top_gene_summary(de),
        "kegg_library": kegg_audit,
        "figure3e_gsea": {
            "parameters": gsea_audit["parameters"],
            "direction_definition": {
                "positive_nes": "Tumor_T2_enriched",
                "negative_nes": "Control_Ctrl_enriched",
            },
            "significant_fdr_0_05_counts": {
                str(cluster): {
                    "tumor": int(
                        (
                            gsea["cluster"].eq(cluster)
                            & gsea["significant_fdr_0_05"].astype(bool)
                            & gsea["nes"].gt(0)
                        ).sum()
                    ),
                    "control": int(
                        (
                            gsea["cluster"].eq(cluster)
                            & gsea["significant_fdr_0_05"].astype(bool)
                            & gsea["nes"].lt(0)
                        ).sum()
                    ),
                }
                for cluster in GSEA_CLUSTERS
            },
            "heatmap_modules": gsea_selection.to_dict(orient="records"),
            "display_rule": (
                "one representative raw KEGG term per de-redundant module; "
                "no NES averaging"
            ),
            "n_significant_raw_terms_accounted": int(
                gsea_selection["n_source_terms"].sum()
            ),
            "reuse_validation": gsea_audit.get("reuse_validation", {}),
            "full_results_table": str(
                tables_dir / f"{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz"
            ),
            "ranking_table": str(
                tables_dir
                / f"{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz"
            ),
        },
        "figure3f": {
            "main_c5split_bm_spleen": {
                "results": {
                    threshold: frame.to_dict(orient="records")
                    for threshold, frame in main_f_summaries.items()
                },
                "strict_min20_all_four_estimable": strict_main_f_estimable,
                "strict_min20_all_four_intersections": strict_main_f_intersections,
                "tissue_de_status": main_f_status.to_dict(orient="records"),
            },
            "supplemental_c5_2_three_tissue": {
                "exploratory_results": {
                    threshold: frame.to_dict(orient="records")
                    for threshold, frame in supplemental_venn_summaries.items()
                },
                "strict_min20_all_three_estimable": strict_three_tissue_estimable,
                "strict_min20_all_three_intersections": strict_supplemental_intersections,
                "tissue_de_status": c5_tissue_status.to_dict(orient="records"),
            },
        },
        "ncam1_raw": {
            "feature_present_samples": int(ncam1_audit["feature_present"].sum()),
            "nonzero_entries": int(ncam1_audit["raw_matrix_nonzero_entries"].sum()),
            "total_counts": float(ncam1_audit["raw_matrix_total_counts"].sum()),
            "status_counts": ncam1_audit["audit_status"]
            .value_counts()
            .to_dict(),
        },
        "h5ad_roundtrip": h5ad_roundtrip,
        "statistical_warning": STATISTICAL_WARNING,
    }
    write_json(
        args.out_dir / f"{RUN_DATE}_summary.json",
        json_ready(summary_payload),
    )
    write_readme(
        args.out_dir / f"{RUN_DATE}_README.md",
        c5=c5,
        metrics=c5_metrics,
        selected_resolution=selected_resolution,
        complementarity=complementarity,
        de_status=de_status,
        de=de,
        ora=ora,
        gsea=gsea,
        gsea_selection=gsea_selection,
        gsea_audit=gsea_audit,
        main_f_summaries=main_f_summaries,
        supplemental_venn_summaries=supplemental_venn_summaries,
        ncam1_audit=ncam1_audit,
        original_delta=original_delta,
        kegg_audit=kegg_audit,
        h5ad_written=not args.skip_h5ad,
    )

    overview_images = [
        figures_dir / f"{RUN_DATE}_C5_subclustering_validation.png",
        figures_dir / f"{RUN_DATE}_full_UMAP_C5_replaced.png",
        figures_dir / f"{RUN_DATE}_IL4_CD94_condition_masks.png",
        figures_dir / f"{RUN_DATE}_IL4_CD94_complementarity.png",
        figures_dir / f"{RUN_DATE}_denominator_reconciliation_original_clusters.png",
        figures_dir / f"{RUN_DATE}_denominator_reconciliation_C5split.png",
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3D_volcano.png",
        figures_dir
        / f"{RUN_DATE}_Fig3E_all_refined_clusters_KEGG_GSEA_signed_NES_heatmap.png",
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3F_C5split_BM_Spleen_UpSet_paper.png",
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3F_C5split_BM_Spleen_UpSet_robust.png",
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3F_C5-2_tissue_Venn_paper.png",
        figures_dir / f"{RUN_DATE}_paper_adapted_Fig3F_C5-2_tissue_Venn_robust.png",
    ]
    overview_pdf = args.out_dir / f"{RUN_DATE}_iNKT_C5_Fig3DEF_results.pdf"
    build_overview_pdf(overview_images, overview_pdf)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_date": RUN_DATE,
        "inputs": {
            "selected_h5ad": {
                "path": str(args.h5ad),
                "size_bytes": args.h5ad.stat().st_size,
                "sha256": sha256_file(args.h5ad),
            },
            "paper": {
                "path": str(args.paper),
                "size_bytes": args.paper.stat().st_size,
                "sha256": sha256_file(args.paper),
                "figure_mapping": "Figure 3D-F; PDF p6-7 / journal p764-765",
            },
            "raw_10x_root": str(args.raw_root),
            "script": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256_file(Path(__file__).resolve()),
            },
        },
        "validation": {
            "input_feature_universe": input_validation,
            "gsea_bundle_reuse": gsea_audit.get("reuse_validation", {}),
            "h5ad_roundtrip": h5ad_roundtrip,
        },
        "parameters": {
            "n_hvg": args.n_hvg,
            "n_pcs": args.n_pcs,
            "n_neighbors": args.n_neighbors,
            "min_dist": args.min_dist,
            "resolutions": args.resolutions,
            "stability_seeds": args.stability_seeds,
            "random_state": args.random_state,
            "kegg_library": args.kegg_library,
            "gsea_clusters": list(GSEA_CLUSTERS),
            "gsea_rank_metric": "scanpy_wilcoxon_score_T2_vs_Ctrl",
            "gsea_positive_direction": "Tumor_T2",
            "gsea_negative_direction": "Control_Ctrl",
            "gsea_permutations": args.gsea_permutations,
            "gsea_seed": args.gsea_seed,
            "gsea_threads": args.gsea_threads,
            "gsea_min_size": GSEA_MIN_SIZE,
            "gsea_max_size": GSEA_MAX_SIZE,
            "gsea_weight": GSEA_WEIGHT,
            "gsea_fdr_threshold": GSEA_FDR_THRESHOLD,
            "gsea_heatmap_n_modules": int(len(gsea_selection)),
            "gsea_heatmap_n_significant_source_terms": int(gsea_selection["n_source_terms"].sum()),
            "reuse_gsea_tables": bool(args.reuse_gsea_tables),
            "paper_venn_threshold": PAPER_VENN_THRESHOLD,
            "robust_threshold": ROBUST_THRESHOLD,
            "figure3f_main_sets": list(FIGURE_F_MAIN_SET_ORDER),
            "figure3f_main_min_cells_per_condition": RECOMMENDED_MIN_CELLS_PER_CONDITION,
            "figure3f_supplemental_sets": [f"C5-2__{tissue}" for tissue in VENN_TISSUES],
        },
        "outputs": {
            "overview_pdf": str(overview_pdf),
            "full_h5ad": str(full_h5ad) if not args.skip_h5ad else None,
            "c5_h5ad": str(c5_h5ad) if not args.skip_h5ad else None,
            "figures_dir": str(figures_dir),
            "tables_dir": str(tables_dir),
            "figure3e_gsea": str(
                figures_dir
                / f"{RUN_DATE}_Fig3E_all_refined_clusters_KEGG_GSEA_signed_NES_heatmap.png"
            ),
            "gsea_full_table": str(
                tables_dir / f"{RUN_DATE}_KEGG_GSEA_all_clusters.csv.gz"
            ),
            "gsea_rankings": str(
                tables_dir
                / f"{RUN_DATE}_KEGG_GSEA_ranked_genes_all_clusters.csv.gz"
            ),
            "gsea_module_selection": str(
                tables_dir
                / f"{RUN_DATE}_KEGG_GSEA_heatmap_module_selection.csv"
            ),
            "gsea_pathway_selection_compatibility_copy": str(
                tables_dir
                / f"{RUN_DATE}_KEGG_GSEA_heatmap_pathway_selection.csv"
            ),
            "gsea_bundle_receipt": str(
                tables_dir / f"{RUN_DATE}_KEGG_GSEA_bundle_receipt.json"
            ),
            "gsea_audit": str(
                tables_dir / f"{RUN_DATE}_KEGG_GSEA_audit.json"
            ),
            "input_validation_json": str(tables_dir / f"{RUN_DATE}_input_feature_universe_validation.json"),
            "h5ad_roundtrip_json": str(args.out_dir / f"{RUN_DATE}_H5AD_roundtrip_validation.json"),
        },
        "statistical_warning": STATISTICAL_WARNING,
    }
    write_json(
        args.out_dir / f"{RUN_DATE}_input_output_manifest.json",
        json_ready(manifest),
    )
    print(json.dumps(json_ready(summary_payload), indent=2)[:8000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
