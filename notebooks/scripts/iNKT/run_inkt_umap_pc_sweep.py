#!/usr/bin/env python3
"""Run the legacy-aligned iNKT UMAP 50/100/200-PC sensitivity analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from io import BytesIO
from pathlib import Path
from typing import Any
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
from sklearn.manifold import trustworthiness
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from inkt_palette import (
    SAMPLE_DISPLAY_ORDER,
    SAMPLE_PALETTE,
    configure_inkt_palettes,
    palette_validation_failures,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_H5AD = (
    ROOT
    / "output/iNKT_legacy_ppt_qc_runs/20260825_095422/preprocess/inkt_scanpy_tutorial_processed.h5ad"
)
DEFAULT_SCORES = (
    ROOT
    / "output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/standardized_legacy/standardized_scores_per_cell.csv.gz"
)
DEFAULT_LEGACY_PPT = ROOT / "input/iNKT/iNKT.pptx"
DEFAULT_OUT_DIR = ROOT / "output/iNKT_reproduction_deck/umap_pc_sweep"
DEFAULT_PCS = (50, 100, 200)
CLUSTER_KEY = "leiden_res_0_5"
PROGRAM_COLUMNS = ("inkt1", "inkt2", "inkt17")
PROGRAM_LABELS = {"inkt1": "iNKT1", "inkt2": "iNKT2", "inkt17": "NKT17"}
PROGRAM_COLORS = {"inkt1": "#1F8A8A", "inkt2": "#ED7D31", "inkt17": "#7A5AF8"}
CLUSTER_COLORS = tuple(plt.get_cmap("tab10")(index) for index in range(10))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_pcs(pcs: Sequence[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in pcs)
    if not values or any(value < 2 for value in values):
        raise ValueError(f"PC values must all be >= 2: {values}")
    if len(set(values)) != len(values):
        raise ValueError(f"PC values must be unique: {values}")
    if values != tuple(sorted(values)):
        raise ValueError(f"PC values must be sorted: {values}")
    return values


def label_centroids(
    embedding: np.ndarray,
    labels: Sequence[str],
    order: Sequence[str],
) -> np.ndarray:
    label_array = np.asarray(labels, dtype=str)
    centroids = []
    for label in order:
        mask = label_array == label
        if not mask.any():
            raise ValueError(f"No cells found for label {label!r}")
        centroids.append(np.asarray(embedding[mask], dtype=float).mean(axis=0))
    return np.asarray(centroids)


def procrustes_concordance(reference: np.ndarray, candidate: np.ndarray) -> float:
    if reference.shape != candidate.shape or reference.ndim != 2 or reference.shape[1] != 2:
        raise ValueError(f"Expected matched n x 2 coordinates, got {reference.shape}, {candidate.shape}")
    _, _, disparity = procrustes(reference, candidate)
    return float(np.clip(1.0 - disparity, 0.0, 1.0))


def distance_concordance(reference: np.ndarray, candidate: np.ndarray) -> float:
    reference_distances = pdist(reference)
    candidate_distances = pdist(candidate)
    if np.allclose(reference_distances, reference_distances[0]) or np.allclose(
        candidate_distances, candidate_distances[0]
    ):
        return 0.0
    return float(pearsonr(reference_distances, candidate_distances).statistic)


def knn_indices(embedding: np.ndarray, n_neighbors: int) -> np.ndarray:
    effective = min(int(n_neighbors) + 1, len(embedding))
    model = NearestNeighbors(n_neighbors=effective, metric="euclidean", n_jobs=1)
    indices = model.fit(np.asarray(embedding)).kneighbors(return_distance=False)
    # sklearn does not guarantee that self is first for tied points. Remove it row-wise.
    rows = np.arange(len(embedding))[:, None]
    without_self = [row[row != index][:n_neighbors] for index, row in enumerate(indices)]
    if any(len(row) != n_neighbors for row in without_self):
        raise RuntimeError("Could not construct the requested UMAP-neighbor index matrix")
    return np.asarray(without_self, dtype=int)


def label_knn_purity(labels: Sequence[str], neighbors: np.ndarray) -> float:
    values = np.asarray(labels, dtype=str)
    return float(np.mean(values[neighbors] == values[:, None]))


def program_knn_smoothness(scores: np.ndarray, neighbors: np.ndarray) -> tuple[float, list[float]]:
    per_program = []
    for index in range(scores.shape[1]):
        values = np.asarray(scores[:, index], dtype=float)
        neighbor_mean = values[neighbors].mean(axis=1)
        if np.std(values) == 0 or np.std(neighbor_mean) == 0:
            correlation = 0.0
        else:
            correlation = float(np.corrcoef(values, neighbor_mean)[0, 1])
        per_program.append(correlation)
    return float(np.mean(per_program)), per_program


def safe_silhouette(
    embedding: np.ndarray,
    labels: Sequence[str],
    *,
    sample_size: int,
    random_state: int,
) -> float:
    values = np.asarray(labels, dtype=str)
    if len(np.unique(values)) < 2:
        return float("nan")
    effective = min(int(sample_size), len(embedding))
    return float(
        silhouette_score(
            np.asarray(embedding),
            values,
            metric="euclidean",
            sample_size=effective,
            random_state=random_state,
        )
    )


def centroid_separation_ratio(embedding: np.ndarray, labels: Sequence[str]) -> float:
    values = np.asarray(labels, dtype=str)
    categories = sorted(np.unique(values))
    centroids = label_centroids(embedding, values, categories)
    between = float(pdist(centroids).mean())
    within = []
    for index, category in enumerate(categories):
        points = embedding[values == category]
        within.append(float(np.sqrt(np.mean(np.sum((points - centroids[index]) ** 2, axis=1)))))
    denominator = float(np.mean(within))
    return between / denominator if denominator > 0 else float("inf")


def legacy_sample_centroids(
    pptx_path: Path,
    *,
    media_path: str = "ppt/media/image8.png",
) -> tuple[np.ndarray, Image.Image]:
    with ZipFile(pptx_path) as archive:
        image = Image.open(BytesIO(archive.read(media_path))).convert("RGB")
        image.load()
    pixels = np.asarray(image)
    plot_width = min(1400, int(image.width * 0.70))
    centroids = []
    for sample in SAMPLE_DISPLAY_ORDER:
        color = tuple(bytes.fromhex(SAMPLE_PALETTE[sample][1:]))
        mask = np.all(pixels[:, :plot_width] == color, axis=2)
        y, x = np.nonzero(mask)
        if len(x) < 100:
            raise RuntimeError(f"Insufficient exact pixels for {sample} in {media_path}: {len(x)}")
        centroids.append([float(x.mean()), float(-y.mean())])
    return np.asarray(centroids), image


def density_mass_threshold(smoothed_histogram: np.ndarray, mass: float = 0.80) -> float:
    if not 0 < mass < 1:
        raise ValueError("mass must be strictly between 0 and 1")
    values = np.asarray(smoothed_histogram, dtype=float).ravel()
    positive = values[values > 0]
    if not len(positive):
        return float("nan")
    ordered = np.sort(positive)[::-1]
    cumulative = np.cumsum(ordered)
    index = int(np.searchsorted(cumulative, mass * cumulative[-1], side="left"))
    return float(ordered[min(index, len(ordered) - 1)])


def add_density_envelopes(
    ax: Any,
    embedding: np.ndarray,
    cluster_labels: Sequence[str],
    *,
    extent: tuple[float, float, float, float],
    bins: int = 120,
    mass: float = 0.80,
) -> None:
    labels = np.asarray(cluster_labels, dtype=str)
    x_min, x_max, y_min, y_max = extent
    x_edges = np.linspace(x_min, x_max, bins + 1)
    y_edges = np.linspace(y_min, y_max, bins + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    for cluster in sorted(np.unique(labels), key=lambda value: int(value)):
        points = embedding[labels == cluster]
        histogram, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=(x_edges, y_edges))
        smoothed = gaussian_filter(histogram, sigma=1.5)
        threshold = density_mass_threshold(smoothed, mass=mass)
        color = CLUSTER_COLORS[int(cluster) % len(CLUSTER_COLORS)]
        if np.isfinite(threshold) and smoothed.max() > threshold:
            ax.contour(
                x_centers,
                y_centers,
                smoothed.T,
                levels=[threshold],
                colors=[color],
                linewidths=1.15,
                zorder=4,
            )
        center = np.median(points, axis=0)
        ax.text(
            center[0],
            center[1],
            f"c{cluster}",
            ha="center",
            va="center",
            fontsize=7,
            fontweight="bold",
            color="white",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": color, "edgecolor": "white", "lw": 0.5},
            zorder=5,
        )


def plot_pc_embedding(
    ax: Any,
    embedding: np.ndarray,
    sample_labels: Sequence[str],
    cluster_labels: Sequence[str],
    pc: int,
    *,
    n_neighbors: int,
    min_dist: float,
    include_legend: bool,
) -> None:
    samples = np.asarray(sample_labels, dtype=str)
    for sample in SAMPLE_DISPLAY_ORDER:
        mask = samples == sample
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            s=1.2,
            alpha=0.72,
            c=SAMPLE_PALETTE[sample],
            edgecolors="none",
            rasterized=True,
            label=sample,
            zorder=2,
        )
    x_min, y_min = embedding.min(axis=0)
    x_max, y_max = embedding.max(axis=0)
    x_pad = max((x_max - x_min) * 0.03, 0.1)
    y_pad = max((y_max - y_min) * 0.03, 0.1)
    extent = (x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    add_density_envelopes(ax, embedding, cluster_labels, extent=extent)
    ax.set_title(f"{pc} PCs", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.grid(color="#D7DEE5", linewidth=0.35, alpha=0.55)
    ax.text(
        0.01,
        0.01,
        f"UMAP1 [{x_min:.2f}, {x_max:.2f}]\nUMAP2 [{y_min:.2f}, {y_max:.2f}]",
        transform=ax.transAxes,
        fontsize=6.8,
        ha="left",
        va="bottom",
        color="#263746",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#CCD6DF"},
        zorder=6,
    )
    ax.text(
        0.99,
        0.01,
        f"neighbors={n_neighbors}; min_dist={min_dist:g}\n80% Leiden density envelopes",
        transform=ax.transAxes,
        fontsize=6.8,
        ha="right",
        va="bottom",
        color="#52616B",
        zorder=6,
    )
    if include_legend:
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markersize=5,
                markerfacecolor=SAMPLE_PALETTE[sample],
                markeredgecolor="none",
                label=sample,
            )
            for sample in SAMPLE_DISPLAY_ORDER
        ]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=7)


def render_candidate_figures(
    out_dir: Path,
    embeddings: Mapping[int, np.ndarray],
    sample_labels: Sequence[str],
    cluster_labels: Sequence[str],
    *,
    n_neighbors: int,
    min_dist: float,
) -> dict[str, str]:
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for pc, embedding in embeddings.items():
        figure, ax = plt.subplots(figsize=(8.8, 6.5))
        plot_pc_embedding(
            ax,
            embedding,
            sample_labels,
            cluster_labels,
            pc,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            include_legend=True,
        )
        figure.suptitle("iNKT UMAP PC sensitivity: tissue hue + condition shade", fontsize=14, fontweight="bold")
        figure.tight_layout(rect=(0, 0, 0.86, 0.96))
        path = figures_dir / f"umap_npcs_{pc:03d}.png"
        figure.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(figure)
        outputs[f"pc_{pc}"] = str(path)

    figure, axes = plt.subplots(1, len(embeddings), figsize=(17.5, 5.5), sharex=False, sharey=False)
    for index, (pc, embedding) in enumerate(embeddings.items()):
        plot_pc_embedding(
            axes[index],
            embedding,
            sample_labels,
            cluster_labels,
            pc,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            include_legend=index == len(embeddings) - 1,
        )
    figure.suptitle(
        "UMAP PC sensitivity — only n_pcs changes; fixed current Leiden labels",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 0.95, 0.95), w_pad=2.0)
    comparison_path = figures_dir / "umap_pc_sweep_tissue_cluster.png"
    figure.savefig(comparison_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    outputs["comparison"] = str(comparison_path)
    return outputs


def render_program_scores(
    path: Path,
    embedding: np.ndarray,
    scores: pd.DataFrame,
    selected_pc: int,
) -> None:
    values = scores.loc[:, PROGRAM_COLUMNS].to_numpy(dtype=float)
    limit = float(np.nanpercentile(np.abs(values), 98))
    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), sharex=True, sharey=True)
    for index, program in enumerate(PROGRAM_COLUMNS):
        scatter = axes[index].scatter(
            embedding[:, 0],
            embedding[:, 1],
            c=values[:, index],
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
            s=1.2,
            edgecolors="none",
            rasterized=True,
        )
        axes[index].set_title(PROGRAM_LABELS[program], color=PROGRAM_COLORS[program], fontweight="bold")
        axes[index].set_xlabel("UMAP1")
        if index == 0:
            axes[index].set_ylabel("UMAP2")
        figure.colorbar(scatter, ax=axes[index], fraction=0.046, pad=0.03, label="standardized score")
    figure.suptitle(
        f"Selected {selected_pc}-PC UMAP: continuous subtype-program scores (not discrete cell labels)",
        fontsize=13,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def normalize_metric_ranks(metrics: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    missing = set(weights) - set(metrics)
    if missing:
        raise ValueError(f"Missing selection metrics: {sorted(missing)}")
    if not math.isclose(sum(weights.values()), 1.0, rel_tol=0, abs_tol=1e-9):
        raise ValueError(f"Selection weights must sum to 1.0: {weights}")
    composite = pd.Series(0.0, index=metrics.index, dtype=float)
    for column, weight in weights.items():
        values = pd.to_numeric(metrics[column], errors="raise")
        if values.nunique(dropna=False) == 1:
            ranks = pd.Series(0.5, index=metrics.index)
        else:
            ranks = values.rank(method="average", pct=True)
        composite = composite + weight * ranks
    return composite


def select_recommended_pc(metrics: pd.DataFrame) -> tuple[int, dict[str, float]]:
    weights = {
        "tissue_silhouette": 0.25,
        "leiden_silhouette": 0.25,
        "trustworthiness": 0.15,
        "mean_stability": 0.15,
        "legacy_centroid_concordance": 0.10,
        "program_knn_smoothness": 0.10,
    }
    composite = normalize_metric_ranks(metrics, weights)
    metrics["selection_score"] = composite
    ordered = metrics.sort_values(["selection_score", "n_pcs"], ascending=[False, True])
    return int(ordered.iloc[0]["n_pcs"]), weights


def render_metrics_figure(path: Path, metrics: pd.DataFrame, selected_pc: int) -> None:
    columns = [
        "n_pcs",
        "tissue_silhouette",
        "leiden_silhouette",
        "tissue_knn_purity",
        "leiden_knn_purity",
        "trustworthiness",
        "mean_stability",
        "legacy_centroid_concordance",
        "program_knn_smoothness",
        "selection_score",
    ]
    labels = [
        "PCs",
        "Tissue\nsilhouette",
        "Leiden\nsilhouette",
        "Tissue\nkNN purity",
        "Leiden\nkNN purity",
        "Trust-\nworthiness",
        "Cross-PC\nstability",
        "Legacy\nconcordance",
        "Program\nsmoothness",
        "Selection\nscore",
    ]
    table_values = []
    row_colors = []
    for _, row in metrics.sort_values("n_pcs").iterrows():
        table_values.append([f"{int(row['n_pcs'])}"] + [f"{row[column]:.3f}" for column in columns[1:]])
        row_colors.append("#E8F6EF" if int(row["n_pcs"]) == selected_pc else "#FFFFFF")
    figure, ax = plt.subplots(figsize=(14.2, 3.1))
    ax.axis("off")
    table = ax.table(
        cellText=table_values,
        colLabels=labels,
        cellLoc="center",
        colLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2)
    table.scale(1, 1.85)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("#D7DEE5")
        if row == 0:
            cell.set_facecolor("#102A43")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor(row_colors[row - 1])
            if column == 0 and int(table_values[row - 1][0]) == selected_pc:
                cell.get_text().set_fontweight("bold")
                cell.get_text().set_color("#2E8B57")
    ax.set_title(
        f"UMAP PC selection metrics — recommended: {selected_pc} PCs",
        fontsize=14,
        fontweight="bold",
        pad=16,
    )
    figure.tight_layout()
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def render_scanpy_selected_figures(
    adata: Any,
    output_dir: Path,
    marker_genes: Sequence[str],
    cluster_key: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(
        adata,
        color=["sample", "condition", "tissue", cluster_key],
        wspace=0.35,
        show=False,
    )
    plt.savefig(output_dir / "umap_selected_sample_condition_tissue_cluster.png", dpi=180, bbox_inches="tight")
    plt.close("all")

    if marker_genes:
        sc.pl.umap(
            adata,
            color=["condition", cluster_key] + list(marker_genes[:8]),
            ncols=4,
            wspace=0.35,
            show=False,
        )
        plt.savefig(output_dir / "umap_selected_condition_cluster_top_markers.png", dpi=180, bbox_inches="tight")
        plt.close("all")

    if "dpt_pseudotime" in adata.obs:
        sc.pl.paga(
            adata,
            color=cluster_key,
            threshold=0.03,
            show=False,
        )
        plt.savefig(
            output_dir / "trajectory_selected_paga_cluster_graph.png",
            dpi=180,
            bbox_inches="tight",
        )
        plt.close("all")
        sc.pl.umap(
            adata,
            color=[cluster_key, "dpt_pseudotime"],
            wspace=0.35,
            show=False,
        )
        plt.savefig(
            output_dir / "trajectory_selected_umap_cluster_dpt.png",
            dpi=180,
            bbox_inches="tight",
        )
        plt.close("all")
        sc.pl.umap(
            adata,
            color=[cluster_key, "dpt_pseudotime", "condition", "tissue"],
            wspace=0.35,
            show=False,
        )
        plt.savefig(output_dir / "trajectory_selected_umap_dpt_pseudotime.png", dpi=180, bbox_inches="tight")
        plt.close("all")
        sc.pl.paga_compare(
            adata,
            basis="umap",
            color="dpt_pseudotime",
            threshold=0.03,
            show=False,
        )
        plt.savefig(output_dir / "trajectory_selected_paga_compare_umap_pseudotime.png", dpi=180, bbox_inches="tight")
        plt.close("all")


def load_marker_genes(h5ad_path: Path) -> list[str]:
    marker_path = h5ad_path.parent / "tables/marker_genes_used.json"
    if not marker_path.exists():
        return []
    groups = json.loads(marker_path.read_text(encoding="utf-8"))
    return [gene for genes in groups.values() for gene in genes]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--standardized-scores", type=Path, default=DEFAULT_SCORES)
    parser.add_argument("--legacy-ppt", type=Path, default=DEFAULT_LEGACY_PPT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--pcs", type=int, nargs="+", default=list(DEFAULT_PCS))
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--min-dist", type=float, default=0.5)
    parser.add_argument("--spread", type=float, default=1.0)
    parser.add_argument("--metric", default="euclidean")
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--silhouette-sample-size", type=int, default=5000)
    parser.add_argument("--trustworthiness-sample-size", type=int, default=3000)
    parser.add_argument("--cluster-key", default=CLUSTER_KEY)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.h5ad = args.h5ad.resolve()
    args.standardized_scores = args.standardized_scores.resolve()
    args.legacy_ppt = args.legacy_ppt.resolve()
    args.out_dir = args.out_dir.resolve()
    pcs = validate_pcs(args.pcs)
    if args.n_neighbors < 2:
        raise ValueError("--n-neighbors must be >= 2")
    if not 0 <= args.min_dist:
        raise ValueError("--min-dist must be >= 0")
    for path in (args.h5ad, args.standardized_scores, args.legacy_ppt):
        if not path.exists():
            raise FileNotFoundError(path)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    if args.cluster_key not in adata.obs:
        raise KeyError(f"Missing cluster key: {args.cluster_key}")
    if "highly_variable" not in adata.var:
        raise KeyError("Missing highly_variable mask")
    n_hvg = int(adata.var["highly_variable"].sum())
    if max(pcs) >= min(adata.n_obs, n_hvg):
        raise ValueError(f"Cannot compute {max(pcs)} PCs from n_obs={adata.n_obs}, n_hvg={n_hvg}")
    configure_inkt_palettes(adata)
    palette_failures = palette_validation_failures(adata)
    if palette_failures:
        raise RuntimeError("Palette validation failed: " + "; ".join(palette_failures))

    scores = pd.read_csv(args.standardized_scores).set_index("cell_id")
    missing_cells = adata.obs_names.difference(scores.index)
    extra_cells = scores.index.difference(adata.obs_names)
    if len(missing_cells) or len(extra_cells):
        raise RuntimeError(
            f"Score/cell mismatch: missing={len(missing_cells)}, extra={len(extra_cells)}"
        )
    scores = scores.reindex(adata.obs_names)
    if scores.loc[:, PROGRAM_COLUMNS].isna().any().any():
        raise RuntimeError("Standardized program scores contain missing values")

    print(f"Computing one nested {max(pcs)}-PC representation from {n_hvg} HVGs")
    sc.tl.pca(
        adata,
        n_comps=max(pcs),
        svd_solver="arpack",
        mask_var="highly_variable",
        random_state=args.random_state,
    )
    embeddings: dict[int, np.ndarray] = {}
    for pc in pcs:
        neighbors_key = f"pc_sweep_neighbors_{pc}"
        print(f"Computing neighbors/UMAP for {pc} PCs")
        sc.pp.neighbors(
            adata,
            n_neighbors=args.n_neighbors,
            n_pcs=pc,
            use_rep="X_pca",
            metric=args.metric,
            random_state=args.random_state,
            key_added=neighbors_key,
        )
        sc.tl.umap(
            adata,
            neighbors_key=neighbors_key,
            min_dist=args.min_dist,
            spread=args.spread,
            random_state=args.random_state,
        )
        embeddings[pc] = np.asarray(adata.obsm["X_umap"], dtype=np.float32).copy()
        neighbor_metadata = dict(adata.uns.pop(neighbors_key))
        for graph_key in (
            neighbor_metadata.get("distances_key"),
            neighbor_metadata.get("connectivities_key"),
        ):
            if graph_key:
                adata.obsp.pop(graph_key, None)

    sample_labels = adata.obs["sample"].astype(str).to_numpy()
    tissue_labels = adata.obs["tissue"].astype(str).to_numpy()
    cluster_labels = adata.obs[args.cluster_key].astype(str).to_numpy()
    score_values = scores.loc[:, PROGRAM_COLUMNS].to_numpy(dtype=float)
    legacy_centroids, legacy_image = legacy_sample_centroids(args.legacy_ppt)
    legacy_reference_path = args.out_dir / "figures/legacy_reference_200pc.png"
    legacy_reference_path.parent.mkdir(parents=True, exist_ok=True)
    legacy_image.save(legacy_reference_path)

    rng = np.random.default_rng(args.random_state)
    trust_indices = np.sort(
        rng.choice(
            adata.n_obs,
            size=min(args.trustworthiness_sample_size, adata.n_obs),
            replace=False,
        )
    )
    metric_rows = []
    for pc, embedding in embeddings.items():
        neighbors = knn_indices(embedding, args.n_neighbors)
        smoothness, per_program = program_knn_smoothness(score_values, neighbors)
        sample_centroids = label_centroids(embedding, sample_labels, SAMPLE_DISPLAY_ORDER)
        row = {
            "n_pcs": pc,
            "umap1_min": float(embedding[:, 0].min()),
            "umap1_max": float(embedding[:, 0].max()),
            "umap2_min": float(embedding[:, 1].min()),
            "umap2_max": float(embedding[:, 1].max()),
            "tissue_silhouette": safe_silhouette(
                embedding,
                tissue_labels,
                sample_size=args.silhouette_sample_size,
                random_state=args.random_state,
            ),
            "leiden_silhouette": safe_silhouette(
                embedding,
                cluster_labels,
                sample_size=args.silhouette_sample_size,
                random_state=args.random_state,
            ),
            "tissue_knn_purity": label_knn_purity(tissue_labels, neighbors),
            "leiden_knn_purity": label_knn_purity(cluster_labels, neighbors),
            "tissue_centroid_separation": centroid_separation_ratio(embedding, tissue_labels),
            "leiden_centroid_separation": centroid_separation_ratio(embedding, cluster_labels),
            "program_knn_smoothness": smoothness,
            "inkt1_knn_smoothness": per_program[0],
            "inkt2_knn_smoothness": per_program[1],
            "inkt17_knn_smoothness": per_program[2],
            "trustworthiness": float(
                trustworthiness(
                    np.asarray(adata.obsm["X_pca"])[trust_indices, :pc],
                    embedding[trust_indices],
                    n_neighbors=args.n_neighbors,
                    metric=args.metric,
                )
            ),
            "legacy_centroid_concordance": procrustes_concordance(
                legacy_centroids, sample_centroids
            ),
            "legacy_distance_concordance": distance_concordance(
                legacy_centroids, sample_centroids
            ),
        }
        metric_rows.append(row)

    pairwise_stability: dict[tuple[int, int], float] = {}
    for left_index, left_pc in enumerate(pcs):
        for right_pc in pcs[left_index + 1 :]:
            pairwise_stability[(left_pc, right_pc)] = procrustes_concordance(
                embeddings[left_pc], embeddings[right_pc]
            )
    metrics = pd.DataFrame(metric_rows)
    metrics["mean_stability"] = metrics["n_pcs"].map(
        lambda pc: float(
            np.mean(
                [
                    value
                    for pair, value in pairwise_stability.items()
                    if int(pc) in pair
                ]
            )
        )
    )
    selected_pc, selection_weights = select_recommended_pc(metrics)
    metrics = metrics.sort_values("n_pcs").reset_index(drop=True)

    tables_dir = args.out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = tables_dir / "umap_pc_sweep_metrics.csv"
    metrics.to_csv(metrics_path, index=False)
    pd.DataFrame(
        [
            {"left_pc": left, "right_pc": right, "procrustes_concordance": value}
            for (left, right), value in pairwise_stability.items()
        ]
    ).to_csv(tables_dir / "umap_pc_pairwise_stability.csv", index=False)

    figure_outputs = render_candidate_figures(
        args.out_dir,
        embeddings,
        sample_labels,
        cluster_labels,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
    )
    metrics_figure = args.out_dir / "figures/umap_pc_sweep_metrics.png"
    render_metrics_figure(metrics_figure, metrics, selected_pc)
    program_figure = args.out_dir / "figures/umap_selected_program_scores.png"
    render_program_scores(program_figure, embeddings[selected_pc], scores, selected_pc)

    np.savez_compressed(
        args.out_dir / "umap_pc_sweep_embeddings.npz",
        cell_ids=np.asarray(adata.obs_names, dtype=str),
        **{f"X_umap_pc{pc}": embedding for pc, embedding in embeddings.items()},
    )

    print(f"Installing selected {selected_pc}-PC graph/embedding")
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=selected_pc,
        use_rep="X_pca",
        metric=args.metric,
        random_state=args.random_state,
    )
    adata.obsm["X_umap"] = embeddings[selected_pc]
    trajectory_summary: dict[str, Any] = {"status": "not_available"}
    # Diffusion components and DPT depend on the neighbor graph.  Recompute
    # them after installing the selected graph so trajectory-colored UMAPs do
    # not silently retain values from the former 30-PC graph.
    if "iroot" in adata.uns:
        root_index = int(adata.uns["iroot"])
        if not 0 <= root_index < adata.n_obs:
            raise ValueError(f"Invalid audited DPT root index: {root_index}")
        sc.tl.diffmap(adata)
        sc.tl.dpt(
            adata,
            n_dcs=min(10, adata.obsm["X_diffmap"].shape[1]),
        )
        cluster_values = adata.obs[args.cluster_key].astype(str)
        dpt_values = pd.to_numeric(adata.obs["dpt_pseudotime"], errors="raise")
        cluster_medians = {
            str(cluster): float(dpt_values.loc[cluster_values.eq(str(cluster))].median())
            for cluster in sorted(cluster_values.unique(), key=int)
        }
        root_cluster = str(cluster_values.iloc[root_index])
        root_cluster_mask = cluster_values.eq(root_cluster)
        root_cluster_tissues = adata.obs.loc[root_cluster_mask, "tissue"].astype(str)
        dominant_tissue = str(root_cluster_tissues.value_counts().index[0])
        trajectory_summary = {
            "status": "recomputed_from_selected_graph",
            "root_index": root_index,
            "root_cell": str(adata.obs_names[root_index]),
            "root_cell_tissue": str(adata.obs["tissue"].iloc[root_index]),
            "root_cluster": root_cluster,
            "root_cluster_n_cells": int(root_cluster_mask.sum()),
            "root_cluster_dominant_tissue": dominant_tissue,
            "root_cluster_dominant_tissue_fraction": float(
                root_cluster_tissues.eq(dominant_tissue).mean()
            ),
            "cluster_median_dpt": cluster_medians,
        }
    sc.tl.paga(adata, groups=args.cluster_key)
    adata.uns["umap_pc_sweep"] = {
        "candidate_pcs": list(pcs),
        "selected_pc": selected_pc,
        "n_neighbors": args.n_neighbors,
        "min_dist": args.min_dist,
        "spread": args.spread,
        "metric": args.metric,
        "random_state": args.random_state,
        "cluster_label_policy": f"fixed {args.cluster_key} labels from audited legacy-QC run",
    }
    selected_h5ad = args.out_dir / "inkt_selected_umap.h5ad"
    adata.write(selected_h5ad, compression="gzip")

    selected_figures_dir = args.out_dir / "selected_figures"
    render_scanpy_selected_figures(
        adata,
        selected_figures_dir,
        load_marker_genes(args.h5ad),
        args.cluster_key,
    )

    selected_row = metrics.loc[metrics["n_pcs"].eq(selected_pc)].iloc[0]
    best_metrics = {
        column: int(metrics.loc[metrics[column].idxmax(), "n_pcs"])
        for column in selection_weights
    }
    summary = {
        "analysis": "iNKT UMAP PC sensitivity",
        "input_h5ad": str(args.h5ad.relative_to(ROOT)),
        "input_h5ad_sha256": sha256_file(args.h5ad),
        "standardized_scores": str(args.standardized_scores.relative_to(ROOT)),
        "legacy_reference": {
            "pptx": str(args.legacy_ppt.relative_to(ROOT)),
            "slide": 4,
            "media": "ppt/media/image8.png",
            "reference_pc": 200,
            "comparison": "six-sample centroid geometry after Procrustes alignment",
        },
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_highly_variable_genes": n_hvg,
        "candidate_pcs": list(pcs),
        "selected_pc": selected_pc,
        "selection_score": float(selected_row["selection_score"]),
        "selection_weights": selection_weights,
        "metric_winners": best_metrics,
        "parameters": {
            "n_neighbors": args.n_neighbors,
            "min_dist": args.min_dist,
            "spread": args.spread,
            "metric": args.metric,
            "random_state": args.random_state,
            "pca_solver": "arpack",
            "pca_mask": "3,000 highly_variable genes",
        },
        "cluster_boundary_policy": (
            f"Fixed audited {args.cluster_key} labels; 80% smoothed density envelopes in each UMAP. "
            "Clusters were not recomputed, so downstream cluster evidence remains comparable."
        ),
        "subtype_policy": (
            "No validated discrete iNKT1/iNKT2/NKT17 cell labels exist. Continuous standardized "
            "program scores are shown and assessed by local kNN smoothness; they are not forced calls."
        ),
        "trajectory_policy": (
            "PAGA, diffusion components, and DPT were recomputed from the selected neighbor graph; "
            "the audited root cell was retained."
        ),
        "trajectory": trajectory_summary,
        "pairwise_stability": {
            f"{left}_vs_{right}": value
            for (left, right), value in pairwise_stability.items()
        },
        "artifacts": {
            "metrics_csv": str(metrics_path.relative_to(ROOT)),
            "metrics_figure": str(metrics_figure.relative_to(ROOT)),
            "program_figure": str(program_figure.relative_to(ROOT)),
            "comparison_figure": str(Path(figure_outputs["comparison"]).relative_to(ROOT)),
            "legacy_reference_figure": str(legacy_reference_path.relative_to(ROOT)),
            "selected_h5ad": str(selected_h5ad.relative_to(ROOT)),
            "selected_h5ad_sha256": sha256_file(selected_h5ad),
            "selected_figures_dir": str(selected_figures_dir.relative_to(ROOT)),
        },
    }
    summary_path = args.out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
