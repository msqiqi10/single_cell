#!/usr/bin/env python3
"""Plot literature-linked iNKT markers on the selected 50-PC UMAP.

The focused audit intentionally distinguishes a feature being present in the
10x annotation from it having measurable counts.  In particular, mouse Ncam1
(the requested CD56 ortholog) is present in every feature table but has no
non-zero entries in any of the six raw matrices.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy import sparse

from inkt_palette import SAMPLE_DISPLAY_ORDER, SAMPLE_LOAD_ORDER, SAMPLE_PALETTE


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_H5AD = (
    ROOT / "output/iNKT_reproduction_deck/umap_pc_sweep/inkt_selected_umap.h5ad"
)
DEFAULT_RAW_ROOT = ROOT / "input/iNKT/data"
DEFAULT_OUT_DIR = (
    ROOT / "output/iNKT_reproduction_deck/20260825_reference_marker_UMAP"
)
CLUSTER_KEY = "leiden_res_0_5"


@dataclass(frozen=True)
class Marker:
    symbol: str
    alias: str
    group: str
    source: str
    note: str = ""

    @property
    def label(self) -> str:
        return f"{self.alias} ({self.symbol})" if self.alias else self.symbol


FOCUSED_MARKERS = (
    Marker(
        "Ncam1",
        "CD56",
        "requested",
        "User-requested protein alias; mouse ortholog Ncam1. Not present in the two local marker workbooks.",
        "A feature row exists in all six raw 10x inputs, but all six matrices have zero Ncam1 entries.",
    ),
    Marker(
        "Klrd1",
        "CD94",
        "iNKT1 / NK-like",
        "iNKT gene list for scRNAseq.xlsx: iNKT1!A4; iNKT1 (Wang et al., 2022)!A2; original PPT slide 6.",
    ),
    Marker(
        "Il4",
        "IL-4",
        "iNKT2 / activation",
        "marker.xlsx: Sheet1!E6; iNKT gene list: iNKT2!A11/A83, My list iNKT2!A8, Wang iNKT2!A10; PPT slides 6-7.",
        "Il4 also occurs in the local iNKT1 list (iNKT1!A25), so it is not interpreted alone as a subtype call.",
    ),
)

REFERENCE_PANEL = (
    Marker(
        "Xcl1",
        "",
        "iNKT1 / NK-like",
        "iNKT gene list for scRNAseq.xlsx: iNKT1!A1; My list iNKT1!A35; Wang iNKT1!A7; PPT slide 6.",
    ),
    Marker(
        "Klrd1",
        "CD94",
        "iNKT1 / NK-like",
        "iNKT gene list for scRNAseq.xlsx: iNKT1!A4; Wang iNKT1!A2; PPT slide 6.",
    ),
    Marker(
        "Tbx21",
        "T-bet",
        "iNKT1 / NK-like",
        "iNKT gene list for scRNAseq.xlsx: iNKT1!A22; My list iNKT1!A6; PPT slide 6.",
    ),
    Marker(
        "Ifng",
        "IFN-gamma",
        "iNKT1 / NK-like",
        "iNKT gene list for scRNAseq.xlsx: iNKT1!A24; My list iNKT1!A7; Wang iNKT1!A47; PPT slide 6.",
    ),
    Marker(
        "Il4",
        "IL-4",
        "iNKT2 / activation",
        "marker.xlsx: Sheet1!E6; iNKT gene list: iNKT2!A11/A83, My list iNKT2!A8, Wang iNKT2!A10; PPT slides 6-7.",
        "Use with the whole iNKT2 panel rather than as a subtype-exclusive marker.",
    ),
    Marker(
        "Il17rb",
        "",
        "iNKT2 / activation",
        "iNKT gene list for scRNAseq.xlsx: iNKT2!A3, My list iNKT2!A1, Wang iNKT2!A2; PPT slide 7.",
    ),
    Marker(
        "Zbtb16",
        "PLZF",
        "iNKT2 / activation",
        "iNKT gene list for scRNAseq.xlsx: iNKT2!A14, My list iNKT2!A7, Wang iNKT2!A13; PPT slide 7.",
    ),
    Marker(
        "Gata3",
        "",
        "iNKT2 / activation",
        "iNKT gene list for scRNAseq.xlsx: iNKT2!A20; Wang iNKT2!A37; PPT slide 7.",
    ),
    Marker(
        "Rorc",
        "ROR-gamma-t",
        "iNKT17",
        "iNKT gene list for scRNAseq.xlsx: iNKT17!A3, My list iNKT17!A8, Wang iNKT17!A12; PPT slide 8.",
    ),
    Marker(
        "Tmem176a",
        "",
        "iNKT17",
        "iNKT gene list for scRNAseq.xlsx: iNKT17!A9; Wang iNKT17!A1; PPT slide 8.",
    ),
    Marker(
        "Il1r1",
        "",
        "iNKT17",
        "iNKT gene list for scRNAseq.xlsx: iNKT17!A37, My list iNKT17!A2, Wang iNKT17!A8; PPT slide 8.",
    ),
    Marker(
        "Il23r",
        "",
        "iNKT17",
        "iNKT gene list for scRNAseq.xlsx: iNKT17!A36, My list iNKT17!A11, Wang iNKT17!A9; PPT slide 8.",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot CD56/CD94/IL-4 and local literature marker panels on the selected iNKT UMAP."
    )
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_ready(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def vector(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    return np.asarray(values).ravel()


def normalized_gene(adata: ad.AnnData, symbol: str) -> np.ndarray:
    source = adata.raw if adata.raw is not None else adata
    if symbol not in source.var_names:
        raise KeyError(symbol)
    return vector(source[:, symbol].X).astype(float, copy=False)


def counted_gene(adata: ad.AnnData, symbol: str) -> np.ndarray:
    if symbol not in adata.var_names:
        raise KeyError(symbol)
    if "counts" not in adata.layers:
        raise KeyError("The selected object has no counts layer")
    return vector(adata[:, symbol].layers["counts"]).astype(float, copy=False)


def feature_rows(features_path: Path, symbol: str) -> list[tuple[int, str, str]]:
    matches: list[tuple[int, str, str]] = []
    with gzip.open(features_path, "rt") as handle:
        for row_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[1] == symbol and fields[2] == "Gene Expression":
                matches.append((row_number, fields[0], fields[2]))
    return matches


def extract_raw_gene_for_selected_cells(
    adata: ad.AnnData,
    raw_root: Path,
    symbol: str,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Recover a filtered-out raw gene and audit all six 10x matrices."""
    selected_position = {str(name): index for index, name in enumerate(adata.obs_names)}
    selected_counts = np.zeros(adata.n_obs, dtype=float)
    audit_rows: list[dict[str, Any]] = []

    for sample in SAMPLE_LOAD_ORDER:
        matrix_dir = raw_root / sample / "sample_feature_bc_matrix"
        features_path = matrix_dir / "features.tsv.gz"
        matrix_path = matrix_dir / "matrix.mtx.gz"
        barcodes_path = matrix_dir / "barcodes.tsv.gz"
        for required in (features_path, matrix_path, barcodes_path):
            if not required.exists():
                raise FileNotFoundError(required)

        matches = feature_rows(features_path, symbol)
        target_rows = [row for row, _, _ in matches]
        raw_entries: list[tuple[int, int, float]] = []
        if target_rows:
            row_pattern = "^(" + "|".join(str(row) for row in target_rows) + ")[[:space:]]"
            completed = subprocess.run(
                ["zgrep", "-E", row_pattern, str(matrix_path)],
                check=False,
                capture_output=True,
                text=True,
            )
            if completed.returncode not in (0, 1):
                raise RuntimeError(
                    f"zgrep failed for {sample}/{symbol}: {completed.stderr.strip()}"
                )
            for line in completed.stdout.splitlines():
                row_number, column_number, count = line.split()
                raw_entries.append((int(row_number), int(column_number), float(count)))

        selected_nonzero = 0
        selected_total = 0.0
        if raw_entries:
            with gzip.open(barcodes_path, "rt") as handle:
                barcodes = [line.strip() for line in handle]
            for _, column_number, count in raw_entries:
                cell_name = f"{sample}_{barcodes[column_number - 1]}"
                position = selected_position.get(cell_name)
                if position is not None:
                    selected_counts[position] += count
                    selected_nonzero += 1
                    selected_total += count

        audit_rows.append(
            {
                "sample": sample,
                "symbol": symbol,
                "feature_present": bool(matches),
                "feature_rows_1_based": ";".join(str(row) for row, _, _ in matches),
                "feature_ids": ";".join(feature_id for _, feature_id, _ in matches),
                "raw_matrix_nonzero_entries": len(raw_entries),
                "raw_matrix_total_counts": float(sum(entry[2] for entry in raw_entries)),
                "selected_cells_nonzero": int(selected_nonzero),
                "selected_cells_total_counts": float(selected_total),
                "features_path": str(features_path.relative_to(ROOT)),
                "matrix_path": str(matrix_path.relative_to(ROOT)),
            }
        )

    denominator = adata.obs["total_counts"].to_numpy(dtype=float)
    if np.any(denominator <= 0):
        raise RuntimeError("Non-positive retained-gene library size encountered")
    normalized = np.log1p(selected_counts * 10_000.0 / denominator)
    return normalized, selected_counts, pd.DataFrame(audit_rows)


def axis_extent(coordinates: np.ndarray) -> tuple[float, float, float, float]:
    minimum = coordinates.min(axis=0)
    maximum = coordinates.max(axis=0)
    padding = np.maximum((maximum - minimum) * 0.035, 0.1)
    return (
        float(minimum[0] - padding[0]),
        float(maximum[0] + padding[0]),
        float(minimum[1] - padding[1]),
        float(maximum[1] + padding[1]),
    )


def shared_vmax(expression: dict[str, np.ndarray], symbols: list[str]) -> float:
    positive = [expression[symbol][expression[symbol] > 0] for symbol in symbols]
    combined = np.concatenate([values for values in positive if len(values)])
    if not len(combined):
        return 1.0
    return max(float(np.quantile(combined, 0.99)), 1.0)


def style_umap_axis(
    axis: Any,
    extent: tuple[float, float, float, float],
    *,
    show_ylabel: bool,
) -> None:
    axis.set_xlim(extent[0], extent[1])
    axis.set_ylim(extent[2], extent[3])
    axis.set_xlabel("UMAP1", fontsize=8)
    axis.set_ylabel("UMAP2" if show_ylabel else "", fontsize=8)
    axis.tick_params(labelsize=6, length=2, colors="#596674")
    axis.grid(color="#D8DEE5", linewidth=0.35, alpha=0.55, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#94A3B8")
    axis.spines["bottom"].set_color("#94A3B8")


def add_cluster_labels(axis: Any, coordinates: np.ndarray, labels: np.ndarray) -> None:
    categories = sorted(set(labels), key=lambda value: int(value) if value.isdigit() else value)
    for category in categories:
        points = coordinates[labels == category]
        center = np.median(points, axis=0)
        axis.text(
            center[0],
            center[1],
            f"c{category}",
            ha="center",
            va="center",
            fontsize=5.5,
            color="#263746",
            bbox={
                "boxstyle": "round,pad=0.12",
                "facecolor": "white",
                "edgecolor": "#CBD5E1",
                "alpha": 0.78,
                "linewidth": 0.4,
            },
            zorder=5,
        )


def plot_feature(
    axis: Any,
    coordinates: np.ndarray,
    values: np.ndarray,
    marker: Marker,
    extent: tuple[float, float, float, float],
    norm: Normalize,
    *,
    show_ylabel: bool,
    cluster_labels: np.ndarray | None = None,
    unavailable_note: str | None = None,
) -> None:
    axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        s=1.8,
        c="#D8DEE5",
        alpha=0.58,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    positive = values > 0
    if positive.any():
        positions = np.flatnonzero(positive)
        order = positions[np.argsort(values[positions])]
        axis.scatter(
            coordinates[order, 0],
            coordinates[order, 1],
            s=2.4,
            c=values[order],
            cmap="viridis",
            norm=norm,
            alpha=0.92,
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
    if unavailable_note:
        axis.text(
            0.5,
            0.52,
            unavailable_note,
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=8,
            color="#7F1D1D",
            linespacing=1.35,
            bbox={
                "boxstyle": "round,pad=0.45",
                "facecolor": "#FFF7ED",
                "edgecolor": "#FDBA74",
                "alpha": 0.94,
            },
            zorder=6,
        )
    detected = int(positive.sum())
    axis.text(
        0.015,
        0.985,
        f"detected {detected:,}/{len(values):,} ({detected / len(values) * 100:.2f}%)",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=6.2,
        color="#334155",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.4},
        zorder=5,
    )
    axis.set_title(marker.label, fontsize=9.5, fontweight="bold", pad=5)
    style_umap_axis(axis, extent, show_ylabel=show_ylabel)
    if cluster_labels is not None:
        add_cluster_labels(axis, coordinates, cluster_labels)


def render_focused_umap(
    path: Path,
    coordinates: np.ndarray,
    expression: dict[str, np.ndarray],
    cluster_labels: np.ndarray,
    selected_pc: int,
) -> None:
    symbols = [marker.symbol for marker in FOCUSED_MARKERS]
    norm = Normalize(vmin=0.0, vmax=shared_vmax(expression, symbols))
    extent = axis_extent(coordinates)
    figure, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), layout="constrained")
    for index, marker in enumerate(FOCUSED_MARKERS):
        note = None
        if marker.symbol == "Ncam1" and not np.any(expression[marker.symbol] > 0):
            note = (
                "0 non-zero Ncam1 reads\n"
                "across all 6 raw 10x matrices\n"
                "feature present; no distribution to infer"
            )
        plot_feature(
            axes[index],
            coordinates,
            expression[marker.symbol],
            marker,
            extent,
            norm,
            show_ylabel=index == 0,
            cluster_labels=cluster_labels,
            unavailable_note=note,
        )
    colorbar = figure.colorbar(
        ScalarMappable(norm=norm, cmap="viridis"),
        ax=axes,
        fraction=0.025,
        pad=0.015,
    )
    colorbar.set_label("log1p normalized expression (shared scale; capped at pooled 99th percentile)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    figure.suptitle(
        f"Requested marker audit on selected {selected_pc}-PC UMAP",
        fontsize=14,
        fontweight="bold",
    )
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def render_reference_panel(
    path: Path,
    coordinates: np.ndarray,
    expression: dict[str, np.ndarray],
    selected_pc: int,
) -> None:
    symbols = [marker.symbol for marker in REFERENCE_PANEL]
    norm = Normalize(vmin=0.0, vmax=shared_vmax(expression, symbols))
    extent = axis_extent(coordinates)
    figure, axes = plt.subplots(3, 4, figsize=(14.8, 13.0), layout="constrained")
    for index, marker in enumerate(REFERENCE_PANEL):
        row, column = divmod(index, 4)
        plot_feature(
            axes[row, column],
            coordinates,
            expression[marker.symbol],
            marker,
            extent,
            norm,
            show_ylabel=column == 0,
        )
        if column == 0:
            axes[row, column].annotate(
                marker.group,
                xy=(-0.23, 0.5),
                xycoords="axes fraction",
                ha="center",
                va="center",
                rotation=90,
                fontsize=9,
                fontweight="bold",
                color="#334155",
            )
    colorbar = figure.colorbar(
        ScalarMappable(norm=norm, cmap="viridis"),
        ax=axes,
        fraction=0.018,
        pad=0.012,
    )
    colorbar.set_label("log1p normalized expression (shared scale; capped at pooled 99th percentile)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    figure.suptitle(
        f"Local literature marker distributions on selected {selected_pc}-PC UMAP",
        fontsize=15,
        fontweight="bold",
    )
    figure.savefig(path, dpi=210, bbox_inches="tight")
    plt.close(figure)


def summarize_by_group(
    adata: ad.AnnData,
    expression: dict[str, np.ndarray],
    counts: dict[str, np.ndarray],
    markers: list[Marker],
    group_key: str,
) -> pd.DataFrame:
    labels = adata.obs[group_key].astype(str).to_numpy()
    if isinstance(adata.obs[group_key].dtype, pd.CategoricalDtype):
        group_order = [str(value) for value in adata.obs[group_key].cat.categories]
    else:
        group_order = sorted(set(labels))
    rows: list[dict[str, Any]] = []
    for marker_order, marker in enumerate(markers):
        marker_expression = expression[marker.symbol]
        marker_counts = counts[marker.symbol]
        for group_order_index, group in enumerate(group_order):
            mask = labels == group
            detected = marker_counts[mask] > 0
            rows.append(
                {
                    "group_key": group_key,
                    "group": group,
                    "group_order": group_order_index,
                    "marker_order": marker_order,
                    "marker_group": marker.group,
                    "symbol": marker.symbol,
                    "display_label": marker.label,
                    "n_cells": int(mask.sum()),
                    "detected_cells": int(detected.sum()),
                    "detected_pct": float(detected.mean() * 100),
                    "mean_log1p_all": float(marker_expression[mask].mean()),
                    "mean_log1p_detected": (
                        float(marker_expression[mask][detected].mean()) if detected.any() else 0.0
                    ),
                    "total_raw_counts": float(marker_counts[mask].sum()),
                }
            )
    return pd.DataFrame(rows)


def contrast_text(hex_color: str) -> str:
    red, green, blue = (int(hex_color[index : index + 2], 16) for index in (1, 3, 5))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#111827" if luminance > 145 else "white"


def render_dotplot(
    path: Path,
    summary: pd.DataFrame,
    markers: list[Marker],
    group_order: list[str],
    title: str,
    *,
    color_group_labels: bool,
) -> None:
    marker_order = [marker.symbol for marker in markers]
    lookup = summary.set_index(["symbol", "group"])
    x_positions: list[int] = []
    y_positions: list[int] = []
    detected: list[float] = []
    means: list[float] = []
    for y_position, symbol in enumerate(marker_order):
        for x_position, group in enumerate(group_order):
            row = lookup.loc[(symbol, group)]
            x_positions.append(x_position)
            y_positions.append(y_position)
            detected.append(float(row["detected_pct"]))
            means.append(float(row["mean_log1p_all"]))

    detected_array = np.asarray(detected)
    mean_array = np.asarray(means)
    sizes = 14.0 + 270.0 * np.sqrt(np.clip(detected_array, 0, 100) / 100.0)
    vmax = max(float(np.quantile(mean_array, 0.98)), 0.25)
    norm = Normalize(vmin=0.0, vmax=vmax)

    figure, axis = plt.subplots(figsize=(12.2, 8.8), layout="constrained")
    scatter = axis.scatter(
        x_positions,
        y_positions,
        s=sizes,
        c=mean_array,
        cmap="viridis",
        norm=norm,
        edgecolors="#475569",
        linewidths=0.35,
    )
    axis.set_xticks(range(len(group_order)), labels=group_order, rotation=28, ha="right")
    labels = [
        "CD56 (Ncam1) — 0 raw reads" if marker.symbol == "Ncam1" else marker.label
        for marker in markers
    ]
    axis.set_yticks(range(len(markers)), labels=labels)
    axis.invert_yaxis()
    axis.set_xlim(-0.65, len(group_order) - 0.35)
    axis.set_ylim(len(markers) - 0.35, -0.65)
    axis.grid(color="#CBD5E1", linewidth=0.5, alpha=0.65)
    axis.set_axisbelow(True)
    axis.set_title(title, fontsize=14, fontweight="bold", pad=12)
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="both", labelsize=8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if color_group_labels:
        for tick, group in zip(axis.get_xticklabels(), group_order, strict=True):
            color = SAMPLE_PALETTE[group]
            tick.set_color(contrast_text(color))
            tick.set_bbox(
                {"facecolor": color, "edgecolor": "#64748B", "linewidth": 0.5, "pad": 2.0}
            )

    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.035, pad=0.025)
    colorbar.set_label("mean log1p normalized expression (all cells)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

    legend_percentages = (0, 10, 25, 50, 90)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=np.sqrt(14.0 + 270.0 * np.sqrt(value / 100.0)) / 1.25,
            markerfacecolor="#CBD5E1",
            markeredgecolor="#475569",
            markeredgewidth=0.4,
            label=f"{value}%",
        )
        for value in legend_percentages
    ]
    axis.legend(
        handles=handles,
        title="detected cells",
        loc="upper left",
        bbox_to_anchor=(1.13, 1.0),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def marker_audit_table(
    adata: ad.AnnData,
    expression: dict[str, np.ndarray],
    counts: dict[str, np.ndarray],
    markers: list[Marker],
    ncam1_raw_audit: pd.DataFrame,
) -> pd.DataFrame:
    h5ad_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    rows: list[dict[str, Any]] = []
    for marker in markers:
        marker_counts = counts[marker.symbol]
        positive = marker_counts > 0
        raw_feature_samples = np.nan
        raw_matrix_nonzero_entries = np.nan
        if marker.symbol == "Ncam1":
            raw_feature_samples = int(ncam1_raw_audit["feature_present"].sum())
            raw_matrix_nonzero_entries = int(
                ncam1_raw_audit["raw_matrix_nonzero_entries"].sum()
            )
        rows.append(
            {
                "requested_focus": marker.symbol in {"Ncam1", "Klrd1", "Il4"},
                "group": marker.group,
                "alias": marker.alias,
                "symbol": marker.symbol,
                "present_in_selected_h5ad": marker.symbol in h5ad_names,
                "raw_10x_feature_samples_of_6": raw_feature_samples,
                "raw_10x_matrix_nonzero_entries": raw_matrix_nonzero_entries,
                "detected_cells": int(positive.sum()),
                "detected_pct": float(positive.mean() * 100),
                "total_raw_counts_selected_cells": float(marker_counts.sum()),
                "min_positive_raw_count": float(marker_counts[positive].min()) if positive.any() else 0.0,
                "max_raw_count": float(marker_counts.max()),
                "mean_log1p_all": float(expression[marker.symbol].mean()),
                "max_log1p": float(expression[marker.symbol].max()),
                "highly_variable": (
                    bool(adata.var.loc[marker.symbol, "highly_variable"])
                    if marker.symbol in adata.var_names and "highly_variable" in adata.var
                    else False
                ),
                "local_source": marker.source,
                "interpretation_note": marker.note,
            }
        )
    return pd.DataFrame(rows)


def write_readme(
    path: Path,
    audit: pd.DataFrame,
    selected_pc: int,
    h5ad_path: Path,
) -> None:
    klrd1 = audit.set_index("symbol").loc["Klrd1"]
    il4 = audit.set_index("symbol").loc["Il4"]
    text = f"""# 20260825 reference-marker UMAP validation

## Input decision

- CD94 (`Klrd1`) and IL-4 (`Il4`) are present in the selected h5ad and can be plotted directly.
- CD56 was mapped to the mouse ortholog `Ncam1`. `Ncam1` is present in all six raw 10x feature tables, but has **zero non-zero matrix entries in all six samples**. It was consequently absent after the legacy `min_cells=100` gene filter. There is no CD56 expression distribution to infer from this experiment; no substitute gene was relabeled as CD56.
- The remaining 12-gene compact panel is present in the selected object. No additional local input was needed to generate these plots.

## Focused result

- `Klrd1` / CD94: {int(klrd1.detected_cells):,} cells ({klrd1.detected_pct:.2f}%) detected.
- `Il4` / IL-4: {int(il4.detected_cells):,} cells ({il4.detected_pct:.2f}%) detected.
- `Ncam1` / CD56: 0 reads; the focused figure shows an explicit unavailable-data annotation rather than a biological feature distribution.

## Files

- `figures/focused_cd56_cd94_il4_umap.png`: requested three-marker audit on the selected {selected_pc}-PC UMAP; fixed Leiden cluster centroids are labeled.
- `figures/reference_marker_panel_umap.png`: 12 local literature markers grouped as iNKT1/NK-like, iNKT2/activation, and iNKT17.
- `figures/reference_marker_dotplot_by_sample.png`: detection fraction and mean expression across six condition-tissue samples. Sample label swatches use the established tissue-hue / condition-shade palette.
- `figures/reference_marker_dotplot_by_cluster.png`: the same markers summarized by fixed `leiden_res_0_5` clusters.
- `tables/marker_gene_audit.csv`: gene availability, detection, expression, and exact local source cells/slides.
- `tables/ncam1_raw_10x_audit.csv`: six-sample feature and raw-matrix audit for CD56/Ncam1.
- `tables/marker_expression_by_sample.csv` and `tables/marker_expression_by_cluster.csv`: tidy numerical summaries behind the dotplots.
- `input_manifest.json` and `summary.json`: hashes, parameters, and machine-readable conclusions.

## Method

- Coordinates: `X_umap` from `{h5ad_path.relative_to(ROOT)}`, selected with {selected_pc} PCs; UMAP was not rerun for this step.
- Expression colors: log1p of library-size-normalized expression (`raw.X`); panels share a scale capped at the pooled 99th percentile so genes remain comparable.
- Detection fractions: integer `layers["counts"] > 0`.
- Zero-expression cells are light gray; positive cells are drawn last so sparse markers remain visible.
- Dot size is detection percentage; dot color is mean log1p expression across all cells in the group.

## Interpretation cautions

- `Il4` appears in both the local iNKT1 and iNKT2 lists. Its UMAP alone is not a discrete iNKT2 classifier; interpret it with `Il17rb`, `Zbtb16`, and `Gata3`.
- Sparse iNKT17 markers (`Rorc`, `Il1r1`, `Il23r`) are best read together with the dotplot rather than from isolated colored points.
- These are expression overlays on fixed clusters, not newly inferred subtype labels.
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    h5ad_path = args.h5ad.resolve()
    raw_root = args.raw_root.resolve()
    out_dir = args.out_dir.resolve()
    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)
    if not raw_root.exists():
        raise FileNotFoundError(raw_root)

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)
    required_obs = {"sample", "condition", "tissue", "total_counts", CLUSTER_KEY}
    missing_obs = sorted(required_obs - set(adata.obs.columns))
    if missing_obs:
        raise KeyError(f"Missing required obs columns: {missing_obs}")
    if "X_umap" not in adata.obsm:
        raise KeyError("Selected object has no X_umap")
    if adata.obsm["X_umap"].shape != (adata.n_obs, 2):
        raise ValueError(f"Unexpected X_umap shape: {adata.obsm['X_umap'].shape}")

    selection_metadata = json_ready(dict(adata.uns.get("umap_pc_sweep", {})))
    selected_pc = int(selection_metadata.get("selected_pc", -1))
    if selected_pc <= 0:
        raise ValueError(f"Missing selected PC metadata: {selection_metadata}")

    panel_missing = [
        marker.symbol
        for marker in REFERENCE_PANEL
        if marker.symbol not in (adata.raw.var_names if adata.raw is not None else adata.var_names)
    ]
    if panel_missing:
        raise RuntimeError(f"Reference-panel genes missing from selected object: {panel_missing}")

    ncam1_expression, ncam1_counts, ncam1_raw_audit = extract_raw_gene_for_selected_cells(
        adata, raw_root, "Ncam1"
    )
    ncam1_raw_audit.to_csv(tables_dir / "ncam1_raw_10x_audit.csv", index=False)

    unique_markers: list[Marker] = []
    seen: set[str] = set()
    for marker in (*FOCUSED_MARKERS, *REFERENCE_PANEL):
        if marker.symbol not in seen:
            unique_markers.append(marker)
            seen.add(marker.symbol)

    expression: dict[str, np.ndarray] = {"Ncam1": ncam1_expression}
    counts: dict[str, np.ndarray] = {"Ncam1": ncam1_counts}
    for marker in unique_markers:
        if marker.symbol == "Ncam1":
            continue
        expression[marker.symbol] = normalized_gene(adata, marker.symbol)
        counts[marker.symbol] = counted_gene(adata, marker.symbol)

    if not bool(ncam1_raw_audit["feature_present"].all()):
        raise RuntimeError("Ncam1 is not present in all six raw feature tables")
    if int(ncam1_raw_audit["raw_matrix_nonzero_entries"].sum()) != 0:
        raise RuntimeError(
            "Ncam1 unexpectedly has non-zero raw entries; review the rescue normalization before reporting"
        )

    coordinates = np.asarray(adata.obsm["X_umap"], dtype=float)
    cluster_labels = adata.obs[CLUSTER_KEY].astype(str).to_numpy()
    focused_path = figures_dir / "focused_cd56_cd94_il4_umap.png"
    panel_path = figures_dir / "reference_marker_panel_umap.png"
    render_focused_umap(
        focused_path,
        coordinates,
        expression,
        cluster_labels,
        selected_pc,
    )
    render_reference_panel(panel_path, coordinates, expression, selected_pc)

    dotplot_markers = [FOCUSED_MARKERS[0], *REFERENCE_PANEL]
    sample_summary = summarize_by_group(
        adata, expression, counts, dotplot_markers, "sample"
    )
    cluster_summary = summarize_by_group(
        adata, expression, counts, dotplot_markers, CLUSTER_KEY
    )
    sample_summary.to_csv(tables_dir / "marker_expression_by_sample.csv", index=False)
    cluster_summary.to_csv(tables_dir / "marker_expression_by_cluster.csv", index=False)

    sample_path = figures_dir / "reference_marker_dotplot_by_sample.png"
    cluster_path = figures_dir / "reference_marker_dotplot_by_cluster.png"
    render_dotplot(
        sample_path,
        sample_summary,
        dotplot_markers,
        list(SAMPLE_DISPLAY_ORDER),
        "Reference markers by condition-tissue sample",
        color_group_labels=True,
    )
    cluster_order = [
        str(value)
        for value in (
            adata.obs[CLUSTER_KEY].cat.categories
            if isinstance(adata.obs[CLUSTER_KEY].dtype, pd.CategoricalDtype)
            else sorted(set(cluster_labels), key=lambda value: int(value))
        )
    ]
    render_dotplot(
        cluster_path,
        cluster_summary,
        dotplot_markers,
        cluster_order,
        f"Reference markers by fixed {CLUSTER_KEY} cluster",
        color_group_labels=False,
    )

    audit = marker_audit_table(
        adata, expression, counts, unique_markers, ncam1_raw_audit
    )
    audit_path = tables_dir / "marker_gene_audit.csv"
    audit.to_csv(audit_path, index=False)

    workbook_path = ROOT / "input/iNKT/iNKT gene list for scRNAseq.xlsx"
    marker_path = ROOT / "input/iNKT/marker.xlsx"
    pptx_path = ROOT / "input/iNKT/iNKT.pptx"
    script_path = Path(__file__).resolve()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "selected_h5ad": {
                "path": str(h5ad_path.relative_to(ROOT)),
                "sha256": sha256_file(h5ad_path),
                "size_bytes": h5ad_path.stat().st_size,
            },
            "local_gene_list": {
                "path": str(workbook_path.relative_to(ROOT)),
                "sha256": sha256_file(workbook_path),
            },
            "local_marker_table": {
                "path": str(marker_path.relative_to(ROOT)),
                "sha256": sha256_file(marker_path),
            },
            "original_pptx": {
                "path": str(pptx_path.relative_to(ROOT)),
                "sha256": sha256_file(pptx_path),
            },
            "raw_10x_root": str(raw_root.relative_to(ROOT)),
            "script": {
                "path": str(script_path.relative_to(ROOT)),
                "sha256": sha256_file(script_path),
            },
        },
        "selected_umap": selection_metadata,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
    }
    (out_dir / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    audit_lookup = audit.set_index("symbol")
    artifacts = {
        "focused_umap": str(focused_path.relative_to(ROOT)),
        "reference_panel_umap": str(panel_path.relative_to(ROOT)),
        "sample_dotplot": str(sample_path.relative_to(ROOT)),
        "cluster_dotplot": str(cluster_path.relative_to(ROOT)),
        "marker_audit": str(audit_path.relative_to(ROOT)),
        "ncam1_raw_audit": str((tables_dir / "ncam1_raw_10x_audit.csv").relative_to(ROOT)),
        "sample_summary": str((tables_dir / "marker_expression_by_sample.csv").relative_to(ROOT)),
        "cluster_summary": str((tables_dir / "marker_expression_by_cluster.csv").relative_to(ROOT)),
    }
    summary = {
        "analysis": "iNKT local-reference marker UMAP validation",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_pc": selected_pc,
        "cluster_key": CLUSTER_KEY,
        "input_complete": {
            "CD94_Klrd1": True,
            "IL4_Il4": True,
            "CD56_Ncam1_meaningful_distribution": False,
            "other_reference_panel_genes": True,
        },
        "focused_detection": {
            symbol: {
                "detected_cells": int(audit_lookup.loc[symbol, "detected_cells"]),
                "detected_pct": float(audit_lookup.loc[symbol, "detected_pct"]),
            }
            for symbol in ("Ncam1", "Klrd1", "Il4")
        },
        "ncam1_conclusion": (
            "Ncam1 is annotated as a Gene Expression feature in all six raw 10x inputs, "
            "but every matrix has zero non-zero Ncam1 entries; no CD56 distribution can be inferred."
        ),
        "reference_panel": [marker.symbol for marker in REFERENCE_PANEL],
        "expression_color": "raw.X log1p normalized expression; shared scale capped at pooled positive-cell 99th percentile",
        "detection_definition": "layers['counts'] > 0",
        "artifacts": artifacts,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_readme(out_dir / "README.md", audit, selected_pc, h5ad_path)

    expected_files = [
        focused_path,
        panel_path,
        sample_path,
        cluster_path,
        audit_path,
        tables_dir / "ncam1_raw_10x_audit.csv",
        tables_dir / "marker_expression_by_sample.csv",
        tables_dir / "marker_expression_by_cluster.csv",
        out_dir / "input_manifest.json",
        out_dir / "summary.json",
        out_dir / "README.md",
    ]
    empty = [str(path) for path in expected_files if not path.exists() or path.stat().st_size == 0]
    if empty:
        raise RuntimeError(f"Missing or empty output files: {empty}")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
