#!/usr/bin/env python3
"""Plot Ctrl versus tumor cluster proportions, stratified by tissue.

Two complementary denominators are exported because they answer different
questions:

1. P(condition | tissue, cluster): literal Ctrl/T2 cell share inside a cluster.
2. P(cluster | tissue, condition): cluster frequency inside a tissue-condition
   sample, which is the safer descriptive comparison when cell yields differ.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from inkt_palette import (
    CONDITION_ORDER,
    SAMPLE_DISPLAY_ORDER,
    SAMPLE_PALETTE,
    TISSUE_ORDER,
    TISSUE_PALETTE,
    palette_validation_failures,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_H5AD = (
    ROOT / "output/iNKT_reproduction_deck/umap_pc_sweep/inkt_selected_umap.h5ad"
)
DEFAULT_OUT_DIR = (
    ROOT
    / "output/iNKT_reproduction_deck/20260825_control_vs_tumor_cluster_proportions_by_tissue"
)
CLUSTER_KEY = "leiden_res_0_5"
CONDITION_LABELS = {"Ctrl": "Control", "T2": "Tumor (T2)"}
TISSUE_LABELS = {
    "thymus": "Thymus",
    "bone_marrow": "Bone marrow",
    "spleen": "Spleen",
}
SAMPLE_BY_TISSUE_CONDITION = {
    ("thymus", "Ctrl"): "Ctrl_Thymus",
    ("thymus", "T2"): "T2_Thymus",
    ("bone_marrow", "Ctrl"): "Ctrl_BM",
    ("bone_marrow", "T2"): "T2_BM",
    ("spleen", "Ctrl"): "Ctrl_Spleen",
    ("spleen", "T2"): "T2_Spleen",
}
DELTA_COLORS = {"Control enriched": "#457B9D", "Tumor enriched": "#B23A48"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Ctrl versus tumor proportions in fixed iNKT clusters, faceted by tissue."
    )
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
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


def natural_cluster_order(series: pd.Series) -> list[str]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        values = [str(value) for value in series.cat.categories]
    else:
        values = sorted({str(value) for value in series})
    return sorted(values, key=lambda value: int(value) if value.isdigit() else value)


def build_proportion_table(
    adata: ad.AnnData,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    required = {"condition", "tissue", "sample", CLUSTER_KEY}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")
    if adata.obs[list(required)].isna().any().any():
        raise ValueError("Missing condition/tissue/sample/cluster metadata")

    clusters = natural_cluster_order(adata.obs[CLUSTER_KEY])
    obs = pd.DataFrame(
        {
            "condition": adata.obs["condition"].astype(str).to_numpy(),
            "tissue": adata.obs["tissue"].astype(str).to_numpy(),
            "sample": adata.obs["sample"].astype(str).to_numpy(),
            "cluster": adata.obs[CLUSTER_KEY].astype(str).to_numpy(),
        },
        index=adata.obs_names,
    )
    unknown_conditions = sorted(set(obs["condition"]) - set(CONDITION_ORDER))
    unknown_tissues = sorted(set(obs["tissue"]) - set(TISSUE_ORDER))
    if unknown_conditions or unknown_tissues:
        raise ValueError(
            f"Unexpected metadata: conditions={unknown_conditions}, tissues={unknown_tissues}"
        )

    observed_samples = set(obs["sample"])
    expected_samples = set(SAMPLE_DISPLAY_ORDER)
    if observed_samples != expected_samples:
        raise ValueError(
            f"Expected six canonical samples; missing={sorted(expected_samples-observed_samples)}, "
            f"extra={sorted(observed_samples-expected_samples)}"
        )

    grouped = obs.groupby(["tissue", "condition", "cluster"], observed=True).size()
    complete_index = pd.MultiIndex.from_product(
        [TISSUE_ORDER, CONDITION_ORDER, clusters],
        names=["tissue", "condition", "cluster"],
    )
    table = grouped.reindex(complete_index, fill_value=0).rename("count").reset_index()
    table["count"] = table["count"].astype(int)
    table["tissue_order"] = table["tissue"].map(
        {value: index for index, value in enumerate(TISSUE_ORDER)}
    )
    table["condition_order"] = table["condition"].map(
        {value: index for index, value in enumerate(CONDITION_ORDER)}
    )
    table["cluster_order"] = table["cluster"].map(
        {value: index for index, value in enumerate(clusters)}
    )
    table["condition_label"] = table["condition"].map(CONDITION_LABELS)
    table["tissue_label"] = table["tissue"].map(TISSUE_LABELS)
    table["sample"] = [
        SAMPLE_BY_TISSUE_CONDITION[(tissue, condition)]
        for tissue, condition in zip(table["tissue"], table["condition"], strict=True)
    ]

    table["total_tissue_condition"] = table.groupby(
        ["tissue", "condition"], observed=True
    )["count"].transform("sum")
    table["pct_cluster_within_tissue_condition"] = (
        table["count"] / table["total_tissue_condition"] * 100.0
    )
    table["total_tissue_cluster"] = table.groupby(
        ["tissue", "cluster"], observed=True
    )["count"].transform("sum")
    table["pct_condition_within_tissue_cluster"] = np.where(
        table["total_tissue_cluster"] > 0,
        table["count"] / table["total_tissue_cluster"] * 100.0,
        np.nan,
    )
    table = table.sort_values(
        ["tissue_order", "condition_order", "cluster_order"]
    ).reset_index(drop=True)

    frequency_wide = table.pivot(
        index=["tissue", "cluster", "tissue_order", "cluster_order"],
        columns="condition",
        values="pct_cluster_within_tissue_condition",
    ).reset_index()
    count_wide = table.pivot(
        index=["tissue", "cluster"], columns="condition", values="count"
    ).reset_index()
    delta = frequency_wide.merge(count_wide, on=["tissue", "cluster"], suffixes=("_pct", "_count"))
    delta = delta.rename(
        columns={
            "Ctrl_pct": "control_pct",
            "T2_pct": "tumor_pct",
            "Ctrl_count": "control_count",
            "T2_count": "tumor_count",
        }
    )
    delta["tumor_minus_control_percentage_points"] = (
        delta["tumor_pct"] - delta["control_pct"]
    )
    delta = delta.sort_values(["tissue_order", "cluster_order"]).reset_index(drop=True)

    within_condition_sums = table.groupby(["tissue", "condition"], observed=True)[
        "pct_cluster_within_tissue_condition"
    ].sum()
    if not np.allclose(within_condition_sums.to_numpy(), 100.0, atol=1e-9):
        raise RuntimeError(f"Within-condition cluster frequencies do not sum to 100: {within_condition_sums}")
    nonempty = table[table["total_tissue_cluster"] > 0]
    within_cluster_sums = nonempty.groupby(["tissue", "cluster"], observed=True)[
        "pct_condition_within_tissue_cluster"
    ].sum()
    if not np.allclose(within_cluster_sums.to_numpy(), 100.0, atol=1e-9):
        raise RuntimeError(f"Within-cluster condition shares do not sum to 100: {within_cluster_sums}")
    if int(table["count"].sum()) != adata.n_obs:
        raise RuntimeError(f"Count table total {table['count'].sum()} != n_obs {adata.n_obs}")
    return table, delta, clusters


def contrast_text(hex_color: str) -> str:
    red, green, blue = (int(hex_color[index : index + 2], 16) for index in (1, 3, 5))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#111827" if luminance > 145 else "white"


def style_tissue_title(axis: Any, tissue: str) -> None:
    color = TISSUE_PALETTE[tissue]
    axis.set_title(
        TISSUE_LABELS[tissue],
        fontsize=12,
        fontweight="bold",
        color=contrast_text(color),
        pad=10,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": color,
            "edgecolor": "#475569",
            "linewidth": 0.7,
        },
    )


def sample_pair(tissue: str) -> tuple[str, str]:
    return (
        SAMPLE_BY_TISSUE_CONDITION[(tissue, "Ctrl")],
        SAMPLE_BY_TISSUE_CONDITION[(tissue, "T2")],
    )


def save_figure(figure: Any, png_path: Path) -> tuple[Path, Path]:
    pdf_path = png_path.with_suffix(".pdf")
    figure.savefig(png_path, dpi=220, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def render_condition_share(
    table: pd.DataFrame,
    clusters: list[str],
    path: Path,
) -> tuple[Path, Path]:
    figure, axes = plt.subplots(1, 3, figsize=(17.2, 5.8), sharey=True, layout="constrained")
    x = np.arange(len(clusters))
    for axis, tissue in zip(axes, TISSUE_ORDER, strict=True):
        tissue_table = table[table["tissue"] == tissue]
        share = tissue_table.pivot(
            index="cluster", columns="condition", values="pct_condition_within_tissue_cluster"
        ).reindex(clusters)
        totals = (
            tissue_table.drop_duplicates("cluster")
            .set_index("cluster")["total_tissue_cluster"]
            .reindex(clusters)
            .to_numpy(dtype=int)
        )
        control_sample, tumor_sample = sample_pair(tissue)
        control = share["Ctrl"].fillna(0).to_numpy(dtype=float)
        tumor = share["T2"].fillna(0).to_numpy(dtype=float)
        axis.bar(
            x,
            control,
            width=0.72,
            color=SAMPLE_PALETTE[control_sample],
            edgecolor="#1F2937",
            linewidth=0.55,
            label="Control",
        )
        axis.bar(
            x,
            tumor,
            bottom=control,
            width=0.72,
            color=SAMPLE_PALETTE[tumor_sample],
            edgecolor="#1F2937",
            linewidth=0.55,
            label="Tumor (T2)",
        )
        for index, (control_value, tumor_value, total) in enumerate(
            zip(control, tumor, totals, strict=True)
        ):
            if total == 0:
                axis.text(index, 50, "no cells", ha="center", va="center", fontsize=7, color="#64748B", rotation=90)
                axis.text(index, 102.2, "n=0", ha="center", va="bottom", fontsize=6.5, color="#64748B")
                continue
            axis.text(index, 102.2, f"n={total:,}", ha="center", va="bottom", fontsize=6.5, color="#334155")
            if control_value >= 8:
                axis.text(
                    index,
                    control_value / 2,
                    f"{control_value:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color=contrast_text(SAMPLE_PALETTE[control_sample]),
                    fontweight="bold",
                )
            if tumor_value >= 8:
                axis.text(
                    index,
                    control_value + tumor_value / 2,
                    f"{tumor_value:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color=contrast_text(SAMPLE_PALETTE[tumor_sample]),
                    fontweight="bold",
                )
        axis.set_xticks(x, [f"c{cluster}" for cluster in clusters])
        axis.set_xlabel("Fixed Leiden cluster", fontsize=9)
        axis.set_ylim(0, 110)
        axis.set_yticks([0, 25, 50, 75, 100])
        axis.grid(axis="y", color="#CBD5E1", linewidth=0.55, alpha=0.75)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        style_tissue_title(axis, tissue)
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=8)
    axes[0].set_ylabel("Condition share within tissue × cluster (%)", fontsize=10)
    figure.suptitle(
        "Control vs tumor share inside each cluster, stratified by tissue",
        fontsize=15,
        fontweight="bold",
        y=1.10,
    )
    figure.text(
        0.5,
        1.035,
        "Each non-empty cluster bar sums to 100%; n above bars is the total captured cell count",
        ha="center",
        va="top",
        fontsize=9,
        color="#475569",
    )
    return save_figure(figure, path)


def render_cluster_frequency(
    table: pd.DataFrame,
    clusters: list[str],
    path: Path,
) -> tuple[Path, Path]:
    figure, axes = plt.subplots(1, 3, figsize=(17.2, 5.9), sharey=True, layout="constrained")
    x = np.arange(len(clusters))
    width = 0.38
    for axis, tissue in zip(axes, TISSUE_ORDER, strict=True):
        tissue_table = table[table["tissue"] == tissue]
        frequency = tissue_table.pivot(
            index="cluster", columns="condition", values="pct_cluster_within_tissue_condition"
        ).reindex(clusters)
        denominators = (
            tissue_table.drop_duplicates("condition")
            .set_index("condition")["total_tissue_condition"]
            .to_dict()
        )
        control_sample, tumor_sample = sample_pair(tissue)
        control = frequency["Ctrl"].to_numpy(dtype=float)
        tumor = frequency["T2"].to_numpy(dtype=float)
        control_bars = axis.bar(
            x - width / 2,
            control,
            width=width,
            color=SAMPLE_PALETTE[control_sample],
            edgecolor="#1F2937",
            linewidth=0.55,
            label="Control",
        )
        tumor_bars = axis.bar(
            x + width / 2,
            tumor,
            width=width,
            color=SAMPLE_PALETTE[tumor_sample],
            edgecolor="#1F2937",
            linewidth=0.55,
            label="Tumor (T2)",
        )
        for bars, values, color in (
            (control_bars, control, SAMPLE_PALETTE[control_sample]),
            (tumor_bars, tumor, SAMPLE_PALETTE[tumor_sample]),
        ):
            for bar, value in zip(bars, values, strict=True):
                if value < 0.5:
                    continue
                if value >= 20:
                    axis.text(
                        bar.get_x() + bar.get_width() / 2,
                        value - 2.0,
                        f"{value:.1f}",
                        ha="center",
                        va="top",
                        fontsize=6.5,
                        color=contrast_text(color),
                        fontweight="bold",
                    )
                else:
                    axis.text(
                        bar.get_x() + bar.get_width() / 2,
                        value + 0.7,
                        f"{value:.1f}",
                        ha="center",
                        va="bottom",
                        fontsize=6.2,
                        color="#334155",
                        rotation=90,
                    )
        axis.text(
            0.02,
            0.985,
            f"Control n={int(denominators['Ctrl']):,}\nTumor n={int(denominators['T2']):,}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            color="#334155",
            bbox={"facecolor": "white", "edgecolor": "#CBD5E1", "alpha": 0.88, "pad": 2.5},
        )
        axis.set_xticks(x, [f"c{cluster}" for cluster in clusters])
        axis.set_xlabel("Fixed Leiden cluster", fontsize=9)
        axis.set_ylim(0, 100)
        axis.set_yticks([0, 25, 50, 75, 100])
        axis.grid(axis="y", color="#CBD5E1", linewidth=0.55, alpha=0.75)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        style_tissue_title(axis, tissue)
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=8)
    axes[0].set_ylabel("Cells in cluster / all cells in tissue × condition (%)", fontsize=10)
    figure.suptitle(
        "Cluster frequencies: control vs tumor, stratified by tissue",
        fontsize=15,
        fontweight="bold",
        y=1.10,
    )
    figure.text(
        0.5,
        1.035,
        "Within each tissue, each condition series sums to 100%; use this view for composition shifts",
        ha="center",
        va="top",
        fontsize=9,
        color="#475569",
    )
    return save_figure(figure, path)


def render_delta(
    delta: pd.DataFrame,
    clusters: list[str],
    path: Path,
) -> tuple[Path, Path]:
    maximum = max(float(delta["tumor_minus_control_percentage_points"].abs().max()), 1.0)
    limit = np.ceil(maximum * 1.28 * 2) / 2
    figure, axes = plt.subplots(1, 3, figsize=(16.4, 4.9), sharey=True, layout="constrained")
    x = np.arange(len(clusters))
    for axis, tissue in zip(axes, TISSUE_ORDER, strict=True):
        values = (
            delta[delta["tissue"] == tissue]
            .set_index("cluster")
            .reindex(clusters)["tumor_minus_control_percentage_points"]
            .to_numpy(dtype=float)
        )
        colors = [
            DELTA_COLORS["Tumor enriched"] if value >= 0 else DELTA_COLORS["Control enriched"]
            for value in values
        ]
        bars = axis.bar(x, values, width=0.68, color=colors, edgecolor="#1F2937", linewidth=0.5)
        axis.axhline(0, color="#111827", linewidth=0.8)
        for bar, value in zip(bars, values, strict=True):
            if abs(value) < 0.2:
                continue
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                value + (0.16 if value >= 0 else -0.16),
                f"{value:+.2f}",
                ha="center",
                va="bottom" if value >= 0 else "top",
                fontsize=6.8,
                color="#334155",
                rotation=90,
            )
        axis.set_xticks(x, [f"c{cluster}" for cluster in clusters])
        axis.set_xlabel("Fixed Leiden cluster", fontsize=9)
        axis.set_ylim(-limit, limit)
        axis.grid(axis="y", color="#CBD5E1", linewidth=0.55, alpha=0.75)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        style_tissue_title(axis, tissue)
    axes[0].set_ylabel("Tumor − Control cluster frequency (percentage points)", fontsize=10)
    figure.suptitle(
        "Descriptive cluster-frequency shifts within each tissue",
        fontsize=15,
        fontweight="bold",
        y=1.10,
    )
    figure.legend(
        handles=[
            Patch(facecolor=DELTA_COLORS["Control enriched"], edgecolor="#1F2937", label="Control enriched"),
            Patch(facecolor=DELTA_COLORS["Tumor enriched"], edgecolor="#1F2937", label="Tumor enriched"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.035),
        ncol=2,
        frameon=False,
        fontsize=8,
    )
    figure.text(
        0.5,
        1.035,
        "Descriptive only: one sample per condition × tissue; no biological-replicate significance test",
        ha="center",
        va="top",
        fontsize=9,
        color="#7F1D1D",
    )
    return save_figure(figure, path)


def result_summary(
    table: pd.DataFrame,
    delta: pd.DataFrame,
) -> tuple[dict[str, Any], dict[str, Any]]:
    dominant: dict[str, Any] = {}
    maximum_shifts: dict[str, Any] = {}
    for tissue in TISSUE_ORDER:
        dominant[tissue] = {}
        for condition in CONDITION_ORDER:
            subset = table[(table["tissue"] == tissue) & (table["condition"] == condition)]
            row = subset.loc[subset["pct_cluster_within_tissue_condition"].idxmax()]
            dominant[tissue][condition] = {
                "cluster": str(row["cluster"]),
                "count": int(row["count"]),
                "pct": float(row["pct_cluster_within_tissue_condition"]),
            }
        tissue_delta = delta[delta["tissue"] == tissue]
        row = tissue_delta.loc[
            tissue_delta["tumor_minus_control_percentage_points"].abs().idxmax()
        ]
        maximum_shifts[tissue] = {
            "cluster": str(row["cluster"]),
            "control_pct": float(row["control_pct"]),
            "tumor_pct": float(row["tumor_pct"]),
            "tumor_minus_control_percentage_points": float(
                row["tumor_minus_control_percentage_points"]
            ),
        }
    return dominant, maximum_shifts


def write_readme(
    path: Path,
    h5ad_path: Path,
    selected_pc: int,
    dominant: dict[str, Any],
    maximum_shifts: dict[str, Any],
) -> None:
    text = f"""# 20260825 Control vs Tumor cluster proportions by tissue

## 结果文件

- `figures/ctrl_vs_tumor_share_within_cluster_by_tissue.png`：最贴近字面要求。每个 tissue×cluster 柱内 Control 与 Tumor 合计 100%，柱顶标注该 cluster 的细胞总数。
- `figures/cluster_frequency_ctrl_vs_tumor_by_tissue.png`：推荐用于比较 composition。每个 tissue×condition 内，各 cluster 比例合计 100%。
- `figures/cluster_frequency_delta_tumor_minus_control_by_tissue.png`：Tumor − Control 的百分点差值。
- 三张图均提供同名 PDF 矢量版本。
- `tables/cluster_counts_and_proportions.csv`：两种分母、原始计数及样本信息。
- `tables/cluster_frequency_delta_tumor_minus_control.csv`：差值图数值。
- `summary.json`、`input_manifest.json`：机器可读结论、参数与输入哈希。

## 两种分母

1. **cluster 内 condition 构成**：`count / all cells in the same tissue × cluster`。回答“这个 cluster 里 Control/Tumor 各占多少”，但会受两个样本捕获细胞总数不同影响。
2. **condition 内 cluster 频率**：`count / all cells in the same tissue × condition`。回答“Control/Tumor 各自有多少比例落在这个 cluster”，更适合描述 composition shift。

## 主要描述性结果

- Thymus：Control 和 Tumor 都由 c{dominant['thymus']['Ctrl']['cluster']} 主导（{dominant['thymus']['Ctrl']['pct']:.2f}% vs {dominant['thymus']['T2']['pct']:.2f}%）。最大变化是 c{maximum_shifts['thymus']['cluster']}：{maximum_shifts['thymus']['tumor_minus_control_percentage_points']:+.2f} percentage points。
- Bone marrow：两组都由 c{dominant['bone_marrow']['Ctrl']['cluster']} 主导（{dominant['bone_marrow']['Ctrl']['pct']:.2f}% vs {dominant['bone_marrow']['T2']['pct']:.2f}%）。最大变化是 c{maximum_shifts['bone_marrow']['cluster']}：{maximum_shifts['bone_marrow']['tumor_minus_control_percentage_points']:+.2f} percentage points。
- Spleen：两组都由 c{dominant['spleen']['Ctrl']['cluster']} 主导（{dominant['spleen']['Ctrl']['pct']:.2f}% vs {dominant['spleen']['T2']['pct']:.2f}%）。最大变化是 c{maximum_shifts['spleen']['cluster']}：{maximum_shifts['spleen']['tumor_minus_control_percentage_points']:+.2f} percentage points。
- 整体模式为 tissue 差异明显大于 Control/Tumor 差异。

## 配色

沿用原 PPT 的逻辑：tissue 决定色相、condition 决定同色相的明暗。

- Ctrl_Thymus `#0000FF`; T2_Thymus `#00008B`
- Ctrl_BM `#FFFF00`; T2_BM `#FFA500`
- Ctrl_Spleen `#FF0000`; T2_Spleen `#8B0000`

## 输入与限制

- 输入：`{h5ad_path.relative_to(ROOT)}`，selected {selected_pc}-PC UMAP 对象。
- Cluster：固定 `{CLUSTER_KEY}`，没有重聚类。
- 每个 condition×tissue 仅有一个样本，condition 与样本效应混杂；细胞数不能当作 biological replicates，因此这里只报告描述性比例，不做显著性检验。
- 这些比例描述 QC 后捕获细胞的 composition，不代表组织中的绝对细胞丰度。
- 小 cluster 的百分比容易不稳定，需同时查看柱顶 `n` 和 CSV 原始计数。
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    h5ad_path = args.h5ad.resolve()
    out_dir = args.out_dir.resolve()
    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)
    palette_failures = palette_validation_failures(adata)
    if palette_failures:
        raise RuntimeError(f"Selected h5ad palette state is invalid: {palette_failures}")
    selection_metadata = json_ready(dict(adata.uns.get("umap_pc_sweep", {})))
    selected_pc = int(selection_metadata.get("selected_pc", -1))
    if selected_pc <= 0:
        raise ValueError(f"Missing selected-PC metadata: {selection_metadata}")

    table, delta, clusters = build_proportion_table(adata)
    table_path = tables_dir / "cluster_counts_and_proportions.csv"
    delta_path = tables_dir / "cluster_frequency_delta_tumor_minus_control.csv"
    table.to_csv(table_path, index=False)
    delta.to_csv(delta_path, index=False)

    condition_share_png = figures_dir / "ctrl_vs_tumor_share_within_cluster_by_tissue.png"
    cluster_frequency_png = figures_dir / "cluster_frequency_ctrl_vs_tumor_by_tissue.png"
    delta_png = figures_dir / "cluster_frequency_delta_tumor_minus_control_by_tissue.png"
    condition_share_files = render_condition_share(table, clusters, condition_share_png)
    cluster_frequency_files = render_cluster_frequency(table, clusters, cluster_frequency_png)
    delta_files = render_delta(delta, clusters, delta_png)

    dominant, maximum_shifts = result_summary(table, delta)
    script_path = Path(__file__).resolve()
    palette_path = script_path.with_name("inkt_palette.py")
    original_pptx = ROOT / "input/iNKT/iNKT.pptx"
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "selected_h5ad": {
                "path": str(h5ad_path.relative_to(ROOT)),
                "sha256": sha256_file(h5ad_path),
                "size_bytes": h5ad_path.stat().st_size,
            },
            "original_pptx_palette_reference": {
                "path": str(original_pptx.relative_to(ROOT)),
                "sha256": sha256_file(original_pptx),
                "slide": 12,
            },
            "palette_code": {
                "path": str(palette_path.relative_to(ROOT)),
                "sha256": sha256_file(palette_path),
            },
            "script": {
                "path": str(script_path.relative_to(ROOT)),
                "sha256": sha256_file(script_path),
            },
        },
        "selected_umap": selection_metadata,
        "cluster_key": CLUSTER_KEY,
        "n_cells": int(adata.n_obs),
        "n_clusters": len(clusters),
    }
    (out_dir / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    sample_sizes = (
        table.drop_duplicates(["tissue", "condition"])
        .set_index("sample")["total_tissue_condition"]
        .astype(int)
        .to_dict()
    )
    artifacts = {
        "condition_share_png": str(condition_share_files[0].relative_to(ROOT)),
        "condition_share_pdf": str(condition_share_files[1].relative_to(ROOT)),
        "cluster_frequency_png": str(cluster_frequency_files[0].relative_to(ROOT)),
        "cluster_frequency_pdf": str(cluster_frequency_files[1].relative_to(ROOT)),
        "delta_png": str(delta_files[0].relative_to(ROOT)),
        "delta_pdf": str(delta_files[1].relative_to(ROOT)),
        "proportion_table": str(table_path.relative_to(ROOT)),
        "delta_table": str(delta_path.relative_to(ROOT)),
    }
    summary = {
        "analysis": "Control versus tumor cluster proportions stratified by tissue",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_pc": selected_pc,
        "cluster_key": CLUSTER_KEY,
        "cluster_order": clusters,
        "tissue_order": list(TISSUE_ORDER),
        "condition_order": list(CONDITION_ORDER),
        "sample_sizes": sample_sizes,
        "denominators": {
            "condition_share_within_cluster": "count / all cells in the same tissue x cluster",
            "cluster_frequency_within_condition": "count / all cells in the same tissue x condition",
        },
        "dominant_clusters": dominant,
        "largest_absolute_frequency_shift_per_tissue": maximum_shifts,
        "sample_palette": SAMPLE_PALETTE,
        "interpretation": "Tissue differences dominate; condition comparisons are descriptive only because there is one sample per condition x tissue.",
        "artifacts": artifacts,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_readme(out_dir / "README.md", h5ad_path, selected_pc, dominant, maximum_shifts)

    expected = [
        *condition_share_files,
        *cluster_frequency_files,
        *delta_files,
        table_path,
        delta_path,
        out_dir / "README.md",
        out_dir / "summary.json",
        out_dir / "input_manifest.json",
    ]
    bad = [str(path) for path in expected if not path.exists() or path.stat().st_size == 0]
    if bad:
        raise RuntimeError(f"Missing or empty artifacts: {bad}")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
