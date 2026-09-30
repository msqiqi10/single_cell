#!/usr/bin/env python3
"""Build a visual audit of two legacy-to-current response-program matches.

The figure triangulates three distinct evidence types without implying a
cell-identity mapping:

1. tissue context;
2. recovery of legacy T2-vs-Ctrl response genes;
3. relative iNKT1 program enrichment in the legacy and current analyses.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "output/iNKT_legacy_ppt_qc_runs/20260818_125231"
H5AD = RUN / "preprocess/inkt_scanpy_tutorial_processed.h5ad"
MATCHES = RUN / "comparison_to_legacy_ppt/legacy_ppt_de_best_matches.csv"
CALLS = RUN / "extended/signatures/standardized_legacy/standardized_cluster_calls.csv"
LEGACY_PPT = ROOT / "input/iNKT/iNKT.pptx"
OUT = ROOT / "output/iNKT_reproduction_deck/figures"

NAVY = "#102A43"
INK = "#172B3A"
MID = "#52616B"
MUTED = "#7B8794"
LIGHT = "#DCE5ED"
PALE = "#F5F8FA"
TEAL = "#1F8A8A"
ORANGE = "#E76F51"
GREEN = "#2E8B57"
BLUE = "#4472C4"
RED = "#C94C4C"


def ppt_image(media_name: str) -> Image.Image:
    with ZipFile(LEGACY_PPT) as archive:
        payload = archive.read(f"ppt/media/{media_name}")
    return Image.open(BytesIO(payload)).convert("RGB")


def card(ax, face: str = "white") -> None:
    ax.set_facecolor(face)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(LIGHT)
        spine.set_linewidth(1.0)


def add_panel_tag(ax, text: str, color: str) -> None:
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.35", facecolor=color, edgecolor=color),
        zorder=20,
    )


def draw_legacy_map(ax, image: Image.Image, pair: dict) -> None:
    ax.imshow(image)
    ax.add_patch(
        Ellipse(
            pair["legacy_ellipse"][:2],
            pair["legacy_ellipse"][2],
            pair["legacy_ellipse"][3],
            fill=False,
            edgecolor=pair["color"],
            linewidth=5.0,
            zorder=10,
        )
    )
    ax.annotate(
        pair["legacy_cluster"],
        xy=pair["legacy_point"],
        xytext=pair["legacy_label_xy"],
        color=pair["color"],
        fontsize=13,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=pair["color"], lw=2.0),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=pair["color"], alpha=0.94),
        zorder=20,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    card(ax)
    add_panel_tag(ax, "ORIGINAL PPT CLUSTER MAP", NAVY)


def draw_legacy_violin(ax, image: Image.Image, pair: dict) -> None:
    ax.imshow(image)
    x = pair["legacy_violin_x"]
    ax.add_patch(
        Rectangle(
            (x - 23, 18),
            46,
            492,
            facecolor=pair["color"],
            edgecolor=pair["color"],
            linewidth=2.5,
            alpha=0.20,
            zorder=10,
        )
    )
    ax.annotate(
        pair["legacy_violin_label"],
        xy=(x, 105),
        xytext=(0.50, 0.87),
        textcoords="axes fraction",
        ha="center",
        va="top",
        color=pair["color"],
        fontsize=10,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=pair["color"], lw=1.8),
        zorder=20,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    card(ax)
    add_panel_tag(ax, "LEGACY iNKT1 RAW SCORE", NAVY)
    ax.text(
        0.5,
        0.015,
        "Relative position only; legacy raw units",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.2,
        color=MUTED,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.5),
        zorder=30,
    )


def evidence_box(ax, y: float, title: str, body: str, accent: str) -> None:
    patch = FancyBboxPatch(
        (0.04, y),
        0.92,
        0.22,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        transform=ax.transAxes,
        facecolor="white",
        edgecolor=LIGHT,
        linewidth=1.0,
    )
    ax.add_patch(patch)
    ax.add_patch(
        Rectangle(
            (0.04, y),
            0.018,
            0.22,
            transform=ax.transAxes,
            facecolor=accent,
            edgecolor="none",
            zorder=3,
        )
    )
    ax.text(0.085, y + 0.178, title, transform=ax.transAxes, ha="left", va="top", fontsize=9.2, color=accent, fontweight="bold")
    ax.text(0.085, y + 0.128, body, transform=ax.transAxes, ha="left", va="top", fontsize=8.0, color=INK, linespacing=1.25)


def draw_evidence(ax, pair: dict) -> None:
    ax.set_axis_off()
    ax.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            transform=ax.transAxes,
            facecolor=PALE,
            edgecolor=LIGHT,
            linewidth=1.0,
        )
    )
    ax.text(0.05, 0.965, "WHY THE PROGRAMS ARE CONCORDANT", transform=ax.transAxes, ha="left", va="top", fontsize=10.2, color=NAVY, fontweight="bold")
    evidence_box(ax, 0.68, "1  Tissue context", pair["tissue_text"], TEAL)
    evidence_box(ax, 0.40, "2  T2-vs-Ctrl response", pair["response_text"], BLUE)
    evidence_box(ax, 0.12, "3  iNKT1 program", pair["subtype_text"], GREEN)
    ax.text(0.05, 0.022, pair["boundary_text"], transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, color=RED, fontweight="bold", linespacing=1.15)


def draw_current_umap(ax, coords: np.ndarray, clusters: np.ndarray, pair: dict) -> None:
    target = clusters == pair["current_cluster"].removeprefix("c")
    assert int(target.sum()) == pair["current_n"]
    ax.scatter(coords[~target, 0], coords[~target, 1], s=1.2, c="#D5DEE6", alpha=0.33, linewidths=0, rasterized=True)
    ax.scatter(coords[target, 0], coords[target, 1], s=3.4, c=pair["color"], alpha=0.90, linewidths=0, rasterized=True)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    card(ax)
    add_panel_tag(ax, "CURRENT UMAP", NAVY)
    ax.text(
        0.97,
        0.96,
        f"current {pair['current_cluster']}\n{pair['current_n']:,} cells",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.5,
        color=pair["color"],
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.30", facecolor="white", edgecolor=pair["color"], alpha=0.94),
    )
    ax.text(0.03, 0.03, pair["current_tissue"], transform=ax.transAxes, ha="left", va="bottom", fontsize=8.0, color=MID)


def draw_score_bars(ax, scores: dict[str, float], pair: dict) -> None:
    labels = [f"c{i}" for i in range(8)]
    values = [scores[label] for label in labels]
    colors = [pair["color"] if label == pair["current_cluster"] else "#C8D2DC" for label in labels]
    ax.bar(np.arange(8), values, color=colors, width=0.72, edgecolor="none")
    ax.axhline(0, color=MUTED, linestyle="--", linewidth=1.0)
    ax.set_ylim(-3.55, 1.15)
    ax.set_xticks(np.arange(8), labels, fontsize=7.5)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("mean standardized iNKT1", fontsize=7.8, color=MID)
    ax.grid(axis="y", color=LIGHT, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    card(ax)
    add_panel_tag(ax, "CURRENT STANDARDIZED SCORE", NAVY)
    value = scores[pair["current_cluster"]]
    ax.text(
        0.97,
        0.91,
        f"{value:+.2f} SD",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        color=pair["color"],
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor=pair["color"], alpha=0.94),
    )
    ax.text(
        0.5,
        0.025,
        "Same standardized SD unit in both rows; zero = dataset-wide mean",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.0,
        color=MUTED,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.5),
    )


def load_inputs() -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    adata = ad.read_h5ad(H5AD, backed="r")
    coords = np.asarray(adata.obsm["X_umap"])
    clusters = np.asarray(adata.obs["leiden_res_0_5"].astype(str))
    calls = pd.read_csv(CALLS)
    scores = {f"c{int(row.cluster)}": float(row.mean_inkt1) for row in calls.itertuples(index=False)}
    return coords, clusters, scores


def validate_matches() -> None:
    matches = pd.read_csv(MATCHES)
    c1 = matches.loc[matches["legacy_unit"].eq("c1_thymus")].iloc[0]
    c2 = matches.loc[matches["legacy_unit"].eq("c2_bone_marrow")].iloc[0]
    assert int(c1.current_cluster) == 6 and int(c1.n_shared) == 26 and int(c1.sign_agree) == 26
    assert int(c2.current_cluster) == 4 and int(c2.n_shared) == 17 and int(c2.sign_agree) == 17


def main() -> None:
    validate_matches()
    coords, clusters, scores = load_inputs()
    legacy_map = ppt_image("image62.png")
    legacy_violin = ppt_image("image19.png")

    pairs = [
        {
            "color": TEAL,
            "legacy_cluster": "legacy c1 = Thymus",
            "legacy_ellipse": (320, 390, 430, 330),
            "legacy_point": (300, 390),
            "legacy_label_xy": (665, 165),
            "legacy_violin_x": 466,
            "legacy_violin_label": "Thymus / legacy c1\nrelative-high iNKT1",
            "current_cluster": "c6",
            "current_n": 2379,
            "current_tissue": "97.3% thymus",
            "tissue_text": "legacy c1 is ~98.9% thymus.\ncurrent c6 is 97.3% thymus.",
            "response_text": "c1-thymus → c6-thymus: 26/30 genes\nrecovered (86.7%); direction 26/26.\nCurrent: 304 DEGs; Jaccard 0.084.\nExamples: Ccl5, Cxcr6, Klrk1, Fos, Hspa1a/b.",
            "subtype_text": "Legacy c1 is relatively high in the old iNKT1\npanel. Current c6 is +0.812 SD\n(95% cell-bootstrap CI +0.777 to +0.846);\nmargin vs iNKT2: +0.631 [+0.576,+0.685].",
            "boundary_text": "Concordant program evidence ≠ direct cell identity.\nNo original per-cell legacy labels are available.",
        },
        {
            "color": ORANGE,
            "legacy_cluster": "legacy c2 = BS_l1",
            "legacy_ellipse": (405, 540, 275, 205),
            "legacy_point": (405, 540),
            "legacy_label_xy": (710, 730),
            "legacy_violin_x": 287,
            "legacy_violin_label": "B_S_l1 / legacy c2\nrelative-high iNKT1",
            "current_cluster": "c4",
            "current_n": 887,
            "current_tissue": "52.6% BM · 45.5% spleen",
            "tissue_text": "legacy c2: ~49.1% BM, 45.9% spleen.\ncurrent c4: 52.6% BM, 45.5% spleen.",
            "response_text": "BM comparison: c2-BM → c4-BM: 17/18 genes\nrecovered (94.4%); direction 17/17.\nCurrent: 83 DEGs; Jaccard 0.202.\nExamples: Cd52, Dnaja1, Hspa8, Cox6c/7c, mt-Nd2/3.",
            "subtype_text": "Legacy c2 is relatively high in the old iNKT1\npanel. Current c4 is +0.725 SD\n(95% cell-bootstrap CI +0.671 to +0.778);\nmargin vs iNKT2: +0.920 [+0.848,+0.991].",
            "boundary_text": "Partial match: the legacy c2 spleen response is split\nacross current c3/c4.",
        },
    ]

    fig = plt.figure(figsize=(18, 10.5), facecolor="white")
    grid = fig.add_gridspec(
        2,
        5,
        left=0.025,
        right=0.985,
        top=0.865,
        bottom=0.105,
        wspace=0.12,
        hspace=0.26,
        width_ratios=(1.30, 0.78, 1.42, 1.12, 1.02),
    )

    fig.text(0.03, 0.955, "Why legacy c1→current c6 and legacy c2→current c4 are program-level concordant", ha="left", va="top", fontsize=24, fontweight="bold", color=NAVY)
    fig.text(0.03, 0.908, "Three evidence axes are shown separately: tissue context, treatment-response recovery, and relative iNKT1 enrichment.", ha="left", va="top", fontsize=11.5, color=MID)

    for row, pair in enumerate(pairs):
        draw_legacy_map(fig.add_subplot(grid[row, 0]), legacy_map, pair)
        draw_legacy_violin(fig.add_subplot(grid[row, 1]), legacy_violin, pair)
        draw_evidence(fig.add_subplot(grid[row, 2]), pair)
        draw_current_umap(fig.add_subplot(grid[row, 3]), coords, clusters, pair)
        draw_score_bars(fig.add_subplot(grid[row, 4]), scores, pair)

    fig.text(
        0.03,
        0.055,
        "Interpretation boundary: the legacy cluster map and iNKT1 violin are raster images from the original PPT. "
        "Response matching uses same-tissue T2-vs-Ctrl DEG overlap, not shared barcodes or whole-transcriptome identity. "
        "A direct one-to-one map requires the original per-cell legacy cluster assignments. "
        "Direction agreement assumes positive legacy logFC means T2−Ctrl; the original group order is unavailable.",
        ha="left",
        va="bottom",
        fontsize=9.0,
        color=MID,
        wrap=True,
    )
    fig.text(0.03, 0.022, "Sources: original iNKT.pptx pages 5/6/12; legacy_ppt_de_best_matches.csv; standardized_cluster_calls.csv. Score intervals are cell-bootstrap stability, not biological-replicate inference.", ha="left", va="bottom", fontsize=8.0, color=MUTED)

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "legacy_current_inkt1_crosswalk_c1_c6_c2_c4.png"
    pdf = OUT / "legacy_current_inkt1_crosswalk_c1_c6_c2_c4.pdf"
    fig.savefig(png, dpi=180, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
