#!/usr/bin/env python3
"""Build a consolidated 16:9 PDF from the dated iNKT results of 2026-08-25."""

from __future__ import annotations

import argparse
import hashlib
import json
import textwrap
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[3]
PC_DIR = ROOT / "output/iNKT_reproduction_deck/20260825_UMAP_PC_sweep"
MARKER_DIR = ROOT / "output/iNKT_reproduction_deck/20260825_reference_marker_UMAP"
PROPORTION_DIR = (
    ROOT
    / "output/iNKT_reproduction_deck/20260825_control_vs_tumor_cluster_proportions_by_tissue"
)
DEFAULT_OUT_DIR = ROOT / "output/iNKT_reproduction_deck/20260825_daily_results_summary"
DEFAULT_PDF_NAME = "20260825_iNKT_daily_results_summary.pdf"

PAGE_SIZE = (13.333333, 7.5)
PAGE_DPI = 180
NAVY = "#102F4A"
TEAL = "#188F91"
TEAL_LIGHT = "#DDF3F2"
BLUE = "#4776C6"
GREEN = "#2E9B61"
ORANGE = "#F07F2F"
RED = "#D34A4A"
PURPLE = "#7654D1"
INK = "#172B3A"
MUTED = "#657789"
PAPER = "#F7F9FC"
BORDER = "#D7E1E8"
DATE_LABEL = "25 August 2026"


@dataclass(frozen=True)
class PageRecord:
    page: int
    section: str
    title: str
    sources: tuple[Path, ...]
    note: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build today's consolidated iNKT results PDF.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--pdf-name", default=DEFAULT_PDF_NAME)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing PDF source files: {missing}")


def new_page(
    page_number: int,
    section: str,
    title: str,
    *,
    subtitle: str | None = None,
) -> plt.Figure:
    figure = plt.figure(figsize=PAGE_SIZE, dpi=PAGE_DPI, facecolor=PAPER)
    canvas = figure.add_axes([0, 0, 1, 1])
    canvas.set_axis_off()
    canvas.add_patch(Rectangle((0, 0.979), 1, 0.021, facecolor=TEAL, edgecolor="none"))
    canvas.add_patch(
        FancyBboxPatch(
            (0.04, 0.894),
            0.16,
            0.052,
            boxstyle="round,pad=0.005,rounding_size=0.016",
            facecolor=TEAL_LIGHT,
            edgecolor="none",
        )
    )
    figure.text(0.12, 0.92, section.upper(), ha="center", va="center", fontsize=9, color=TEAL, fontweight="bold")
    figure.text(0.04, 0.858, title, ha="left", va="top", fontsize=20, color=NAVY, fontweight="bold")
    if subtitle:
        figure.text(0.04, 0.812, subtitle, ha="left", va="top", fontsize=8.8, color=MUTED)
    canvas.plot([0.04, 0.96], [0.785, 0.785], color=BORDER, linewidth=0.8)
    canvas.plot([0.04, 0.96], [0.047, 0.047], color=BORDER, linewidth=0.7)
    figure.text(0.04, 0.022, "iNKT daily results summary · 2026-08-25", ha="left", va="center", fontsize=6.8, color="#82909E")
    figure.text(0.96, 0.022, f"{page_number:02d}", ha="right", va="center", fontsize=7.2, color="#82909E", fontweight="bold")
    return figure


def full_slide_page(source: Path, page_number: int) -> plt.Figure:
    figure = plt.figure(figsize=PAGE_SIZE, dpi=PAGE_DPI, facecolor="white")
    axis = figure.add_axes([0, 0, 1, 1])
    axis.imshow(Image.open(source).convert("RGB"))
    axis.set_axis_off()
    axis.add_patch(
        Rectangle(
            (0, 0),
            1,
            0.042,
            transform=axis.transAxes,
            facecolor="white",
            edgecolor="none",
            zorder=5,
        )
    )
    axis.plot(
        [0.04, 0.96],
        [0.041, 0.041],
        transform=axis.transAxes,
        color=BORDER,
        linewidth=0.7,
        zorder=6,
    )
    figure.text(0.04, 0.018, "iNKT daily results summary · 2026-08-25", ha="left", va="center", fontsize=6.8, color="#82909E", zorder=7)
    figure.text(0.96, 0.018, f"{page_number:02d}", ha="right", va="center", fontsize=7.2, color="#82909E", fontweight="bold", zorder=7)
    return figure


def add_image(
    figure: plt.Figure,
    source: Path,
    bounds: tuple[float, float, float, float],
    *,
    border: bool = True,
) -> None:
    axis = figure.add_axes(bounds)
    axis.set_facecolor("white")
    axis.imshow(Image.open(source).convert("RGB"))
    axis.set_axis_off()
    if border:
        axis.add_patch(
            Rectangle(
                (0, 0),
                1,
                1,
                transform=axis.transAxes,
                facecolor="none",
                edgecolor=BORDER,
                linewidth=0.7,
                zorder=10,
            )
        )


def add_card(
    figure: plt.Figure,
    bounds: tuple[float, float, float, float],
    heading: str,
    body: str,
    *,
    accent: str = TEAL,
    facecolor: str = "white",
    heading_size: float = 10.5,
    body_size: float = 8.5,
) -> None:
    axis = figure.add_axes(bounds)
    axis.set_axis_off()
    axis.add_patch(
        FancyBboxPatch(
            (0.01, 0.01),
            0.98,
            0.98,
            boxstyle="round,pad=0.012,rounding_size=0.03",
            facecolor=facecolor,
            edgecolor=BORDER,
            linewidth=0.8,
        )
    )
    axis.add_patch(Rectangle((0.01, 0.01), 0.018, 0.98, facecolor=accent, edgecolor="none"))
    axis.text(0.07, 0.82, heading, ha="left", va="top", fontsize=heading_size, fontweight="bold", color=accent)
    width_chars = max(20, int(PAGE_SIZE[0] * bounds[2] * 72 / (body_size * 0.58)) - 6)
    wrapped_body = "\n".join(textwrap.fill(line, width=width_chars) for line in body.splitlines())
    axis.text(0.07, 0.64, wrapped_body, ha="left", va="top", fontsize=body_size, color=INK, linespacing=1.35)


def add_bullet_panel(
    figure: plt.Figure,
    bounds: tuple[float, float, float, float],
    heading: str,
    bullets: list[str],
    *,
    accent: str = TEAL,
    width_chars: int = 36,
    font_size: float = 8.2,
) -> None:
    axis = figure.add_axes(bounds)
    axis.set_axis_off()
    axis.add_patch(
        FancyBboxPatch(
            (0.005, 0.005),
            0.99,
            0.99,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor="white",
            edgecolor=BORDER,
            linewidth=0.8,
        )
    )
    axis.add_patch(Rectangle((0.005, 0.005), 0.018, 0.99, facecolor=accent, edgecolor="none"))
    axis.text(0.07, 0.94, heading, ha="left", va="top", fontsize=10.2, color=accent, fontweight="bold")
    y = 0.82
    line_height = font_size / 72 / (PAGE_SIZE[1] * bounds[3]) * 1.3
    for bullet in bullets:
        wrapped = textwrap.fill(bullet, width=width_chars)
        axis.text(0.075, y, "• " + wrapped, ha="left", va="top", fontsize=font_size, color=INK, linespacing=1.3)
        y -= line_height * (wrapped.count("\n") + 1) + 0.055


def add_takeaway_strip(
    figure: plt.Figure,
    text: str,
    *,
    accent: str = TEAL,
    y: float = 0.075,
    height: float = 0.085,
) -> None:
    axis = figure.add_axes([0.04, y, 0.92, height])
    axis.set_axis_off()
    axis.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor="#EDF7F6" if accent == TEAL else "white",
            edgecolor=BORDER,
            linewidth=0.8,
        )
    )
    axis.add_patch(Rectangle((0, 0), 0.012, 1, facecolor=accent, edgecolor="none"))
    width_chars = max(40, int(PAGE_SIZE[0] * 0.92 * 72 / (8.7 * 0.58)) - 12)
    wrapped = textwrap.fill(text, width=width_chars)
    axis.text(0.035, 0.5, wrapped, ha="left", va="center", fontsize=8.7, color=INK, linespacing=1.25)


def draw_cover(
    pc_summary: dict[str, Any],
    marker_summary: dict[str, Any],
    proportion_summary: dict[str, Any],
) -> plt.Figure:
    figure = plt.figure(figsize=PAGE_SIZE, dpi=PAGE_DPI, facecolor=NAVY)
    canvas = figure.add_axes([0, 0, 1, 1])
    canvas.set_axis_off()
    canvas.add_patch(Rectangle((0, 0.975), 1, 0.025, facecolor=TEAL, edgecolor="none"))
    figure.text(0.055, 0.885, "iNKT ANALYSIS", fontsize=11, color="#8FD7D4", fontweight="bold", ha="left")
    figure.text(0.055, 0.74, "Daily results summary", fontsize=34, color="white", fontweight="bold", ha="left")
    figure.text(0.055, 0.675, DATE_LABEL, fontsize=14, color="#B8C8D6", ha="left")
    figure.text(
        0.055,
        0.61,
        "UMAP PC sensitivity · literature-marker validation · tissue-stratified Control vs Tumor composition",
        fontsize=10,
        color="#D7E1E8",
        ha="left",
    )

    cards = [
        (
            "01  UMAP decision",
            f"{pc_summary['selected_pc']} PCs selected\nBest balance of cluster separation, trustworthiness, and legacy geometry.",
            TEAL,
        ),
        (
            "02  Marker audit",
            f"CD94: {marker_summary['focused_detection']['Klrd1']['detected_pct']:.2f}%\nIL-4: {marker_summary['focused_detection']['Il4']['detected_pct']:.2f}%\nCD56/Ncam1: 0 RNA reads.",
            ORANGE,
        ),
        (
            "03  Composition",
            "Thymus → c6\nBone marrow → c0\nSpleen → c3\nTissue effects dominate Ctrl/T2 shifts.",
            GREEN,
        ),
    ]
    for index, (heading, body, accent) in enumerate(cards):
        x = 0.055 + index * 0.307
        axis = figure.add_axes([x, 0.255, 0.275, 0.255])
        axis.set_axis_off()
        axis.add_patch(
            FancyBboxPatch(
                (0, 0),
                1,
                1,
                boxstyle="round,pad=0.015,rounding_size=0.04",
                facecolor="#173B58",
                edgecolor="#31526B",
                linewidth=0.9,
            )
        )
        axis.add_patch(Rectangle((0, 0), 0.018, 1, facecolor=accent, edgecolor="none"))
        axis.text(0.07, 0.84, heading, ha="left", va="top", fontsize=10, color=accent, fontweight="bold")
        axis.text(0.07, 0.65, body, ha="left", va="top", fontsize=9, color="white", linespacing=1.45)

    palette = [
        ("Ctrl thymus", "#0000FF"),
        ("T2 thymus", "#00008B"),
        ("Ctrl BM", "#FFFF00"),
        ("T2 BM", "#FFA500"),
        ("Ctrl spleen", "#FF0000"),
        ("T2 spleen", "#8B0000"),
    ]
    for index, (label, color) in enumerate(palette):
        x = 0.055 + index * 0.151
        canvas.add_patch(Rectangle((x, 0.13), 0.14, 0.022, facecolor=color, edgecolor="#FFFFFF", linewidth=0.4))
        figure.text(x, 0.105, label, ha="left", va="top", fontsize=6.5, color="#B8C8D6")
    figure.text(0.055, 0.055, "15,532 QC-retained cells · 10,670 genes · 6 condition–tissue samples", fontsize=9, color="#D7E1E8", ha="left")
    figure.text(0.96, 0.055, "01", fontsize=8, color="#91A6B6", ha="right", fontweight="bold")
    return figure


def build_contact_sheet(page_paths: list[Path], output_path: Path) -> None:
    columns = 3
    rows = 4
    thumb_width = 480
    thumb_height = 270
    label_height = 24
    sheet = Image.new("RGB", (columns * thumb_width, rows * (thumb_height + label_height)), "#E9EEF3")
    draw = ImageDraw.Draw(sheet)
    for index, path in enumerate(page_paths):
        image = Image.open(path).convert("RGB")
        image.thumbnail((thumb_width - 12, thumb_height - 12), Image.Resampling.LANCZOS)
        column = index % columns
        row = index // columns
        x = column * thumb_width + (thumb_width - image.width) // 2
        y = row * (thumb_height + label_height) + 6
        sheet.paste(image, (x, y))
        draw.text((column * thumb_width + 8, row * (thumb_height + label_height) + thumb_height + 2), f"Page {index + 1:02d}", fill="#17324D")
    sheet.save(output_path)


def write_readme(path: Path, pdf_path: Path, pages: list[PageRecord]) -> None:
    lines = [
        "# 20260825 iNKT daily results summary",
        "",
        f"- PDF: `{relative(pdf_path)}`",
        f"- Page count: {len(pages)}",
        "- Format: 16:9 landscape; all analytical figures are contained without cropping.",
        "- Language: English, matching the existing figure labels and avoiding unavailable CJK-font substitution.",
        "",
        "## Page index",
        "",
    ]
    for page in pages:
        lines.append(f"{page.page}. **{page.title}** — {page.note}")
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- Tumor means T2.",
            "- The 50/100/200-PC comparison changes only `n_pcs`; fixed audited `leiden_res_0_5` labels are overlaid.",
            "- Program scores and marker overlays are continuous evidence, not validated discrete subtype calls.",
            "- `Ncam1` has no RNA counts in all six raw matrices; this does not prove absence of CD56 protein biology.",
            "- Composition results are descriptive: one sample per condition×tissue, no biological-replicate significance test.",
            "- Percentages describe post-QC captured cells, not absolute tissue abundance.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    pages_dir = out_dir / "pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / args.pdf_name

    slide_pc = PC_DIR / "deck_selected_slides/slide_04.png"
    slide_selection = PC_DIR / "deck_selected_slides/slide_05.png"
    program_scores = PC_DIR / "umap_pc_sweep/figures/umap_selected_program_scores.png"
    focused_markers = MARKER_DIR / "figures/focused_cd56_cd94_il4_umap.png"
    marker_panel = MARKER_DIR / "figures/reference_marker_panel_umap.png"
    marker_cluster = MARKER_DIR / "figures/reference_marker_dotplot_by_cluster.png"
    marker_sample = MARKER_DIR / "figures/reference_marker_dotplot_by_sample.png"
    cluster_frequency = PROPORTION_DIR / "figures/cluster_frequency_ctrl_vs_tumor_by_tissue.png"
    cluster_delta = PROPORTION_DIR / "figures/cluster_frequency_delta_tumor_minus_control_by_tissue.png"
    condition_share = PROPORTION_DIR / "figures/ctrl_vs_tumor_share_within_cluster_by_tissue.png"
    pc_summary_path = PC_DIR / "umap_pc_sweep/summary.json"
    marker_summary_path = MARKER_DIR / "summary.json"
    proportion_summary_path = PROPORTION_DIR / "summary.json"
    marker_cluster_table = MARKER_DIR / "tables/marker_expression_by_cluster.csv"
    marker_sample_table = MARKER_DIR / "tables/marker_expression_by_sample.csv"

    source_files = [
        slide_pc,
        slide_selection,
        program_scores,
        focused_markers,
        marker_panel,
        marker_cluster,
        marker_sample,
        cluster_frequency,
        cluster_delta,
        condition_share,
        pc_summary_path,
        marker_summary_path,
        proportion_summary_path,
        marker_cluster_table,
        marker_sample_table,
    ]
    require_files(source_files)

    pc_summary = json.loads(pc_summary_path.read_text())
    marker_summary = json.loads(marker_summary_path.read_text())
    proportion_summary = json.loads(proportion_summary_path.read_text())
    cluster_table = pd.read_csv(marker_cluster_table)
    sample_table = pd.read_csv(marker_sample_table)

    klrd1_cluster = cluster_table[cluster_table["symbol"] == "Klrd1"].set_index("group")
    il4_cluster = cluster_table[cluster_table["symbol"] == "Il4"].set_index("group")
    klrd1_sample = sample_table[sample_table["symbol"] == "Klrd1"].set_index("group")
    il4_sample = sample_table[sample_table["symbol"] == "Il4"].set_index("group")

    pages: list[tuple[plt.Figure, PageRecord]] = []
    pages.append(
        (
            draw_cover(pc_summary, marker_summary, proportion_summary),
            PageRecord(1, "Cover", "Daily results summary", tuple(), "Three completed analysis blocks and the headline decisions."),
        )
    )
    pages.append(
        (
            full_slide_page(slide_pc, 2),
            PageRecord(2, "UMAP", "UMAP sensitivity across 50, 100, and 200 PCs", (slide_pc,), "Controlled sensitivity experiment; only the number of PCs changes."),
        )
    )
    pages.append(
        (
            full_slide_page(slide_selection, 3),
            PageRecord(3, "UMAP", "Selection result: 50 PCs", (slide_selection,), "Quantitative selection evidence and comparison with legacy geometry."),
        )
    )

    figure = new_page(
        4,
        "Biological sanity check",
        "Continuous subtype programs remain spatially coherent on the selected map",
        subtitle="Selected 50-PC UMAP · shared standardized scale · program gradients are not discrete cell labels",
    )
    add_image(figure, program_scores, (0.04, 0.23, 0.92, 0.50))
    add_takeaway_strip(
        figure,
        "iNKT1, iNKT2, and NKT17 evidence forms coherent but overlapping gradients. Use these scores as continuous support; do not force per-cell subtype calls.",
        accent=PURPLE,
        y=0.08,
        height=0.095,
    )
    pages.append(
        (
            figure,
            PageRecord(4, "Programs", "Continuous subtype-program distributions", (program_scores,), "Biological sanity check on the selected embedding."),
        )
    )

    focus = marker_summary["focused_detection"]
    figure = new_page(
        5,
        "Requested marker audit",
        "CD56, CD94, and IL-4 on the selected 50-PC UMAP",
        subtitle="Mouse RNA symbols: Ncam1, Klrd1, and Il4 · shared log1p expression scale",
    )
    add_image(figure, focused_markers, (0.035, 0.24, 0.93, 0.49))
    add_card(
        figure,
        (0.04, 0.075, 0.285, 0.12),
        "CD56 / Ncam1",
        "0 non-zero RNA reads across all six raw 10x matrices; no distribution can be inferred.",
        accent=RED,
        heading_size=9.2,
        body_size=7.4,
    )
    add_card(
        figure,
        (0.355, 0.075, 0.285, 0.12),
        "CD94 / Klrd1",
        f"{focus['Klrd1']['detected_cells']:,} cells ({focus['Klrd1']['detected_pct']:.2f}%) detected.",
        accent=GREEN,
        heading_size=9.2,
        body_size=7.4,
    )
    add_card(
        figure,
        (0.67, 0.075, 0.285, 0.12),
        "IL-4 / Il4",
        f"{focus['Il4']['detected_cells']:,} cells ({focus['Il4']['detected_pct']:.2f}%) detected; not subtype-exclusive alone.",
        accent=ORANGE,
        heading_size=9.2,
        body_size=7.4,
    )
    pages.append(
        (
            figure,
            PageRecord(5, "Markers", "Focused CD56/CD94/IL-4 audit", (focused_markers, marker_summary_path), "Requested marker distributions and explicit CD56 data limitation."),
        )
    )

    figure = new_page(
        6,
        "Literature marker atlas",
        "Twelve local-reference markers on the same selected UMAP",
        subtitle="iNKT1/NK-like · iNKT2/activation · iNKT17 context",
    )
    add_image(figure, marker_panel, (0.04, 0.07, 0.66, 0.69))
    add_bullet_panel(
        figure,
        (0.73, 0.13, 0.23, 0.56),
        "How to read this page",
        [
            "The shared color scale supports cross-gene comparison; zero-expression cells are gray.",
            "Il4 occurs in both local iNKT1 and iNKT2 lists and is not an isolated classifier.",
            "Sparse Rorc, Il1r1, and Il23r signals require joint interpretation with the dotplots.",
            "All panels use the same fixed 50-PC coordinates; no marker-driven reclustering was performed.",
        ],
        accent=ORANGE,
        width_chars=34,
        font_size=7.7,
    )
    pages.append(
        (
            figure,
            PageRecord(6, "Markers", "Twelve-gene literature marker atlas", (marker_panel,), "Compact multi-gene context for iNKT programs."),
        )
    )

    figure = new_page(
        7,
        "Marker summary by cluster",
        "Fixed clusters carry distinct but overlapping marker programs",
        subtitle="Dot size = detected-cell fraction · color = mean log1p expression across all cells in the cluster",
    )
    add_image(figure, marker_cluster, (0.04, 0.07, 0.66, 0.70))
    add_bullet_panel(
        figure,
        (0.73, 0.13, 0.23, 0.56),
        "Cluster-level findings",
        [
            f"CD94/Klrd1 detection is strongest in c6 ({klrd1_cluster.loc[6, 'detected_pct']:.1f}%) and c4 ({klrd1_cluster.loc[4, 'detected_pct']:.1f}%).",
            f"IL-4/Il4 detection is highest in c3 ({il4_cluster.loc[3, 'detected_pct']:.1f}%).",
            "c5 concentrates Il17rb, Zbtb16, Tmem176a, and sparse iNKT17-context markers.",
            "c7 contains only 92 cells; extreme percentages in this cluster should not be over-interpreted.",
        ],
        accent=GREEN,
        width_chars=34,
        font_size=7.7,
    )
    pages.append(
        (
            figure,
            PageRecord(7, "Markers", "Marker programs by fixed cluster", (marker_cluster, marker_cluster_table), "Detection fraction and mean expression by c0–c7."),
        )
    )

    figure = new_page(
        8,
        "Marker summary by sample",
        "Tissue differences exceed Control–Tumor differences",
        subtitle="Six columns are individual condition–tissue samples, not replicate means",
    )
    add_image(figure, marker_sample, (0.04, 0.07, 0.66, 0.70))
    add_bullet_panel(
        figure,
        (0.73, 0.13, 0.23, 0.56),
        "Sample-level findings",
        [
            f"CD94/Klrd1 is highest in thymus: {klrd1_sample.loc['Ctrl_Thymus', 'detected_pct']:.1f}% Ctrl and {klrd1_sample.loc['T2_Thymus', 'detected_pct']:.1f}% Tumor.",
            f"IL-4/Il4 is highest in spleen: {il4_sample.loc['Ctrl_Spleen', 'detected_pct']:.1f}% Ctrl and {il4_sample.loc['T2_Spleen', 'detected_pct']:.1f}% Tumor.",
            "Control and Tumor generally remain more similar within a tissue than across tissues.",
            "With one sample per condition×tissue, this page is descriptive and cannot separate condition from sample effects.",
        ],
        accent=BLUE,
        width_chars=34,
        font_size=7.7,
    )
    pages.append(
        (
            figure,
            PageRecord(8, "Markers", "Marker patterns by condition–tissue sample", (marker_sample, marker_sample_table), "Tissue-stratified marker detection and expression."),
        )
    )

    dominant = proportion_summary["dominant_clusters"]
    figure = new_page(
        9,
        "Primary composition comparison",
        "Control vs Tumor cluster frequencies within each tissue",
        subtitle="Denominator: all cells in the same tissue×condition; each condition series sums to 100%",
    )
    add_image(figure, cluster_frequency, (0.03, 0.16, 0.94, 0.60))
    add_takeaway_strip(
        figure,
        f"Thymus is c6-dominant ({dominant['thymus']['Ctrl']['pct']:.1f}% Ctrl; {dominant['thymus']['T2']['pct']:.1f}% Tumor), bone marrow is c0-dominant ({dominant['bone_marrow']['Ctrl']['pct']:.1f}%; {dominant['bone_marrow']['T2']['pct']:.1f}%), and spleen is c3-dominant ({dominant['spleen']['Ctrl']['pct']:.1f}%; {dominant['spleen']['T2']['pct']:.1f}%).",
        accent=TEAL,
        y=0.07,
        height=0.095,
    )
    pages.append(
        (
            figure,
            PageRecord(9, "Composition", "Cluster frequency within tissue×condition", (cluster_frequency, proportion_summary_path), "Recommended denominator for descriptive composition shifts."),
        )
    )

    shifts = proportion_summary["largest_absolute_frequency_shift_per_tissue"]
    figure = new_page(
        10,
        "Descriptive condition shifts",
        "Tumor minus Control cluster-frequency changes",
        subtitle="Percentage-point differences within each tissue; no biological-replicate significance test",
    )
    add_image(figure, cluster_delta, (0.03, 0.16, 0.94, 0.60))
    add_takeaway_strip(
        figure,
        f"Largest absolute shifts: spleen c{shifts['spleen']['cluster']} {shifts['spleen']['tumor_minus_control_percentage_points']:+.2f} pp, thymus c{shifts['thymus']['cluster']} {shifts['thymus']['tumor_minus_control_percentage_points']:+.2f} pp, and bone marrow c{shifts['bone_marrow']['cluster']} {shifts['bone_marrow']['tumor_minus_control_percentage_points']:+.2f} pp. Tissue structure is much larger than these within-tissue shifts.",
        accent=RED,
        y=0.07,
        height=0.095,
    )
    pages.append(
        (
            figure,
            PageRecord(10, "Composition", "Tumor − Control frequency deltas", (cluster_delta, proportion_summary_path), "Locations and magnitudes of descriptive composition shifts."),
        )
    )

    figure = new_page(
        11,
        "Appendix · alternative denominator",
        "Control vs Tumor share inside each tissue×cluster",
        subtitle="Denominator: all cells in the same tissue×cluster; every non-empty stacked bar sums to 100%",
    )
    add_image(figure, condition_share, (0.03, 0.16, 0.94, 0.60))
    add_takeaway_strip(
        figure,
        "This view answers a different question and is sensitive to unequal captured-cell totals between Control and Tumor. Read it with the n above each bar; tiny or empty clusters can create extreme percentages.",
        accent=ORANGE,
        y=0.07,
        height=0.095,
    )
    pages.append(
        (
            figure,
            PageRecord(11, "Appendix", "Condition share within tissue×cluster", (condition_share, proportion_summary_path), "Legacy-style denominator audit; not the primary composition comparison."),
        )
    )

    figure = new_page(
        12,
        "Decision register",
        "Conclusions, limits, and output locations",
        subtitle="Today's three dated result packages remain unchanged; this PDF is a presentation layer over them",
    )
    add_card(
        figure,
        (0.04, 0.58, 0.285, 0.17),
        "Use 50 PCs",
        "Best overall separation and fidelity. Keep n_neighbors=20, min_dist=0.5, and fixed audited c0–c7 labels.",
        accent=TEAL,
    )
    add_card(
        figure,
        (0.357, 0.58, 0.285, 0.17),
        "Use multi-gene context",
        "CD94 and IL-4 are measurable; CD56 RNA is unavailable. Program scores and markers do not create discrete subtype calls.",
        accent=ORANGE,
    )
    add_card(
        figure,
        (0.674, 0.58, 0.285, 0.17),
        "Stratify by tissue",
        "Tissue composition dominates. Ctrl/Tumor differences are descriptive because each condition×tissue has one sample.",
        accent=GREEN,
    )
    add_bullet_panel(
        figure,
        (0.04, 0.12, 0.44, 0.39),
        "Interpretation guardrails",
        [
            "Tumor = T2; percentages describe post-QC captured cells, not absolute tissue abundance.",
            "The weighted PC score is a ranking aid, not an effect size or probability.",
            "Leiden density contours are descriptive boundaries; fixed clusters were not recomputed per PC candidate.",
            "Small clusters and sparse genes require raw n and multi-marker context.",
        ],
        accent=RED,
        width_chars=58,
        font_size=7.8,
    )
    add_bullet_panel(
        figure,
        (0.52, 0.12, 0.44, 0.39),
        "Dated source packages",
        [
            "20260825_UMAP_PC_sweep — PC candidates, metrics, selected h5ad, and regenerated deck.",
            "20260825_reference_marker_UMAP — marker UMAPs, dotplots, gene audit, and raw Ncam1 audit.",
            "20260825_control_vs_tumor_cluster_proportions_by_tissue — both denominators, frequency deltas, and tables.",
            "Every package includes a README and machine-readable JSON/CSV provenance.",
        ],
        accent=BLUE,
        width_chars=58,
        font_size=7.8,
    )
    pages.append(
        (
            figure,
            PageRecord(
                12,
                "Summary",
                "Conclusions and artifact index",
                (pc_summary_path, marker_summary_path, proportion_summary_path),
                "Decision register and interpretation limits.",
            ),
        )
    )

    page_paths: list[Path] = []
    records: list[PageRecord] = []
    with PdfPages(
        pdf_path,
        metadata={
            "Title": "iNKT analysis — daily results summary — 2026-08-25",
            "Author": "Codex",
            "Subject": "UMAP PC sensitivity, marker validation, and tissue-stratified Control vs Tumor composition",
            "Keywords": "iNKT, UMAP, single-cell RNA-seq, CD94, IL-4, cluster composition",
            "CreationDate": datetime.now(timezone.utc),
        },
    ) as pdf:
        for figure, record in pages:
            if record.page != len(records) + 1:
                raise RuntimeError(f"Non-sequential page record: {record}")
            page_path = pages_dir / f"page_{record.page:02d}.png"
            figure.savefig(page_path, dpi=PAGE_DPI, facecolor=figure.get_facecolor())
            pdf.savefig(figure, dpi=PAGE_DPI, facecolor=figure.get_facecolor())
            plt.close(figure)
            page_paths.append(page_path)
            records.append(record)

    contact_sheet_path = out_dir / "contact_sheet.png"
    build_contact_sheet(page_paths, contact_sheet_path)

    index = pd.DataFrame(
        [
            {
                "page": record.page,
                "section": record.section,
                "title": record.title,
                "source_paths": ";".join(relative(path) for path in record.sources),
                "note": record.note,
            }
            for record in records
        ]
    )
    index_path = out_dir / "page_index.csv"
    index.to_csv(index_path, index=False)
    readme_path = out_dir / "README.md"
    write_readme(readme_path, pdf_path, records)

    unique_sources = sorted({path.resolve() for record in records for path in record.sources})
    script_path = Path(__file__).resolve()
    manifest = {
        "analysis": "iNKT daily results consolidated PDF",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "page_count": len(records),
        "page_size_inches": list(PAGE_SIZE),
        "pdf": {
            "path": relative(pdf_path),
            "sha256": sha256_file(pdf_path),
            "size_bytes": pdf_path.stat().st_size,
        },
        "contact_sheet": {
            "path": relative(contact_sheet_path),
            "sha256": sha256_file(contact_sheet_path),
        },
        "page_previews": [
            {"path": relative(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in page_paths
        ],
        "sources": [
            {"path": relative(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in unique_sources
        ],
        "script": {"path": relative(script_path), "sha256": sha256_file(script_path)},
        "source_packages": [relative(PC_DIR), relative(MARKER_DIR), relative(PROPORTION_DIR)],
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    expected = [pdf_path, contact_sheet_path, index_path, readme_path, manifest_path, *page_paths]
    bad = [str(path) for path in expected if not path.exists() or path.stat().st_size == 0]
    if bad:
        raise RuntimeError(f"Missing or empty PDF artifacts: {bad}")
    print(json.dumps({"pdf": relative(pdf_path), "page_count": len(records), "contact_sheet": relative(contact_sheet_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
