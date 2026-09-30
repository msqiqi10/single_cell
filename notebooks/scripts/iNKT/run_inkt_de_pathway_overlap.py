#!/usr/bin/env python3
"""Extended iNKT differential-expression and enrichment analyses.

This script deliberately keeps the old PowerPoint/XLSX results and the current
Scanpy clustering in separate namespaces.  It supports three independently
runnable tasks so that they can be placed in separate tmux jobs::

    python run_inkt_de_pathway_overlap.py --task legacy
    python run_inkt_de_pathway_overlap.py --task current-de --n-jobs 12
    python run_inkt_de_pathway_overlap.py --task overlap-pathway

The current DE p-values are *cell-level exploratory* statistics.  The six input
matrices do not currently expose independent animal identifiers, so cells must
not be treated as biological replicates in confirmatory inference.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

# Each DE comparison is already parallelized at the process level.  Prevent
# BLAS/OpenMP oversubscription and make it impossible for libraries to select a
# GPU (GPU 6 on this host is known to be faulty).
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_H5AD = ROOT / "output/iNKT_scanpy_tutorial_run/inkt_scanpy_tutorial_processed.h5ad"
DEFAULT_GENE_LIST = ROOT / "input/iNKT/iNKT gene list for scRNAseq.xlsx"
DEFAULT_MARKER = ROOT / "input/iNKT/marker.xlsx"
DEFAULT_OUT = ROOT / "output/iNKT_xlsx_extended_analysis"

CURRENT_CLUSTER_NAMESPACE = "current_leiden_res_0_5_c0_c8"
LEGACY_CLUSTER_NAMESPACE = "legacy_ppt_c1_c11"
STATISTICAL_UNIT_WARNING = (
    "cell_level_exploratory_no_independent_animal_ids; "
    "p-values_do_not_establish_biological_replication"
)

# Conservative corrections for clear historical symbols/typos documented in
# the supplied workbook review.  Ambiguous values such as Slam4 and Opct are
# intentionally left untouched.
GENE_ALIASES = {
    "sepp1": "Selenop",
    "cd103": "Itgae",
    "cd244": "Cd244a",
    "slam7": "Slamf7",
    "ckd6": "Cdk6",
    "hnmpa1": "Hnrnpa1",
    "dfna5": "Gsdme",
    "itgab1": "Itgb1",
    "itgab4": "Itgb4",
    "itgab5": "Itgb5",
}


@dataclass(frozen=True)
class DEUnit:
    unit_id: str
    scope: str
    tissue: str
    cluster_current: str | None = None


def slug(value: Any) -> str:
    text = re.sub(r"[^0-9A-Za-z._-]+", "-", str(value)).strip("-._")
    return text or "unnamed"


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _json_safe_number(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    return value


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg correction with NaN preservation."""

    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    finite = np.isfinite(arr)
    if not finite.any():
        return out
    p = arr[finite]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[finite] = restored
    return out


def normalize_tissue(value: Any) -> str | None:
    text = str(value).strip().casefold()
    if not text or text == "nan":
        return None
    if text in {"bm", "bone marrow", "bone_marrow", "bonemarrow"}:
        return "bone_marrow"
    if "thym" in text:
        return "thymus"
    if "spleen" in text or text == "spl":
        return "spleen"
    return None


def infer_legacy_scope(label: str, sheet: str) -> dict[str, Any]:
    text = f"{label} {sheet}".casefold()
    tissue = None
    for token in ("bone marrow", "bone_marrow", "bm", "spleen", "thymus", "thym"):
        if re.search(rf"(?<![a-z]){re.escape(token)}(?![a-z])", text):
            tissue = normalize_tissue(token)
            break
    cluster_match = re.search(r"(?<![a-z0-9])c(?:luster)?\s*[-_ ]?\s*(1[01]|[0-9])(?![0-9])", text)
    cluster = cluster_match.group(1) if cluster_match else None
    is_global = bool(re.search(r"(?<![a-z])global(?![a-z])", text))
    if cluster is not None and tissue is not None:
        scope = "legacy_cluster_tissue"
    elif cluster is not None:
        scope = "legacy_cluster"
    elif tissue is not None:
        scope = "legacy_tissue"
    elif is_global:
        scope = "legacy_global"
    else:
        scope = "legacy_unresolved"
    return {
        "legacy_scope": scope,
        "tissue": tissue,
        "cluster_legacy": cluster,
        "cluster_namespace": LEGACY_CLUSTER_NAMESPACE if cluster is not None else "legacy_no_cluster",
    }


def normalize_header(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().casefold())


def canonical_header(value: Any, *, pathway_context: bool = False) -> str | None:
    token = normalize_header(value)
    if not token:
        return None
    exact = {
        "gene": "gene",
        "genes": "gene",
        "genename": "gene",
        "genesymbol": "gene",
        "symbol": "gene",
        "officialgenesymbol": "gene",
        "feature": "gene",
        "avglog2fc": "logfoldchange",
        "avglogfc": "logfoldchange",
        "log2fc": "logfoldchange",
        "logfc": "logfoldchange",
        "logfoldchange": "logfoldchange",
        "foldchange": "foldchange",
        "pval": "pvalue",
        "pvalue": "pvalue",
        "pvalues": "pvalue",
        "rawp": "pvalue",
        "pvaladj": "fdr",
        "padj": "fdr",
        "adjustedpvalue": "fdr",
        "adjpvalue": "fdr",
        "qvalue": "fdr",
        "fdr": "fdr",
        "pct1": "pct_t2",
        "pct2": "pct_ctrl",
        "pathway": "term",
        "pathwayname": "term",
        "term": "term",
        "description": "term",
        "geneset": "term",
        "genesetname": "term",
        "pagname": "term",
        "gsid": "pathway_id",
        "pagid": "pathway_id",
        "pathwayid": "pathway_id",
        "source": "pathway_source",
        "overlap": "overlap",
        "overlapgenes": "overlap_genes",
        "overlappinggenes": "overlap_genes",
        "gssize": "gene_set_size",
        "genesetsize": "gene_set_size",
        "oddsratio": "odds_ratio",
        "score": "score",
        "zscore": "zscore",
    }
    if token in exact:
        return exact[token]
    if token in {"names", "name"}:
        return "term" if pathway_context else "gene"
    if "log2fold" in token or ("log" in token and "foldchange" in token):
        return "logfoldchange"
    if "adjust" in token and ("pval" in token or token.endswith("p")):
        return "fdr"
    if token.startswith("pval"):
        return "pvalue"
    return None


def _nonempty(value: Any) -> bool:
    return not (value is None or (isinstance(value, float) and np.isnan(value)) or str(value).strip() == "")


def _sheet_context_is_pathway(sheet_name: str) -> bool:
    text = sheet_name.casefold()
    return any(term in text for term in ("pathway", "pager", "enrich", "gsea", "ontology"))


def read_excel_raw(path: Path) -> dict[str, pd.DataFrame]:
    """Read every XLSX sheet without trusting any particular header row."""

    return pd.read_excel(path, sheet_name=None, header=None, dtype=object, engine="openpyxl")


def _label_above(df: pd.DataFrame, header_row: int, left: int, right: int) -> str:
    values: list[str] = []
    for row in range(max(0, header_row - 4), header_row):
        for col in range(left, right + 1):
            value = df.iat[row, col]
            if _nonempty(value):
                text = str(value).strip()
                if text not in values:
                    values.append(text)
    return " | ".join(values[-4:])


def detect_tables(
    df: pd.DataFrame,
    *,
    workbook: str,
    sheet: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Detect side-by-side or vertically stacked DE/pathway tables.

    Returns row records plus a table-level audit.  The supplied workbooks use
    both layouts, so table boundaries are inferred from repeated gene/pathway
    header anchors rather than from fixed Excel cell ranges.
    """

    records: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    pathway_sheet = _sheet_context_is_pathway(sheet)
    table_counter = 0

    for header_row in range(df.shape[0]):
        mappings: dict[int, str] = {}
        for col in range(df.shape[1]):
            key = canonical_header(df.iat[header_row, col], pathway_context=pathway_sheet)
            if key is not None:
                mappings[col] = key

        gene_anchors = [col for col, key in mappings.items() if key == "gene"]
        term_anchors = [col for col, key in mappings.items() if key == "term"]
        id_anchors = [col for col, key in mappings.items() if key == "pathway_id"]
        if pathway_sheet or term_anchors or id_anchors:
            anchors = term_anchors or id_anchors
            kind = "pathway"
        else:
            anchors = gene_anchors
            kind = "de"
        if not anchors:
            continue

        # Midpoints between anchors isolate parallel tables on the same row.
        for anchor_index, anchor in enumerate(anchors):
            left = 0 if anchor_index == 0 else (anchors[anchor_index - 1] + anchor) // 2 + 1
            right = (
                df.shape[1] - 1
                if anchor_index == len(anchors) - 1
                else (anchor + anchors[anchor_index + 1]) // 2
            )
            local = {col: key for col, key in mappings.items() if left <= col <= right}
            if kind == "de" and not ({"gene", "logfoldchange", "pvalue", "fdr"} & set(local.values())):
                continue
            if kind == "pathway" and not ({"term", "pathway_id"} & set(local.values())):
                continue

            table_counter += 1
            table_id = f"{slug(workbook)}__{slug(sheet)}__{kind}_{table_counter:02d}"
            label = _label_above(df, header_row, left, right)
            scope = infer_legacy_scope(label, sheet)
            n_rows = 0
            blank_run = 0

            for row in range(header_row + 1, df.shape[0]):
                anchor_value = df.iat[row, anchor]
                local_values = [df.iat[row, col] for col in local]
                if not any(_nonempty(value) for value in local_values):
                    blank_run += 1
                    if blank_run >= 2:
                        break
                    continue
                blank_run = 0

                # Stop when a new vertical block begins in the same columns.
                if canonical_header(anchor_value, pathway_context=(kind == "pathway")) in {
                    "gene",
                    "term",
                    "pathway_id",
                }:
                    break
                if not _nonempty(anchor_value):
                    continue

                record: dict[str, Any] = {
                    "source_workbook": workbook,
                    "source_sheet": sheet,
                    "source_table_id": table_id,
                    "source_excel_row": row + 1,
                    "legacy_table_label": label,
                    **scope,
                }
                for col, key in local.items():
                    value = df.iat[row, col]
                    if key not in record or not _nonempty(record.get(key)):
                        record[key] = value
                records.append(record)
                n_rows += 1

            audit.append(
                {
                    "source_workbook": workbook,
                    "source_sheet": sheet,
                    "source_table_id": table_id,
                    "detected_kind": kind,
                    "header_excel_row": header_row + 1,
                    "anchor_excel_column": anchor + 1,
                    "left_excel_column": left + 1,
                    "right_excel_column": right + 1,
                    "legacy_table_label": label,
                    "n_rows_detected": n_rows,
                    **scope,
                }
            )
    return records, audit


_GENE_PATTERN = re.compile(r"^(?:[A-Za-z][A-Za-z0-9._-]{1,39}|[0-9]+[A-Za-z][A-Za-z0-9._-]{1,39})$")
_GENE_STOPWORDS = {
    "gene",
    "genes",
    "symbol",
    "marker",
    "markers",
    "list",
    "inkt1",
    "inkt2",
    "inkt17",
    "control",
    "ctrl",
    "t2",
    "case",
}


def looks_like_gene(value: Any) -> bool:
    if not _nonempty(value):
        return False
    text = str(value).strip()
    if text.casefold() in _GENE_STOPWORDS:
        return False
    if any(char.isspace() for char in text):
        return False
    if text.replace(".", "", 1).isdigit():
        return False
    return bool(_GENE_PATTERN.fullmatch(text))


def normalize_curated_gene(value: Any) -> tuple[str, bool]:
    source = str(value).strip()
    replacement = GENE_ALIASES.get(source.casefold())
    return (replacement or source, replacement is not None)


def extract_signature_sheet(
    df: pd.DataFrame,
    *,
    workbook: str,
    sheet: str,
) -> list[dict[str, Any]]:
    """Extract one or more gene columns from a curated signature sheet."""

    records: list[dict[str, Any]] = []
    for col in range(df.shape[1]):
        populated = [(row, df.iat[row, col]) for row in range(df.shape[0]) if _nonempty(df.iat[row, col])]
        gene_values = [(row, value) for row, value in populated if looks_like_gene(value)]
        if len(gene_values) < 2:
            continue
        first_row, first_value = populated[0]
        # A descriptive first cell becomes the column label.  With genuinely
        # headerless lists the sheet name remains the label and no gene is lost.
        if not looks_like_gene(first_value) or str(first_value).strip().casefold() in _GENE_STOPWORDS:
            gene_set = str(first_value).strip()
        else:
            gene_set = sheet if len([c for c in range(df.shape[1]) if df.iloc[:, c].notna().any()]) == 1 else f"{sheet}__col_{col + 1}"
        for row, value in gene_values:
            if row == first_row and str(value).strip().casefold() == str(gene_set).strip().casefold():
                continue
            gene, alias_applied = normalize_curated_gene(value)
            records.append(
                {
                    "library": "iNKT_gene_list_xlsx",
                    "gene_set": gene_set,
                    "gene": gene,
                    "source_gene": str(value).strip(),
                    "alias_applied": alias_applied,
                    "source_workbook": workbook,
                    "source_sheet": sheet,
                    "source_excel_row": row + 1,
                    "source_excel_column": col + 1,
                }
            )
    return records


def _coerce_numeric_columns(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = [
        "logfoldchange",
        "foldchange",
        "pvalue",
        "fdr",
        "pct_t2",
        "pct_ctrl",
        "overlap",
        "gene_set_size",
        "odds_ratio",
        "score",
        "zscore",
    ]
    for column in numeric:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def run_legacy(args: argparse.Namespace) -> None:
    out = args.outdir / "legacy"
    raw_dir = out / "raw_sheets"
    raw_dir.mkdir(parents=True, exist_ok=True)

    inventory: list[dict[str, Any]] = []
    all_detected: list[dict[str, Any]] = []
    all_audit: list[dict[str, Any]] = []
    all_signatures: list[dict[str, Any]] = []
    workbooks = [(args.gene_list_xlsx, "gene_list"), (args.marker_xlsx, "marker")]

    for path, workbook_tag in workbooks:
        if not path.exists():
            raise FileNotFoundError(path)
        sheets = read_excel_raw(path)
        for sheet, frame in sheets.items():
            raw_name = f"{slug(workbook_tag)}__{slug(sheet)}.csv"
            frame.to_csv(raw_dir / raw_name, index=False, header=False)
            nonempty = int(frame.notna().sum().sum())
            records, table_audit = detect_tables(frame, workbook=path.name, sheet=sheet)
            all_detected.extend(records)
            all_audit.extend(table_audit)

            sheet_lower = sheet.strip().casefold()
            is_known_result_sheet = (
                workbook_tag == "marker"
                or sheet_lower in {"ctrl_t2", "ctrlt2", "pathway", "c1"}
                or _sheet_context_is_pathway(sheet)
            )
            signatures = [] if is_known_result_sheet else extract_signature_sheet(
                frame, workbook=path.name, sheet=sheet
            )
            all_signatures.extend(signatures)
            inventory.append(
                {
                    "workbook_tag": workbook_tag,
                    "source_workbook": path.name,
                    "source_path": str(path),
                    "source_sheet": sheet,
                    "n_rows_raw": int(frame.shape[0]),
                    "n_columns_raw": int(frame.shape[1]),
                    "n_nonempty_cells": nonempty,
                    "n_detected_tables": len(table_audit),
                    "n_signature_gene_rows": len(signatures),
                    "raw_export": str((raw_dir / raw_name).relative_to(args.outdir)),
                }
            )

    detected = pd.DataFrame(all_detected)
    if detected.empty:
        detected = pd.DataFrame(columns=["source_workbook", "source_sheet", "source_table_id"])
    detected = _coerce_numeric_columns(detected)

    de = detected[detected.get("gene", pd.Series(index=detected.index, dtype=object)).notna()].copy()
    if not de.empty:
        de["gene"] = de["gene"].astype(str).str.strip()
        de = de[de["gene"].map(looks_like_gene)].copy()
        if "logfoldchange" in de:
            de["direction_as_stored"] = np.select(
                [de["logfoldchange"] > 0, de["logfoldchange"] < 0],
                ["positive", "negative"],
                default="zero_or_unknown",
            )
        else:
            de["direction_as_stored"] = "unknown"
        de["contrast"] = "legacy_ctrl_T2_workbook_direction_requires_manual_confirmation"
        de["result_namespace"] = "legacy_xlsx_ppt"

    pathway_anchor = detected.get("term", pd.Series(index=detected.index, dtype=object)).notna()
    pathway_anchor |= detected.get("pathway_id", pd.Series(index=detected.index, dtype=object)).notna()
    pathways = detected[pathway_anchor].copy()
    if not pathways.empty:
        pathways["result_namespace"] = "legacy_xlsx_ppt"
        pathways["rerun_status"] = "recovered_result_only_not_recomputed"

    signatures = pd.DataFrame(all_signatures)
    if signatures.empty:
        signatures = pd.DataFrame(
            columns=["library", "gene_set", "gene", "source_gene", "alias_applied", "source_workbook", "source_sheet"]
        )
    else:
        signatures = signatures.drop_duplicates(["library", "gene_set", "gene"]).sort_values(
            ["gene_set", "gene"], kind="stable"
        )

    pd.DataFrame(inventory).to_csv(out / "workbook_inventory.csv", index=False)
    pd.DataFrame(all_audit).to_csv(out / "detected_table_audit.csv", index=False)
    de.to_csv(out / "legacy_de_candidates.csv", index=False)
    pathways.to_csv(out / "legacy_pathway_candidates.csv", index=False)
    signatures.to_csv(out / "curated_gene_sets.csv", index=False)

    unresolved = pd.DataFrame(all_audit)
    unresolved_count = (
        int((unresolved["legacy_scope"] == "legacy_unresolved").sum())
        if not unresolved.empty and "legacy_scope" in unresolved
        else 0
    )
    write_json(
        out / "legacy_manifest.json",
        {
            "task": "legacy",
            "workbooks": [str(path) for path, _ in workbooks],
            "n_sheets": len(inventory),
            "n_detected_tables": len(all_audit),
            "n_legacy_de_rows": int(len(de)),
            "n_legacy_pathway_rows": int(len(pathways)),
            "n_curated_gene_sets": int(signatures["gene_set"].nunique()) if not signatures.empty else 0,
            "n_curated_gene_rows": int(len(signatures)),
            "n_unresolved_table_scopes": unresolved_count,
            "legacy_cluster_namespace": LEGACY_CLUSTER_NAMESPACE,
            "important_caveats": [
                "Raw sheet CSV exports are the source-of-truth audit trail for the adaptive parser.",
                (
                    "Legacy c1-c11 identifiers are never numerically mapped to the rerun "
                    f"cluster namespace {args.cluster_namespace}."
                ),
                "The sign convention of stored legacy fold changes requires confirmation from the original analysis notes.",
                "Unlabelled marker.xlsx pathway blocks are retained with unresolved scope rather than guessed.",
            ],
        },
    )


_WORKER_ADATA: Any = None


def _init_de_worker(h5ad_path: str) -> None:
    global _WORKER_ADATA
    import scanpy as sc

    sc.settings.verbosity = 0
    _WORKER_ADATA = sc.read_h5ad(h5ad_path)


def _matrix_summary(matrix: Any, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    selected = matrix[mask]
    means = np.asarray(selected.mean(axis=0)).ravel()
    if sparse.issparse(selected):
        fractions = np.asarray((selected > 0).mean(axis=0)).ravel()
    else:
        fractions = np.mean(np.asarray(selected) > 0, axis=0)
    return means, fractions


def _run_de_unit(
    unit_payload: Mapping[str, Any],
    cluster_key: str,
    cluster_namespace: str,
    min_cells: int,
    unit_dir: str,
) -> dict[str, Any]:
    import scanpy as sc

    unit = DEUnit(**unit_payload)
    adata = _WORKER_ADATA
    obs = adata.obs
    mask = obs["tissue"].astype(str).to_numpy() == unit.tissue
    if unit.cluster_current is not None:
        mask &= obs[cluster_key].astype(str).to_numpy() == unit.cluster_current
    sub = adata[mask].copy()
    condition = sub.obs["condition"].astype(str).to_numpy()
    n_t2 = int(np.sum(condition == "T2"))
    n_ctrl = int(np.sum(condition == "Ctrl"))
    result_path = Path(unit_dir) / f"{slug(unit.unit_id)}.csv"

    base = {
        **asdict(unit),
        "n_t2": n_t2,
        "n_ctrl": n_ctrl,
        "output_csv": str(result_path),
    }
    if min(n_t2, n_ctrl) < min_cells:
        return {**base, "status": "skipped_insufficient_cells"}
    if sub.raw is None:
        raise RuntimeError("Current H5AD has no .raw; full-gene log-normalized expression is required")

    sub.obs["_de_condition"] = pd.Categorical(condition, categories=["Ctrl", "T2"])
    key = "_rank_t2_vs_ctrl"
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
    frame = sc.get.rank_genes_groups_df(sub, group="T2", key=key)
    frame = frame.rename(columns={"names": "gene"})

    raw_matrix = sub.raw.X
    raw_names = pd.Index(sub.raw.var_names.astype(str))
    mean_t2, pct_t2 = _matrix_summary(raw_matrix, condition == "T2")
    mean_ctrl, pct_ctrl = _matrix_summary(raw_matrix, condition == "Ctrl")
    summaries = pd.DataFrame(
        {
            "gene": raw_names,
            "mean_log1p_expression_t2": mean_t2,
            "mean_log1p_expression_ctrl": mean_ctrl,
            "pct_expressing_t2": pct_t2,
            "pct_expressing_ctrl": pct_ctrl,
        }
    )
    frame = frame.merge(summaries, on="gene", how="left", validate="one_to_one")
    frame.insert(0, "unit_id", unit.unit_id)
    frame.insert(1, "scope", unit.scope)
    frame.insert(2, "tissue", unit.tissue)
    frame.insert(3, "cluster_current", unit.cluster_current)
    frame.insert(4, "cluster_namespace", cluster_namespace if unit.cluster_current is not None else "current_no_cluster")
    frame.insert(5, "contrast", "T2_vs_Ctrl")
    frame.insert(6, "n_t2", n_t2)
    frame.insert(7, "n_ctrl", n_ctrl)
    frame["direction"] = np.select(
        [frame["logfoldchanges"] > 0, frame["logfoldchanges"] < 0],
        ["up_in_T2", "down_in_T2"],
        default="zero_or_unknown",
    )
    frame["statistical_unit_warning"] = STATISTICAL_UNIT_WARNING
    frame["biological_replicates_per_condition_within_tissue"] = 1
    frame = frame.sort_values(["pvals_adj", "pvals"], na_position="last", kind="stable")
    result_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(result_path, index=False)
    return {**base, "status": "completed", "n_genes": int(len(frame))}


def discover_de_units(obs: pd.DataFrame, cluster_key: str) -> list[DEUnit]:
    tissues = sorted({str(value) for value in obs["tissue"].dropna().astype(str)})
    units = [DEUnit(unit_id=f"tissue__{tissue}", scope="current_tissue", tissue=tissue) for tissue in tissues]
    pairs = (
        obs[["tissue", cluster_key]]
        .astype(str)
        .drop_duplicates()
        .sort_values(["tissue", cluster_key], kind="stable")
    )
    for tissue, cluster in pairs.itertuples(index=False, name=None):
        units.append(
            DEUnit(
                unit_id=f"cluster_tissue__c{cluster}__{tissue}",
                scope="current_cluster_tissue",
                tissue=str(tissue),
                cluster_current=str(cluster),
            )
        )
    return units


def run_current_de(args: argparse.Namespace) -> None:
    import anndata as ad

    if not args.h5ad.exists():
        raise FileNotFoundError(args.h5ad)
    out = args.outdir / "current_de"
    unit_dir = out / "unit_csv"
    unit_dir.mkdir(parents=True, exist_ok=True)

    # Read only obs to plan work.  Each worker then loads one complete copy and
    # reuses it for multiple comparisons.
    backed = ad.read_h5ad(args.h5ad, backed="r")
    required = {"condition", "tissue", args.cluster_key}
    missing = required - set(backed.obs.columns)
    if missing:
        backed.file.close()
        raise KeyError(f"Missing required obs columns: {sorted(missing)}")
    obs = backed.obs[["condition", "tissue", args.cluster_key]].copy()
    has_raw = backed.raw is not None
    has_counts = "counts" in backed.layers
    shape = [int(backed.n_obs), int(backed.n_vars)]
    backed.file.close()
    if not has_raw:
        raise RuntimeError("H5AD .raw is required for full-gene DE")

    all_units = discover_de_units(obs, args.cluster_key)
    if args.shard_count > 1:
        units = [unit for index, unit in enumerate(all_units) if index % args.shard_count == args.shard_index]
    else:
        units = all_units
    if not units:
        raise RuntimeError("No DE units selected by this shard")

    workers = max(1, min(args.n_jobs, len(units)))
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_de_worker,
        initargs=(str(args.h5ad),),
    ) as pool:
        futures = {
            pool.submit(
                _run_de_unit,
                asdict(unit),
                args.cluster_key,
                args.cluster_namespace,
                args.min_cells_per_group,
                str(unit_dir),
            ): unit
            for unit in units
        }
        for future in as_completed(futures):
            unit = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    **asdict(unit),
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            results.append(result)
            print(json.dumps({key: _json_safe_number(value) for key, value in result.items()}, sort_keys=True), flush=True)

    results.sort(key=lambda row: row["unit_id"])
    pd.DataFrame(results).to_csv(
        out / f"current_de_unit_status_shard-{args.shard_index}-of-{args.shard_count}.csv",
        index=False,
    )
    completed_paths = [Path(row["output_csv"]) for row in results if row["status"] == "completed"]
    frames = [pd.read_csv(path, low_memory=False) for path in completed_paths]
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not combined.empty:
        combined[combined["scope"] == "current_tissue"].to_csv(out / f"current_de_tissue_shard-{args.shard_index}-of-{args.shard_count}.csv", index=False)
        combined[combined["scope"] == "current_cluster_tissue"].to_csv(
            out / f"current_de_cluster_tissue_shard-{args.shard_index}-of-{args.shard_count}.csv", index=False
        )

    manifest = {
        "task": "current-de",
        "h5ad": str(args.h5ad),
        "h5ad_shape": shape,
        "raw_available": has_raw,
        "counts_layer_available": has_counts,
        "cluster_key": args.cluster_key,
        "cluster_namespace": args.cluster_namespace,
        "n_units_all": len(all_units),
        "n_units_this_shard": len(units),
        "n_workers": workers,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "min_cells_per_group": args.min_cells_per_group,
        "completed": int(sum(row["status"] == "completed" for row in results)),
        "skipped": int(sum(row["status"].startswith("skipped") for row in results)),
        "failed": int(sum(row["status"] == "failed" for row in results)),
        "statistical_unit_warning": STATISTICAL_UNIT_WARNING,
    }
    write_json(out / f"current_de_manifest_shard-{args.shard_index}-of-{args.shard_count}.json", manifest)
    if manifest["failed"]:
        raise RuntimeError(f"{manifest['failed']} DE units failed; inspect unit status CSV")


def load_current_de(outdir: Path) -> pd.DataFrame:
    unit_dir = outdir / "current_de/unit_csv"
    paths = sorted(unit_dir.glob("*.csv"))
    if not paths:
        raise FileNotFoundError(f"No current DE unit CSVs found in {unit_dir}")
    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    combined = pd.concat(frames, ignore_index=True)
    if "cluster_current" in combined:
        combined["cluster_current"] = combined["cluster_current"].map(
            lambda value: None if pd.isna(value) else str(int(float(value)))
        )
    return combined


def significant_gene_sets(
    de: pd.DataFrame,
    *,
    fdr: float,
    min_abs_logfc: float,
) -> dict[tuple[str, str], set[str]]:
    selected = de[
        (pd.to_numeric(de["pvals_adj"], errors="coerce") <= fdr)
        & (pd.to_numeric(de["logfoldchanges"], errors="coerce").abs() >= min_abs_logfc)
    ].copy()
    selected["direction_simple"] = np.where(selected["logfoldchanges"] > 0, "up", "down")
    return {
        (str(unit), str(direction)): set(group["gene"].dropna().astype(str))
        for (unit, direction), group in selected.groupby(["unit_id", "direction_simple"], observed=True)
    }


def make_overlap_outputs(args: argparse.Namespace, de: pd.DataFrame, output: Path) -> None:
    sets = significant_gene_sets(de, fdr=args.de_fdr, min_abs_logfc=args.de_min_abs_logfc)
    tissue_rows = de[de["scope"] == "current_tissue"][["unit_id", "tissue"]].drop_duplicates()
    cluster_rows = de[de["scope"] == "current_cluster_tissue"][
        ["unit_id", "tissue", "cluster_current"]
    ].drop_duplicates()
    summaries: list[dict[str, Any]] = []
    memberships: list[dict[str, Any]] = []

    tissue_unit_by_tissue = dict(tissue_rows[["tissue", "unit_id"]].itertuples(index=False, name=None))
    for cluster_unit, tissue, cluster in cluster_rows.itertuples(index=False, name=None):
        tissue_unit = tissue_unit_by_tissue.get(tissue)
        if tissue_unit is None:
            continue
        for direction in ("up", "down"):
            tissue_set = sets.get((str(tissue_unit), direction), set())
            cluster_set = sets.get((str(cluster_unit), direction), set())
            shared = tissue_set & cluster_set
            union = tissue_set | cluster_set
            denom = min(len(tissue_set), len(cluster_set))
            summaries.append(
                {
                    "tissue": tissue,
                    "cluster_current": cluster,
                    "cluster_namespace": args.cluster_namespace,
                    "direction": direction,
                    "n_tissue_significant": len(tissue_set),
                    "n_cluster_tissue_significant": len(cluster_set),
                    "n_shared": len(shared),
                    "n_tissue_only": len(tissue_set - cluster_set),
                    "n_cluster_tissue_only": len(cluster_set - tissue_set),
                    "jaccard": len(shared) / len(union) if union else np.nan,
                    "overlap_coefficient": len(shared) / denom if denom else np.nan,
                }
            )
            for gene in sorted(union):
                memberships.append(
                    {
                        "gene": gene,
                        "tissue": tissue,
                        "cluster_current": cluster,
                        "cluster_namespace": args.cluster_namespace,
                        "direction": direction,
                        "membership": (
                            "shared"
                            if gene in shared
                            else "tissue_only"
                            if gene in tissue_set
                            else "cluster_tissue_only"
                        ),
                    }
                )

    summary = pd.DataFrame(summaries)
    membership = pd.DataFrame(memberships)
    summary.to_csv(output / "current_tissue_vs_cluster_overlap_summary.csv", index=False)
    membership.to_csv(output / "current_tissue_vs_cluster_overlap_membership.csv", index=False)

    if not summary.empty:
        summary["row"] = "c" + summary["cluster_current"].astype(str) + "__" + summary["direction"]
        heat = summary.pivot(index="row", columns="tissue", values="jaccard")
        fig, ax = plt.subplots(figsize=(max(6, 1.7 * heat.shape[1]), max(5, 0.33 * heat.shape[0])))
        image = ax.imshow(heat.fillna(0).to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(np.arange(heat.shape[1]), heat.columns, rotation=35, ha="right")
        ax.set_yticks(np.arange(heat.shape[0]), heat.index)
        ax.set_title("Current tissue DEG vs cluster×tissue DEG: Jaccard")
        for row in range(heat.shape[0]):
            for col in range(heat.shape[1]):
                value = heat.iloc[row, col]
                if np.isfinite(value):
                    ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=7, color="white" if value > 0.45 else "black")
        fig.colorbar(image, ax=ax, label="Jaccard")
        fig.tight_layout()
        fig.savefig(output / "current_tissue_vs_cluster_overlap_jaccard.png", dpi=180)
        plt.close(fig)

    # UpSet-like Boolean membership matrices and top intersection-pattern plots
    # are generated without adding a dependency on upsetplot.
    for tissue in sorted(tissue_unit_by_tissue):
        comparisons: dict[str, set[str]] = {}
        tissue_unit = str(tissue_unit_by_tissue[tissue])
        for direction in ("up", "down"):
            comparisons[f"tissue_{direction}"] = sets.get((tissue_unit, direction), set())
        subset = cluster_rows[cluster_rows["tissue"] == tissue]
        for unit, _, cluster in subset.itertuples(index=False, name=None):
            for direction in ("up", "down"):
                comparisons[f"c{cluster}_{direction}"] = sets.get((str(unit), direction), set())
        universe = sorted(set().union(*comparisons.values())) if comparisons else []
        matrix = pd.DataFrame(
            {name: [gene in genes for gene in universe] for name, genes in comparisons.items()},
            index=pd.Index(universe, name="gene"),
        )
        matrix.to_csv(output / f"current_upset_membership__{slug(tissue)}.csv")
        if matrix.empty:
            continue
        pattern = matrix.astype(int).astype(str).agg("".join, axis=1)
        counts = pattern.value_counts().head(args.top_patterns).rename_axis("binary_pattern").reset_index(name="n_genes")
        labels = []
        columns = list(matrix.columns)
        for binary in counts["binary_pattern"]:
            members = [columns[index] for index, flag in enumerate(binary) if flag == "1"]
            labels.append(" & ".join(members) if members else "none")
        counts["members"] = labels
        counts.to_csv(output / f"current_upset_top_patterns__{slug(tissue)}.csv", index=False)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.36 * len(counts))))
        ax.barh(np.arange(len(counts)), counts["n_genes"].to_numpy(), color="#4C78A8")
        ax.set_yticks(np.arange(len(counts)), counts["members"], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("genes")
        ax.set_title(f"Top current DEG intersection patterns: {tissue}")
        fig.tight_layout()
        fig.savefig(output / f"current_upset_top_patterns__{slug(tissue)}.png", dpi=180)
        plt.close(fig)


def read_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = {gene for gene in fields[2:] if gene}
    return gene_sets


def load_gene_set_libraries(args: argparse.Namespace) -> tuple[dict[str, dict[str, set[str]]], list[dict[str, Any]]]:
    libraries: dict[str, dict[str, set[str]]] = {}
    audit: list[dict[str, Any]] = []
    curated_path = args.outdir / "legacy/curated_gene_sets.csv"
    if curated_path.exists():
        curated = pd.read_csv(curated_path)
        sets = {
            str(name): set(group["gene"].dropna().astype(str))
            for name, group in curated.groupby("gene_set", observed=True)
        }
        libraries["curated_xlsx_signatures"] = sets
        audit.append({"library": "curated_xlsx_signatures", "path": str(curated_path), "n_gene_sets": len(sets), "status": "available"})
    else:
        audit.append({"library": "curated_xlsx_signatures", "path": str(curated_path), "n_gene_sets": 0, "status": "missing_run_legacy_first"})

    for gmt in args.gmt:
        if gmt.exists():
            name = f"gmt__{gmt.stem}"
            sets = read_gmt(gmt)
            libraries[name] = sets
            audit.append({"library": name, "path": str(gmt), "n_gene_sets": len(sets), "status": "available"})
        else:
            audit.append({"library": f"gmt__{gmt.stem}", "path": str(gmt), "n_gene_sets": 0, "status": "missing"})
    return libraries, audit


def ora_for_current_de(
    args: argparse.Namespace,
    de: pd.DataFrame,
    libraries: Mapping[str, Mapping[str, set[str]]],
) -> pd.DataFrame:
    universe_original = set(de["gene"].dropna().astype(str))
    universe_case = {gene.casefold(): gene for gene in universe_original}
    sets = significant_gene_sets(de, fdr=args.de_fdr, min_abs_logfc=args.de_min_abs_logfc)
    unit_meta = de[["unit_id", "scope", "tissue", "cluster_current", "cluster_namespace"]].drop_duplicates("unit_id")
    rows: list[dict[str, Any]] = []

    for meta in unit_meta.to_dict(orient="records"):
        unit_id = str(meta["unit_id"])
        for direction in ("up", "down"):
            selected = {gene.casefold() for gene in sets.get((unit_id, direction), set()) if gene.casefold() in universe_case}
            for library, gene_sets in libraries.items():
                for term, genes in gene_sets.items():
                    members = {gene.casefold() for gene in genes if gene.casefold() in universe_case}
                    if len(members) < args.min_gene_set_size or len(members) > args.max_gene_set_size:
                        continue
                    overlap = selected & members
                    pvalue = hypergeom.sf(
                        len(overlap) - 1,
                        len(universe_case),
                        len(members),
                        len(selected),
                    ) if selected else 1.0
                    rows.append(
                        {
                            **meta,
                            "direction": direction,
                            "library": library,
                            "term": term,
                            "universe_size": len(universe_case),
                            "selected_de_size": len(selected),
                            "gene_set_size_in_universe": len(members),
                            "overlap_size": len(overlap),
                            "overlap_genes": ";".join(sorted(universe_case[gene] for gene in overlap)),
                            "pvalue": float(pvalue),
                        }
                    )
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(
            columns=["unit_id", "scope", "tissue", "cluster_current", "direction", "library", "term", "pvalue", "fdr"]
        )
    result["fdr"] = result.groupby(["unit_id", "direction", "library"], observed=True)["pvalue"].transform(
        lambda values: bh_adjust(values.to_numpy())
    )
    return result.sort_values(["unit_id", "direction", "library", "fdr", "pvalue"], kind="stable")


def make_legacy_current_tissue_overlap(args: argparse.Namespace, current: pd.DataFrame, output: Path) -> None:
    legacy_path = args.outdir / "legacy/legacy_de_candidates.csv"
    columns = [
        "tissue",
        "n_legacy_genes_any_direction",
        "n_current_genes_significant",
        "n_shared_any_direction",
        "jaccard_any_direction",
        "important_note",
    ]
    if not legacy_path.exists():
        pd.DataFrame(columns=columns).to_csv(output / "legacy_vs_current_tissue_gene_overlap.csv", index=False)
        return
    legacy = pd.read_csv(legacy_path, low_memory=False)
    if legacy.empty or "gene" not in legacy or "tissue" not in legacy:
        pd.DataFrame(columns=columns).to_csv(output / "legacy_vs_current_tissue_gene_overlap.csv", index=False)
        return
    current_sets = significant_gene_sets(current, fdr=args.de_fdr, min_abs_logfc=args.de_min_abs_logfc)
    tissue_units = current[current["scope"] == "current_tissue"][["tissue", "unit_id"]].drop_duplicates()
    rows = []
    membership = []
    for tissue, unit in tissue_units.itertuples(index=False, name=None):
        old = set(legacy.loc[legacy["tissue"] == tissue, "gene"].dropna().astype(str))
        new = current_sets.get((str(unit), "up"), set()) | current_sets.get((str(unit), "down"), set())
        shared = old & new
        union = old | new
        rows.append(
            {
                "tissue": tissue,
                "n_legacy_genes_any_direction": len(old),
                "n_current_genes_significant": len(new),
                "n_shared_any_direction": len(shared),
                "jaccard_any_direction": len(shared) / len(union) if union else np.nan,
                "important_note": "direction_agnostic_because_legacy_fold_change_sign_convention_requires_confirmation",
            }
        )
        for gene in sorted(union):
            membership.append(
                {
                    "tissue": tissue,
                    "gene": gene,
                    "membership": "shared" if gene in shared else "legacy_only" if gene in old else "current_only",
                }
            )
    pd.DataFrame(rows, columns=columns).to_csv(output / "legacy_vs_current_tissue_gene_overlap.csv", index=False)
    pd.DataFrame(membership).to_csv(output / "legacy_vs_current_tissue_gene_overlap_membership.csv", index=False)


def run_overlap_pathway(args: argparse.Namespace) -> None:
    output = args.outdir / "overlap_pathway"
    output.mkdir(parents=True, exist_ok=True)
    de = load_current_de(args.outdir)
    make_overlap_outputs(args, de, output)
    make_legacy_current_tissue_overlap(args, de, output)

    libraries, library_audit = load_gene_set_libraries(args)
    ora = ora_for_current_de(args, de, libraries)
    ora.to_csv(output / "current_offline_gene_set_enrichment.csv", index=False)
    pd.DataFrame(library_audit).to_csv(output / "offline_gene_set_library_audit.csv", index=False)

    pathway_libraries = [row for row in library_audit if str(row["library"]).startswith("gmt__") and row["status"] == "available"]
    write_json(
        output / "overlap_pathway_manifest.json",
        {
            "task": "overlap-pathway",
            "de_fdr_threshold": args.de_fdr,
            "de_min_abs_logfc": args.de_min_abs_logfc,
            "n_current_de_rows": int(len(de)),
            "n_enrichment_rows": int(len(ora)),
            "gene_set_library_audit": library_audit,
            "current_pathway_rerun_status": (
                "recomputed_with_supplied_offline_gmt"
                if pathway_libraries
                else "no_offline_pathway_gmt_available; curated_XLSX_signature_ORA_only"
            ),
            "legacy_pathway_status": "recovered_from_XLSX_by_legacy_task; original_unlabelled_blocks_not_reassigned",
            "cluster_namespace_rule": (
                f"{LEGACY_CLUSTER_NAMESPACE} is never numerically matched to {args.cluster_namespace}"
            ),
            "statistical_unit_warning": STATISTICAL_UNIT_WARNING,
        },
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task",
        choices=["legacy", "current-de", "overlap-pathway", "all"],
        default="all",
    )
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--gene-list-xlsx", type=Path, default=DEFAULT_GENE_LIST)
    parser.add_argument("--marker-xlsx", type=Path, default=DEFAULT_MARKER)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cluster-key", default="leiden_res_0_5")
    parser.add_argument(
        "--cluster-namespace",
        default=CURRENT_CLUSTER_NAMESPACE,
        help="Provenance label written beside rerun cluster identifiers.",
    )
    parser.add_argument("--n-jobs", type=int, default=min(12, os.cpu_count() or 1))
    parser.add_argument("--min-cells-per-group", type=int, default=20)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--de-fdr", type=float, default=0.05)
    parser.add_argument("--de-min-abs-logfc", type=float, default=0.25)
    parser.add_argument("--min-gene-set-size", type=int, default=3)
    parser.add_argument("--max-gene-set-size", type=int, default=5000)
    parser.add_argument("--top-patterns", type=int, default=20)
    parser.add_argument(
        "--gmt",
        action="append",
        type=Path,
        default=[],
        help="Optional offline GMT; repeat for multiple pathway libraries.",
    )
    args = parser.parse_args(argv)
    if args.n_jobs < 1:
        parser.error("--n-jobs must be >= 1")
    if args.shard_count < 1:
        parser.error("--shard-count must be >= 1")
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must satisfy 0 <= index < count")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.task in {"legacy", "all"}:
        run_legacy(args)
    if args.task in {"current-de", "all"}:
        run_current_de(args)
    if args.task in {"overlap-pathway", "all"}:
        run_overlap_pathway(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
