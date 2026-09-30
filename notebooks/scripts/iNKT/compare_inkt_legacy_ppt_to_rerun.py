from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[3]
A_NS = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main"}

DE_SPECS = [
    {"page": 10, "table": 0, "legacy_unit": "global", "scope": "global", "tissue": None},
    {"page": 11, "table": 0, "legacy_unit": "tissue_thymus", "scope": "tissue", "tissue": "thymus"},
    {"page": 11, "table": 1, "legacy_unit": "tissue_bone_marrow", "scope": "tissue", "tissue": "bone_marrow"},
    {"page": 11, "table": 2, "legacy_unit": "tissue_spleen", "scope": "tissue", "tissue": "spleen"},
    {"page": 13, "table": 0, "legacy_unit": "c1_thymus", "scope": "cluster_tissue", "tissue": "thymus"},
    {"page": 15, "table": 0, "legacy_unit": "c2_bone_marrow", "scope": "cluster_tissue", "tissue": "bone_marrow"},
    {"page": 15, "table": 1, "legacy_unit": "c2_spleen", "scope": "cluster_tissue", "tissue": "spleen"},
    {"page": 17, "table": 0, "legacy_unit": "c5_bone_marrow", "scope": "cluster_tissue", "tissue": "bone_marrow"},
    {"page": 19, "table": 0, "legacy_unit": "c6_spleen", "scope": "cluster_tissue", "tissue": "spleen"},
    {"page": 21, "table": 0, "legacy_unit": "c7_bone_marrow", "scope": "cluster_tissue", "tissue": "bone_marrow"},
    {"page": 22, "table": 0, "legacy_unit": "c8_spleen", "scope": "cluster_tissue", "tissue": "spleen"},
    # Slide XML table order differs from visual order: table 0 is spleen, table 1 is BM.
    {"page": 24, "table": 0, "legacy_unit": "c9_spleen", "scope": "cluster_tissue", "tissue": "spleen"},
    {"page": 24, "table": 1, "legacy_unit": "c9_bone_marrow", "scope": "cluster_tissue", "tissue": "bone_marrow"},
]

PATHWAY_SPECS = [
    {"page": 14, "table": 0, "legacy_unit": "c1_thymus"},
    {"page": 16, "table": 0, "legacy_unit": "c2_bone_marrow"},
    {"page": 16, "table": 1, "legacy_unit": "c2_spleen"},
    {"page": 18, "table": 0, "legacy_unit": "c5_bone_marrow"},
    {"page": 20, "table": 0, "legacy_unit": "c6_spleen"},
    {"page": 23, "table": 0, "legacy_unit": "c8_spleen"},
    # Slide 25 table 0 is spleen and table 1 is BM despite the incomplete page title.
    {"page": 25, "table": 0, "legacy_unit": "c9_spleen"},
    {"page": 25, "table": 1, "legacy_unit": "c9_bone_marrow"},
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare legacy PPT tables with a completed iNKT rerun.")
    parser.add_argument("--pptx", type=Path, default=ROOT / "input/iNKT/iNKT.pptx")
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-abs-logfc", type=float, default=0.25)
    return parser.parse_args()


def slide_tables(archive: zipfile.ZipFile, page: int) -> list[pd.DataFrame]:
    root = ET.fromstring(archive.read(f"ppt/slides/slide{page}.xml"))
    tables: list[pd.DataFrame] = []
    for table in root.findall(".//a:tbl", A_NS):
        rows: list[list[str]] = []
        for row in table.findall("./a:tr", A_NS):
            rows.append(
                [
                    "".join(node.text or "" for node in cell.findall(".//a:t", A_NS)).strip()
                    for cell in row.findall("./a:tc", A_NS)
                ]
            )
        if not rows:
            continue
        width = max(map(len, rows))
        rows = [row + [""] * (width - len(row)) for row in rows]
        header = [value or f"column_{index}" for index, value in enumerate(rows[0])]
        tables.append(pd.DataFrame(rows[1:], columns=header))
    return tables


def significant(frame: pd.DataFrame, fdr: float, min_abs_logfc: float) -> pd.DataFrame:
    return frame.loc[
        (pd.to_numeric(frame["pvals_adj"], errors="coerce") <= fdr)
        & (pd.to_numeric(frame["logfoldchanges"], errors="coerce").abs() >= min_abs_logfc)
    ].copy()


def comparison_row(spec: dict, legacy: pd.DataFrame, current: pd.DataFrame, unit_id: str) -> dict:
    legacy = legacy.copy()
    legacy["gene_key"] = legacy["names"].astype(str).str.casefold()
    current = current.copy()
    current["gene_key"] = current["gene"].astype(str).str.casefold()
    legacy_by_key = legacy.drop_duplicates("gene_key").set_index("gene_key")
    current_by_key = current.drop_duplicates("gene_key").set_index("gene_key")
    shared = sorted(set(legacy_by_key.index) & set(current_by_key.index))
    union = set(legacy_by_key.index) | set(current_by_key.index)
    sign_total = 0
    sign_agree = 0
    for key in shared:
        old_sign = np.sign(float(legacy_by_key.loc[key, "logFC"]))
        new_sign = np.sign(float(current_by_key.loc[key, "logfoldchanges"]))
        if old_sign and new_sign:
            sign_total += 1
            sign_agree += int(old_sign == new_sign)
    first = current.iloc[0] if len(current) else pd.Series(dtype=object)
    return {
        **spec,
        "current_unit_id": unit_id,
        "current_cluster": first.get("cluster_current"),
        "n_t2": int(first.get("n_t2")) if pd.notna(first.get("n_t2")) else None,
        "n_ctrl": int(first.get("n_ctrl")) if pd.notna(first.get("n_ctrl")) else None,
        "n_legacy_genes": int(legacy_by_key.shape[0]),
        "n_current_significant": int(current_by_key.shape[0]),
        "n_shared": len(shared),
        "legacy_recovery": len(shared) / len(legacy_by_key) if len(legacy_by_key) else np.nan,
        "jaccard": len(shared) / len(union) if union else np.nan,
        "sign_agree": sign_agree,
        "sign_total": sign_total,
        "sign_agreement": sign_agree / sign_total if sign_total else np.nan,
        "shared_genes": ";".join(legacy_by_key.loc[shared, "names"].astype(str)) if shared else "",
        "legacy_only_genes": ";".join(
            legacy_by_key.loc[sorted(set(legacy_by_key.index) - set(current_by_key.index)), "names"].astype(str)
        ),
    }


def main() -> int:
    args = parse_args()
    run_root = args.run_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    unit_dir = run_root / "extended/de_pathway/current_de/unit_csv"

    current_units: dict[str, pd.DataFrame] = {}
    for path in sorted(unit_dir.glob("*.csv")):
        frame = pd.read_csv(path, low_memory=False)
        current_units[str(frame.iloc[0]["unit_id"])] = significant(frame, args.fdr, args.min_abs_logfc)

    h5ad = run_root / "preprocess/inkt_scanpy_tutorial_processed.h5ad"
    adata = sc.read_h5ad(h5ad)
    global_de = sc.get.rank_genes_groups_df(
        adata,
        group="T2",
        key="rank_genes_condition_t2_vs_ctrl",
    ).rename(columns={"names": "gene"})
    global_de["unit_id"] = "global"
    global_de["cluster_current"] = np.nan
    global_de["n_t2"] = int((adata.obs["condition"].astype(str) == "T2").sum())
    global_de["n_ctrl"] = int((adata.obs["condition"].astype(str) == "Ctrl").sum())
    current_units["global"] = significant(global_de, args.fdr, args.min_abs_logfc)

    with zipfile.ZipFile(args.pptx) as archive:
        ppt_tables = {
            page: slide_tables(archive, page)
            for page in sorted({spec["page"] for spec in DE_SPECS + PATHWAY_SPECS})
        }

    all_rows: list[dict] = []
    extracted_de: list[pd.DataFrame] = []
    for spec in DE_SPECS:
        legacy = ppt_tables[spec["page"]][spec["table"]].copy()
        legacy["slide"] = spec["page"]
        legacy["legacy_unit"] = spec["legacy_unit"]
        extracted_de.append(legacy)
        if spec["scope"] == "global":
            candidates = ["global"]
        elif spec["scope"] == "tissue":
            candidates = [f"tissue__{spec['tissue']}"]
        else:
            suffix = f"__{spec['tissue']}"
            candidates = sorted(
                unit for unit in current_units if unit.startswith("cluster_tissue__") and unit.endswith(suffix)
            )
        for unit_id in candidates:
            all_rows.append(comparison_row(spec, legacy, current_units[unit_id], unit_id))

    comparisons = pd.DataFrame(all_rows)
    comparisons.to_csv(output / "legacy_ppt_de_all_matches.csv", index=False)
    best_rows: list[pd.Series] = []
    for _, group in comparisons.groupby("legacy_unit", sort=False):
        best_rows.append(
            group.sort_values(
                ["n_shared", "legacy_recovery", "jaccard", "n_current_significant"],
                ascending=[False, False, False, True],
                kind="stable",
            ).iloc[0]
        )
    best = pd.DataFrame(best_rows).reset_index(drop=True)
    best.to_csv(output / "legacy_ppt_de_best_matches.csv", index=False)
    pd.concat(extracted_de, ignore_index=True).to_csv(output / "legacy_ppt_de_tables_extracted.csv", index=False)

    best_by_unit = best.set_index("legacy_unit")
    pathway_rows: list[dict] = []
    extracted_pathway: list[pd.DataFrame] = []
    measured = {str(gene).casefold(): str(gene) for gene in adata.var_names}
    for spec in PATHWAY_SPECS:
        table = ppt_tables[spec["page"]][spec["table"]].copy()
        table["slide"] = spec["page"]
        table["legacy_unit"] = spec["legacy_unit"]
        extracted_pathway.append(table)
        drivers = sorted(
            {
                token.strip()
                for value in table["GENE_SYM"].dropna().astype(str)
                for token in value.split(",")
                if token.strip()
            },
            key=str.casefold,
        )
        match = best_by_unit.loc[spec["legacy_unit"]]
        unit_id = str(match["current_unit_id"])
        current = current_units[unit_id]
        current_keys = set(current["gene"].astype(str).str.casefold())
        driver_keys = {gene.casefold() for gene in drivers}
        retained_keys = sorted(driver_keys & current_keys)
        missing_keys = sorted(driver_keys - current_keys)
        pathway_rows.append(
            {
                **spec,
                "current_unit_id": unit_id,
                "current_cluster": match["current_cluster"],
                "n_legacy_terms": len(table),
                "n_unique_drivers": len(driver_keys),
                "n_drivers_measured": len(driver_keys & set(measured)),
                "n_drivers_significant": len(retained_keys),
                "driver_retention": len(retained_keys) / len(driver_keys) if driver_keys else np.nan,
                "significant_drivers": ";".join(measured.get(key, key) for key in retained_keys),
                "not_significant_or_unmeasured_drivers": ";".join(
                    measured.get(key, key) for key in missing_keys
                ),
            }
        )
    pd.DataFrame(pathway_rows).to_csv(output / "legacy_ppt_pathway_driver_retention.csv", index=False)
    pd.concat(extracted_pathway, ignore_index=True).to_csv(
        output / "legacy_ppt_pathway_tables_extracted.csv", index=False
    )

    cluster_key = "leiden_res_0_5"
    composition_rows: list[dict] = []
    for cluster, obs in adata.obs.groupby(cluster_key, observed=True):
        tissue = obs["tissue"].astype(str).value_counts()
        condition = obs["condition"].astype(str).value_counts()
        row = {"cluster": str(cluster), "n_cells": len(obs)}
        for value in ("bone_marrow", "spleen", "thymus"):
            row[f"n_{value}"] = int(tissue.get(value, 0))
            row[f"pct_{value}"] = float(100 * tissue.get(value, 0) / len(obs))
        for value in ("Ctrl", "T2"):
            row[f"n_{value.lower()}"] = int(condition.get(value, 0))
            row[f"pct_{value.lower()}"] = float(100 * condition.get(value, 0) / len(obs))
        composition_rows.append(row)
    pd.DataFrame(composition_rows).to_csv(output / "rerun_cluster_composition.csv", index=False)

    summary = {
        "pptx": str(args.pptx.resolve()),
        "run_root": str(run_root),
        "h5ad": str(h5ad),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "cluster_namespace": "legacy_ppt_qc_leiden_res_0_5",
        "n_clusters": int(adata.obs[cluster_key].nunique()),
        "de_thresholds": {"fdr_lte": args.fdr, "min_abs_logfc": args.min_abs_logfc},
        "legacy_sign_caveat": (
            "Sign agreement assumes the PPT slide title T2 vs Ctrl defines positive as up in T2; "
            "the original code/group order is unavailable."
        ),
        "cluster_matching_rule": (
            "Within the same tissue, maximize recovered legacy DE genes; cluster IDs are never mapped numerically."
        ),
        "outputs": sorted(path.name for path in output.glob("*")),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
