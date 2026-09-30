#!/usr/bin/env python3
"""Reinterpret an iNKT DPT trajectory with curated biomarker programs.

This is intentionally a separate, read-only consumer of the processed H5AD.  It
does not overwrite the notebook output or pretend that the legacy c1-c11 PPT
clusters are identical to clusters produced by a rerun.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_H5AD = ROOT / "output/iNKT_scanpy_tutorial_run/inkt_scanpy_tutorial_processed.h5ad"
DEFAULT_GENE_LIST = ROOT / "input/iNKT/iNKT gene list for scRNAseq.xlsx"
DEFAULT_MARKER = ROOT / "input/iNKT/marker.xlsx"
DEFAULT_OUT = ROOT / "output/iNKT_extended_analysis/trajectory_biomarkers"
CLUSTER_KEY = "leiden_res_0_5"


ALIASES = {
    "sepp1": "Selenop",
    "cd103": "Itgae",
    "cd244": "Cd244a",
    "slam7": "Slamf7",
    "ckd6": "Cdk6",
    "hnmpa1": "Hnrnpa1",
    "if27l2a": "Ifi27l2a",
    "gnb2l1": "Rack1",
    "dfna5": "Gsdme",
    "izumol1r": "Izumo1r",
    "dnala1": "Dnaja1",
    "dnala4": "Dnaja4",
    "itgab1": "Itgb1",
    "itgab4": "Itgb4",
    "itgab5": "Itgb5",
    "il23rb": "Il23r",
    "rpgrip": "Rpgrip1",
    "cd24": "Cd24a",
}


STAGE_PROGRAMS = {
    "stage0_immature": ["Cd24a"],
    "stage1_early_effector": ["Tbx21", "Ifng", "Cd27"],
    "stage2_plzf_il4_inkt17": ["Cd44", "Zbtb16", "Il4", "Rorc", "Il17a", "Ccr7"],
    "stage3_cytotoxic_mature": ["Cd44", "Klrb1c", "Prf1", "Gzma", "Gzmb", "Fasl", "Nkg7", "Ccl5"],
    "emigration": ["Klf2", "S1pr1", "Sell"],
}


SIGNATURE_SHEETS = {
    "inkt1_literature": "iNKT1",
    "inkt2_literature": "iNKT2",
    "inkt17_literature": "iNKT17",
    "inkt1_inhouse": "My list iNKT1",
    "inkt2_inhouse": "My list iNKT2",
    "inkt17_inhouse": "My list iNKT17",
    "inkt17_wang2022": "Wang 2022 iNKT17",
    "inkt1_wang2022": "Wang 2022 iNKT1",
    "inkt2_wang2022": "Wang 2022 iNKT2",
    "circulatory": "Circulatory",
    "direct_tcr_activation": "Direct TCR activation",
    "residency": "Residency",
}


def clean_gene(value: object) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    gene = re.sub(r"\s+", "", str(value).strip())
    if not gene or gene.lower() in {"gene", "genes", "marker", "markers", "nan"}:
        return None
    return ALIASES.get(gene.casefold(), gene)


def load_single_column_signatures(path: Path) -> dict[str, list[str]]:
    book = pd.ExcelFile(path, engine="openpyxl")
    available = {name.casefold(): name for name in book.sheet_names}
    signatures: dict[str, list[str]] = {}
    for program, wanted in SIGNATURE_SHEETS.items():
        actual = available.get(wanted.casefold())
        if actual is None:
            if "wang2022" in program:
                subtype = program.split("_", 1)[0]
                subtype_pattern = re.compile(re.escape(subtype) + r"(?!\d)")
                candidates = [
                    name for name in book.sheet_names
                    if subtype_pattern.search(name.casefold()) and "wang" in name.casefold()
                ]
            else:
                candidates = [
                    name for name in book.sheet_names
                    if name.casefold().startswith(wanted.casefold())
                ]
            if len(candidates) != 1:
                raise KeyError(f"Could not uniquely resolve sheet {wanted!r}; found {candidates}")
            actual = candidates[0]
        frame = pd.read_excel(path, sheet_name=actual, header=None, engine="openpyxl")
        genes: list[str] = []
        for value in frame.iloc[:, 0].tolist():
            gene = clean_gene(value)
            if gene and gene.casefold() not in {g.casefold() for g in genes}:
                genes.append(gene)
        signatures[program] = genes
    return signatures


def resolve_programs(
    programs: dict[str, list[str]], var_names: pd.Index
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    present = {str(gene).casefold(): str(gene) for gene in var_names}
    resolved: dict[str, list[str]] = {}
    audit: list[dict[str, object]] = []
    for program, genes in programs.items():
        kept: list[str] = []
        for original in genes:
            canonical = ALIASES.get(original.casefold(), original)
            matched = present.get(canonical.casefold())
            audit.append(
                {
                    "program": program,
                    "input_gene": original,
                    "canonical_gene": canonical,
                    "matched_gene": matched,
                    "present": matched is not None,
                }
            )
            if matched and matched not in kept:
                kept.append(matched)
        if kept:
            resolved[program] = kept
    return resolved, pd.DataFrame(audit)


def expression_matrix(adata: sc.AnnData, genes: list[str]) -> np.ndarray:
    source = adata.raw if adata.raw is not None else adata
    values = source[:, genes].X
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=np.float32)


def zscore_columns(values: np.ndarray) -> np.ndarray:
    means = np.nanmean(values, axis=0, keepdims=True)
    stds = np.nanstd(values, axis=0, keepdims=True)
    stds[~np.isfinite(stds) | (stds < 1e-8)] = 1.0
    return (values - means) / stds


def bh_fdr(frame: pd.DataFrame, p_col: str = "pvalue") -> pd.DataFrame:
    out = frame.copy()
    valid = np.isfinite(out[p_col].to_numpy(dtype=float))
    out["fdr"] = np.nan
    if valid.any():
        out.loc[valid, "fdr"] = multipletests(out.loc[valid, p_col], method="fdr_bh")[1]
    return out


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10 or np.nanstd(x[mask]) < 1e-10 or np.nanstd(y[mask]) < 1e-10:
        return np.nan, np.nan, int(mask.sum())
    rho, pvalue = spearmanr(x[mask], y[mask])
    return float(rho), float(pvalue), int(mask.sum())


def save_heatmap(frame: pd.DataFrame, path: Path, title: str, center: float | None = 0.0) -> None:
    width = max(7.0, 0.55 * frame.shape[1] + 3)
    height = max(4.5, 0.33 * frame.shape[0] + 2)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(frame, cmap="vlag" if center is not None else "viridis", center=center, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--gene-list-xlsx", type=Path, default=DEFAULT_GENE_LIST)
    parser.add_argument("--marker-xlsx", type=Path, default=DEFAULT_MARKER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cluster-namespace", default="rerun_leiden_res_0_5")
    parser.add_argument("--n-bins", type=int, default=20)
    parser.add_argument("--top-dynamic", type=int, default=40)
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for required in (args.input_h5ad, args.gene_list_xlsx, args.marker_xlsx):
        if not required.exists():
            raise FileNotFoundError(required)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_dir = args.output_dir / "tables"
    figure_dir = args.output_dir / "figures"
    table_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    signatures = load_single_column_signatures(args.gene_list_xlsx)
    if args.validate_only:
        adata_backed = sc.read_h5ad(args.input_h5ad, backed="r")
        missing_obs = [key for key in (CLUSTER_KEY, "condition", "tissue", "sample", "dpt_pseudotime") if key not in adata_backed.obs]
        if missing_obs:
            raise KeyError(f"Missing obs columns: {missing_obs}")
        print(json.dumps({"shape": list(adata_backed.shape), "signatures": {k: len(v) for k, v in signatures.items()}}, indent=2))
        return

    adata = sc.read_h5ad(args.input_h5ad)
    source_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    all_programs = {**STAGE_PROGRAMS, **signatures}
    programs, audit = resolve_programs(all_programs, source_names)
    audit.to_csv(table_dir / "trajectory_program_gene_audit.csv", index=False)

    union_genes = sorted({gene for genes in programs.values() for gene in genes})
    values = expression_matrix(adata, union_genes)
    zvalues = zscore_columns(values)
    gene_index = {gene: idx for idx, gene in enumerate(union_genes)}

    score_frame = pd.DataFrame(index=adata.obs_names)
    for program, genes in programs.items():
        cols = [gene_index[gene] for gene in genes]
        score_frame[program] = np.nanmean(zvalues[:, cols], axis=1)

    metadata = adata.obs[["sample", "condition", "tissue", CLUSTER_KEY, "dpt_pseudotime"]].copy()
    cell_scores = pd.concat([metadata, score_frame], axis=1)
    cell_scores.to_csv(table_dir / "cell_trajectory_program_scores.csv.gz", compression="gzip")

    program_names = list(programs)
    grouped_rows: list[pd.DataFrame] = []
    for grouping in ([CLUSTER_KEY], ["condition"], ["tissue"], ["condition", "tissue"], [CLUSTER_KEY, "condition", "tissue"]):
        summary = cell_scores.groupby(grouping, observed=True)[program_names].agg(["mean", "median", "std", "count"])
        summary.columns = [f"{program}__{stat}" for program, stat in summary.columns]
        summary = summary.reset_index()
        summary.insert(0, "grouping", "+".join(grouping))
        grouped_rows.append(summary)
    pd.concat(grouped_rows, ignore_index=True, sort=False).to_csv(
        table_dir / "trajectory_program_scores_by_group.csv", index=False
    )

    cluster_means = cell_scores.groupby(CLUSTER_KEY, observed=True)[program_names].mean()
    save_heatmap(cluster_means.T, figure_dir / "program_scores_by_cluster.png", "Curated programs by rerun cluster")
    condition_tissue_means = cell_scores.groupby(["condition", "tissue"], observed=True)[program_names].mean()
    condition_tissue_means.index = [" | ".join(map(str, idx)) for idx in condition_tissue_means.index]
    save_heatmap(
        condition_tissue_means.T,
        figure_dir / "program_scores_by_condition_tissue.png",
        "Curated programs by condition and tissue",
    )

    pt = cell_scores["dpt_pseudotime"].to_numpy(dtype=float)
    correlation_rows: list[dict[str, object]] = []
    strata = [("all", np.ones(adata.n_obs, dtype=bool))]
    strata += [(f"condition={value}", cell_scores["condition"].astype(str).to_numpy() == str(value)) for value in sorted(cell_scores["condition"].astype(str).unique())]
    strata += [(f"tissue={value}", cell_scores["tissue"].astype(str).to_numpy() == str(value)) for value in sorted(cell_scores["tissue"].astype(str).unique())]
    for program in program_names:
        score = score_frame[program].to_numpy(dtype=float)
        for stratum, mask in strata:
            rho, pvalue, n = safe_spearman(score[mask], pt[mask])
            correlation_rows.append({"program": program, "stratum": stratum, "rho": rho, "pvalue": pvalue, "n": n})
    program_corr = bh_fdr(pd.DataFrame(correlation_rows))
    program_corr.to_csv(table_dir / "program_spearman_vs_dpt.csv", index=False)

    finite_pt = np.isfinite(pt)
    cell_scores["dpt_bin"] = pd.NA
    if finite_pt.sum() >= args.n_bins:
        bin_values = pd.qcut(
            pd.Series(pt[finite_pt]), q=args.n_bins, labels=False, duplicates="drop"
        )
        cell_scores.loc[finite_pt, "dpt_bin"] = bin_values.astype("Int64").to_numpy()
    binned = (
        cell_scores.dropna(subset=["dpt_bin"])
        .groupby(["condition", "tissue", "dpt_bin"], observed=True)[program_names]
        .agg(["mean", "count"])
    )
    binned.columns = [f"{program}__{stat}" for program, stat in binned.columns]
    binned.reset_index().to_csv(table_dir / "program_scores_by_dpt_bin.csv", index=False)

    occupancy = (
        cell_scores.dropna(subset=["dpt_bin"])
        .groupby(["condition", "tissue", "dpt_bin"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    occupancy["fraction_within_condition_tissue"] = occupancy["n_cells"] / occupancy.groupby(
        ["condition", "tissue"], observed=True
    )["n_cells"].transform("sum")
    occupancy.to_csv(table_dir / "dpt_occupancy_by_condition_tissue.csv", index=False)

    dynamic_rows: list[dict[str, object]] = []
    for gene, col in gene_index.items():
        vector = values[:, col].astype(float)
        for stratum, mask in strata:
            rho, pvalue, n = safe_spearman(vector[mask], pt[mask])
            dynamic_rows.append({"gene": gene, "stratum": stratum, "rho": rho, "pvalue": pvalue, "n": n})
    dynamic = bh_fdr(pd.DataFrame(dynamic_rows))
    dynamic.to_csv(table_dir / "curated_gene_spearman_vs_dpt.csv", index=False)

    all_dynamic = dynamic.loc[dynamic["stratum"] == "all"].copy()
    all_dynamic["abs_rho"] = all_dynamic["rho"].abs()
    top_genes = all_dynamic.sort_values(["abs_rho", "fdr"], ascending=[False, True]).head(args.top_dynamic)["gene"].tolist()
    if top_genes and cell_scores["dpt_bin"].notna().any():
        top_cols = [gene_index[gene] for gene in top_genes]
        trend_frame = pd.DataFrame(zvalues[:, top_cols], columns=top_genes, index=adata.obs_names)
        trend_frame["dpt_bin"] = cell_scores["dpt_bin"]
        trend = trend_frame.dropna(subset=["dpt_bin"]).groupby("dpt_bin", observed=True)[top_genes].mean().T
        trend.to_csv(table_dir / "top_dynamic_gene_binned_zscores.csv")
        save_heatmap(trend, figure_dir / "top_curated_dynamic_genes_along_dpt.png", "Curated genes along rerun DPT")

    root_cluster_current: str | None = None
    if "iroot" in adata.uns:
        root_idx = int(adata.uns["iroot"])
        if 0 <= root_idx < adata.n_obs:
            root_cluster_current = str(adata.obs.iloc[root_idx][CLUSTER_KEY])
    root_columns = [name for name in ("stage0_immature", "stage1_early_effector", "emigration", "stage3_cytotoxic_mature") if name in cluster_means]
    root_table = cluster_means[root_columns].copy()
    for column in root_columns:
        std = root_table[column].std(ddof=0)
        root_table[f"{column}__cluster_z"] = (root_table[column] - root_table[column].mean()) / (std if std > 1e-8 else 1.0)
    early_terms = [c for c in root_table if c.endswith("__cluster_z") and not c.startswith("stage3")]
    mature_term = "stage3_cytotoxic_mature__cluster_z"
    root_table["biomarker_root_score"] = root_table[early_terms].mean(axis=1)
    if mature_term in root_table:
        root_table["biomarker_root_score"] -= root_table[mature_term]
    root_table = root_table.sort_values("biomarker_root_score", ascending=False)
    root_table.to_csv(table_dir / "biomarker_root_cluster_scores.csv")
    candidate_root = str(root_table.index[0])

    alternative_summary: dict[str, object] = {
        "current_root_cluster": root_cluster_current,
        "biomarker_candidate_root_cluster": candidate_root,
        "alternative_dpt_computed": False,
    }
    if "X_diffmap" in adata.obsm and "neighbors" in adata.uns:
        candidates = np.flatnonzero(adata.obs[CLUSTER_KEY].astype(str).to_numpy() == candidate_root)
        if len(candidates):
            original_pt = adata.obs["dpt_pseudotime"].copy()
            adata.uns["iroot"] = int(candidates[0])
            sc.tl.dpt(adata, n_dcs=min(10, adata.obsm["X_diffmap"].shape[1]))
            alternative_pt = adata.obs["dpt_pseudotime"].to_numpy(dtype=float).copy()
            rho, pvalue, n = safe_spearman(original_pt.to_numpy(dtype=float), alternative_pt)
            pd.DataFrame(
                {
                    "dpt_pseudotime_current": original_pt.to_numpy(dtype=float),
                    "dpt_pseudotime_biomarker_root": alternative_pt,
                },
                index=adata.obs_names,
            ).to_csv(table_dir / "alternative_biomarker_root_dpt.csv.gz", compression="gzip")
            alternative_summary.update(
                {
                    "alternative_dpt_computed": True,
                    "current_vs_alternative_spearman_rho": rho,
                    "current_vs_alternative_pvalue": pvalue,
                    "n_finite_comparison": n,
                }
            )
            adata.obs["dpt_pseudotime"] = original_pt

    summary = {
        "input_h5ad": str(args.input_h5ad),
        "gene_list_xlsx": str(args.gene_list_xlsx),
        "marker_xlsx": str(args.marker_xlsx),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "cluster_key": CLUSTER_KEY,
        "cluster_namespace": args.cluster_namespace,
        "n_clusters": int(adata.obs[CLUSTER_KEY].nunique()),
        "n_programs": len(programs),
        "n_union_genes": len(union_genes),
        "n_finite_dpt": int(np.isfinite(pt).sum()),
        "n_nonfinite_dpt": int((~np.isfinite(pt)).sum()),
        "root_validation": alternative_summary,
        "limitations": [
            (
                f"This {adata.obs[CLUSTER_KEY].nunique()}-cluster rerun in namespace "
                f"{args.cluster_namespace} is not numerically identical to the legacy c1-c11 PPT run."
            ),
            "DPT is an expression-similarity ordering, not lineage tracing.",
            "No branch labels were available; curated dynamics are fitted against global DPT.",
            "Biological replicate IDs must be recovered before inferential condition tests.",
        ],
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
