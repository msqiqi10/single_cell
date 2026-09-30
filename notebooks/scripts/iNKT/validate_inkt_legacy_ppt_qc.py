from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse

from inkt_palette import palette_state, palette_validation_failures


EXPECTED_SHAPE = (15532, 10670)
EXPECTED_SAMPLE_COUNTS = {
    "Ctrl_BM": 3379,
    "Ctrl_Spleen": 3422,
    "Ctrl_Thymus": 1372,
    "T2_BM": 2731,
    "T2_Spleen": 3371,
    "T2_Thymus": 1257,
}
EXPECTED_CONDITION_COUNTS = {"Ctrl": 8173, "T2": 7359}
EXPECTED_TISSUE_COUNTS = {"bone_marrow": 6110, "spleen": 6793, "thymus": 2629}
EXPECTED_CELL_ID_SHA256 = "641109be6ece372c9ef5f17f7e8dd4c12e512692bcb34ad733998eba3ff7d1f5"
EXPECTED_GENE_ID_SHA256 = "6ea3ce0ccf3bda34d72ccb64eaa4a14139dbdc9603cb5cd44c81754f4307592f"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate an iNKT legacy-PPT-QC pipeline output.")
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--filter-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    adata = sc.read_h5ad(args.h5ad, backed="r")
    observed_shape = tuple(map(int, adata.shape))
    sample_counts = (
        adata.obs["sample"].astype(str).value_counts().sort_index().astype(int).to_dict()
    )
    condition_counts = (
        adata.obs["condition"].astype(str).value_counts().sort_index().astype(int).to_dict()
    )
    tissue_counts = (
        adata.obs["tissue"].astype(str).value_counts().sort_index().astype(int).to_dict()
    )
    cell_id_sha256 = hashlib.sha256(
        ("\n".join(sorted(adata.obs_names.astype(str))) + "\n").encode()
    ).hexdigest()
    gene_id_sha256 = hashlib.sha256(
        ("\n".join(sorted(adata.var_names.astype(str))) + "\n").encode()
    ).hexdigest()
    failures: list[str] = []
    if observed_shape != EXPECTED_SHAPE:
        failures.append(f"shape {observed_shape} != {EXPECTED_SHAPE}")
    if sample_counts != EXPECTED_SAMPLE_COUNTS:
        failures.append(f"sample counts {sample_counts} != {EXPECTED_SAMPLE_COUNTS}")
    if condition_counts != EXPECTED_CONDITION_COUNTS:
        failures.append(f"condition counts {condition_counts} != {EXPECTED_CONDITION_COUNTS}")
    if tissue_counts != EXPECTED_TISSUE_COUNTS:
        failures.append(f"tissue counts {tissue_counts} != {EXPECTED_TISSUE_COUNTS}")
    if cell_id_sha256 != EXPECTED_CELL_ID_SHA256:
        failures.append(f"cell ID SHA256 {cell_id_sha256} != {EXPECTED_CELL_ID_SHA256}")
    if gene_id_sha256 != EXPECTED_GENE_ID_SHA256:
        failures.append(f"gene ID SHA256 {gene_id_sha256} != {EXPECTED_GENE_ID_SHA256}")
    if "counts" not in adata.layers:
        failures.append("counts layer is missing")
        counts_nonnegative_integer = False
    else:
        counts = adata.layers["counts"]
        if hasattr(counts, "to_memory"):
            counts = counts.to_memory()
        values = counts.data if sparse.issparse(counts) else np.asarray(counts)
        counts_nonnegative_integer = bool(
            np.isfinite(values).all()
            and (values >= 0).all()
            and np.equal(values, np.rint(values)).all()
        )
        if tuple(counts.shape) != EXPECTED_SHAPE:
            failures.append(f"counts layer shape {counts.shape} != {EXPECTED_SHAPE}")
        if not counts_nonnegative_integer:
            failures.append("counts layer contains nonfinite, negative, or non-integer values")
    if adata.raw is None:
        failures.append("raw log-normalized expression is missing")
        raw_shape = None
    else:
        raw_shape = tuple(map(int, adata.raw.shape))
        if raw_shape != EXPECTED_SHAPE:
            failures.append(f"raw shape {raw_shape} != {EXPECTED_SHAPE}")
    if not adata.obs_names.is_unique:
        failures.append("cell IDs are not unique")
    if not adata.var_names.is_unique:
        failures.append("gene IDs are not unique")
    for column in ("n_genes_by_counts", "pct_counts_mt", "leiden_res_0_5", "dpt_pseudotime"):
        if column not in adata.obs:
            failures.append(f"obs column is missing: {column}")
    for column in ("sample", "condition", "tissue"):
        if column in adata.obs and adata.obs[column].isna().any():
            failures.append(f"obs column contains missing values: {column}")
    failures.extend(palette_validation_failures(adata))

    hvg_count = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else None
    if hvg_count != 3000:
        failures.append(f"HVG count {hvg_count} != 3000")
    expected_embeddings = {"X_pca": 50, "X_umap": 2, "X_diffmap": 15}
    embedding_shapes: dict[str, list[int] | None] = {}
    for key, n_components in expected_embeddings.items():
        if key not in adata.obsm:
            failures.append(f"embedding is missing: {key}")
            embedding_shapes[key] = None
            continue
        embedding = np.asarray(adata.obsm[key])
        embedding_shapes[key] = list(map(int, embedding.shape))
        if embedding.shape != (EXPECTED_SHAPE[0], n_components):
            failures.append(
                f"embedding shape {key}={embedding.shape} != {(EXPECTED_SHAPE[0], n_components)}"
            )
        if not np.isfinite(embedding).all():
            failures.append(f"embedding contains nonfinite values: {key}")
    for key in ("neighbors", "paga"):
        if key not in adata.uns:
            failures.append(f"uns entry is missing: {key}")
    expected_cluster_keys = {"leiden_res_0_2", "leiden_res_0_5", "leiden_res_1_0"}
    missing_cluster_keys = expected_cluster_keys - set(adata.obs)
    if missing_cluster_keys:
        failures.append(f"Leiden columns are missing: {sorted(missing_cluster_keys)}")
    neighbor_params = dict(adata.uns.get("neighbors", {}).get("params", {}))
    if neighbor_params.get("n_neighbors") != 15 or neighbor_params.get("n_pcs") != 30:
        failures.append(f"unexpected neighbor parameters: {neighbor_params}")

    if "n_genes_by_counts" in adata.obs:
        min_genes = int(adata.obs["n_genes_by_counts"].min())
        max_genes = int(adata.obs["n_genes_by_counts"].max())
        if min_genes < 200 or max_genes >= 2500:
            failures.append(f"cell gene-count bounds violated: min={min_genes}, max={max_genes}")
    else:
        min_genes = None
        max_genes = None
    if "pct_counts_mt" in adata.obs:
        max_pct_mt = float(adata.obs["pct_counts_mt"].max())
        if max_pct_mt >= 5.0:
            failures.append(f"mitochondrial threshold violated: max={max_pct_mt}")
    else:
        max_pct_mt = None
    if "dpt_pseudotime" in adata.obs:
        dpt = adata.obs["dpt_pseudotime"].to_numpy(dtype=float)
        n_finite_dpt = int(np.isfinite(dpt).sum())
        dpt_min = float(np.nanmin(dpt))
        dpt_max = float(np.nanmax(dpt))
        if n_finite_dpt != EXPECTED_SHAPE[0] or dpt_min < 0.0 or dpt_max > 1.0:
            failures.append(
                f"invalid DPT values: finite={n_finite_dpt}, min={dpt_min}, max={dpt_max}"
            )
    else:
        n_finite_dpt = 0
        dpt_min = None
        dpt_max = None

    filter_summary = json.loads(args.filter_summary.read_text())
    if filter_summary.get("profile") != "legacy_ppt_exact_reconstruction":
        failures.append("filter summary profile is not legacy_ppt_exact_reconstruction")
    if filter_summary.get("post_gene_prefilter") != {"n_obs": 18458, "n_vars": 10670}:
        failures.append("post-gene-prefilter shape is not 18458 x 10670")
    if filter_summary.get("post_filter") != {"n_obs": 15532, "n_vars": 10670}:
        failures.append("post-filter shape is not 15532 x 10670")

    result = {
        "status": "passed" if not failures else "failed",
        "h5ad": str(args.h5ad.resolve()),
        "shape": list(observed_shape),
        "raw_shape": list(raw_shape) if raw_shape is not None else None,
        "counts_nonnegative_integer": counts_nonnegative_integer,
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
        "highly_variable_genes": hvg_count,
        "embedding_shapes": embedding_shapes,
        "neighbor_params": neighbor_params,
        "n_finite_dpt": n_finite_dpt,
        "dpt_min": dpt_min,
        "dpt_max": dpt_max,
        "sample_counts": sample_counts,
        "condition_counts": condition_counts,
        "tissue_counts": tissue_counts,
        "palettes": palette_state(adata),
        "cell_id_sha256": cell_id_sha256,
        "gene_id_sha256": gene_id_sha256,
        "n_clusters_leiden_res_0_5": (
            int(adata.obs["leiden_res_0_5"].nunique())
            if "leiden_res_0_5" in adata.obs
            else None
        ),
        "min_n_genes_by_counts": min_genes,
        "max_n_genes_by_counts": max_genes,
        "max_pct_counts_mt": max_pct_mt,
        "failures": failures,
    }
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end="")
    if failures:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
