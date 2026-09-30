from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXPECTED_SHAPE = [15532, 10670]
EXPECTED_NAMESPACE = "legacy_ppt_qc_leiden_res_0_5"
EXPECTED_SIGNATURE_MATCHES = [25, 76, 30, 35, 16, 15, 31, 44, 36, 24, 11, 36]
PIPELINE_STAGES = [
    "qc_validation",
    "gene_sets",
    "signature_scores",
    "marker_validation",
    "legacy",
    "de0",
    "de1",
    "trajectory",
    "overlap_pathway",
    "overall",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit a completed iNKT legacy-PPT-QC run.")
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def parse_pipeline_status(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def inventory(paths: list[Path]) -> list[dict[str, object]]:
    return [
        {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(paths)
    ]


def main() -> int:
    args = parse_args()
    run_root = args.run_root.resolve()
    output = (args.output or (run_root / "run_audit.json")).resolve()
    h5ad = run_root / "preprocess/inkt_scanpy_tutorial_processed.h5ad"
    expected_h5ad = str(h5ad.resolve())
    failures: list[str] = []

    required = [
        h5ad,
        run_root / "notebooks/scanpy_iNKT_legacy_ppt_qc_pipeline.ipynb",
        run_root / "preprocess/summary.json",
        run_root / "preprocess/tables/filter_summary.json",
        run_root / "qc_validation.json",
        run_root / "pipeline_status.txt",
        run_root / "extended/signatures/manifest_gene-sets.json",
        run_root / "extended/signatures/manifest_signature-scores.json",
        run_root / "extended/signatures/manifest_marker-validation.json",
        run_root / "extended/trajectory/summary.json",
        run_root / "extended/de_pathway/legacy/legacy_manifest.json",
        run_root / "extended/de_pathway/current_de/current_de_manifest_shard-0-of-2.json",
        run_root / "extended/de_pathway/current_de/current_de_manifest_shard-1-of-2.json",
        run_root / "extended/de_pathway/overlap_pathway/overlap_pathway_manifest.json",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        failures.append(f"missing required artifacts: {missing}")
        raise SystemExit("Cannot audit incomplete run: " + ", ".join(missing))

    pipeline_status = parse_pipeline_status(run_root / "pipeline_status.txt")
    bad_stages = {
        stage: pipeline_status.get(stage)
        for stage in PIPELINE_STAGES
        if pipeline_status.get(stage) != "0"
    }
    if bad_stages:
        failures.append(f"nonzero or missing pipeline stages: {bad_stages}")

    qc = read_json(run_root / "qc_validation.json")
    if qc.get("status") != "passed" or qc.get("shape") != EXPECTED_SHAPE:
        failures.append("QC validation did not pass the expected shape")

    preprocess = read_json(run_root / "preprocess/summary.json")
    if preprocess.get("shape") != EXPECTED_SHAPE:
        failures.append("preprocess summary shape mismatch")
    if preprocess.get("qc_profile") != "legacy_ppt_exact_reconstruction":
        failures.append("preprocess summary QC profile mismatch")
    n_clusters = int(preprocess.get("n_clusters", -1))

    signature_manifest_paths = sorted((run_root / "extended/signatures").glob("manifest_*.json"))
    signature_manifests = [read_json(path) for path in signature_manifest_paths]
    for manifest in signature_manifests:
        if manifest.get("input_h5ad") != expected_h5ad:
            failures.append(f"signature manifest used another H5AD: {manifest.get('input_h5ad')}")
        if manifest.get("validation", {}).get("shape") != EXPECTED_SHAPE:
            failures.append(f"signature validation shape mismatch: {manifest.get('task_requested')}")
    gene_set_manifest = next(
        manifest for manifest in signature_manifests if manifest.get("task_requested") == "gene-sets"
    )
    matched = [
        int(row["n_matched_unique"])
        for row in gene_set_manifest["validation"]["signature_sets"]
    ]
    if matched != EXPECTED_SIGNATURE_MATCHES:
        failures.append(f"signature matched counts mismatch: {matched}")

    score_summary = read_json(run_root / "extended/signatures/signature_scores/summary.json")
    if score_summary.get("n_cells") != 15532 or score_summary.get("n_sets_scored") != 12:
        failures.append("signature score coverage mismatch")
    marker_summary = read_json(run_root / "extended/signatures/marker_validation/summary.json")
    if len(marker_summary.get("markers_measured", [])) != 10:
        failures.append("marker validation did not measure exactly ten retained markers")
    if set(marker_summary.get("markers_missing", [])) != {"Il17a", "Il17f"}:
        failures.append(f"unexpected missing markers: {marker_summary.get('markers_missing')}")

    trajectory = read_json(run_root / "extended/trajectory/summary.json")
    if trajectory.get("input_h5ad") != expected_h5ad:
        failures.append(f"trajectory used another H5AD: {trajectory.get('input_h5ad')}")
    if trajectory.get("shape") != EXPECTED_SHAPE:
        failures.append("trajectory shape mismatch")
    if trajectory.get("cluster_namespace") != EXPECTED_NAMESPACE:
        failures.append("trajectory cluster namespace mismatch")
    if trajectory.get("n_clusters") != n_clusters:
        failures.append("trajectory cluster count mismatch")
    if trajectory.get("n_finite_dpt") != 15532 or trajectory.get("n_nonfinite_dpt") != 0:
        failures.append("trajectory contains nonfinite DPT")
    if trajectory.get("n_programs") != 16 or trajectory.get("n_union_genes") != 275:
        failures.append("trajectory program coverage mismatch")

    de_dir = run_root / "extended/de_pathway/current_de"
    de_manifests = [
        read_json(de_dir / f"current_de_manifest_shard-{index}-of-2.json")
        for index in (0, 1)
    ]
    for manifest in de_manifests:
        if manifest.get("h5ad") != expected_h5ad:
            failures.append(f"DE shard used another H5AD: {manifest.get('h5ad')}")
        if manifest.get("h5ad_shape") != EXPECTED_SHAPE:
            failures.append("DE shard shape mismatch")
        if manifest.get("cluster_namespace") != EXPECTED_NAMESPACE:
            failures.append("DE shard namespace mismatch")
        if manifest.get("failed") != 0:
            failures.append("DE shard contains failed units")

    shard_status = [
        pd.read_csv(de_dir / f"current_de_unit_status_shard-{index}-of-2.csv")
        for index in (0, 1)
    ]
    shard_unit_sets = [set(frame["unit_id"].astype(str)) for frame in shard_status]
    if shard_unit_sets[0] & shard_unit_sets[1]:
        failures.append("DE shard unit sets overlap")
    statuses = pd.concat(shard_status, ignore_index=True)
    if len(statuses) != 26 or statuses["unit_id"].nunique() != 26:
        failures.append("DE shards do not cover 26 unique units")
    status_counts = statuses["status"].value_counts().sort_index().to_dict()
    if status_counts != {"completed": 19, "skipped_insufficient_cells": 7}:
        failures.append(f"unexpected DE unit statuses: {status_counts}")

    unit_paths = sorted((de_dir / "unit_csv").glob("*.csv"))
    unit_row_counts: dict[str, int] = {}
    for path in unit_paths:
        frame = pd.read_csv(path, low_memory=False)
        unit_row_counts[path.stem] = int(len(frame))
        if len(frame) != 10670 or frame["gene"].nunique() != 10670:
            failures.append(f"DE unit is not one row per retained gene: {path.name}")
        if frame.duplicated(["unit_id", "gene"]).any():
            failures.append(f"duplicate unit/gene DE rows: {path.name}")
        for column in ("pvals", "pvals_adj", "pct_expressing_t2", "pct_expressing_ctrl"):
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            # Sparse means can overshoot 1 by a few ulps (observed max 6.6e-14).
            if ((values < -1e-12) | (values > 1 + 1e-12)).any():
                failures.append(f"out-of-range DE values in {path.name}: {column}")
    if len(unit_paths) != 19 or sum(unit_row_counts.values()) != 202730:
        failures.append("completed DE unit files do not total 19 x 10,670 rows")

    legacy = read_json(run_root / "extended/de_pathway/legacy/legacy_manifest.json")
    legacy_expected = {
        "n_sheets": 18,
        "n_detected_tables": 27,
        "n_legacy_de_rows": 298,
        "n_legacy_pathway_rows": 625,
        "n_curated_gene_sets": 12,
        "n_curated_gene_rows": 437,
        "n_unresolved_table_scopes": 12,
    }
    for key, value in legacy_expected.items():
        if legacy.get(key) != value:
            failures.append(f"legacy manifest {key}={legacy.get(key)} != {value}")

    overlap = read_json(run_root / "extended/de_pathway/overlap_pathway/overlap_pathway_manifest.json")
    if overlap.get("n_current_de_rows") != 202730 or overlap.get("n_enrichment_rows") != 456:
        failures.append("overlap/pathway row counts mismatch")
    if "curated_XLSX_signature_ORA_only" not in overlap.get("current_pathway_rerun_status", ""):
        failures.append("overlap/pathway status incorrectly claims an external pathway rerun")
    enrichment = pd.read_csv(
        run_root / "extended/de_pathway/overlap_pathway/current_offline_gene_set_enrichment.csv"
    )
    significant_enrichment = int((enrichment["fdr"] < 0.05).sum())

    raw_input_files = sorted((ROOT / "input/iNKT/data").glob("*/sample_feature_bc_matrix/*"))
    provenance_inputs = raw_input_files + [
        ROOT / "input/iNKT/iNKT gene list for scRNAseq.xlsx",
        ROOT / "input/iNKT/marker.xlsx",
        ROOT / "input/iNKT/iNKT.pptx",
    ]
    if len(raw_input_files) != 18:
        failures.append(f"expected 18 raw 10x files, found {len(raw_input_files)}")

    code_files = [
        ROOT / "run_inkt_legacy_ppt_qc_pipeline.sh",
        ROOT / "notebooks/scripts/iNKT/build_inkt_legacy_ppt_qc_notebook.py",
        ROOT / "notebooks/scripts/iNKT/execute_inkt_notebook.py",
        ROOT / "notebooks/scripts/iNKT/validate_inkt_legacy_ppt_qc.py",
        ROOT / "notebooks/scripts/iNKT/audit_inkt_legacy_ppt_qc_run.py",
        ROOT / "notebooks/scripts/iNKT/analyze_inkt_gene_signatures.py",
        ROOT / "notebooks/scripts/iNKT/test_analyze_inkt_gene_signatures.py",
        ROOT / "notebooks/scripts/iNKT/analyze_inkt_trajectory_biomarkers.py",
        ROOT / "notebooks/scripts/iNKT/run_inkt_de_pathway_overlap.py",
    ]
    artifact_paths = [
        path for path in run_root.rglob("*") if path.is_file() and path.resolve() != output
    ]

    result = {
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "run_root": str(run_root),
        "pipeline_status": pipeline_status,
        "qc_profile": preprocess.get("qc_profile"),
        "shape": EXPECTED_SHAPE,
        "n_clusters_leiden_res_0_5": n_clusters,
        "cluster_namespace": EXPECTED_NAMESPACE,
        "h5ad": {
            "path": expected_h5ad,
            "size_bytes": h5ad.stat().st_size,
            "sha256": sha256_file(h5ad),
        },
        "signatures": {
            "n_sets": 12,
            "matched_gene_counts": matched,
            "markers_measured": marker_summary.get("markers_measured"),
            "markers_missing": marker_summary.get("markers_missing"),
        },
        "trajectory": {
            "n_programs": trajectory.get("n_programs"),
            "n_union_genes": trajectory.get("n_union_genes"),
            "n_finite_dpt": trajectory.get("n_finite_dpt"),
            "root_validation": trajectory.get("root_validation"),
        },
        "de": {
            "unit_status_counts": status_counts,
            "completed_unit_csvs": len(unit_paths),
            "total_gene_rows": sum(unit_row_counts.values()),
        },
        "enrichment": {
            "scope": overlap.get("current_pathway_rerun_status"),
            "rows": int(len(enrichment)),
            "fdr_lt_0_05_rows": significant_enrichment,
        },
        "artifact_counts_before_audit": {
            "files": len(artifact_paths),
            "png": sum(path.suffix == ".png" for path in artifact_paths),
            "csv": sum(path.suffix == ".csv" for path in artifact_paths),
            "json": sum(path.suffix == ".json" for path in artifact_paths),
        },
        "input_inventory": inventory(provenance_inputs),
        "code_inventory": inventory(code_files),
        "artifact_inventory": inventory(artifact_paths),
        "environment": {
            "python": platform.python_version(),
            "packages": {
                name: importlib.metadata.version(name)
                for name in ("scanpy", "anndata", "numpy", "pandas", "scipy", "matplotlib")
            },
            "cuda_visible_devices": "",
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({key: value for key, value in result.items() if not key.endswith("inventory")}, indent=2))
    if failures:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
