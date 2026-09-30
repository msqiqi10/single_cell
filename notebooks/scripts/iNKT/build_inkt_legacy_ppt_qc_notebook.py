from __future__ import annotations

import argparse
import json
from pathlib import Path

import nbformat


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TEMPLATE = ROOT / "notebooks/iNKT/scanpy_iNKT_preprocess_plotting_trajectory.ipynb"

EXPECTED_SAMPLE_COUNTS = {
    "Ctrl_BM": 3379,
    "Ctrl_Spleen": 3422,
    "Ctrl_Thymus": 1372,
    "T2_BM": 2731,
    "T2_Spleen": 3371,
    "T2_Thymus": 1257,
}

LEGACY_QC_SOURCE = r'''
# Legacy-PPT QC, reconstructed exactly from the six reported sample counts.
# Order is important: gene filtering precedes QC metric calculation.
raw_shape = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}
raw_sample_counts = (
    adata.obs["sample"].astype(str).value_counts().sort_index().astype(int).to_dict()
)

sc.pp.filter_genes(adata, min_cells=100)
post_gene_prefilter = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)}

adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl", "RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.match(r"^(Hb[ab]|HBA|HBB)")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

cell_mask = (
    (adata.obs["n_genes_by_counts"] >= 200)
    & (adata.obs["n_genes_by_counts"] < 2500)
    & (adata.obs["pct_counts_mt"] < 5.0)
)
adata = adata[cell_mask].copy()

post_sample_counts = (
    adata.obs["sample"].astype(str).value_counts().sort_index().astype(int).to_dict()
)
expected_sample_counts = {
    "Ctrl_BM": 3379,
    "Ctrl_Spleen": 3422,
    "Ctrl_Thymus": 1372,
    "T2_BM": 2731,
    "T2_Spleen": 3371,
    "T2_Thymus": 1257,
}
if adata.shape != (15532, 10670):
    raise RuntimeError(f"Legacy-PPT QC shape mismatch: {adata.shape} != (15532, 10670)")
if post_sample_counts != expected_sample_counts:
    raise RuntimeError(
        "Legacy-PPT QC sample counts mismatch: "
        f"{post_sample_counts} != {expected_sample_counts}"
    )

qc_cols = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "pct_counts_hb",
]
qc_summary = adata.obs.groupby(
    ["condition", "tissue", "sample"], observed=True
)[qc_cols].describe()
qc_summary.to_csv(TABLE_DIR / "qc_summary_by_condition_tissue_sample.csv")

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample",
    jitter=0.2,
    multi_panel=True,
    show=False,
)
save_and_show(FIG_DIR / "qc_violin_by_sample.png", dpi=160, bbox_inches="tight")

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="sample", show=False)
save_and_show(FIG_DIR / "qc_scatter_counts_genes_sample.png", dpi=160, bbox_inches="tight")

filter_summary = {
    "profile": "legacy_ppt_exact_reconstruction",
    "inference_basis": "exact total and six-of-six per-sample count reconstruction",
    "filter_order": [
        "genes: min_cells >= 100 on all raw cells",
        "recompute cell QC metrics on the retained 10,670 genes",
        "cells: 200 <= n_genes_by_counts < 2500",
        "cells: pct_counts_mt < 5.0",
    ],
    "raw": raw_shape,
    "post_gene_prefilter": post_gene_prefilter,
    "post_filter": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
    "raw_sample_counts": raw_sample_counts,
    "post_filter_sample_counts": post_sample_counts,
    "thresholds": {
        "min_cells_per_gene": 100,
        "min_genes_by_counts": 200,
        "max_genes_by_counts_exclusive": 2500,
        "max_pct_counts_mt_exclusive": 5.0,
    },
}
(TABLE_DIR / "filter_summary.json").write_text(
    json.dumps(filter_summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
)

print(json.dumps(filter_summary, indent=2))
print(qc_summary.head())
'''.strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a clean copy of the iNKT pipeline with the exactly reconstructed legacy-PPT QC."
    )
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output-notebook", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    return parser.parse_args()


def relative_to_root(path: Path) -> Path:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT)
    except ValueError as exc:
        raise ValueError(f"Path must be inside repository root {ROOT}: {resolved}") from exc


def main() -> int:
    args = parse_args()
    template = args.template.resolve()
    output_notebook = args.output_notebook.resolve()
    run_dir_relative = relative_to_root(args.run_dir)
    relative_to_root(output_notebook)

    notebook = nbformat.read(template, as_version=4)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.execution_count = None
            cell.outputs = []

    setup_matches = [
        cell
        for cell in notebook.cells
        if cell.cell_type == "code"
        and 'RUN_DIR = ROOT / "output/iNKT_scanpy_tutorial_run"' in cell.source
    ]
    if len(setup_matches) != 1:
        raise RuntimeError(f"Expected one pipeline setup cell; found {len(setup_matches)}")
    setup_matches[0].source = setup_matches[0].source.replace(
        'RUN_DIR = ROOT / "output/iNKT_scanpy_tutorial_run"',
        f"RUN_DIR = ROOT / {json.dumps(run_dir_relative.as_posix())}",
        1,
    )

    qc_matches = [
        cell
        for cell in notebook.cells
        if cell.cell_type == "code"
        and 'min_genes_by_counts": 200' in cell.source
        and "sc.pp.calculate_qc_metrics" in cell.source
    ]
    if len(qc_matches) != 1:
        raise RuntimeError(f"Expected one current-QC cell; found {len(qc_matches)}")
    qc_matches[0].source = LEGACY_QC_SOURCE

    summary_matches = [
        cell
        for cell in notebook.cells
        if cell.cell_type == "code" and '"notebook": "scanpy_iNKT_preprocess_plotting_trajectory"' in cell.source
    ]
    if len(summary_matches) != 1:
        raise RuntimeError(f"Expected one summary cell; found {len(summary_matches)}")
    summary_matches[0].source = summary_matches[0].source.replace(
        '"notebook": "scanpy_iNKT_preprocess_plotting_trajectory",',
        '"notebook": "scanpy_iNKT_legacy_ppt_qc_pipeline",\n'
        '    "qc_profile": "legacy_ppt_exact_reconstruction",',
        1,
    )
    summary_matches[0].source = summary_matches[0].source.replace(
        '"source_archive": "input.zip",',
        '"source_data": "input/iNKT/data (six extracted 10x matrices)",',
        1,
    )

    title = notebook.cells[0]
    if title.cell_type != "markdown":
        raise RuntimeError("Expected the first notebook cell to be markdown")
    title.source = title.source.replace(
        "# iNKT Scanpy preprocessing, plotting, and trajectory workflow",
        "# iNKT pipeline using the reconstructed legacy-PPT QC",
        1,
    )
    title.source += (
        "\n\nQC order: gene `min_cells=100`, recompute QC on the retained genes, then "
        "retain `200 <= n_genes_by_counts < 2500` and `pct_counts_mt < 5`."
    )

    notebook.metadata["inkt_qc_profile"] = "legacy_ppt_exact_reconstruction"
    notebook.metadata["expected_shape"] = [15532, 10670]
    notebook.metadata["expected_sample_counts"] = EXPECTED_SAMPLE_COUNTS

    output_notebook.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, output_notebook)
    print(output_notebook)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
