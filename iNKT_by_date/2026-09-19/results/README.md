# Results

- `tables/cytotoxicity_contrasts.csv`: all core/E1/prior-effector contrasts; low-count groups explicitly untested.
- `tables/cytotoxicity_QC_sensitivity.csv`: five balanced QC resamples, effect direction and balance metrics.
- `tables/cytotoxicity_integrated_summary.csv`: original-cluster core results with QC sensitivity.
- `tables/focus_cytotoxic_gene_DE.csv`: five component genes in BM C4 and spleen C3.
- `tables/focus_stable_cluster_crosswalk.csv`: overlap with stable tissue clusters and corresponding effects.
- `tables/GO_ORA_all.csv.gz`: every tested GO BP/MF/CC term, including zero overlaps; family and cross-contrast BH.
- `tables/GO_cytotoxicity_targeted_audit.csv`: killing/degranulation-related GO evidence, including nonsignificant terms.
- `tables/GO_contrast_coverage.csv`: query, expressed, annotated background and tested term counts.
- `tables/trajectory_root_audit.csv`: actual root genes, root cells, connectivity and root sensitivity.
- `tables/trajectory_binned_curves.csv`: shared-bin means and cell counts underlying curves.
- `de/`: new cluster marker DE and frozen copied condition DE.
- `figures/`: PNG and editable SVG for each presentation page.
- `objects/`: lightweight local trajectory AnnData files; expression stays in source objects.
- `verification.json`: independent numerical and output checks.
- `logs/`: analysis, rendering, tests and verification. Initial orchestration stopped before the presentation builder existed; `render_verify` resumed the completed analysis and generated the final deliverables. The initial exit file is historical; use `workflow_status.json` for final status.
