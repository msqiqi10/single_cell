# iNKT follow-up for Rob and Yue

September 25, 2026 — analysis of the existing dataset following the September 24 meeting. This is a discussion brief; no new wet-lab validation is claimed.

## What was done

We retained all 26 existing full differential-expression tables, including global and low-cell-count comparisons with their eligibility flags. The main functional analysis covers 17 eligible tissue-by-cluster comparisons and three whole-tissue comparisons. The rule is gene BH FDR ≤ 0.05 and |log2FC| ≥ 0.25, with T2-up and control-up analyzed separately.

For every main comparison, we audited complete candidate coverage across GO BP, MF and CC. We also reconstructed the original 15,532 retained cells from the six raw expression matrices, retaining all 32,285 shared gene features for a marker-detection audit.

## Main findings

**AP-1 merits priority in bone marrow C0 and spleen C3.** GO:0035976, “transcription factor AP-1 complex,” is enriched among T2-up genes in both populations. Family/global GO FDR values are approximately 0.000261/0.00286 for BM C0 and 0.000193/0.00229 for spleen C3. The overlapping genes are Fos/Fosl2/Jun/Junb in BM C0 and Fos/Jun/Junb/Jund in spleen C3. This is distinct from “AP-1 adaptor complex.” It supports a transcriptional program, not direct measurement of protein activation or DNA binding.

**MAPK activation remains unproven.** BM C0 MAPK cascade GO FDR is 1; the previous KEGG result also did not pass FDR under either retained selection rule. The AP-1 result must not be used to relabel MAPK as significant. Targeted QC matching preserved the principal expression directions, but cannot distinguish environmental biology from processing effects.

**Ccr6 was missed by the original feature filter, not absent from the data.** It was detected in 92 input cell barcodes, below the original min_cells=100 filter. Among retained cells, 60 expressed Ccr6 and 54 of these belonged to C5-2. This supplements the existing Rorc/Il23r identity evidence. In contrast, Il17a was not detected in retained cells and Il17f was detected in only two cells, both outside C5-2. These observations cannot determine stimulated IL-17 secretion capacity. Thymic control C5-2 contains only 11 cells.

**Many GO terms reuse the same genes.** BM C0 has 274 T2-up candidates, but its 52 significant BP terms collectively cover only 66 distinct genes. Spleen C3 has 275 candidates with 96 covered by significant BP terms. Uncovered genes and non-significant annotations remain available in the full tables.

## Functional grouping

We built a mouse GO network with 97 significant BP nodes and 577 edges based on complete measured gene-set overlap. We ran an overlap baseline, a local adaptation of GOLDEN description/adjacency fusion, and three GraphSAGE runs followed by semantic/GNN fusion. The GNN-based reference partition has seven exploratory groups covering translation/ribosome, respiration/ATP, respiratory-chain assembly, nucleotide metabolism, transport, and heterogeneous shared annotations.

Some terms contain muscle, synaptic or reproductive labels driven by shared genes. They should not be interpreted literally as iNKT cell functions. One group remains heterogeneous. Group-to-term-to-gene-to-condition evidence is exported.

This is explicitly a mouse GO adaptation: Jaccard edges replace unavailable original m-type/PAGER inputs. Held-out graph reconstruction was tested against semantic-only and graph-neighborhood baselines. These small, internal tests are not independent biological validation or an external method benchmark.

## Proposed experimental discussion

1. Map BM C0 and spleen C3 to experimentally identifiable phenotypes, then assess whether AP-1-related transcription corresponds to protein/signaling or standardized stimulation responses.
2. For C5-2, separate identity, abundance and per-cell secretion measurements. Match these to the conditions of the existing IL-17 experiments.
3. Retain the prior BM C4 and spleen C3 cytotoxicity/CCT observations. Test whether activation-associated transcription and direct effector readouts diverge; do not infer causality from a short gene score.

Each sample label pools three mice according to the meeting, but animal/pool identities, independent replicate counts and cross-tissue pairing are still missing. These remain exploratory cell-level results. Independent biological replication, processing metadata, the existing IL-17 assay details, and the actual AML data package are required for the next experimental or external-validation stage.

See [candidate discussion](EXPERIMENT_CANDIDATES.md), [gene evidence](tables/Rob_focus_all_20_comparisons.csv), [AP-1 GO results](tables/AP1_complex_GO_exact_ID.csv), and [module evidence](tables/module_gene_evidence.csv.gz).
