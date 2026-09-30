// Builds the 30 Sept drill-down deck (slides 32-43, continuing the 31-slide deck of 25 Sept).
// Usage: NODE_PATH=<dir with pptxgenjs> node build_deck_0930.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const ROOT = path.resolve(__dirname, "..");
const DD = path.join(ROOT, "drilldown");
const OUT = path.join(ROOT, "presentation", "iNKT_pathway_network_drilldown_20260930.pptx");
const C = { ink: "1F2A30", teal: "137D8A", tealLight: "D7ECEE", rust: "C8694A", rustLight: "F6E2DA", grey: "5F6B72", line: "D5DADD", bg: "FFFFFF", dark: "12343A" };
const H = "Cambria", B = "Calibri";
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "iNKT pathway network drill-down (30 Sept)";
let pageNo = 31;
function base(title, subtitle) {
  const s = pres.addSlide(); s.background = { color: C.bg }; pageNo += 1;
  s.addText(title, { x: 0.5, y: 0.3, w: 12.3, h: 0.7, fontFace: H, fontSize: 28, bold: true, color: C.ink, margin: 0, isTextBox: true });
  if (subtitle) s.addText(subtitle, { x: 0.5, y: 0.98, w: 12.3, h: 0.4, fontFace: B, fontSize: 14, color: C.grey, margin: 0, isTextBox: true });
  s.addText(`${pageNo}`, { x: 12.3, y: 7.0, w: 0.6, h: 0.3, fontFace: B, fontSize: 10, color: C.grey, align: "right", margin: 0, isTextBox: true });
  s.addText("Exploratory cell-level results · 3 mice pooled per tissue × condition label (n = 1 pool)", { x: 0.5, y: 7.0, w: 9, h: 0.3, fontFace: B, fontSize: 10, color: C.grey, margin: 0, isTextBox: true });
  return s;
}
function card(s, x, y, w, h, fill) { s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: fill }, rectRadius: 0.08 }); }
function bullets(s, items, opt) {
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1 } })),
    Object.assign({ fontFace: B, fontSize: 15, color: C.ink, paraSpaceAfter: 6, valign: "top", isTextBox: true, margin: 0 }, opt));
}
function img(s, file, x, y, w, h, pw, ph) { // contain-fit
  const r = Math.min(w / pw, h / ph), iw = pw * r, ih = ph * r;
  s.addImage({ path: file, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih });
}
const hdr = (a) => a.map((t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: C.teal } } }));
const tbl = (s, rows, o) => s.addTable(rows, Object.assign({ fontFace: B, fontSize: 13, color: C.ink, border: { type: "solid", pt: 0.75, color: C.line }, valign: "middle" }, o));

// 32 section
{
  const s = pres.addSlide(); pageNo += 1; s.background = { color: C.dark };
  s.addText("Pathway network drill-down", { x: 0.8, y: 2.0, w: 11.5, h: 0.9, fontFace: H, fontSize: 40, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
  s.addText("GO MF + full KEGG network, anchor grouping, and gene-level bridges between AP-1, MAPK, IL-17 and TCR", { x: 0.8, y: 3.0, w: 11.5, h: 0.9, fontFace: B, fontSize: 20, color: "CFE6E8", margin: 0, isTextBox: true });
  s.addText("Answering Yue's question of 29 Sept · analysis of 30 Sept 2026 · no new data, no recomputation of DE", { x: 0.8, y: 4.3, w: 11.5, h: 0.5, fontFace: B, fontSize: 15, color: "9FC3C7", margin: 0, isTextBox: true });
  s.addNotes("Continues the numbering after the 25 Sept deck (slides 20-31). Source: 2026-09-30/README.md. This section answers Yue's question of 29 Sept: do the anchor pathways (AP-1, MAPK, IL-17, TCR, TNF) group together in a functional network, and what is the gene-level basis if they do?");
}
// 33 asked vs done
{
  const s = base("What Yue asked on 29 Sept, and what we did", "Same DE tables (20 eligible comparisons) and frozen GO; KEGG_2019_Mouse added");
  const rows = [
    ["Add GO MF and the full KEGG", "ORA on GO BP, GO MF, KEGG; term size 5–500; BH within comparison × direction × source"],
    ["Keep every significant term", "All terms with q ≤ 0.05 and ≥ 3 hit genes enter the network (no top-N cut)"],
    ["Mixed-source network", "GO BP / MF and KEGG nodes in one graph; edges = Jaccard ≥ 0.25 and ≥ 3 shared genes"],
    ["Encode significance and group", "Node size = best −log10 q; colour = group (GNN + semantic)"],
    ["Do anchors group together?", "Louvain, semantic, GOLDEN fusion, GNN + semantic; k by silhouette (2–15)"],
    ["Go below pathway level", "Gene-level drill-down and an interactive HTML explorer (pathway ↔ gene)"],
  ];
  tbl(s, [hdr(["Request", "Done"]), ...rows.map((r, i) => r.map((t, j) => ({ text: t, options: { bold: j === 0 } })))], { x: 0.5, y: 1.6, w: 12.3, colW: [4.2, 8.1], rowH: 0.62 });
  s.addText("Not done: no new datasets, no ORA/GNN/DE recomputation for the drill-down; no STRING/PPI edges.", { x: 0.5, y: 6.3, w: 12.3, h: 0.4, fontFace: B, fontSize: 14, italic: true, color: C.rust, margin: 0, isTextBox: true });
  s.addNotes("Sources: 2026-09-30/README.md ('What was done', 'Drill-down'). ORA reuse check vs the 19 Sept GO results agreed term by term (results/tables/ORA_reuse_check_vs_0919.csv; 312,550 rows, max q difference 1e-16). Yue's request is paraphrased from the 29 Sept discussion.");
}
// 34 design choices
{
  const s = base("Three design choices to agree on", "All are visible in the figure legends and can be toggled in the HTML explorer");
  const cols = [
    ["1. Primary version excludes ribo / mito / Hsp", "Regex ^(Rp[ls]\\d|Rplp\\d|Rpsa$|Mrp[ls]\\d|mt-|Hsp(?!g)|Dnaj) removes 242 of 10,670 genes from both query and background. Rps6k* kinases and Hspg2 are kept. Version (a) with all genes is kept for comparison."],
    ["2. Fixed anchor nodes", "12 anchors are always drawn, even if not significant: KEGG MAPK, IL-17, TCR, Th17, TNF; GO AP-1 complex, MAPK / ERK / p38 / JNK cascades, Th17 response, IL-17 production. Anchors are never merged into other nodes."],
    ["3. Redundancy collapse", "Non-anchor terms with Jaccard ≥ 0.75 to a more significant term are merged, so near-identical GO terms do not inflate the graph. Merged members stay listed in node_attributes."],
  ];
  cols.forEach((c, i) => {
    const x = 0.5 + i * 4.15; card(s, x, 1.6, 3.95, 4.9, i === 1 ? C.rustLight : C.tealLight);
    s.addText(c[0], { x: x + 0.2, y: 1.75, w: 3.55, h: 0.9, fontFace: H, fontSize: 17, bold: true, color: C.ink, margin: 0, isTextBox: true, valign: "top" });
    s.addText(c[1], { x: x + 0.2, y: 2.75, w: 3.55, h: 3.6, fontFace: B, fontSize: 14, color: C.ink, margin: 0, isTextBox: true, valign: "top" });
  });
  s.addNotes("Sources: README.md 'What was done'; code/common.py (EXCL_REGEX, ANCHOR_KEGG, ANCHOR_GO); results/network/network_summary.json (242/10,670 excluded genes is from README; anchors n_anchor_nodes = 12); results/tables/redundancy_merges_excl.csv. The exclusion is the primary version because ribosomal, mitochondrial-encoded and heat-shock genes are the classic dissociation / QC signature; Rob raised this on 24 Sept.");
}
// 35 network
{
  const s = base("The network (primary version b)", "66 nodes, 157 edges after excluding ribo / mito / Hsp genes");
  img(s, path.join(DD, "figures", "network_full_b_primary_excl_ribo_mito_hsp_v2.png"), 0.4, 1.5, 8.6, 5.4, 3022, 1614);
  card(s, 9.2, 1.6, 3.65, 5.2, C.tealLight);
  s.addText("Counts", { x: 9.4, y: 1.7, w: 3.3, h: 0.4, fontFace: H, fontSize: 18, bold: true, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, ["Version b: 94 → 66 nodes after collapse, 157 edges, 21 isolated nodes", "Version a (all genes): 159 → 112 nodes, 231 edges, 32 isolated", "Groups b (Louvain / semantic / GOLDEN / GNN + semantic): 31 / 15 / 3 / 4", "One large connected block is energy metabolism; anchors are mostly isolated nodes"], { x: 9.4, y: 2.2, w: 3.3, h: 4.4, fontSize: 14 });
  s.addNotes("Sources: results/network/network_summary.json (excl: n_nodes_before_collapse 94, n_nodes 66, n_edges 157, n_isolates 21; all: 159, 112, 231, 32); README.md key-number table for group counts (a: 49/15/2/15; b: 31/15/3/4). Figure: drilldown/figures/network_full_b_primary_excl_ribo_mito_hsp_v2.png (label-avoidance version of the 09-30 figure; labels no longer overlap each other, a few still sit on the dense energy cluster).");
}
// 36 energy core
{
  const s = base("The significant core is energy metabolism", "It survives the exclusion: six strongest terms, version b (ribo / mito-encoded / Hsp genes already removed)");
  const labels = ["Oxidative phosphorylation (KEGG)", "oxidative phosphorylation (GO)", "Thermogenesis (KEGG)", "aerobic respiration (GO)", "Parkinson disease (KEGG)", "energy derivation by oxidation (GO)"];
  const vals = [25.77, 23.88, 22.55, 19.96, 19.40, 18.17];
  s.addChart(pres.charts.BAR, [{ name: "−log10 q", labels, values: vals }], { x: 0.5, y: 1.5, w: 7.4, h: 5.2, barDir: "bar", chartColors: [C.teal], catAxisOrientation: "maxMin", catAxisLabelFontSize: 12, catAxisLabelFontFace: B, valAxisLabelFontSize: 11, valAxisLabelFontFace: B, showValue: true, dataLabelFontSize: 11, dataLabelFormatCode: "0.0", showLegend: false, valGridLine: { color: "E5E8EA", size: 0.5 }, valAxisTitle: "best −log10 q", showValAxisTitle: true, valAxisTitleFontSize: 11 });
  card(s, 8.2, 1.6, 4.65, 5.1, C.rustLight);
  bullets(s, ["OXPHOS / ATP synthesis / respiratory chain terms dominate; the driver genes are nuclear-encoded (Ndu*, Cox*, Atp5*, Uqcr*)", "Strongest in BM C0 and Spleen C3 (T2-up); BM and spleen C5-2 have no significant term", "Honest framing: energy metabolism is not new for the bench, but it is a broader signal than the G01 stress program (ribo / Hsp) from 25 Sept", "Enrichment of a large, well-annotated pathway is not evidence of increased activity"], { x: 8.4, y: 1.75, w: 4.25, h: 4.8, fontSize: 14 });
  s.addNotes("Sources: results/tables/node_attributes_excl.csv (best_q: 1.68e-26, 1.32e-24, 2.85e-23, 1.09e-20, 3.94e-20, 6.82e-19; chart shows -log10 of these). Group L01/G01 in results/tables/group_summary_excl.csv. Comparison statements about C0/C3 and C5-2: README.md ('热图显示...') and drilldown/tables/pathway_drill_summary.csv (Energy module DEG counts BM C0 31, Spleen C3 36, BM C5-2 1). Note that Parkinson disease is a KEGG term that is essentially a respiratory-chain gene list here.");
}
// 37 anchors table
{
  const s = base("Anchors: only AP-1 is significant", "Best q across comparisons, version b (BH within comparison × direction × source)");
  const R = [["GO AP-1 complex (GO:0035976)", "1.7e-4", "6", "Yes: BM C0, Spleen C3"], ["KEGG IL-17 signaling", "0.071", "51", "Trend only"], ["KEGG TNF signaling", "0.56", "78", "No"], ["KEGG MAPK signaling", "0.64", "159", "No"], ["KEGG TCR, KEGG Th17, GO MAPK / ERK / p38 / JNK, GO Th17 response", "1.0", "16–114", "No"], ["GO interleukin-17 production", "n/a", "0", "No measured genes in frozen annotation"]];
  tbl(s, [hdr(["Anchor", "Best q", "Measured genes", "Significant?"]), ...R.map((r, i) => r.map((t, j) => ({ text: t, options: { bold: i === 0, fill: i === 0 ? { color: C.tealLight } : undefined } })))], { x: 0.5, y: 1.6, w: 12.3, colW: [5.6, 1.5, 2.0, 3.2], rowH: 0.6 });
  s.addText("AP-1 is the only anchor with statistical support, and it rests on five genes: Fos, Fosl2, Jun, Junb, Jund. MAPK is not enriched among T2-up genes even though 23 of its members are DE (next slides).", { x: 0.5, y: 5.6, w: 12.3, h: 0.9, fontFace: B, fontSize: 15, color: C.ink, margin: 0, isTextBox: true, valign: "top" });
  s.addNotes("Source: results/tables/anchor_by_method_excl.csv (best_q, n_measured_genes: AP-1 0.000170, 6 genes; IL-17 0.0712, 51; TNF 0.561, 78; MAPK 0.636, 159; GO MAPK cascade 1.0, 114; JNK 30; p38 22; ERK 22; Th17 response 16; TCR 92; Th17 KEGG 79; IL-17 production no q, 0 measured genes). The AP-1 members in the hit list are from README.md. The '16-114' range covers the GO anchors, TCR (92) and Th17 (79).");
}
// 38 co-grouping
{
  const s = base("Do the anchors group together? At pathway level: it depends", "Group of each anchor by method, version b (letter + number = group id)");
  const R = [
    ["AP-1 complex", "L05", "S11", "F03", "G03"], ["IL-17 signaling", "L25", "S14", "F03", "G03"], ["MAPK signaling", "L26", "S15", "F03", "G03"], ["TCR signaling", "L17", "S09", "F03", "G03"], ["Th17 differentiation", "L17", "S14", "F03", "G03"], ["TNF signaling", "L25", "S09", "F03", "G03"],
  ];
  tbl(s, [hdr(["Anchor", "Louvain", "Semantic", "GOLDEN fusion", "GNN + semantic"]), ...R.map((r) => r.map((t, j) => ({ text: t, options: { bold: j === 0, align: j ? "center" : "left" } })))], { x: 0.5, y: 1.55, w: 7.3, colW: [2.3, 1.1, 1.2, 1.4, 1.3], fontSize: 12, rowH: 0.42 });
  s.addText("Louvain and semantic scatter them; GOLDEN and GNN at the silhouette k put all in G03, a 36-node catch-all of isolates and 2-node fragments.", { x: 0.5, y: 4.7, w: 7.3, h: 0.9, fontFace: B, fontSize: 14, color: C.ink, margin: 0, isTextBox: true, valign: "top" });
  card(s, 8.1, 1.55, 4.75, 5.1, C.rustLight);
  s.addText("k-sensitivity (GNN + semantic)", { x: 8.3, y: 1.65, w: 4.35, h: 0.4, fontFace: H, fontSize: 16, bold: true, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, ["k ≤ 6: all anchors in one group (k chosen by silhouette = 4)", "k ≥ 7: MAPK separates from AP-1 / IL-17 / TCR", "Semantic-only separates them already at k ≥ 5", "Version a at k = 15: AP-1, MAPK, and IL-17 / TCR are three separate places"], { x: 8.3, y: 2.15, w: 4.35, h: 3.0, fontSize: 14 });
  s.addText("Honest message: the question cannot be answered at pathway level. Co-grouping reflects k and method, not biology.", { x: 8.3, y: 5.3, w: 4.35, h: 1.2, fontFace: B, fontSize: 15, bold: true, color: C.rust, margin: 0, isTextBox: true, valign: "top" });
  s.addNotes("Sources: results/tables/anchor_by_method_excl.csv (overlap_louvain, semantic_only, golden_fusion, gnn_semantic columns); results/tables/anchor_k_sensitivity_excl.csv (gnn_semantic rows: AP1_with_MAPK True for k=2..6, False for k>=7; silhouette peak at k=4, 0.261; semantic_only False from k=5); version a statement from README.md. G03 composition: results/tables/group_summary_excl.csv. Therefore the co-grouping in the primary version is not evidence of a functional relationship.");
}
// 39 bridge heatmap
{
  const s = base("Gene level: only Fos and Jun bridge the anchors", "Genes in ≥ 2 anchor pathways and DE in ≥ 1 key comparison, plus AP-1 members (* = DEG)");
  img(s, path.join(DD, "figures", "bridge_genes_heatmap_v2.png"), 0.4, 1.5, 6.9, 5.4, 1779, 1215);
  bullets(s, ["Fos and Jun (plus Junb / Jund) are the only DE-driven bridge between AP-1 and all KEGG anchors; up in BM C0, Spleen C3, BM C4", "Other shared genes are non-DE scaffolds: Mapk1 / 3 / 14, Nfkb1, Rela, Ikbkb", "Chuk, Traf2, Akt2, Map2k7, Grb2 go down in Spleen C3", "Fosb, Dusp1, Nr4a1, Egr1 are strongly DE but do not bridge (Fosb only IL-17; Dusp1 / Nr4a1 only MAPK)", "The energy module shares 0 genes with AP-1 and 0–3 non-DE genes with each KEGG anchor"], { x: 7.5, y: 1.6, w: 5.4, h: 5.2, fontSize: 14 });
  s.addNotes("Sources: drilldown/tables/anchor_bridge_genes.csv (AP-1 vs each KEGG pathway share 3 genes: Fos/Jun + Jund, Junb or Nfatc2; energy module vs AP-1 0 shared, vs KEGG 0-3, none DEG), drilldown/tables/bridge_candidate_genes.csv (per-gene log2FC/FDR in six key comparisons; Fos BM C0 +0.93, Spleen C3 +1.78; Fosb and Dusp1 / Nr4a1 pathway membership), drilldown/figures/bridge_genes_heatmap_v2.png (v2 adds a log2FC colour bar and moves the footer; original kept). Spleen C5-2 columns show large log2FC without significance; BM and spleen C5-2 have no DEG in any drilled set. Direction of Chuk / Traf2 / Akt2 / Map2k7 / Grb2 down in Spleen C3 is from README.md and the heatmap.");
}
// 40 MAPK example
{
  const s = base("Drill-down example: KEGG MAPK signaling", "294 annotated members · 166 detected in the data · 23 DEG in at least one key comparison");
  img(s, path.join(DD, "figures", "gene_network_KEGG_MAPK.png"), 0.4, 1.45, 5.6, 5.45, 2672, 2910);
  const R = [["BM C0", "14"], ["Spleen C3", "14"], ["BM C4", "6"], ["BM C5-2", "0"], ["Spleen C5-2", "0"], ["Thymus C6", "7"]];
  tbl(s, [hdr(["Comparison", "DEG among MAPK members"]), ...R.map((r) => r.map((t, j) => ({ text: t, options: { align: j ? "center" : "left" } })))], { x: 6.4, y: 1.6, w: 6.4, colW: [3.2, 3.2], fontSize: 13, rowH: 0.4 });
  bullets(s, ["Red-boxed stars = DEG; grey hollow = not detected (128 members); panels show BM C0 and Spleen C3", "The 23 DEG mix immediate-early genes (Fos, Jun, Jund, Dusp1/2, Nr4a1, Gadd45g) with a few signalling genes (Chuk, Traf2, Map2k7, Grb2, Akt2); the core kinases (Mapk1/3/14) are not DE", "Pathway-level q = 0.64: a few DEG among many unchanged members"], { x: 6.4, y: 4.4, w: 6.4, h: 2.4, fontSize: 14 });
  s.addNotes("Sources: drilldown/tables/pathway_drill_summary.csv (KEGG MAPK: n_members 294, n_detected 166, DEG BM C0 14, Spleen C3 14, BM C4 6, BM C5-2 0, Spleen C5-2 0, Thymus C6 7, DEG_any_key 23); figure drilldown/figures/gene_network_KEGG_MAPK.png; q = 0.636 from results/tables/anchor_by_method_excl.csv. 294 - 166 = 128 undetected. Note that the pathway-level q is computed per comparison and direction on the enrichment universe, so it is not directly derived from the 23 DEG count; treat that last bullet as an interpretation. Per-pathway figures for IL-17, TCR, Th17, TNF, AP-1 and the energy module are in drilldown/figures/.");
}
// 41 GNN ablation
{
  const s = base("GNN ablation: it adds little on a graph this small", "Held-out edge prediction, mean of seeds 42 / 43 / 44");
  const labels = ["GNN", "Common neighbours", "Adamic-Adar", "Semantic cosine"];
  s.addChart(pres.charts.BAR, [{ name: "(a) all genes", labels, values: [0.973, 0.922, 0.928, 0.840] }, { name: "(b) primary", labels, values: [0.981, 0.922, 0.922, 0.723] }], { x: 0.5, y: 1.5, w: 6.6, h: 4.0, barDir: "col", barGrouping: "clustered", chartColors: [C.teal, C.rust], showLegend: true, legendPos: "b", legendFontSize: 12, legendFontFace: B, valAxisMinVal: 0.5, valAxisMaxVal: 1.0, showValue: true, dataLabelFontSize: 10, dataLabelFormatCode: "0.00", catAxisLabelFontSize: 12, catAxisLabelFontFace: B, valAxisLabelFontSize: 11, valAxisTitle: "AUC", showValAxisTitle: true, valAxisTitleFontSize: 11, valGridLine: { color: "E5E8EA", size: 0.5 } });
  s.addText("Test edges: 23 positive + 23 negative (a), 15 + 15 (b). AUC differences of ~0.05 on so few edges are within seed-to-seed noise.", { x: 0.5, y: 5.6, w: 6.6, h: 0.9, fontFace: B, fontSize: 13, italic: true, color: C.grey, margin: 0, isTextBox: true, valign: "top" });
  const R = [["Seed pair", "(a) all genes", "(b) primary"], ["42 vs 43", "0.71", "0.32"], ["42 vs 44", "0.62", "0.55"], ["43 vs 44", "0.61", "0.36"]];
  tbl(s, R.map((r, i) => i === 0 ? hdr(r) : r.map((t, j) => ({ text: t, options: { align: j ? "center" : "left" } }))), { x: 7.5, y: 1.6, w: 5.3, colW: [1.7, 1.8, 1.8], fontSize: 13, rowH: 0.42 });
  s.addText("Cluster stability across seeds (ARI)", { x: 7.5, y: 1.2, w: 5.3, h: 0.35, fontFace: B, fontSize: 13, bold: true, color: C.grey, margin: 0, isTextBox: true });
  bullets(s, ["GNN AUC beats graph baselines by ~0.05 but on tiny test sets", "In (b) the GNN grouping agrees with GOLDEN fusion (ARI 0.70) and is unstable across seeds (0.32–0.55)", "Conclusion: on ~100 nodes the GNN adds little beyond common neighbours and is less stable; groups are exploratory"], { x: 7.5, y: 3.5, w: 5.3, h: 3.2, fontSize: 14 });
  s.addNotes("Sources: results/tables/gnn_vs_baselines_all.csv and gnn_vs_baselines_excl.csv (AUC means over seeds 42/43/44 - all: GNN 0.973, semantic 0.840, common neighbours 0.922, Adamic-Adar 0.928; excl: 0.981, 0.723, 0.922, 0.922; n_test_pos 23 / 15); results/tables/ARI_across_gnn_seeds_all.csv (0.710, 0.623, 0.614) and ARI_across_gnn_seeds_excl.csv (0.325, 0.549, 0.364); GOLDEN-GNN ARI 0.70 in version b from README.md (results/tables/ARI_between_methods_excl.csv). Transductive link reconstruction on a small graph, not biological validation.");
}
// 42 HTML explorer
{
  const s = base("Interactive explorer", "One HTML file: pathway network ↔ gene network");
  img(s, path.join(DD, "html", "qa", "terms_after_click_hub.png"), 0.4, 1.5, 6.3, 3.5, 1500, 813);
  img(s, path.join(DD, "html", "qa", "genes_from_AP1.png"), 6.7, 1.5, 6.3, 3.5, 1500, 813);
  s.addText("(i) Pathway view: click a node to highlight neighbours; anchor labels visible at load", { x: 0.5, y: 5.05, w: 6.1, h: 0.5, fontFace: B, fontSize: 12, color: C.grey, margin: 0, isTextBox: true, valign: "top" });
  s.addText("(ii) Gene view: opened from the AP-1 node; stars = DEG in the chosen comparison", { x: 6.8, y: 5.05, w: 6.1, h: 0.5, fontFace: B, fontSize: 12, color: C.grey, margin: 0, isTextBox: true, valign: "top" });
  bullets(s, ["Open drilldown/html/inkt_network_explorer.html in a browser (needs internet: vis-network loads from a CDN); Esc or an empty click resets", "Version a / b switch, group filters; button opens the gene view filtered to that pathway; comparison dropdown recolours genes by log2FC; slider sets minimum shared terms"], { x: 0.5, y: 5.7, w: 12.3, h: 1.2, fontSize: 14 });
  s.addNotes("Source: drilldown/html/inkt_network_explorer.html (built by drilldown/code/08_html.py from explorer_template.html); screenshots from headless-Chrome QA in drilldown/html/qa/ (terms_after_click_hub.png, genes_from_AP1.png). Tested: node click highlight and side panel, AP-1 and KEGG MAPK buttons open the filtered gene view (6 and 294 genes), comparison dropdown recolours, edge slider, Esc / empty click reset, no console errors. Fixes on 30 Sept: anchor labels were invisible at load (zoom too small), Esc did not reset. Not tested on Safari/Firefox or offline.");
}
// 43 interpretation & next steps
{
  const s = base("Interpretation and next steps", "What the T2 signature is, and what would test it");
  card(s, 0.5, 1.6, 6.0, 4.9, C.tealLight); card(s, 6.8, 1.6, 6.0, 4.9, C.rustLight);
  s.addText("Interpretation", { x: 0.7, y: 1.7, w: 5.6, h: 0.4, fontFace: H, fontSize: 18, bold: true, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, ["The T2 signature in BM C0 and Spleen C3 is an immediate-early / AP-1 (Fos–Jun) response", "It could reflect niche activation in vivo or tissue dissociation stress; these data cannot separate the two", "Pathway networks add no support for a linked MAPK–IL-17–TCR program; gene-level evidence is Fos / Jun only", "Energy metabolism is the statistically strongest but least specific signal"], { x: 0.7, y: 2.2, w: 5.6, h: 4.2, fontSize: 14 });
  s.addText("Next steps", { x: 7.0, y: 1.7, w: 5.6, h: 0.4, fontFace: H, fontSize: 18, bold: true, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, ["Bench test (Rob): p-ERK, p-p38 and AP-1 (c-Fos / c-Jun) protein in sorted BM C0 and Spleen C3 iNKT cells, with dissociation controls (e.g. transcription-inhibitor during digestion)", "Needed from the collaborators: pool metadata (mouse / batch per pool), the AML dataset, FASTQ for velocity", "Until then all statements stay exploratory (n = 1 pool per label)"], { x: 7.0, y: 2.2, w: 5.6, h: 4.2, fontSize: 14 });
  s.addNotes("Sources: synthesis of slides 37-40; Fos/Jun direction and magnitude from drilldown/tables/bridge_candidate_genes.csv. The dissociation-stress vs niche-activation distinction and the proposed bench design are proposals for discussion, not results. No numbers on this slide are new.");
}
pres.writeFile({ fileName: OUT }).then((f) => console.log("wrote", f));
