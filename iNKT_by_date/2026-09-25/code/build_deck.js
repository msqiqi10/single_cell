// Builds the Sept 25 follow-up deck (continues after the 19-page cytotoxicity/GO deck).
// Usage: NODE_PATH=<dir with pptxgenjs> node build_deck.js
const path = require("path");
const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const FIG = (f) => path.join(ROOT, "results", "figures", f);
const OUT = path.join(ROOT, "presentation", "iNKT_followup_AP1_iNKT17_GO_network_20260925.pptx");

const C = {
  ink: "1F2A30", teal: "137D8A", tealLight: "D7ECEE", rust: "C8694A",
  rustLight: "F6E2DA", grey: "5F6B72", line: "D5DADD", bg: "FFFFFF", dark: "12343A",
};
const H = "Cambria", B = "Calibri";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.title = "iNKT follow-up: AP-1, iNKT17 and GO functional grouping";

let pageNo = 19; // continues the 19-page cytotoxicity / GO deck
function base(title, subtitle) {
  const s = pres.addSlide();
  s.background = { color: C.bg };
  pageNo += 1;
  s.addText(title, { x: 0.5, y: 0.3, w: 12.3, h: 0.7, fontFace: H, fontSize: 28, bold: true, color: C.ink, margin: 0, isTextBox: true });
  if (subtitle) s.addText(subtitle, { x: 0.5, y: 0.98, w: 12.3, h: 0.4, fontFace: B, fontSize: 14, color: C.grey, margin: 0, isTextBox: true });
  s.addText(`${pageNo}`, { x: 12.3, y: 7.0, w: 0.6, h: 0.3, fontFace: B, fontSize: 10, color: C.grey, align: "right", margin: 0, isTextBox: true });
  s.addText("Exploratory cell-level results · 3 mice pooled per tissue × condition label (n = 1 pool)", { x: 0.5, y: 7.0, w: 9, h: 0.3, fontFace: B, fontSize: 10, color: C.grey, margin: 0, isTextBox: true });
  return s;
}
function card(s, x, y, w, h, fill) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: fill }, rectRadius: 0.08 });
}
function bullets(s, items, opt) {
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1 } })),
    Object.assign({ fontFace: B, fontSize: 15, color: C.ink, paraSpaceAfter: 6, valign: "top", isTextBox: true }, opt));
}

// 20 — Section title
{
  const s = pres.addSlide();
  pageNo += 1;
  s.background = { color: C.dark };
  s.addText("Follow-up to the 24 Sept meeting", { x: 0.8, y: 2.0, w: 11.5, h: 0.9, fontFace: H, fontSize: 40, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
  s.addText("AP-1 / MAPK, iNKT17 and GO functional grouping — continuing after the cytotoxicity and ontology analysis", { x: 0.8, y: 3.0, w: 11.5, h: 0.9, fontFace: B, fontSize: 20, color: "CFE6E8", margin: 0, isTextBox: true });
  s.addText("Same data: 15,532 iNKT cells · BM / spleen / thymus · Ctrl vs T2 · analysis of 25 Sept 2026", { x: 0.8, y: 4.3, w: 11.5, h: 0.5, fontFace: B, fontSize: 15, color: "9FC3C7", margin: 0, isTextBox: true });
  s.addNotes("This section continues directly after the 19 cytotoxicity and GO slides. Everything here uses the same dataset; what changed is the question: we followed Rob's advice to focus on MAPK, AP-1 and the iNKT17 signature, and Yue's advice to use the full candidate lists and group GO terms into functional modules.");
}

// 21 — What was asked vs done
{
  const s = base("What the meeting asked for, and what we did", "Rob: focus on biology seen at the bench · Yue: use full candidate lists and group ontology terms");
  const rows = [
    ["Rob", "Focus on MAPK, AP-1 (Fos/Jun family), iNKT17", "Gene panel across 17 clusters; exact-ID GO for AP-1 complex; MAPK re-tested"],
    ["Rob", "Ribosome / stress may reflect an inflamed niche", "20 QC-matched resamples; GO re-run excluding Rpl/Rps/Hsp/Dnaj genes"],
    ["Rob", "Samples: 3 mice pooled per label", "Recorded as design; still n = 1 pool per tissue × condition"],
    ["Yue", "Look beyond top genes; count pathway coverage", "All 26 DE tables kept; coverage of every T2-up candidate"],
    ["Yue", "Group GO terms; GOLDEN / Fusion / GNN", "97 BP nodes, 577 edges → 7 exploratory functional groups"],
    ["Yue", "Keep ontology + cytotoxicity direction", "Integrated with BM C4 / spleen C3 cytotoxic signals"],
  ];
  const hdr = ["From", "Request", "Done"].map((t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: C.teal } } }));
  s.addTable([hdr, ...rows.map((r) => r.map((t, i) => ({ text: t, options: { bold: i === 0, color: C.ink } })))], {
    x: 0.5, y: 1.6, w: 12.3, colW: [1.0, 4.9, 6.4], fontFace: B, fontSize: 14, border: { type: "solid", pt: 0.75, color: C.line }, rowH: 0.62, valign: "middle",
  });
  s.addText("Pending (external): AML NK data via Box, FASTQ for RNA velocity, wet-lab validation.", { x: 0.5, y: 6.3, w: 12.3, h: 0.4, fontFace: B, fontSize: 14, italic: true, color: C.rust, margin: 0, isTextBox: true });
  s.addNotes("Left column is who raised it, middle is the request, right is what was actually computed. Three items depend on data or experiments we do not have yet: the AML dataset Rob offered, FASTQ files needed for RNA velocity, and any wet-lab validation.");
}

// 22 — Rob gene panel
{
  const s = base("Rob's candidate genes across every eligible cluster", "log2FC T2 / Ctrl; * = gene BH FDR ≤ 0.05 and |log2FC| ≥ 0.25");
  s.addImage({ path: FIG("04_Rob_gene_panel.png"), x: 0.5, y: 1.5, w: 8.0, h: 5.4, sizing: { type: "contain", w: 8.0, h: 5.4 } });
  card(s, 8.9, 1.6, 3.9, 5.1, C.tealLight);
  s.addText("Reading", { x: 9.15, y: 1.8, w: 3.5, h: 0.4, fontFace: H, fontSize: 18, bold: true, color: C.teal, margin: 0, isTextBox: true });
  bullets(s, [
    "AP-1 family (Fos, Jun, Junb, Fosb, Fosl2) and Dusp1 rise in T2 in BM C0 and spleen C3",
    "Same direction, fewer passing genes, in BM C4, spleen C4 / C5-1, thymus C6",
    "Il1r1, Il6ra, Cd8a, Ciita: large colours but no FDR support — sparse detection",
  ], { x: 9.15, y: 2.3, w: 3.5, h: 4.2, fontSize: 14 });
  s.addNotes("Rows are tissue by cluster, columns are the genes Rob named. Red means higher in T2. Stars mark genes that pass the same DE rule as before. The AP-1 block is consistent in BM C0 and spleen C3. Some of the non-AP-1 genes Rob mentioned show strong colours but these come from very few expressing cells and do not pass FDR, so we do not list them as significant.");
}

// 23 — AP-1 vs MAPK
{
  const s = base("AP-1 complex is enriched; MAPK cascade is not", "GO:0035976 'transcription factor AP-1 complex' — T2-up genes, exact GO ID (not 'AP-1 adaptor complex')");
  const stats = [
    { x: 0.5, t: "BM C0", q: "q = 0.0029", sub: "global FDR (family 0.00026)", g: "Fos · Fosl2 · Jun · Junb" },
    { x: 4.7, t: "Spleen C3", q: "q = 0.0023", sub: "global FDR (family 0.00019)", g: "Fos · Jun · Junb · Jund" },
  ];
  stats.forEach((d) => {
    card(s, d.x, 1.7, 3.9, 3.2, C.tealLight);
    s.addText(d.t, { x: d.x + 0.3, y: 1.9, w: 3.3, h: 0.45, fontFace: H, fontSize: 20, bold: true, color: C.teal, margin: 0, isTextBox: true });
    s.addText(d.q, { x: d.x + 0.3, y: 2.45, w: 3.3, h: 0.9, fontFace: H, fontSize: 40, bold: true, color: C.ink, margin: 0, isTextBox: true });
    s.addText(d.sub, { x: d.x + 0.3, y: 3.35, w: 3.3, h: 0.4, fontFace: B, fontSize: 13, color: C.grey, margin: 0, isTextBox: true });
    s.addText(d.g, { x: d.x + 0.3, y: 3.95, w: 3.3, h: 0.6, fontFace: B, fontSize: 15, italic: true, color: C.ink, margin: 0, isTextBox: true });
  });
  card(s, 8.9, 1.7, 3.9, 3.2, C.rustLight);
  s.addText("MAPK", { x: 9.2, y: 1.9, w: 3.3, h: 0.45, fontFace: H, fontSize: 20, bold: true, color: C.rust, margin: 0, isTextBox: true });
  s.addText("q = 1", { x: 9.2, y: 2.45, w: 3.3, h: 0.9, fontFace: H, fontSize: 40, bold: true, color: C.ink, margin: 0, isTextBox: true });
  s.addText("GO 'MAPK cascade', BM C0", { x: 9.2, y: 3.35, w: 3.3, h: 0.4, fontFace: B, fontSize: 13, color: C.grey, margin: 0, isTextBox: true });
  s.addText("KEGG MAPK: q = 0.078 (legacy rule), not significant under the robust rule", { x: 9.2, y: 3.85, w: 3.3, h: 0.9, fontFace: B, fontSize: 13, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, [
    "Supports an AP-1 transcriptional programme in T2 — not protein activation or DNA binding",
    "AP-1 enrichment should not be used to relabel MAPK as significant; p-ERK / p-p38 would be the bench test",
    "Only 6 genes annotate the AP-1 complex, so this is a focused, small-term result",
  ], { x: 0.5, y: 5.2, w: 12.3, h: 1.6, fontSize: 15 });
  s.addNotes("This is the strongest new result. The AP-1 complex term is enriched among T2-up genes in both BM C0 and spleen C3, and survives the global correction across all comparisons. However, the MAPK cascade term itself is not significant. So the honest statement is: an AP-1 transcriptional response in T2, consistent with Rob's bench experience, while MAPK activation remains to be tested at the protein level.");
}

// 24 — QC matching chart
{
  const s = base("AP-1 genes keep their direction after QC matching", "Mean log1p expression difference (T2 − Ctrl) in the full data; all 20 QC-matched resamples agree in direction");
  const genes = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Dusp1"];
  s.addChart(pres.charts.BAR, [
    { name: "BM C0", labels: genes, values: [0.226, 0.203, 0.282, 0.376, 0.101, 0.249] },
    { name: "Spleen C3", labels: genes, values: [0.168, 0.202, 0.498, 0.391, 0.186, 0.112] },
  ], {
    x: 0.5, y: 1.5, w: 8.2, h: 5.3, barDir: "col", barGrouping: "clustered", chartColors: [C.teal, C.rust],
    showLegend: true, legendPos: "t", legendFontSize: 13, legendFontFace: B,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 10, dataLabelFormatCode: "0.00",
    catAxisLabelFontSize: 13, valAxisLabelFontSize: 11, catAxisLabelColor: C.ink, valAxisLabelColor: C.grey,
    valGridLine: { color: "E6E9EB", size: 0.5 }, catGridLine: { style: "none" },
    showValAxisTitle: true, valAxisTitle: "Δ mean log1p (T2 − Ctrl)", valAxisTitleFontSize: 12,
  });
  card(s, 9.1, 1.7, 3.7, 4.9, C.tealLight);
  s.addText("20 / 20", { x: 9.35, y: 1.95, w: 3.2, h: 0.9, fontFace: H, fontSize: 44, bold: true, color: C.teal, margin: 0, isTextBox: true });
  s.addText("resamples same direction for every gene shown", { x: 9.35, y: 2.85, w: 3.2, h: 0.7, fontFace: B, fontSize: 14, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, [
    "Ctrl and T2 matched on QC metrics (max SMD ≈ 0.06)",
    "Rules out a simple QC-depth artefact",
    "Cannot separate niche biology from dissociation / processing time",
  ], { x: 9.35, y: 3.75, w: 3.2, h: 2.7, fontSize: 14 });
  s.addNotes("Rob pointed out that some signals may just reflect the stressful tumour environment. Immediate-early genes like Fos and Jun are also known to rise with dissociation. Matching cells on QC metrics does not change the direction, which rules out a sequencing-depth artefact, but it cannot tell niche biology from processing effects. Processing records per sample would help.");
}

// 25 — iNKT17
{
  const s = base("iNKT17: identity markers in C5-2, but no IL-17 transcript", "Re-read all 32,285 features from the six raw matrices to recover genes dropped by the min_cells filter");
  s.addImage({ path: FIG("05_full_feature_iNKT17.png"), x: 0.5, y: 1.5, w: 8.6, h: 3.9, sizing: { type: "contain", w: 8.6, h: 3.9 } });
  const st = [["60", "retained cells express Ccr6"], ["54", "of these are in C5-2"], ["0", "cells with Il17a"]];
  st.forEach((d, i) => {
    const y = 1.6 + i * 1.3;
    card(s, 9.5, y, 3.3, 1.1, i === 2 ? C.rustLight : C.tealLight);
    s.addText(d[0], { x: 9.7, y: y + 0.1, w: 1.1, h: 0.9, fontFace: H, fontSize: 34, bold: true, color: i === 2 ? C.rust : C.teal, valign: "middle", margin: 0, isTextBox: true });
    s.addText(d[1], { x: 10.8, y: y + 0.1, w: 1.9, h: 0.9, fontFace: B, fontSize: 13, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
  });
  bullets(s, [
    "Ccr6 was missing only because of the min_cells = 100 filter (92 raw cells), not absent from the data",
    "Rorc / Il23r / Ccr6 support an iNKT17 identity for C5-2; Il17f in 2 cells, both outside C5-2",
    "Transcript detection ≠ secretion: IL-17 capacity needs stimulation assays; thymus Ctrl C5-2 has only 11 cells",
  ], { x: 0.5, y: 5.55, w: 12.3, h: 1.4, fontSize: 14 });
  s.addNotes("Rob highlighted the Th17-like signature. Rorc and Il23r were already there; Ccr6 had been filtered out and is now recovered, with almost all expressing cells in C5-2. Il17a and Il17f are essentially undetected, which is typical for unstimulated cells, so we cannot comment on IL-17 secretion. That is exactly the kind of thing Rob's stimulation experiments can answer.");
}

// 26 — Coverage
{
  const s = base("Most T2-up genes are not explained by significant GO terms", "Yue: look past the top genes — coverage of every selected candidate (gene BH FDR ≤ 0.05, log2FC ≥ 0.25)");
  s.addImage({ path: FIG("01_full_candidate_coverage.png"), x: 0.5, y: 1.5, w: 7.6, h: 5.4, sizing: { type: "contain", w: 7.6, h: 5.4 } });
  card(s, 8.5, 1.6, 4.3, 5.1, C.tealLight);
  s.addText("66 / 274", { x: 8.75, y: 1.8, w: 3.8, h: 0.8, fontFace: H, fontSize: 36, bold: true, color: C.teal, margin: 0, isTextBox: true });
  s.addText("BM C0 T2-up genes inside 52 significant BP terms", { x: 8.75, y: 2.6, w: 3.8, h: 0.6, fontFace: B, fontSize: 13, color: C.ink, margin: 0, isTextBox: true });
  s.addText("96 / 275", { x: 8.75, y: 3.35, w: 3.8, h: 0.8, fontFace: H, fontSize: 36, bold: true, color: C.teal, margin: 0, isTextBox: true });
  s.addText("Spleen C3", { x: 8.75, y: 4.15, w: 3.8, h: 0.4, fontFace: B, fontSize: 13, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, [
    "Many terms reuse the same genes",
    "Uncovered genes stay in the full tables for follow-up",
  ], { x: 8.75, y: 4.75, w: 3.8, h: 1.8, fontSize: 14 });
  s.addNotes("For each comparison, the bar is all selected T2-up genes. Dark teal is genes inside at least one significant biological-process term. In BM C0 only 66 of 274 genes are covered: the 52 significant terms reuse the same genes. Light grey genes have no BP annotation at all. This is why grouping terms, rather than listing them, is useful.");
}

// 27 — Network
{
  const s = base("Grouping GO terms into functional modules", "97 significant BP terms · edges = shared measured genes (Jaccard) · GOLDEN-style fusion + 3 GraphSAGE seeds");
  s.addImage({ path: FIG("03_function_network.png"), x: 0.5, y: 1.5, w: 8.4, h: 4.6, sizing: { type: "contain", w: 8.4, h: 4.6 } });
  const g = [["G02 / G03", "Translation, ribosome"], ["G04 / G05 / G06", "Respiration, ATP, OXPHOS"], ["G07", "Ion / proton transport"], ["G01", "Mixed — inspect drivers"]];
  g.forEach((d, i) => {
    const y = 1.6 + i * 1.05;
    card(s, 9.3, y, 3.5, 0.9, i === 3 ? "EEF0F1" : C.tealLight);
    s.addText([{ text: d[0], options: { bold: true, color: C.teal, breakLine: true } }, { text: d[1], options: { color: C.ink } }],
      { x: 9.45, y: y + 0.05, w: 3.25, h: 0.8, fontFace: B, fontSize: 14, valign: "middle", margin: 0, isTextBox: true });
  });
  s.addText("Semantic / GNN agreement ARI 0.91 / 0.80. Mouse GO adaptation of GOLDEN — not a full reproduction of the published model; edges are shared genes, not regulation.",
    { x: 0.5, y: 6.2, w: 12.3, h: 0.6, fontFace: B, fontSize: 13, italic: true, color: C.grey, margin: 0, isTextBox: true });
  s.addNotes("Following Yue's suggestion we built a network of the significant GO terms, linked when they share genes, then grouped them with a GOLDEN-style fusion of term descriptions and graph structure, plus three GraphSAGE runs. Seven groups emerge. Labels like 'synapse' or 'muscle' appear inside ribosome groups because they share ribosomal genes; they are not iNKT functions.");
}

// 28 — Groups by cluster + sensitivity
{
  const s = base("BM C0 and spleen C3 carry energy-metabolism groups", "Unique T2-up driver genes per functional group; dash = no significant evidence");
  s.addImage({ path: FIG("02_function_groups_by_cluster.png"), x: 0.5, y: 1.5, w: 6.9, h: 5.4, sizing: { type: "contain", w: 6.9, h: 5.4 } });
  card(s, 7.8, 1.6, 5.0, 2.3, C.tealLight);
  s.addText("Ribosome / chaperone exclusion test", { x: 8.05, y: 1.75, w: 4.6, h: 0.4, fontFace: H, fontSize: 17, bold: true, color: C.teal, margin: 0, isTextBox: true });
  s.addText("After removing Rpl / Rps / Hsp / Dnaj genes, OXPHOS stays enriched in BM C0 (29 genes: Atp5, Cox, Nduf, Uqcr). Cytoplasmic translation and protein folding drop out.",
    { x: 8.05, y: 2.2, w: 4.6, h: 1.6, fontFace: B, fontSize: 14, color: C.ink, margin: 0, isTextBox: true });
  bullets(s, [
    "BM C0 and spleen C3 load on respiration / OXPHOS groups (G04–G06)",
    "Thymus C6, BM C5-2, spleen C4 are dominated by ribosome groups — consistent with Rob's 'stressed cell' reading",
    "OXPHOS is a supplementary lead, below AP-1 and iNKT17 in priority",
  ], { x: 7.8, y: 4.15, w: 5.0, h: 2.7, fontSize: 14 });
  s.addNotes("Mapping groups back to clusters: BM C0 and spleen C3, the same populations with the AP-1 signal, also carry respiration and ATP-production groups. When we remove ribosomal and chaperone genes, the OXPHOS enrichment remains while translation and folding disappear, so the energy signal is not an artefact of the stress genes Rob warned about. We keep it as a supplementary lead.");
}

// 29 — Integration with cytotoxicity
{
  const s = base("Putting it together with the cytotoxicity results", "Same populations, different readouts — candidates for the bench");
  const hdr = ["Population", "AP-1 (this section)", "Cytotoxic core (Prf1/Gzma/Gzmb)", "Other"].map((t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: C.teal } } }));
  const rows = [
    ["BM C0", "Enriched, q = 0.003", "n.s.", "OXPHOS up; MAPK n.s."],
    ["Spleen C3", "Enriched, q = 0.002", "Down, Δ −0.023, q = 0.0015", "CCT / folding down in T2"],
    ["BM C4", "Fos / Junb / Fosb / Dusp1 up", "Up, Δ +0.086, q = 0.025", "Tcf7 / Itga4 / Emb down"],
    ["C5-2", "—", "n.s.", "iNKT17 identity; IL-17 undetected"],
  ];
  s.addTable([hdr, ...rows.map((r) => r.map((t, i) => ({ text: t, options: { bold: i === 0, color: C.ink } })))], {
    x: 0.5, y: 1.6, w: 12.3, colW: [2.0, 3.3, 3.8, 3.2], fontFace: B, fontSize: 14, border: { type: "solid", pt: 0.75, color: C.line }, rowH: 0.7, valign: "middle",
  });
  card(s, 0.5, 5.4, 12.3, 1.3, C.rustLight);
  s.addText("Spleen C3 is the interesting case: AP-1 transcription up while the cytotoxic core goes down. Whether activation and effector function diverge needs protein / killing readouts in matched cells — no causal claim from RNA.",
    { x: 0.75, y: 5.5, w: 11.8, h: 1.1, fontFace: B, fontSize: 15, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
  s.addNotes("Combining this section with the cytotoxicity slides: spleen C3 shows an AP-1 increase together with a small decrease in the Prf1/Gzma/Gzmb core. BM C4 shows the opposite pattern for cytotoxicity. These are candidates for experiments; the RNA alone does not show that AP-1 changes killing.");
}

// 30 — Proposed experiments
{
  const s = base("Proposed questions for the bench", "Discussion list — not experiments already done; unit of analysis = independent animal or pool");
  const items = [
    ["A", "AP-1 in BM C0 / spleen C3", "Sortable phenotype → AP-1 protein, p-ERK / p-p38 with and without standard stimulation"],
    ["B", "iNKT17 (C5-2)", "Separate identity, abundance and per-cell IL-17 secretion; align with existing IL-17 assays"],
    ["C", "AP-1 vs effector function", "Spleen C3 and BM C4: degranulation / killing readout in the same phenotype"],
    ["D", "Energy metabolism (supplementary)", "Mitochondrial / metabolic readouts if cell numbers allow"],
  ];
  items.forEach((d, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = 0.5 + col * 6.25, y = 1.6 + row * 2.55;
    card(s, x, y, 6.05, 2.3, i === 3 ? "EEF0F1" : C.tealLight);
    s.addShape(pres.shapes.OVAL, { x: x + 0.25, y: y + 0.25, w: 0.7, h: 0.7, fill: { color: i === 3 ? C.grey : C.teal }, line: { color: i === 3 ? C.grey : C.teal } });
    s.addText(d[0], { x: x + 0.25, y: y + 0.25, w: 0.7, h: 0.7, fontFace: H, fontSize: 22, bold: true, color: "FFFFFF", align: "center", valign: "middle", margin: 0, isTextBox: true });
    s.addText(d[1], { x: x + 1.15, y: y + 0.25, w: 4.7, h: 0.7, fontFace: H, fontSize: 18, bold: true, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
    s.addText(d[2], { x: x + 0.3, y: y + 1.1, w: 5.5, h: 1.0, fontFace: B, fontSize: 14, color: C.ink, margin: 0, isTextBox: true });
  });
  s.addNotes("These are questions for Rob's team to choose from, ordered by priority. A and B follow directly from Rob's comments; C connects to the cytotoxicity analysis; D is a supplementary lead from the GO network.");
}

// 31 — Needs / next steps
{
  const s = base("What we need next", "Items that block the next analysis or validation stage");
  const needs = [
    ["Sample metadata", "Which 3 mice form each label; any independent pools; T2 time point, dissociation time, sequencing batch"],
    ["IL-17 assay details", "Cell source, gating, stimulation, time point and replicates of existing experiments"],
    ["AML NK / patient data", "Box access; species, disease / control definitions, NK-only vs whole marrow, FASTQ"],
    ["RNA velocity", "FASTQ → spliced / unspliced counts (current DPT is ordering only, not velocity)"],
  ];
  needs.forEach((d, i) => {
    const y = 1.6 + i * 1.3;
    card(s, 0.5, y, 12.3, 1.1, i % 2 ? "FFFFFF" : C.tealLight);
    s.addText(d[0], { x: 0.8, y, w: 3.4, h: 1.1, fontFace: H, fontSize: 18, bold: true, color: C.teal, valign: "middle", margin: 0, isTextBox: true });
    s.addText(d[1], { x: 4.3, y, w: 8.3, h: 1.1, fontFace: B, fontSize: 15, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
  });
  s.addNotes("To move from exploratory cell-level results to animal-level conclusions we need the pool composition and processing records. For the new directions Rob offered, we need Box access to the AML datasets and FASTQ files, which also enable RNA velocity.");
}

pres.writeFile({ fileName: OUT }).then(() => console.log(OUT));
