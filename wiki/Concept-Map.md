# 概念地图（Concept Map）

Mermaid 图中每个节点都对应一个概念页，下表列出映射。

## 1. 分析主线：基因 → 功能组

```mermaid
flowchart TD
  gene["基因 gene<br/>Gene-mRNA-Protein"] --> DE["差异表达 DE<br/>Differential-Expression"]
  DE --> thr["筛选阈值<br/>Thresholds-Legacy-vs-Robust"]
  thr --> ORA["过表征 ORA<br/>Over-Representation-Analysis"]
  DE --> GSEA["GSEA / NES"]
  ORA --> GO["GO / KEGG 术语"]
  GO --> red["术语冗余 Jaccard"]
  red --> net["术语网络 Term-Network"]
  net --> emb["MiniLM 语义 + 邻接<br/>GOLDEN"]
  net --> gnn["GraphSAGE GNN"]
  emb --> grp["功能组 G01-G07"]
  gnn --> grp
  grp --> exp["实验候选<br/>EXPERIMENT_CANDIDATES"]
```

## 2. 生物学线：Fos/Jun → AP-1 → MAPK → 台面实验

```mermaid
flowchart LR
  fj["Fos / Jun 家族基因<br/>Fos Fosb Fosl2 Jun Junb Jund"] --> ap1["AP-1 转录因子复合体<br/>GO:0035976"]
  mapk["MAPK 级联<br/>ERK / p38 / JNK"] -.激活.-> ap1
  ap1 --> tgt["下游基因：免疫激活、增殖、应激"]
  ieg["即刻早期基因 / 解离应激"] -.替代解释.-> fj
  ap1 --> bench["Rob 的台面读出：<br/>蛋白信号、刺激后 IL-17 分泌"]
  i17["iNKT17：Rorc Il23r Ccr6"] --> bench
  note1["本项目：AP-1 富集显著<br/>MAPK 通路 FDR=1（不显著）"]:::warn
  ap1 -.-> note1
  classDef warn fill:#fff3cd,stroke:#b58900
```

## 3. 统计与证据边界

```mermaid
flowchart TD
  cells["15,532 个细胞"] --> pv["细胞层面 p 值"]
  pv --> fdr["BH FDR / q"]
  fdr --> claim["‘显著’"]
  pool["每组 n=1 个 pool<br/>3 只小鼠混合"] --> pseudo["伪重复：不是动物层面证据"]
  claim --> pseudo
  pseudo --> qc["QC 匹配 / 敏感性 / 排除分析"]
  qc --> plaus["方向稳健 ≠ 生物学重复"]
  plaus --> wet["需要独立动物 + 台面验证"]
```

## 4. iNKT 亚型与细胞因子

```mermaid
flowchart LR
  inkt["iNKT 细胞"] --> s1["iNKT1<br/>Tbx21 / Ifng"]
  inkt --> s2["iNKT2<br/>Gata3 / Il4"]
  inkt --> s17["iNKT17<br/>Rorc / Il23r / Ccr6 / Il17"]
  s1 --> c1["IFN-γ"]
  s2 --> c2["IL-4"]
  s17 --> c17["IL-17"]
  s17 --> c52["本项目：C5-2 有 iNKT17 相关特征<br/>Il17a 保留细胞中 0 个检出"]
```

## 5. 09-30：从通路到基因

```mermaid
flowchart TD
  ora["ORA：GO BP / MF + KEGG"] --> nodes["显著节点 + 12 个锚点"]
  nodes --> net["混合网络 a / b<br/>Mixed-GO-KEGG-Network-0930"]
  net --> grp["分组：Louvain / 语义 / GOLDEN / GNN"]
  grp --> catch["锚点同组 = 杂项组<br/>不是功能联系"]
  nodes --> anc["锚点节点<br/>Anchor-Nodes"]
  anc --> drill["基因层面钻取<br/>Gene-Level-Drilldown"]
  drill --> bridge["桥接基因 Fos + Jun<br/>Bridge-Genes"]
  drill --> html["交互 HTML<br/>Pathway-Network-Explorer-HTML"]
  bridge --> ieg["即刻早期 / 解离应激？<br/>需台面实验"]
```

## 节点 → 页面

| 图中节点 | 页面 |
|---|---|
| 基因 | [[基因、mRNA 与蛋白质（Gene / mRNA / Protein, Central Dogma）|Gene-mRNA-Protein]] |
| DE / 阈值 | [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] [[两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）|Thresholds-Legacy-vs-Robust]] |
| ORA / GSEA | [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] [[基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）|GSEA-and-NES]] |
| GO / KEGG | [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] |
| Jaccard 冗余 | [[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]] |
| 术语网络 / GOLDEN / GNN | [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] |
| 功能组 | [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] |
| Fos/Jun / AP-1 / MAPK | [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] |
| 解离应激 | [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] |
| iNKT 亚型 / 细胞因子 | [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] |
| 伪重复 / QC 匹配 / 敏感性 | [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] [[QC 匹配与标准化均值差（QC Matching, SMD）|QC-Matching-and-SMD]] [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] |
| 台面验证 | [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]]（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md`） |
| 混合网络 / 锚点 | [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] [[锚点节点（Anchor Nodes）|Anchor-Nodes]] |
| 桥接基因 / 基因钻取 | [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] |
| 交互网络浏览器 | [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] |
