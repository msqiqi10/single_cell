# 候选基因被 GO 覆盖了多少（Candidate Coverage by GO）

> Yue 的问题：输入了多少基因，图上只展示了多少，还有多少没被任何显著条目覆盖？

## 定义

“覆盖”（coverage）指完整候选基因列表中，有多少在任何显著术语里出现。未覆盖的原因有：没有注释、有注释但没形成通过阈值的富集条目。“没被覆盖”≠“没有功能”。

## CS 类比

≈ 测试覆盖率：多少代码行被至少一个测试执行过；未覆盖不代表没有用。

## 常见误读与注意

- 图中列出的 driver 只是候选的一部分。

## 在本项目中

- Yue 1:40:50–1:42:06 提出该问题（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 2 节）。
- 旧规则（nominal p、FC≥1.5）：骨髓 C0 输入 166，Top 10 图仅 22 个不同基因（13.3%），在任一已检验 KEGG 条目中出现 43，未覆盖 123；脾脏 C5-1 286/24/89/197；胸腺 C6 188/16/59/129（同上 第 3 节）。
- 新规则（BH FDR≤0.05，|log2FC|≥0.25）BP：骨髓 C0 274 个候选、213 有 BP 注释、66 被显著 BP 覆盖、208 未覆盖；脾脏 C3 275/201/96/179；其中骨髓 C0 有 61 个基因缺当前冻结 BP 背景注释（`iNKT_by_date/2026-09-25/results/README.md` 第 2 节；`iNKT_by_date/2026-09-25/results/tables/full_candidate_GO_coverage.csv`）。
- 图：`iNKT_by_date/2026-09-25/results/figures/01_full_candidate_coverage.png`。
- 9/29：Yue 要求把每个 PAG 内的基因信息展开进网络（gene-level drill-down），标出与 DE 基因重叠的成员（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:30:27–16:35:22）。

## 相关概念

- [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] — GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。
- [[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]] — 很多 GO/KEGG 术语共享同一批基因；Jaccard = 交集/并集，用来量化两个术语有多像。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] — Yue 想要的最终交付：可点击的网络，节点大小 = 显著性，颜色 = 功能组，点开能看到基因。
