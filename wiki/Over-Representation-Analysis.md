# 过表征分析（Over-Representation Analysis, ORA; hypergeometric test）

> ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。

## 定义

过表征分析（over-representation analysis, ORA）用超几何检验（hypergeometric test，等价于单侧 Fisher 精确检验）：背景基因共 M 个，其中属于某术语的有 K 个；从背景里抽 n 个入选基因，其中 k 个属于该术语。p = P(命中数 ≥ k)。背景宇宙（background universe）应是“实际可检出且有注释”的基因。

## CS 类比

≈ 从一个装了 M 个球（K 个红色）的罐子里不放回抽 n 个，问抽到 ≥k 个红球有多罕见。

## 常见误读与注意

- 背景选错会得到假富集；本项目背景 = 实际检出基因 ∩ 该 GO 分支已注释基因。
- 富集不等于激活、不等于因果，且同一批基因可以撑起很多术语。

## 在本项目中

- 实现：`scipy.stats.hypergeom.sf(k-1, M, K, n)`，词条大小 5–500，上调/下调分开，基因筛选 gene BH≤0.05 且 |log2FC|≥0.25（`iNKT_by_date/2026-09-19/code/analyze.py`；`iNKT_by_date/2026-09-19/README.md`）。
- 示例：AP-1 复合体，骨髓 C0：M=9562，K=6，n=225，k=4，p=4.3e-6，家族 q=0.000261（`iNKT_by_date/2026-09-25/results/tables/AP1_complex_GO_exact_ID.csv`）。
- 早期 ORA：19 单元×2 方向×12 signature=456 次检验（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 10.2）。
- 09-30：同一规则再加 KEGG，并复核了对 09-19 结果的复用（312,550 行 q 值最大差 1.1e-16）：`iNKT_by_date/2026-09-30/results/tables/ORA_reuse_check_vs_0919.csv`。

## 相关概念

- [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] — GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
- [[基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）|GSEA-and-NES]] — GSEA 不需要先选“显著基因”，而是把全部基因按变化排序，看基因集是否集中在顶端或底端；NES 带方向。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
- [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] — Yue 的问题：输入了多少基因，图上只展示了多少，还有多少没被任何显著条目覆盖？
