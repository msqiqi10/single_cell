# 术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）

> 很多 GO/KEGG 术语共享同一批基因；Jaccard = 交集/并集，用来量化两个术语有多像。

## 定义

Jaccard 指数 J(A,B) = |A∩B| / |A∪B|，取值 0–1。两个术语的成员基因集合 Jaccard 高，说明它们本质是同一批基因的不同名字。

## CS 类比

≈ 两个代码文件的 token 集合相似度（查重）。

## 常见误读与注意

- 条目数量 ≠ 独立机制数量：骨髓 C0 的 52 个显著 BP 名称只覆盖 66 个不同基因，反复复用 Cox、Nduf、Atp 基因。

## 在本项目中

- 显示用：Jaccard≥0.75 去冗余，完整结果保留（`iNKT_by_date/2026-09-19/README.md`）。
- 网络边：Jaccard≥0.25 且共享基因≥3，得到 97 节点 577 条边（`iNKT_by_date/2026-09-25/results/README.md` 第 5 节；`iNKT_by_date/2026-09-25/results/network/mouse_GO_edges.tsv`）。
- GO:0002181 与 GO:0006412 的 Jaccard=0.507（194 个共享基因，主要是 Rpl/Rps）（`iNKT_by_date/2026-09-25/results/network/mouse_GO_edges.tsv`）。
- legacy 与 current DEG 集合的 Recovery、Jaccard（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 9.2）也是同一个公式。

## 相关概念

- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] — GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] — Yue 的问题：输入了多少基因，图上只展示了多少，还有多少没被任何显著条目覆盖？
