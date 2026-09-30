# 术语网络（Term Network: nodes = terms, edges = shared genes）

> 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。

## 定义

术语网络（term network）中节点是 GO/KEGG 术语（不是基因、也不是细胞），边表示两术语的成员基因足够重叠（本项目：Jaccard≥0.25 且共享基因≥3）。它不是基因调控网络，边也不表示因果。

## CS 类比

≈ 模块依赖图：两个模块如果共用了很多同一批底层库，就画一条边。

## 常见误读与注意

- 一张图里不要把节点含义混用（术语、基因、功能组）再赋予因果箭头。
- 边由“测得的基因成员”决定，因此对基因宇宙（10,670 个基因）敏感。

## 在本项目中

- 97 个 BP 节点（至少在一个主比较中通过家族 FDR 的 GO BP 条目），577 条边（`iNKT_by_date/2026-09-25/results/network/mouse_GO_edges.tsv` 共 577 行数据；节点表 `iNKT_by_date/2026-09-25/results/network/mouse_GO_nodes.csv`）。
- 隔离点不丢弃；只对 BP 分组以避免混合“过程/功能/组分”层级；AP-1 的 CC 结果单独呈现（`iNKT_by_date/2026-09-25/results/README.md` 第 5 节）。
- 图：`iNKT_by_date/2026-09-25/results/figures/03_function_network.png`；组间共享边：`iNKT_by_date/2026-09-25/results/network/module_edges.csv`（G04–G06 共 155 条术语边）。
- 9/24 Yue：GO 不是终点，要把许多条目聚成少数功能组（转录稿 01:42:06–01:44:59，见 `iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 2 节）。
- 9/29 Yue：97 个节点、577 条边“会不会太小”——关键不在大小，而在这些节点是否重要；还要加入 GO MF 和 KEGG（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:20:15–16:20:35、16:36:29–16:36:57）。

## 相关概念

- [[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]] — 很多 GO/KEGG 术语共享同一批基因；Jaccard = 交集/并集，用来量化两个术语有多像。
- [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] — GOLDEN 是 Yue 团队的“功能条目语义 + 关系图”聚类方法；本项目做的是小鼠 GO 本地适配版。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] — Yue 想要的最终交付：可点击的网络，节点大小 = 显著性，颜色 = 功能组，点开能看到基因。
- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 09-30 已完成：GO BP/MF + KEGG 混合网络、锚点与分组结果。
- [[锚点节点（Anchor Nodes）|Anchor-Nodes]] — 网络里固定加入、不论显著与否都保留的节点。
