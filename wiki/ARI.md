# 调整兰德指数（Adjusted Rand Index, ARI）

> ARI 衡量两种分组的一致程度：1 = 完全相同，约 0 = 随机水平。

## 定义

调整兰德指数（adjusted Rand index, ARI）比较两种聚类结果中“成对的点是否被分在一起”的一致性，并按随机期望校正。

## CS 类比

≈ 两次聚类输出的 diff 相似度，不依赖标签编号。

## 常见误读与注意

- ARI 高只说明结果可重复，不说明分组对应真实生物类别。

## 在本项目中

- C5 拆分：五个种子平均 ARI=0.9868（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 5）。
- GNN 功能分组：三个种子在同一 k 下与参考种子 42 的 ARI 为 0.906、0.800（`iNKT_by_date/2026-09-25/results/README.md` 第 5 节）。
- 09-30 最终：GNN+语义不同种子之间 ARI 版本 (a) 0.61–0.71、主版本 (b) 0.32–0.55；方法间 ARI (a) 最高 0.58、(b) GOLDEN–GNN 0.70（`iNKT_by_date/2026-09-30/results/tables/ARI_across_gnn_seeds_all.csv`、`iNKT_by_date/2026-09-30/results/tables/ARI_between_methods_all.csv`、`ARI_across_gnn_seeds_excl.csv`、`ARI_between_methods_excl.csv`）。

## 相关概念

- [[Leiden 聚类（Leiden Clustering）|Leiden-Clustering]] — Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 09-30 已完成：GO BP/MF + KEGG 混合网络、锚点与分组结果。
