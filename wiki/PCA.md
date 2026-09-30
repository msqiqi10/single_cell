# 主成分分析（Principal Component Analysis, PCA）

> PCA 把 3,000 维表达压缩成 50 个主轴，供后续邻居图使用。

## 定义

主成分分析（principal component analysis, PCA）找出数据方差最大的正交方向。每个细胞在前 50 个主成分（PC）上的坐标就是它的低维表示。

## CS 类比

≈ 用 SVD 做有损压缩：保留最大的奇异值方向。

## 常见误读与注意

- 单个 PC 只解释 1–2% 方差很常见；没有做 scaling 或整合，前几个 PC 可能同时携带组织、样本、线粒体/核糖体变异。不要把 PC1/PC2 直接称为“iNKT1→iNKT2 分化轴”。

## 在本项目中

- `n_comps=50, svd_solver="arpack", mask_var="highly_variable"`；PC1 2.169%，前 50 PC 共 15.105%（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 4.2）。
- 8/25 比较 50/100/200 PCs 的 UMAP，选定 50；该比较固定 Leiden 标签，只检验显示结构，不是独立重聚类稳定性（`docs/inkt_history_and_supervisor_requirements_20260915.md` 阶段 4；37 页版第 4 页）。

## 相关概念

- [[高变基因（Highly Variable Genes, HVG）|Highly-Variable-Genes]] — 选 3,000 个在细胞间波动最大的基因给 PCA 用，降低噪声和计算量。
- [[近邻图（Neighbour Graph, k-NN）|Neighbour-Graph]] — 在 PCA 空间给每个细胞找 15 个最相似的邻居，UMAP 和聚类都建立在这张图上。
- [[UMAP 降维可视化（UMAP: why position is not physical）|UMAP]] — UMAP 把近邻图摊平成二维图，位置没有物理意义。
- [[批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）|Harmony-and-Batch-Integration]] — 批次效应是技术带来的系统偏差；Harmony 把不同批次的细胞在低维空间“对齐”。
