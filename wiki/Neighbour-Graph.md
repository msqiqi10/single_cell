# 近邻图（Neighbour Graph, k-NN）

> 在 PCA 空间给每个细胞找 15 个最相似的邻居，UMAP 和聚类都建立在这张图上。

## 定义

近邻图（neighbour graph）把每个细胞与表达最相似的 k 个细胞相连（k-nearest neighbours）。它是 UMAP、Leiden 聚类、PAGA、diffusion map 共用的输入。

## CS 类比

≈ 相似度图 / 最近邻索引。

## 常见误读与注意

- 图的质量取决于前面的 PCA、HVG 选择和有没有整合。

## 在本项目中

- `n_neighbors=15, n_pcs=30, use_rep="X_pca"`，欧氏距离，random_state 0（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 5.1）。注意：PCA 计算了 50 个成分，邻居图只用前 30 个，而 UMAP PC 比较用“50 PCs”标签——`docs/` 中未见对此的进一步说明，若要引用请先核对（不确定）。
- PAGA 使用同一个 15-neighbour 图（同上 11.1）。

## 相关概念

- [[主成分分析（Principal Component Analysis, PCA）|PCA]] — PCA 把 3,000 维表达压缩成 50 个主轴，供后续邻居图使用。
- [[UMAP 降维可视化（UMAP: why position is not physical）|UMAP]] — UMAP 把近邻图摊平成二维图，位置没有物理意义。
- [[Leiden 聚类（Leiden Clustering）|Leiden-Clustering]] — Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。
- [[拟时序（Diffusion Map, DPT, PAGA）与 RNA 速度（Pseudotime vs RNA Velocity）|Pseudotime-DPT-vs-RNA-Velocity]] — 拟时序按表达相似度给细胞排序；RNA velocity 用剪接/未剪接 RNA 估计变化方向，需要额外输入。
