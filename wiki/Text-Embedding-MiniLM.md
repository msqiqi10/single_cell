# 文本嵌入与 MiniLM（Sentence Embedding, all-MiniLM-L6-v2）

> 把 GO 术语的文字定义变成 384 维向量，语义相近的术语向量也相近。

## 定义

句向量模型（sentence embedding model）把一段文字映射到固定维度的向量。本项目使用 sentence-transformers/all-MiniLM-L6-v2（384 维）对 GO 术语定义嵌入，语义相似度用余弦相似度衡量。

## CS 类比

≈ 把每个条目的描述文本做成搜索引擎的向量索引。

## 常见误读与注意

- 语义相似 ≠ 基因重叠；语义余弦作为链接预测基线 AUC 只有 0.894–0.916，弱于共同邻居。
- 嵌入的是公开的 GO 定义文字，不含研究数据。

## 在本项目中

- 模型在本地运行，公开 GO 定义，没有上传研究数据（`iNKT_by_date/2026-09-25/results/fusion_manifest.json` `semantic_model.purpose`）。
- 嵌入文件：`iNKT_by_date/2026-09-25/results/network/semantic_embeddings.npz`。
- 09-30 脚本对 KEGG 术语只用名称作嵌入文本，对 GO 用“名称+定义”（`iNKT_by_date/2026-09-30/code/02_network.py` 中 `semtext`）。

## 相关概念

- [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] — GOLDEN 是 Yue 团队的“功能条目语义 + 关系图”聚类方法；本项目做的是小鼠 GO 本地适配版。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] — 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。
- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
