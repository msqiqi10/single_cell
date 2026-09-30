# UMAP 降维可视化（UMAP: why position is not physical）

> UMAP 把近邻图摊平成二维图，位置没有物理意义。

## 定义

UMAP（uniform manifold approximation and projection）把高维邻居关系投影到二维用于可视化。相邻的点通常表达相似；但轴没有单位，岛间距离不是发育时间，岛面积不代表细胞数。

## CS 类比

≈ 把社交网络图力导向布局后截图：邻接关系可信，绝对坐标和岛屿大小不可信。

## 常见误读与注意

- 参数（PC 数、邻居数、随机种子）改变，局部形状会变（转录稿 00:39:45–00:40:03：ZERU 也强调 PC 对比不能确定哪种“最好”）。
- 不要读“UMAP 上的路径”为发育轨迹。

## 在本项目中

- `sc.tl.umap(random_state=0)`（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.2）。
- 37 页版第 3、4 页：50/100/200 PCs 与样本/群/条件着色（`iNKT_by_date/2026-09-20/package/slides.json`）。
- 9/9 会议主要讨论 R01 UMAP 的含义（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 9）。

## 相关概念

- [[近邻图（Neighbour Graph, k-NN）|Neighbour-Graph]] — 在 PCA 空间给每个细胞找 15 个最相似的邻居，UMAP 和聚类都建立在这张图上。
- [[Leiden 聚类（Leiden Clustering）|Leiden-Clustering]] — Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。
- [[批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）|Harmony-and-Batch-Integration]] — 批次效应是技术带来的系统偏差；Harmony 把不同批次的细胞在低维空间“对齐”。
- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
