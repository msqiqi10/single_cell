# Leiden 聚类（Leiden Clustering）

> Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。

## 定义

Leiden 算法在图上寻找连接紧密的社区（community）。分辨率（resolution）参数越大，社区越小越多。结果是算法标签，不是生物学分类。

## CS 类比

≈ 图上的社区检测 / 无监督分区，resolution 相当于粒度参数。

## 常见误读与注意

- 换分辨率或种子群会变；C5 拆分为 C5-1/C5-2 要靠多个随机种子检验其稳定性。

## 在本项目中

- 三个分辨率：0.2→6 群，0.5→8 群，1.0→16 群；主分析用 0.5（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.3）。
- C5 共 1,129 细胞拆为 C5-1=812、C5-2=317，五种子平均 ARI=0.9868（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 5）。
- 单组织稳定分群：骨髓 5 群、脾脏 4 群、胸腺 6 群（同上 阶段 7）。

## 相关概念

- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
- [[近邻图（Neighbour Graph, k-NN）|Neighbour-Graph]] — 在 PCA 空间给每个细胞找 15 个最相似的邻居，UMAP 和聚类都建立在这张图上。
- [[调整兰德指数（Adjusted Rand Index, ARI）|ARI]] — ARI 衡量两种分组的一致程度：1 = 完全相同，约 0 = 随机水平。
- [[细胞类型与细胞状态（Cell Type vs Cell State）|Cell-Type-vs-Cell-State]] — “类型”是相对稳定的身份，“状态”是同一类细胞当前的活动模式；本项目的 C0…C7 是算法分出的表达状态簇，不是已验证的细胞类型。
