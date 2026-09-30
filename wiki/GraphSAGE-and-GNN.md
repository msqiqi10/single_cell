# 图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）

> GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。

## 定义

图神经网络（graph neural network, GNN）在图上逐层做“消息传递”：每个节点把邻居的向量聚合后与自身合并。GraphSAGE 通过采样并聚合邻居来学习节点表示，这里以 384 维 MiniLM 向量为节点输入特征。

## CS 类比

≈ 图上的卷积/PageRank 式迭代：每一轮每个节点读取邻居的状态更新自己。

## 常见误读与注意

- GNN 的表示由“图 + 节点特征”决定；在一张由共享基因构成的小图上，它很可能只是把已经在图里的信息重新编码。

## 在本项目中

- 配置：官方 GraphSageEncoder 不变；256 隐藏维、2 层、dropout 0.1、Adam 1e-3、120 epochs；种子 42/43/44；输入 97 节点×384 维，边权 = Jaccard，训练/验证/测试边 463/57/57（`iNKT_by_date/2026-09-30/provenance/gnn_adapter.json`；`iNKT_by_date/2026-09-25/results/models/gnn_seed42.summary.json`）。
- 语义 + GNN 融合：各行 L2 归一化后等权拼接，参考种子 42 的 k=2–12 → 7 组（G01–G07）（`iNKT_by_date/2026-09-25/results/network_validation.json`）。
- 本项目改动：训练负样本排除全部真实边及留出负样本；验证/测试负样本互斥；测试只在验证集选定 checkpoint 后评价（同 gnn_adapter.json）。
- 9/29：ZERU 已确认论文中的 “edge drop” 与 dropout 作用相同、不参与 loss（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:10:18–16:10:44，转写不完全清晰）。

## 相关概念

- [[链接预测评估（Link Prediction: held-out edges, negative sampling）|Link-Prediction-Evaluation]] — 把一部分真实边藏起来，再看模型能否把它们排在随机“非边”之前。
- [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] — 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。
- [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] — 在很小、由共享基因直接定义边的图上，简单基线已经很强，GNN 只多出 0.02–0.03 AUC，且测试集只有 57 对——这是解读，不是文件里的结论。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] — GOLDEN 是 Yue 团队的“功能条目语义 + 关系图”聚类方法；本项目做的是小鼠 GO 本地适配版。
