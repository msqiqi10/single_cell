# 链接预测评估（Link Prediction: held-out edges, negative sampling）

> 把一部分真实边藏起来，再看模型能否把它们排在随机“非边”之前。

## 定义

链接预测（link prediction）把边随机分成训练/验证/测试；测试正样本是被藏起来的真实边，负样本是随机抽取的非边（negative sampling）。用 AUC、AP 衡量排序质量。训练时的负样本必须排除所有真实边（包括被留出的），否则会泄漏。

## CS 类比

≈ 训练/验证/测试拆分 + 数据泄漏检查；测试集只能在选定 checkpoint 后碰一次。

## 常见误读与注意

- 这是“内部 transductive”重建：图本身由当前 GO 共享成员构成，只能检验图重建，不是新数据上的方法优越性，也不是生物学验证。

## 在本项目中

- 每次 57 条测试正边 + 57 条负边；GNN AUC 0.982–0.993，语义余弦 0.894–0.916，训练图上的共同邻居 0.960–0.968（`iNKT_by_date/2026-09-25/results/network_validation.json`；`iNKT_by_date/2026-09-25/results/network/gnn_baseline_comparison.csv`）。
- `external_validation: false`、“not biological validation or a method benchmark”（同 network_validation.json）。
- 09-30：测试边只有 23+23（全部）/15+15（排除版），AUC 波动更大（`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_all.csv`）。

## 相关概念

- [[曲线下面积（Area Under the ROC Curve, AUC）|AUC]] — AUC 衡量分类器把正样本排在负样本前的概率：0.5 = 随机，1 = 完美。
- [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] — 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] — 在很小、由共享基因直接定义边的图上，简单基线已经很强，GNN 只多出 0.02–0.03 AUC，且测试集只有 57 对——这是解读，不是文件里的结论。
