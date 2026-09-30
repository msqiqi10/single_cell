# 基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）

> 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。

## 定义

共同邻居（common neighbours）给一对节点打分 = 它们共有的邻居个数；Adamic-Adar 对共有邻居按 1/log(度) 加权，稀有邻居更重要。这两个是链接预测里的经典无监督基线，深度模型必须明显超过它们才有价值。

## CS 类比

≈ 不用机器学习的 baseline heuristics：先看简单规则能做到多少。

## 常见误读与注意

- 评价新模型必须与这类基线比较；只报告自己 AUC 高是不够的。

## 在本项目中

- 共同邻居 AUC 0.9595–0.9685，Adamic-Adar 相近（0.9599–0.9683）（`iNKT_by_date/2026-09-25/results/network_validation.json`）。
- 09-30 最终（3 种子均值）：GNN 0.973（a）/0.981（b），共同邻居 0.922 / 0.922，Adamic-Adar 0.928 / 0.922；GNN 仅高约 0.05，测试边只有 23 / 15 对（`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_all.csv`、`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_excl.csv`）。

## 相关概念

- [[链接预测评估（Link Prediction: held-out edges, negative sampling）|Link-Prediction-Evaluation]] — 把一部分真实边藏起来，再看模型能否把它们排在随机“非边”之前。
- [[曲线下面积（Area Under the ROC Curve, AUC）|AUC]] — AUC 衡量分类器把正样本排在负样本前的概率：0.5 = 随机，1 = 完美。
- [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] — 在很小、由共享基因直接定义边的图上，简单基线已经很强，GNN 只多出 0.02–0.03 AUC，且测试集只有 57 对——这是解读，不是文件里的结论。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
