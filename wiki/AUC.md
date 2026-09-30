# 曲线下面积（Area Under the ROC Curve, AUC）

> AUC 衡量分类器把正样本排在负样本前的概率：0.5 = 随机，1 = 完美。

## 定义

ROC 曲线以不同阈值下的假阳性率为横轴、真阳性率为纵轴；曲线下面积（AUC）等价于“随机取一对正、负样本，正样本得分更高的概率”。AP（average precision）是类似的排序指标。

## CS 类比

≈ 排序质量指标：搜索结果里相关文档排在不相关文档之前的概率。

## 常见误读与注意

- 测试集很小（这里 57 对正负边）时，AUC 差 0.02–0.03 可能在噪声之内。
- 图由共享基因构成，测试的是“重建自身”，不能当作新数据上的优势。

## 在本项目中

- 09-25 内部留边测试：GNN AUC 0.982–0.993；语义余弦基线 0.894–0.916；共同邻居 0.960–0.968；每次 57 条测试正边与 57 条负边（`iNKT_by_date/2026-09-25/results/network_validation.json`；`iNKT_by_date/2026-09-25/results/network/gnn_baseline_comparison.csv`）。
- Borra et al. 2026 的基因 panel 也报告 accuracy 与 AUC（`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`）。

## 相关概念

- [[链接预测评估（Link Prediction: held-out edges, negative sampling）|Link-Prediction-Evaluation]] — 把一部分真实边藏起来，再看模型能否把它们排在随机“非边”之前。
- [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] — 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。
- [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] — 在很小、由共享基因直接定义边的图上，简单基线已经很强，GNN 只多出 0.02–0.03 AUC，且测试集只有 57 对——这是解读，不是文件里的结论。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
