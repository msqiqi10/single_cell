# 为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）

> 在很小、由共享基因直接定义边的图上，简单基线已经很强，GNN 只多出 0.02–0.03 AUC，且测试集只有几十对——这是解读，不是文件里的结论。

## 定义

这一页是基于项目数据的解读，标注为“推断”。事实：GNN AUC 0.982–0.993，共同邻居 0.960–0.968；图有 97 节点、577 条边；测试正负各 57；验证与测试都在同一张内部图内。推断：因为边本身由“共享基因 Jaccard≥0.25”定义，共享邻居信息已经隐含在图结构里，GNN 主要是重新编码它；在如此小的测试集上 0.02–0.03 的差可能在噪声范围内（推断，未做显著性检验）。

## CS 类比

≈ 在一个小到能装进内存的数据集上，用 Transformer 去对比一个哈希表：结果可能略好，但工程价值有限。

## 常见误读与注意

- 不要把 AUC 0.99 写成“GNN 优于传统方法”。

## 在本项目中

- 事实来源：`iNKT_by_date/2026-09-25/results/network_validation.json`，`scope` 字段：内部 transductive 重建，不是生物学验证或方法基准。
- 什么能让 GNN 增益更有说服力（推断，需要 Yue 确认）：更大的图（混合 GO+KEGG，见 [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]]）；能提供“共享基因之外”信息的特征（表达效应、语义、知识库关系）；外部/未见数据上的验证；比较分组与生物学解释的增益（例如 GNN 是否把 AP-1/MAPK/免疫功能放进了比重叠聚类更有意义的组）；多种子稳定性（ARI 0.80–0.91）。
- Yue 9/29：GNN 不是做一次就好，要反复推敲，并加入验证信息（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:31:40–16:31:55）。
- **09-30 最终数字**（`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_all.csv`、`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_excl.csv`，3 种子均值）：留出边 AUC GNN 0.973（版本 a）/ 0.981（版本 b 主版本），共同邻居 0.922 / 0.922，Adamic-Adar 0.928 / 0.922，语义余弦 0.840 / 0.723。GNN 比图基线高约 0.05，但测试集只有 23 / 15 对正边（各配等量负边），单个种子的 AUC 在 0.94–1.00 之间波动。
- 稳定性：GNN+语义分组的跨种子 ARI 为 (a) 0.61–0.71、(b) 0.32–0.55（`iNKT_by_date/2026-09-30/results/tables/ARI_across_gnn_seeds_all.csv`、`iNKT_by_date/2026-09-30/results/tables/ARI_across_gnn_seeds_excl.csv`）；(b) 中 GNN 与 GOLDEN 融合 ARI 0.70，分组层面相对 GOLDEN 融合没有明显额外增益。
- 推断：图更大（混合 GO/KEGG）后，GNN 仍只比共同邻居好一点，而且分组不稳定——支持“增益有限”的判断，但没做显著性检验。

## 相关概念

- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] — 不用任何学习，仅按“共享多少邻居”就能预测两节点是否相连。
- [[链接预测评估（Link Prediction: held-out edges, negative sampling）|Link-Prediction-Evaluation]] — 把一部分真实边藏起来，再看模型能否把它们排在随机“非边”之前。
- [[曲线下面积（Area Under the ROC Curve, AUC）|AUC]] — AUC 衡量分类器把正样本排在负样本前的概率：0.5 = 随机，1 = 完美。
- [[调整兰德指数（Adjusted Rand Index, ARI）|ARI]] — ARI 衡量两种分组的一致程度：1 = 完全相同，约 0 = 随机水平。
