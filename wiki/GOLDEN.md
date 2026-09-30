# GOLDEN 与融合聚类（GOLDEN, Fusion）

> GOLDEN 是 Yue 团队的“功能条目语义 + 关系图”聚类方法；本项目做的是小鼠 GO 本地适配版。

## 定义

GOLDEN（及 GOLDEN-GNN）是 ai-pharm-AU 团队开源的方法：把术语的文字描述嵌入（sentence embedding）与术语关系邻接矩阵拼接后聚类（fusion），可再接 GraphSAGE 学到的图表示。

## CS 类比

≈ 多模态 embedding 拼接：文本向量 + 图结构向量，再做聚类。

## 常见误读与注意

- 官方输入中的 m-type/PAGER 关系未随代码提供，本项目换成冻结小鼠 GO 的 Jaccard 网络，所以不能称完全复现（`iNKT_by_date/2026-09-25/results/fusion_manifest.json` “not_claimed”）。
- α 或 k 改变，分组会变；最佳 silhouette 只给出 2 个粗组，不能当作稳定的精细生物模块。

## 在本项目中

- 融合公式：concat(α·L2(description embedding), (1−α)·raw weighted adjacency)，α=0.5，Ward 聚类，k=2–12 取最大 silhouette → k=2（F01=73 个节点、F02=24 个节点；`iNKT_by_date/2026-09-25/results/fusion_manifest.json`；`iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`）。
- α 稳定性：`iNKT_by_date/2026-09-25/results/network/fusion_alpha_stability.csv`（α=0 时与 α=0.5 的 ARI≈0.755）。
- 上游代码：https://github.com/ai-pharm-AU/GOLDEN 与 https://github.com/ai-pharm-AU/GOLDEN-GNN（`iNKT_by_date/2026-09-25/results/README.md` 第 5 节）。
- 9/29 前 10 分钟：ZERU 向 Yue 汇报论文稿修改（notation、edge drop、最终模型配置未写清等，推测为 GOLDEN/知识图谱稿件）（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:09:58–16:12:16）；这是稿件修改，不是 iNKT 结果。

## 相关概念

- [[文本嵌入与 MiniLM（Sentence Embedding, all-MiniLM-L6-v2）|Text-Embedding-MiniLM]] — 把 GO 术语的文字定义变成 384 维向量，语义相近的术语向量也相近。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] — PAG = 通路/注释基因列表/基因签名；PAGER 定义 PAG 之间的关系；Yue 用 super-PAG 指“把许多 PAG 聚成的功能大组”。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
