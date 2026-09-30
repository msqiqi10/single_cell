# 功能组 G01–G07 及含义（Functional Groups G01–G07）

> 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。

## 定义

功能组（G01–G07）是对 97 个 GO BP 节点的探索性聚类（融合 MiniLM 语义与 GraphSAGE 表示，k=7）。它们是“功能条目的分组”，不是细胞亚型。下表的含义取自 09-25 报告的解释，节点数按 `annotated_nodes.csv` 统计。

## CS 类比

≈ 对模块做的目录分类，不是新的运行时对象。

## 常见误读与注意

- 组名是我们为讨论起的“便于讨论的解释”，边界为 09-25 README 明列。
- 名称如“肌肉/突触”不能当作细胞身份或新功能。

## 在本项目中


| 组 | 节点数 | 便于讨论的解释 | 边界（源：`iNKT_by_date/2026-09-25/results/README.md` 第 5 节） |
|---|---:|---|---|
| G01 | 26 | 混合注释：热应激、蛋白折叠、多种共享基因条目（名称含 spermatogenesis 等） | 异质性高，不能当作一个确定功能 |
| G02 | 9 | 翻译、核糖体生成、rRNA 相关 | 不等于直接测量蛋白合成速率 |
| G03 | 10 | 共享核糖体成员带出的条目（translation at synapse 等） | “肌肉/突触”名称不是新功能 |
| G04 | 15 | 呼吸链、氧化磷酸化、ATP 生成 | 表达富集不等于代谢通量 |
| G05 | 9 | 呼吸链复合体组装 | 与 G04 有成员重叠 |
| G06 | 24 | 核苷酸代谢/合成，与能量相关条目重叠 | 与 G04 的 155 条条目边反映共享注释，不是上下游关系 |
| G07 | 4 | 离子/质子跨膜运输 | driver 包括线粒体/ATP 相关基因 |

- 节点数来自 `iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`（`gnn_module` 列）。
- 骨髓 C0 的 52 个显著 BP 条目落在 G03/G04/G05/G06，分别涉及 12/36/18/34 个 driver（重叠，合并仍只有 66 个不同基因）（`iNKT_by_date/2026-09-25/results/README.md` 第 5 节）。
- 9/29：Yue 认为 G04–G06 至少给出“不只是 G01 应激”的额外信息，而核糖体/应激相关组是已知的、不新颖；他希望进一步用 GO MF 和 KEGG 再做（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:39:06–16:40:56；ASR 文本，G 编号为推测）。
- 图：`iNKT_by_date/2026-09-25/results/figures/02_function_groups_by_cluster.png`；追溯链：`iNKT_by_date/2026-09-25/results/tables/module_gene_evidence.csv.gz`。
- 注意：`annotated_nodes.csv` 里还有一套 GOLDEN 语义+邻接融合分组（F01/F02，k=2），与这里的 G01–G07（GNN 融合）是不同方法的输出。

## 相关概念

- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量；GraphSAGE 是常用实现。
- [[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] — 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。
- [[线粒体、氧化磷酸化与 ATP（Mitochondria, Oxidative Phosphorylation, ATP: OXPHOS）|Mitochondria-and-OXPHOS]] — 线粒体通过氧化磷酸化（OXPHOS）造 ATP；线粒体基因占比也是 QC 指标。
- [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] — 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。
