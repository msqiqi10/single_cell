# 混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）

> 把 GO BP、GO MF、KEGG 放进同一张功能网络，再用 12 个锚点节点问“AP-1 / MAPK / IL-17 / TCR / TNF 到底有没有功能联系”；09-30 已完成。

## 定义

按 9/29 会议 Yue 的要求，09-30 分析在同一 ORA 规则下加入 GO MF 与 KEGG_2019_Mouse，与 GO BP 一起建**混合术语网络**（[[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]]）。ORA：基因 BH FDR≤0.05 且 |log2FC|≥0.25，T2-up / Ctrl-up 分开，词条大小 5–500，BH 在“比较×方向×来源”内。节点：q≤0.05 且命中≥3，另加 12 个固定[[锚点节点（Anchor Nodes）|Anchor-Nodes]]；Jaccard≥0.75 的冗余节点合并（锚点不会被并入别人）；边：Jaccard≥0.25 且共享≥3 个基因。分组用五种方法：重叠 Louvain、语义（MiniLM）、GOLDEN 融合（α=0.5）、GNN+语义（GraphSAGE，seed 42/43/44）、GNN 单独。

两个版本：
- (a) 全部基因；
- (b) **主版本**：query 与背景同时剔除核糖体/线粒体核糖体/线粒体编码/Hsp/Dnaj 基因（正则 `^(Rp[ls]\d|Rplp\d|Rpsa$|Mrp[ls]\d|mt-|Hsp(?!g)|Dnaj)`，242/10,670 个基因；保留 Rps6k* 激酶与 Hspg2）。

## CS 类比

≈ 已合并进 main 的 feature branch：ORA 是“单元测试”，网络分组是“重构后的目录结构”，锚点是我们手写的几个固定用例——不管测试是否通过都要跑它们。

## 常见误读与注意

- **“所有锚点在同一组”不是功能联系。** 版本 (b) 中 GOLDEN 与 GNN+语义（按轮廓系数 k）把全部 AP-1/MAPK/IL-17/TCR 锚点放进同一个大组 G03（36 节点），但 G03 是孤立点与 2 节点碎片组成的“杂项组”；锚点彼此和与能量核心之间几乎没有共享基因边。基因层面的回答见 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。
- 分组答案依赖方法和 k：GNN+语义 k≥7 时 MAPK 与 AP-1/IL-17/TCR 分开；语义与 Louvain 本来就分散；版本 (a) 中 k=15 时 AP-1、MAPK、IL-17/TCR 三处分离（`iNKT_by_date/2026-09-30/results/tables/anchor_k_sensitivity_excl.csv`、`iNKT_by_date/2026-09-30/results/tables/anchor_k_sensitivity_all.csv`）。
- 显著核心是能量代谢，不是信号通路；只有 GO AP-1 复合体这一个锚点显著。
- 富集≠激活；细胞层面探索性统计，每个组织×条件是 3 只小鼠 pooled（n=1 pool）。

## 在本项目中

**规模**（`iNKT_by_date/2026-09-30/results/network/network_summary.json`）

| | (a) 全部基因 | (b) 主版本 |
|---|---|---|
| 显著术语（合并前） | 147 → 159 节点（含锚点） | 82 → 94 节点 |
| 合并后节点 / 边 / 孤立点 | 112 / 231 / 32 | 66 / 157 / 21 |
| 组数 Louvain / 语义 / GOLDEN / GNN+语义 | 49 / 15 / 2 / 15 | 31 / 15 / 3 / 4 |

**显著核心 = 能量代谢。** 氧化磷酸化、ATP/嘌呤核苷酸合成、呼吸链在 BM C0 与 Spleen C3 的 T2-up 驱动基因最多；剔除核糖体/线粒体/Hsp/Dnaj 后仍然存在（版本 (b) 的 G01 氧化磷酸化 15 个词条、G02 嘌呤核苷酸合成 9 个，`iNKT_by_date/2026-09-30/results/tables/group_summary_excl.csv`）。BM/脾 C5-2 无显著词条。见 [[线粒体、氧化磷酸化与 ATP（Mitochondria, Oxidative Phosphorylation, ATP: OXPHOS）|Mitochondria-and-OXPHOS]]。

**锚点**（`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv`）：只有 GO:0035976 AP-1 复合体显著（q≈1.7e-4）；KEGG IL-17 q=0.071，MAPK 0.64，TNF 0.56，其余 1；“interleukin-17 production”在冻结注释中无实测基因，是孤立空节点。

**方法一致性与 GNN**：ARI(a) Louvain–GNN 0.58、语义–GNN 0.53，GOLDEN–其他约 0.03–0.07；(b) GOLDEN–GNN 0.70，其余 0.10–0.35（`iNKT_by_date/2026-09-30/results/tables/ARI_between_methods_excl.csv`）。GNN 跨种子 ARI：(a) 0.61–0.71，(b) 0.32–0.55。留出边 AUC：GNN 0.973（a）/0.981（b），共同邻居 0.922 / 0.922，语义 0.840 / 0.723，但测试边仅 23 / 15 对。详见 [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]]。

**与 09-19 一致性**：GO 部分 312,550 行 q 值最大差 1.1e-16（`iNKT_by_date/2026-09-30/results/tables/ORA_reuse_check_vs_0919.csv`）。

**图**（`iNKT_by_date/2026-09-30/results/figures/`，PNG + SVG）：`network_full_a_all_genes`、`network_full_b_primary_excl_ribo_mito_hsp`（节点大小 = 显著性，颜色 = 功能组；标签避让版在 `iNKT_by_date/2026-09-30/drilldown/figures/network_full_b_primary_excl_ribo_mito_hsp_v2.png`）、`superPAG_panels_*`（大图 + 小图）、`heatmap_group_driver_genes_*`。

**交互 HTML**：`iNKT_by_date/2026-09-30/drilldown/html/inkt_network_explorer.html`，打开方式与三个视图见 [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]]（需联网加载 CDN）。

**基因层面**：[[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]]；结论页 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。

**幻灯片**：`iNKT_by_date/2026-09-30/presentation/iNKT_pathway_network_drilldown_20260930.pptx`（第 32–43 页；构建代码 `iNKT_by_date/2026-09-30/code/build_deck_0930.js`）。

**偏离规格（README 已记）**：锚点在合并时受保护；Louvain 用 networkx 代替 Leiden；额外加了 GNN-only 与 k 敏感性表；只用 CPU，未联系任何人、未上传或下载数据。

## 相关概念

- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[锚点节点（Anchor Nodes）|Anchor-Nodes]] — 为回答“这些通路是否相关”而固定加入的节点，即使不显著也保留。
- [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] — 基因层面的回答：AP-1 与 KEGG 锚点之间唯一有差异表达支持的桥是 Fos + Jun。
- [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] — 每个通路的成员基因、注释频次、共享 vs 特异基因、共成员边。
- [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] — 09-30 交付的单文件 HTML：点击高亮邻居、跳转到基因层。
- [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] — Yue 团队的“语义 + 关系图”聚类方法；本项目是小鼠 GO 本地适配版。
- [[图神经网络与 GraphSAGE（Graph Neural Network, GraphSAGE）|GraphSAGE-and-GNN]] — GNN 让每个节点通过聚合邻居特征来更新自己的向量。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的通路图集合；条目名常是疾病名。
- [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] — Yue 想要的可点击网络；本页 HTML 与 GOLDEN 论文 KG 是两条线。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
