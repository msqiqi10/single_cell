# 基因层面钻取（Gene-level Drill-down）

> 把网络里的“通路节点”展开成基因：每个通路有哪些基因、哪些被检出且差异表达、每个基因被注释了多少次、哪些是共享基因、哪些是特异基因。

## 定义

Yue 在 9/29 要求：每个 PAG 里有多少基因、哪些与 DE 重叠、每个基因被注释多少次（共享 vs 特异）。09-30 的 drill-down 复用已有 ORA 和网络（未重算 ORA/GNN，无新数据，仅 CPU），钻取 7 个对象：KEGG MAPK、IL-17、TCR、Th17、TNF，GO:0035976 AP-1 复合体，以及版本 (b) GNN+语义 G01（氧化磷酸化）作为“能量”对照。

- **成员**：冻结注释的全部基因（KEGG 用冻结 gmt；GO 用冻结 GAF 加 is_a/part_of 传播），包括未检测到的基因。能量模块的成员 = 在 G01 的 15 个词条中出现≥3 次的基因，共 181 个（并集 773 个含疾病类词条，过宽）。
- **“检出”** = 出现在 DE 表（10,670 个基因）。DEG 判定用 6 个关键比较：BM C0、Spleen C3、BM C4、BM C5-2、Spleen C5-2、Thymus C6。
- **注释频次**：一个基因落在多少个显著节点里；**shared** = 落在≥2 个被钻取对象里，**specific** = 只在 1 个里。成员表共 885 行、681 个不同基因，其中 560 行 specific、325 行 shared。
- **共成员边**：两个基因共同出现在≥2 个显著词条/锚点里就连一条边，权重 = 共享词条数（18,416 条，`iNKT_by_date/2026-09-30/drilldown/tables/gene_comembership_edges.csv`）。**这不是 PPI**（没用 STRING），不表示物理相互作用。

## CS 类比

≈ 从“包级依赖图”下钻到“函数级调用表”：先看每个包里有哪些函数、被几个包引用（shared/specific），再问哪些函数在这次运行里真的被执行（DEG）。

## 常见误读与注意

- 大通路里多数成员基因没有差异表达（例如 KEGG MAPK 294 个成员，166 个被检出，任一关键比较里为 DEG 的只有 23 个）。
- 共享基因多不等于功能相关；见 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。
- 未检出的基因不是“没表达”，而是不在 10,670 基因的 DE 表里（受过滤阈值影响，见 [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]]）。
- 细胞层面探索性统计，每个组织×条件是 3 只小鼠 pooled（n=1 pool）；注释不等于通路激活。

## 在本项目中

- `iNKT_by_date/2026-09-30/drilldown/tables/pathway_gene_membership.csv`：每基因是否检出、6 个比较的 log2FC/FDR/DEG、注释频次、shared/specific。
- `iNKT_by_date/2026-09-30/drilldown/tables/pathway_drill_summary.csv`：每个对象的成员数/检出数/DEG 数（KEGG MAPK 294/166/23；IL-17 91/54/13；TCR 101/92/11；Th17 102/81/6；TNF 110/79/14；AP-1 6/6/5；能量 181/154/39，其中 DEG 数为任一关键比较）。
- `iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv` 与 `iNKT_by_date/2026-09-30/drilldown/tables/bridge_candidate_genes.csv`。
- 静态图：`iNKT_by_date/2026-09-30/drilldown/figures/gene_network_*.png`（每个通路一幅，BM C0 与 Spleen C3 双面板，DEG 红框加星，未检出灰色空心）。
- 交互版：[[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] 的“基因网络”标签页。
- 结果：AP-1 与 KEGG 的桥只有 Fos+Jun 加一个；能量模块与 AP-1 无共享基因。

## 相关概念

- [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] — 钻取的核心结论。
- [[锚点节点（Anchor Nodes）|Anchor-Nodes]] — 钻取对象大多是锚点。
- [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] — 可点击版本。
- [[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] — Yue 所说的 PAG 及其基因。
- [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] — 与“检出/被覆盖”相关的另一个问题。
- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 术语层面的网络；这里下钻到基因。
