# 锚点节点（Anchor Nodes）

> 锚点 = 我们出于生物学问题手动固定加入网络的节点；它们即使不显著也保留，因为“不显著”本身就是答案的一部分。

## 定义

在 09-30 混合 GO/KEGG 网络里，常规节点只来自 ORA 显著的术语（q≤0.05 且命中≥3）。**锚点**是另外固定加入的 12 个术语，对应 Rob 与 Yue 关心的问题：MAPK（GO MAPK cascade、ERK1/ERK2 cascade、p38MAPK cascade、JNK cascade、KEGG MAPK signaling）、AP-1（GO:0035976 转录因子 AP-1 复合体，CC）、IL-17/Th17（KEGG IL-17 signaling、KEGG Th17 cell differentiation、GO T-helper 17 type immune response、GO interleukin-17 production）、TCR（KEGG T cell receptor signaling）、TNF（KEGG TNF signaling）（`iNKT_by_date/2026-09-30/code/common.py`；`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv`）。合并冗余节点时锚点受保护，不会被并入别人。

## CS 类比

≈ 单元测试里的固定用例（fixture）：无论覆盖率工具认为它“有没有被触发”，都要把它放进测试集里，才能回答“这条路径到底通不通”。

## 常见误读与注意

- **保留不等于显著。** 12 个锚点中只有 GO AP-1 复合体显著（q≈1.7e-4，Fos/Fosl2/Jun/Junb/Jund；骨髓 C0、脾 C3）。KEGG IL-17 q=0.071，MAPK 0.64，TNF 0.56，其余 1（`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv`）。
- “interleukin-17 production” 在冻结注释中没有任何实测基因（`n_measured_genes=0`），只是孤立的空节点。
- 如果只画显著节点，就看不到“MAPK、IL-17 没有出现”这件事，读者容易误以为它们被漏掉了。保留锚点，正是为了把阴性结果摆在图上（对应 Rob 关心的 MAPK/AP-1/IL-17，见 [[人物与各自关心的点（Who Is Who）|Who-is-who]]）。
- 锚点被分到同一组不代表相关：版本 (b) 中 GOLDEN 与 GNN+语义把它们放进同一个杂项组 G03，主要因为它们大多是孤立点或碎片；分组结果随方法与 k 变化（`iNKT_by_date/2026-09-30/results/tables/anchor_k_sensitivity_excl.csv`）。

## 在本项目中

- 锚点 vs 显著核心：显著核心是能量代谢（氧化磷酸化/ATP/呼吸链），剔除核糖体/线粒体/Hsp/Dnaj 后仍在；锚点与它几乎没有共享基因边。
- 基因层面：锚点两两共享的基因见 `iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv`，结论汇总在 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。
- 网络构建位置：[[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]]。

## 相关概念

- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 锚点所在的网络与分组结果。
- [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] — 锚点之间在基因层面靠什么相连。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — 决定节点是否显著的检验。
- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — q 值的含义。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — 不显著 ≠ 不激活，显著也 ≠ 激活。
- [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] 与 [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] — 两个最重要的锚点家族。
