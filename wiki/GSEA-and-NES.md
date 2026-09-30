# 基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）

> GSEA 不需要先选“显著基因”，而是把全部基因按变化排序，看基因集是否集中在顶端或底端；NES 带方向。

## 定义

GSEA（gene set enrichment analysis）对全部基因的排序列表（如 T2 vs Ctrl 的统计量）做走动求和，检验某基因集成员是否集中在列表两端；通过置换得到显著性。归一化富集分数（normalized enrichment score, NES）是按基因集大小归一化后的富集分，正 = 靠近 T2 上调端，负 = 靠近下调端。

## CS 类比

≈ 不设阈值的“排序质量”检验：看某个模块的成员是否倾向于排在榜单前列。

## 常见误读与注意

- ORA 与 GSEA 回答不同问题：ORA 只用过线基因，GSEA 用全排序。骨髓 C0 的剪接子集 ORA 有支持，但 GSEA 不显著。
- 基因集越大越容易被核糖体/热休克这类成员多的基因“撑起来”。

## 在本项目中

- 8/30 起 Figure E 扩展为 9 个 refined cluster 的 signed GSEA，用带方向 NES；后来要求收窄到 tissue×cluster（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 5）。
- 基因集大小 15–500，默认 5000 次置换（`notebooks/scripts/iNKT/run_inkt_c5_paper_followup.py` 常量 `GSEA_MIN_SIZE`、`GSEA_MAX_SIZE`、`GSEA_DEFAULT_PERMUTATIONS`）。
- 37 页版第 31 页 “signed GSEA and heat-shock sensitivity”：红/蓝表示带符号 NES，符号使用 GSEA q（`iNKT_by_date/2026-09-20/package/notes/presenter_guide.zh.md`）。
- TNF、NF-κB、JAK–STAT、IL-17/Th17 在组织×簇 GSEA 中未出现 q≤0.05（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 7）；第 31 页已有数值的格子仍无 q≤0.05（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。

## 相关概念

- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[图中的富集分数：−ln p 与 −log10 q（Enrichment Scores in Figures）|Enrichment-Scores-in-Figures]] — 不同图用了不同的“分数”：旧图 −ln(名义 p)，新 GO 气泡图 −log10(BH q)，GSEA 图用带符号 NES。
- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
