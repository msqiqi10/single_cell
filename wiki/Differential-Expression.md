# 差异表达（Differential Expression, DE; Wilcoxon）

> DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。

## 定义

差异表达（differential expression, DE）对每个基因问：两组细胞里它的表达有没有系统性差别。本项目用 Wilcoxon 秩和检验（scanpy `rank_genes_groups`，`tie_correct=True`），逐基因给出 log2FC、p 值，再对 10,670 个基因做 BH 校正。

## CS 类比

≈ 对 10,670 个指标各做一次 A/B 测试，再统一做多重比较修正。

## 常见误读与注意

- 细胞被当成独立观测，而每个 tissue×condition 只有一个 pool，所以 p 值是探索性的（见 [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]]）。
- “0 个 DEG” ≠ “被跳过”：没做检验（<20 细胞）与检验了但没基因过线是两回事。

## 在本项目中

- 全局：7,359 T2 vs 8,173 Ctrl，仅 3,000 HVG，345 个基因 FDR≤0.05，132 个同时 |logFC|≥0.25（57 上调、75 下调）（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 8.1）。
- 组织内：骨髓 718（350↑/368↓）、脾脏 864（356/508）、胸腺 352（195/157）；主要簇：C0 骨髓 592（274↑）、C3 脾脏 735（275↑）、C6 胸腺 304（140↑）、C4 骨髓 83（56↑）（同上 8.3–8.4）。
- 每组至少 20 个细胞，否则跳过；规划 26 个单元，19 完成、7 跳过（同上 8.2；这是 8 个簇版本的记录）。
- 09-25 保留 26 份既有完整 DE（277,420 行），主功能分析使用 17 个组织×簇比较和 3 个组织整体比较：`iNKT_by_date/2026-09-25/results/tables/all_26_existing_DE_with_eligibility.csv.gz`、`iNKT_by_date/2026-09-25/results/tables/all_comparison_status.csv`。
- 37 页版第 11–21 页展示 T2 vs Ctrl DEG（`iNKT_by_date/2026-09-20/package/slides.json`）。

## 相关概念

- [[对数倍数变化（Log2 Fold Change, log2FC）|Log2-Fold-Change]] — log2FC 衡量效应大小：+1 表示 T2 约为 Ctrl 的 2 倍。
- [[p 值（p-value）|P-value]] — p 值是“如果两组其实没差别，看到这么大差异的概率”，不是效应大小，也不是结论正确的概率。
- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
- [[两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）|Thresholds-Legacy-vs-Robust]] — legacy = 名义 p≤0.05 且 FC≥1.5；robust = FDR≤0.05 且 |log2FC|≥0.25。两套同时改变校正和效应门槛，不能称“更严/更松”。
- [[火山图（Volcano Plot）|Volcano-Plot]] — 火山图以 log2FC 为横轴、−log10 p 为纵轴，同时看效应大小和显著性。
