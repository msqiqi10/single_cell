# 多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）

> 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。

## 定义

多重检验（multiple testing）问题：检验的次数越多，偶然出现小 p 的次数越多。Benjamini-Hochberg（BH）程序把 p 值转成错误发现率（false discovery rate, FDR）——即被判显著的结果中期望有多少比例是假阳性——调整后的值称 q 值。

## CS 类比

≈ 跑一万个 A/B 测试时，把“不显著”的阈值按测试数量放宽/收紧。

## 常见误读与注意

- FDR 依赖“校正范围”（family）：不同范围得到不同 q 值，不能随意比较。

## 在本项目中

- 基因层面：每个 DE 比较内对 10,670 个基因做 BH（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 8.2）。
- GO 富集：每个 对照×方向×GO 分支 在全部可检验条目内 BH（含零命中条目），另给 `q_analysis_global`（跨全部对照）（`iNKT_by_date/2026-09-19/README.md` GO 方法）。
- AP-1 骨髓 C0：家族 q=0.000261，全局 q=0.00286；骨髓 C4：0.0351 / 0.159（家族内通过，全局不通过）（`iNKT_by_date/2026-09-25/results/tables/AP1_complex_GO_exact_ID.csv`）。
- ORA 表列名 `q_family`/`q_analysis_global` 即这两种范围。

## 相关概念

- [[p 值（p-value）|P-value]] — p 值是“如果两组其实没差别，看到这么大差异的概率”，不是效应大小，也不是结论正确的概率。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）|Thresholds-Legacy-vs-Robust]] — legacy = 名义 p≤0.05 且 FC≥1.5；robust = FDR≤0.05 且 |log2FC|≥0.25。两套同时改变校正和效应门槛，不能称“更严/更松”。
- [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] — DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。
