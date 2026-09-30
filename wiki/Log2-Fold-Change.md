# 对数倍数变化（Log2 Fold Change, log2FC）

> log2FC 衡量效应大小：+1 表示 T2 约为 Ctrl 的 2 倍。

## 定义

倍数变化（fold change, FC）是两组平均表达之比；取以 2 为底的对数得 log2FC。0 表示无变化，+1 表示 2 倍，−1 表示一半。log2(1.5)≈0.585，所以“FC≥1.5”等价于“log2FC≥0.585”。

## CS 类比

≈ 相对性能提升比（speedup ratio）的对数。

## 常见误读与注意

- 稀疏基因（阳性细胞很少）的 fold change 很容易被极少数细胞放大，需要与阳性比例一起读。
- 不同图里的“差值”含义不同：模块分数的 T2−Ctrl 差是 log1p 表达均值之差，不是 log2FC。

## 在本项目中

- DE 里正 logFC = T2 上调（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 8.2）。
- 描述性 marker log2FC = log2((T2 均值+0.001)/(Ctrl 均值+0.001))，不是 DE 模型估计（同上 Stage 7）。
- Il1r1 的 fold change 很大但检出比例低，不能替代不确定性检查（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。

## 相关概念

- [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] — DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。
- [[两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）|Thresholds-Legacy-vs-Robust]] — legacy = 名义 p≤0.05 且 FC≥1.5；robust = FDR≤0.05 且 |log2FC|≥0.25。两套同时改变校正和效应门槛，不能称“更严/更松”。
- [[火山图（Volcano Plot）|Volcano-Plot]] — 火山图以 log2FC 为横轴、−log10 p 为纵轴，同时看效应大小和显著性。
- [[归一化与 log1p（Normalisation, log1p）|Normalisation-and-log1p]] — 把每个细胞缩放到 1 万计数再取 log(1+x)，消除测序深度差异。
