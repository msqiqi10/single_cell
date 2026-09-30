# 火山图（Volcano Plot）

> 火山图以 log2FC 为横轴、−log10 p 为纵轴，同时看效应大小和显著性。

## 定义

火山图（volcano plot）每个点是一个基因：横轴 log2FC（效应大小），纵轴 −log10(p)（显著性）。右上和左上角是既显著又变化大的基因。

## CS 类比

≈ 散点图：x = 影响幅度，y = 置信度。

## 常见误读与注意

- 细胞多时纵轴会被拉得很高；不能只挑最显著的点当结论。
- 图上的点是否标名有筛选规则，“没标”不等于“没变化”。

## 在本项目中

- 37 页版第 22 页 “D Differential expression: Control vs. Tumor”：四个面板（C0/BM、C6/Thy、C5-1/Spl、C5-2/BM），沿用 Blood Fig. 3D 的轴与蓝/红方向；着色规则 基因 BH FDR≤0.05 且 |log2FC|≥0.25，竖线阈值 ±0.25；公共轴 |log2FC|≤5，−log10(P) 截顶到 100，未截顶的值另有导出（`iNKT_by_date/2026-09-20/package/slides.json`；`iNKT_by_date/2026-09-20/package/notes/presenter_guide.zh.md` 第 10 节）。
- 9/1 会议讨论了火山图轴与点的含义（同上 阶段 6）。

## 相关概念

- [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] — DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。
- [[对数倍数变化（Log2 Fold Change, log2FC）|Log2-Fold-Change]] — log2FC 衡量效应大小：+1 表示 T2 约为 Ctrl 的 2 倍。
- [[p 值（p-value）|P-value]] — p 值是“如果两组其实没差别，看到这么大差异的概率”，不是效应大小，也不是结论正确的概率。
