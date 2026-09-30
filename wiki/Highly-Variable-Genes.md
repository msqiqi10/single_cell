# 高变基因（Highly Variable Genes, HVG）

> 选 3,000 个在细胞间波动最大的基因给 PCA 用，降低噪声和计算量。

## 定义

高变基因（highly variable genes, HVG）指细胞之间表达波动明显超出随机预期的基因，是降维前的特征选择。

## CS 类比

≈ 特征选择（feature selection）：只保留方差大的列再做 PCA。

## 常见误读与注意

- `batch_key="sample"` 只让 HVG 选择对样本敏感，并不是批次校正。
- 初始 marker 排名只在 3,000 个 HVG 里检验，所以 Il17a 等不在 HVG 的基因“没被测试”而不是“没差异”。

## 在本项目中

- `n_top_genes=3000, batch_key="sample", flavor="seurat"`；355 个基因在六个样本中都是 HVG（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 4.1）。
- global 组织合并 DE 只测了 3,000 个 HVG；后来的 tissue 与 cluster×tissue DE 用全部 10,670 基因（同上 Stage 8）——两者“基因宇宙”不同，图不可直接比较。

## 相关概念

- [[主成分分析（Principal Component Analysis, PCA）|PCA]] — PCA 把 3,000 维表达压缩成 50 个主轴，供后续邻居图使用。
- [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] — DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。
- [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] — 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。
