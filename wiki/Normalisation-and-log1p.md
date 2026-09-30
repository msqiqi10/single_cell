# 归一化与 log1p（Normalisation, log1p）

> 把每个细胞缩放到 1 万计数再取 log(1+x)，消除测序深度差异。

## 定义

归一化（normalisation）把每个细胞的总计数缩放到相同值（这里 `target_sum=1e4`），再做 `log1p`（log(1+x)）压缩高表达基因的动态范围。

## CS 类比

≈ 把每条请求按总字节数归一化后再取对数，防止大请求淹没小请求。

## 常见误读与注意

- 归一化不会消除 batch effect，也不会使不同样本成为生物学重复，也不能消除组织构成差异。

## 在本项目中

- `sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)`（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 3）。
- 细胞毒性模块的“T2−Ctrl 差”是 log1p 表达均值之差，不是 log2FC，也不是杀伤率（`iNKT_by_date/2026-09-19/README.md`）。

## 相关概念

- [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]] — 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。
- [[高变基因（Highly Variable Genes, HVG）|Highly-Variable-Genes]] — 选 3,000 个在细胞间波动最大的基因给 PCA 用，降低噪声和计算量。
- [[对数倍数变化（Log2 Fold Change, log2FC）|Log2-Fold-Change]] — log2FC 衡量效应大小：+1 表示 T2 约为 Ctrl 的 2 倍。
- [[标志基因与模块分数（Marker Genes, Module Scores）|Marker-Genes-and-Module-Scores]] — marker 用来给簇命名，模块分数把一组基因合成一个数。
