# 计数矩阵（Count Matrix, AnnData）

> 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。

## 定义

计数矩阵（count matrix）的每行是一个细胞、每列是一个基因、值是 UMI 计数。Scanpy 把它存为 AnnData（.h5ad）：`X` 主矩阵，`layers` 存其它版本，`obs` 存细胞元数据，`var` 存基因元数据，`raw` 存备份。

## CS 类比

一个 15,532 × 10,670 的稀疏 DataFrame + 两张元数据表；`layers` ≈ 同一张表的多个视图。

## 常见误读与注意

- Scanpy 的 `adata.raw` 不一定是原始计数：本项目的 `.raw` 是 log-normalised 表达，真正的整数在 `layers["counts"]`。

## 在本项目中

- `layers["counts"]` = QC 后整数计数；`X` = normalize+log1p；`raw` = 全部 10,670 基因的 log-normalised（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 3）。
- 09-25 又从六份原始矩阵重建同一批 15,532 个细胞的 32,285 基因计数对象，用于补查 Ccr6/Il17a/Il17f（`iNKT_by_date/2026-09-25/results/README.md` 第 4 节）。
- 细胞元数据（sample、condition、tissue、leiden、C5 亚群、各类分数）：`iNKT_by_date/2026-09-19/results/tables/cell_metadata_scores.csv.gz`。

## 相关概念

- [[单细胞 RNA 测序与 10x 矩阵（scRNA-seq, 10x Matrices: barcodes, features）|scRNA-seq-and-10x-Matrices]] — scRNA-seq 给每个细胞打一个条形码，数每个基因的 mRNA 分子数；10x 输出三个文件：matrix、barcodes、features。
- [[归一化与 log1p（Normalisation, log1p）|Normalisation-and-log1p]] — 把每个细胞缩放到 1 万计数再取 log(1+x)，消除测序深度差异。
- [[质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）|QC-Metrics]] — QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。
- [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] — 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。
