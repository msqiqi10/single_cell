# 单细胞 RNA 测序与 10x 矩阵（scRNA-seq, 10x Matrices: barcodes, features）

> scRNA-seq 给每个细胞打一个条形码，数每个基因的 mRNA 分子数；10x 输出三个文件：matrix、barcodes、features。

## 定义

单细胞 RNA 测序（single-cell RNA sequencing, scRNA-seq）在液滴中给每个细胞的 mRNA 加上细胞条形码（barcode）和分子标签（UMI），测序后统计每个细胞、每个基因的 UMI 数。10x Genomics 的标准输出包括稀疏矩阵 `matrix.mtx`、细胞条形码 `barcodes.tsv`、基因/特征 `features.tsv`。

## CS 类比

barcodes = 行索引（主键），features = 列名，matrix.mtx = 稀疏矩阵的 COO 存储；UMI ≈ 去重后的计数器。

## 常见误读与注意

- 一个“细胞”是一个条形码，可能是空液滴、死细胞或双细胞（doublet），要靠 QC 过滤。
- UMI 数依赖测序深度，不同细胞不可直接比较（需归一化）。

## 在本项目中

- 每个样本用 `sc.read_10x_mtx` 读取；barcode 加样本名前缀避免重名；六个 RNA 矩阵 outer join；每个输入含 32,285 个 Gene Expression features 与 3 个 Multiplexing Capture features（BM、Spleen、Thymus），capture features 只被记录，未用于重新拆分或过滤（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 1）。
- 原始 18,458 个细胞（Ctrl_BM 3,998；Ctrl_Spleen 3,829；Ctrl_Thymus 2,101；T2_BM 3,131；T2_Spleen 3,760；T2_Thymus 1,639）→ QC 后 15,532（同上）。
- 数据压缩包 `input.zip`（约 7 GB，未解压）：`README.md`。

## 相关概念

- [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]] — 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。
- [[质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）|QC-Metrics]] — QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。
- [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] — 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。
- [[拟时序（Diffusion Map, DPT, PAGA）与 RNA 速度（Pseudotime vs RNA Velocity）|Pseudotime-DPT-vs-RNA-Velocity]] — 拟时序按表达相似度给细胞排序；RNA velocity 用剪接/未剪接 RNA 估计变化方向，需要额外输入。
