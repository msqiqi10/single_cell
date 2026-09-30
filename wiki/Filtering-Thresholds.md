# 过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）

> 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。

## 定义

过滤（filtering）同时发生在基因和细胞两个方向：基因必须在至少 100 个细胞中被检出（`min_cells=100`）；细胞必须 200 ≤ n_genes < 2500 且线粒体比例 < 5%。

## CS 类比

类似日志采样：先丢掉出现次数少于 100 的 key（稀有事件），再丢掉异常请求。稀有但重要的 key（Ccr6）会一起被丢掉。

## 常见误读与注意

- “不在矩阵里”不等于“不表达”：Ccr6 在原始矩阵中有 92 个细胞检出（低于 100），因此被过滤；不能据此写成 iNKT17 不存在。
- 不同顺序会得到不同的细胞集合。

## 在本项目中

- 基因过滤：32,285 → 10,670，移除 21,615 个低频基因；细胞 18,458 → 15,532（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 2）。
- Ccr6：原始 92 个细胞检出；保留细胞中仍有 60 个，其中 54 个在 C5-2；Il17a 保留细胞 0 个（原始输入条码中 3 个），Il17f 保留细胞 2 个（原始 3 个）：`iNKT_by_date/2026-09-25/results/README.md` 第 4 节；`iNKT_by_date/2026-09-25/results/tables/full_feature_marker_status.csv`。
- 旧稿曾把 Il17a/Il17f/Ccr6 “不在 gene universe” 解释为无法评价（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 6.2）；9/25 用原始矩阵补齐。
- 9/24 会议上 ZERU 说 Ccr6 “was not retained in the future gene set”（转录稿 00:42:37–00:42:44）。

## 相关概念

- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）|QC-Metrics]] — QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。
- [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]] — 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。
- [[高变基因（Highly Variable Genes, HVG）|Highly-Variable-Genes]] — 选 3,000 个在细胞间波动最大的基因给 PCA 用，降低噪声和计算量。
