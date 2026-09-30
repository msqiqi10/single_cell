# 质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）

> QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。

## 定义

质量控制（quality control, QC）通过每个细胞的指标筛掉空液滴、受损细胞和双细胞：`n_genes_by_counts`（检出的基因数）、`total_counts`（UMI 总数）、`pct_counts_mt`（线粒体基因占比，高说明细胞破裂）。

## CS 类比

QC ≈ 数据清洗/输入校验：丢弃格式异常的记录。

## 常见误读与注意

- 本流程没有做 doublet 检测、环境 RNA 校正、cell-cycle 回归、批次整合，所以“通过 QC”不等于满足现代标准的全部 QC（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 2.4 节）。
- `n_genes_by_counts` 要在过滤基因之后重新计算，顺序不同得到的细胞集合不同。

## 在本项目中

- mt 基因：`mt-`/`MT-` 开头；ribosomal：`Rps/Rpl`；hemoglobin：`^(Hb[ab]|HBA|HBB)`（同上 2.2 节）。
- 保留细胞最低检出基因 284，最高 2,499；最高线粒体比例 4.99866%（同上 2.3 节）。
- 37 页版第 2 页展示样本、检出基因、UMI、线粒体比例（`iNKT_by_date/2026-09-20/package/slides.json`，old_02）。
- QC 精确重建证据：`docs/inkt_qc_exact_reconstruction.csv`。

## 相关概念

- [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] — 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。
- [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]] — 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。
- [[线粒体、氧化磷酸化与 ATP（Mitochondria, Oxidative Phosphorylation, ATP: OXPHOS）|Mitochondria-and-OXPHOS]] — 线粒体通过氧化磷酸化（OXPHOS）造 ATP；线粒体基因占比也是 QC 指标。
- [[QC 匹配与标准化均值差（QC Matching, SMD）|QC-Matching-and-SMD]] — 让 Ctrl 与 T2 细胞的测序深度、检出基因数、线粒体比例匹配后再比，检查结果是否由 QC 差异造成。
