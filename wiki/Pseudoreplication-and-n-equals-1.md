# 伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）

> 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。

## 定义

伪重复（pseudoreplication）指把非独立的观测当作独立重复。scRNA-seq 中同一动物/文库的细胞高度相关，真正的重复单位是动物（或独立 pool）。

## CS 类比

≈ 对同一台机器测一万次延迟，然后宣称“这类机器”的延迟有显著差异：样本量其实是机器数，不是测量次数。

## 常见误读与注意

- 随机拆细胞、QC 重抽样、bootstrap 都不能变成动物重复。
- 因此本项目所有基因/富集 FDR 都是细胞层面的探索性统计，“可重复方向、效应大小、跨分析一致性”才是重点。

## 在本项目中

- Rob 确认每个标签由 3 只小鼠混合（转录稿 00:37:00），所以每组 n=1 个 pool（`iNKT_by_date/2026-09-25/README.md`）。
- `docs/inkt_pipeline_stage_by_stage_walkthrough.md` “统计设计只有六个样本”与 8.6；`iNKT_by_date/2026-09-25/results/README.md` 第 7 节引 Squair 等 2021。
- bootstrap 区间只说明细胞重抽样稳定性，不是生物学重复的置信区间（同 W 6.4）。

## 相关概念

- [[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]] — 每个样本标签由 3 只小鼠混合而成：每组只有 n=1 个 pool。
- [[QC 匹配与标准化均值差（QC Matching, SMD）|QC-Matching-and-SMD]] — 让 Ctrl 与 T2 细胞的测序深度、检出基因数、线粒体比例匹配后再比，检查结果是否由 QC 差异造成。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
- [[p 值（p-value）|P-value]] — p 值是“如果两组其实没差别，看到这么大差异的概率”，不是效应大小，也不是结论正确的概率。
