# QC 匹配与标准化均值差（QC Matching, SMD）

> 让 Ctrl 与 T2 细胞的测序深度、检出基因数、线粒体比例匹配后再比，检查结果是否由 QC 差异造成。

## 定义

QC 匹配（QC matching）在每个组织×簇内，按 total_counts、检出基因数、线粒体比例的联合四分位分箱，抽取两组等量的细胞，重复多次。标准化均值差（standardized mean difference, SMD）= |均值差| / 合并标准差，用来检查匹配后两组的 QC 指标是否接近。

## CS 类比

≈ 在对照实验中做 propensity matching / 分层抽样，排除混杂变量。

## 常见误读与注意

- “匹配后方向一致”只说明不是 QC 造成，不能排除批次、解离应激，也不是动物重复。
- NA 表示匹配后每组不足 20 个细胞，不是“0 次同向”。

## 在本项目中

- 9/25：20 个固定种子（0–19），骨髓 C0 每条件保留 1,742 个细胞，最大 QC SMD≈0.057；脾脏 C3 每条件 2,061 个，最大≈0.070（`iNKT_by_date/2026-09-25/results/README.md` 第 3 节；`iNKT_by_date/2026-09-25/results/sensitivity_manifest.json`）。
- 9/19：五次 QC 匹配；骨髓 C4 5/5 同向、脾脏 C3 5/5 同向（`iNKT_by_date/2026-09-19/README.md`）。
- SMD 的计算代码：`iNKT_by_date/2026-09-25/code_local/targeted_sensitivity.py`。

## 相关概念

- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
- [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] — 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。
- [[质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）|QC-Metrics]] — QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。
- [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] — Fos/Jun 等即刻早期基因几分钟内就会被诱导，组织解离/处理本身也能诱导它们——这是 AP-1 结果的主要替代解释。
