# 细胞毒性（Cytotoxicity: Perforin, Granzymes）

> 细胞毒性指杀伤靶细胞的能力，转录层面常用 Prf1/Gzma/Gzmb 三个基因的均值近似。

## 定义

细胞毒性（cytotoxicity）指免疫细胞杀死靶细胞。穿孔素（perforin，基因 Prf1）在靶细胞膜上打孔，颗粒酶（granzyme，Gzma、Gzmb）经孔进入并诱导凋亡。

## CS 类比

穿孔素 ≈ 开门的钥匙，颗粒酶 ≈ 从门里投放的载荷。

## 常见误读与注意

- 三基因均值只是“转录评分”，不是杀伤率；签名依赖强：加 Ctla2a 后 BM C4/Spl C3 均不显著（q=0.93、0.53）。
- 杀伤相关 GO 条目没有通过 q≤0.05，最小 q=0.632。

## 在本项目中

- 核心模块 Prf1/Gzma/Gzmb：骨髓 C4 T2−Ctrl 差 +0.08620（q=0.02513）；脾脏 C3 −0.02297（q=0.001538）；加 Nkg7 后为 +0.07593/−0.03367（`docs/audits/inkt_cytotoxicity_source_audit_20260922.md`，数据来自 `iNKT_by_date/2026-09-19/results/tables/cytotoxicity_contrasts.csv`）。
- 93 条合格模块对照中 11 条通过 BH q≤0.05，不能算 11 个独立发现（`iNKT_by_date/2026-09-19/README.md`）。
- 19 页报告第 3–6、16–17 页（`iNKT_by_date/2026-09-19/presentation/slide_index.csv`）。
- 矩阵覆盖：Prf1、Gzma、Gzmb、Nkg7、Gzmm、Ctsw、Fasl、Ccl5、Ifng 等存在；Gzmk、Gnly、Gzmh、Ncr3 不在过滤后矩阵（同 A 文档）。

## 相关概念

- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[标志基因与模块分数（Marker Genes, Module Scores）|Marker-Genes-and-Module-Scores]] — marker 用来给簇命名，模块分数把一组基因合成一个数。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
