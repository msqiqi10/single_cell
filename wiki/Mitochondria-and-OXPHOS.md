# 线粒体、氧化磷酸化与 ATP（Mitochondria, Oxidative Phosphorylation, ATP: OXPHOS）

> 线粒体通过氧化磷酸化（OXPHOS）造 ATP；线粒体基因占比也是 QC 指标。

## 定义

线粒体（mitochondrion）是细胞的“发电厂”；氧化磷酸化（oxidative phosphorylation, OXPHOS）由呼吸链复合体（Nduf*、Sdh*、Uqcr*、Cox*）和 ATP 合酶（Atp5*）组成，产生能量货币 ATP（adenosine triphosphate）。

## CS 类比

OXPHOS ≈ 数据中心的供电系统；ATP ≈ 电。

## 常见误读与注意

- 线粒体基因（mt-）比例高也是细胞受损的 QC 信号，所以“OXPHOS 上调”要区分是真实代谢状态还是质控/捕获偏差（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md` D 节）。
- 表达富集不等于代谢通量。

## 在本项目中

- QC 阈值：线粒体比例 <5%（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 2.3 节）。
- 骨髓 C0 与脾脏 C3 的 T2 上调基因富集氧化磷酸化（`iNKT_by_date/2026-09-19/README.md`）；剔除 Rpl/Rps/Hsp/Dnaj 后 OXPHOS 仍通过家族 FDR（`iNKT_by_date/2026-09-25/results/README.md` 第 3 节）。
- 功能组 G04（呼吸链、OXPHOS、ATP 生成）、G05（复合体组装）、G06（核苷酸代谢，与 G04 有 155 条条目边共享基因）：`iNKT_by_date/2026-09-25/results/network/module_edges.csv`。
- Yue（9/29）：G04–G06 提供了“额外信息”，不只是 G01 的应激（ASR 文本，约 16:40:25–16:40:47）。

## 相关概念

- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[质量控制指标（QC Metrics: n_genes, total counts, % mitochondrial）|QC-Metrics]] — QC 用三个数字识别坏细胞：检出基因数、总 UMI 数、线粒体比例。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
