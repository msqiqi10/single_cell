# 标志基因与模块分数（Marker Genes, Module Scores）

> marker 用来给簇命名，模块分数把一组基因合成一个数。

## 定义

标志基因（marker gene）是在某类细胞中特异高表达的基因。模块分数（module score / signature score）把一套基因的表达合成每个细胞一个数：本项目主分析先对每个基因做 z-score，再对成员等权平均，再跨细胞标准化到均值 0、标准差 1。

## CS 类比

≈ 特征工程：把多个相关特征聚合成一个复合特征。

## 常见误读与注意

- 分数不是概率，也不是表达比例；同一亚型不同来源的基因表重叠很低。
- 稀疏基因必须同时看 fold change、均值和阳性细胞比例：Rorc global log2FC +0.350，但阳性细胞只有约 1–1.5%（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 7）。

## 在本项目中

- 簇 marker：每簇 vs rest，Wilcoxon，仅 3,000 HVG，BH，保存 top 100（同上 6.1）。
- 12 套 signature 覆盖度、评分公式、bootstrap 判定规则（同上 6.2–6.4）；c5 → iNKT17-enriched，但保守表述为 mixed iNKT2/iNKT17。
- PAGER 模块：`PAGER_core`（Prf1/Gzma/Gzmb）与 `PAGER_E1`（加 Ctla2a）；`iNKT_by_date/2026-09-19/notes/method_adaptation.md`。
- 9/24 会议：ZERU 说明 score 是相对背景的表达，不是身份概率（转录稿 00:41:16）。

## 相关概念

- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] — 细胞毒性指杀伤靶细胞的能力，转录层面常用 Prf1/Gzma/Gzmb 三个基因的均值近似。
- [[差异表达（Differential Expression, DE; Wilcoxon）|Differential-Expression]] — DE 逐基因比较两组细胞的表达，本项目用 Wilcoxon 秩和检验比较 T2 与 Ctrl。
- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
