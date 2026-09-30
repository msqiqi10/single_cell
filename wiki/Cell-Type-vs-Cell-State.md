# 细胞类型与细胞状态（Cell Type vs Cell State）

> “类型”是相对稳定的身份，“状态”是同一类细胞当前的活动模式；本项目的 C0…C7 是算法分出的表达状态簇，不是已验证的细胞类型。

## 定义

细胞类型（cell type）指在发育中确定、相对稳定的身份，如 T 细胞、NK 细胞；细胞状态（cell state）指同一类型的细胞因环境、激活、应激而暂时呈现的表达模式。

## CS 类比

类型 ≈ 类（class）；状态 ≈ 对象在某一时刻的字段取值。聚类只能按“字段取值相似”分组，并不知道两个对象是否属于同一个类。

## 常见误读与注意

- 聚类编号（C0、C1…）由算法输出，没有固定生物学含义，不能把旧 PPT 的 c1 直接当作新分析的 c6（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.5 节）。
- “C5-2 更像 iNKT17”是身份特征；“T2 是否改变 C5-2 的 IL-17 功能”是另一个问题（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。

## 在本项目中

- 组织（BM/脾/胸腺）是数据中最强的结构：C0 是骨髓锚点，C3 是脾脏锚点，C6 是胸腺锚点（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.4 节）。
- Rob 强调很多“组织差异”反映的是细胞所处邻域（环境）而不是单个细胞的生物学（转录稿约 01:22:42–01:23:10）。
- 旧 PPT 11 群 vs 当前 8 群（+C5 拆分为 9 个 refined cluster）：`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 3 与 `docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md`。

## 相关概念

- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
- [[Leiden 聚类（Leiden Clustering）|Leiden-Clustering]] — Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。
- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[标志基因与模块分数（Marker Genes, Module Scores）|Marker-Genes-and-Module-Scores]] — marker 用来给簇命名，模块分数把一组基因合成一个数。
