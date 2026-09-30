# 小鼠模型与混样（Mouse Model and Pooled Samples）

> 每个样本标签由 3 只小鼠混合而成：每组只有 n=1 个 pool。

## 定义

“混样”（pooling）指把多只动物的细胞混在一起测序。这节省成本，但个体差异无法被拆分：一个 pool 只算一个生物学观测。

## CS 类比

相当于把 3 台机器的日志合并成一份再做统计：你无法再估计机器之间的方差。

## 常见误读与注意

- 随机拆细胞或 QC 重抽样不能产生动物重复（`iNKT_by_date/2026-09-25/results/README.md` 第 7 节，引 Squair 等 2021）。

## 在本项目中

- Rob 在 9/24 会议确认：“3 animals went into each label”（转录稿 00:37:00–00:37:12）。转录稿 00:36:45 中 ZERU 那句“can now treat these 6 labels as 6 independent animal replicates”疑似 ASR 把 “cannot” 识别成 “can”（推测），以 `iNKT_by_date/2026-09-25/README.md` 的结论为准。
- 每个 tissue×condition 只有一个样本 → 没有独立动物重复（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 开篇）。
- 样本设计登记表：`iNKT_by_date/2026-09-25/results/existing_evidence/sample_design_registry.csv`。
- 仍待确认：pool 与动物对应、跨组织是否取自同一批动物、批次、时间点（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 1 节）。

## 相关概念

- [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] — 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。
- [[白血病 CML / AML 与 “T2” 条件（Leukemia, CML, AML and T2）|Leukemia-CML-AML-and-T2]] — T2 = 肿瘤条件；本项目文件没有说明 T2 具体对应哪种白血病，只能说是 tumour。
- [[批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）|Harmony-and-Batch-Integration]] — 批次效应是技术带来的系统偏差；Harmony 把不同批次的细胞在低维空间“对齐”。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
