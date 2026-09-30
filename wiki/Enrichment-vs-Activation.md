# 富集不等于激活（Enrichment ≠ Activation）

> “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。

## 定义

富集（enrichment）是统计重叠；激活（activation）是生物学状态（蛋白磷酸化、DNA 结合、通量）。前者可由 mRNA 表达提供线索，后者需要专门实验。

## CS 类比

≈ “日志里出现了 error 关键词很多次” ≠ “服务宕机了”。

## 常见误读与注意

- 本项目多处显式警示：AP-1 富集 ≠ MAPK 激活；OXPHOS 富集 ≠ 代谢通量；细胞毒性评分 ≠ 杀伤实验；G02 ≠ 蛋白合成速率。

## 在本项目中

- `iNKT_by_date/2026-09-25/results/README.md` 第 3、5 节；`iNKT_by_date/2026-09-19/README.md` “解释边界”；`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 10.3。
- 实验讨论清单为每个候选写明“会削弱该解释的情形”（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md`）。

## 相关概念

- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] — MAPK 是一条三层激酶级联，末端 ERK/p38/JNK 会激活 AP-1 等转录因子。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[目前的主要发现与证据等级（Key Findings so Far）|Key-Findings]] — 结论与证据等级汇总。
