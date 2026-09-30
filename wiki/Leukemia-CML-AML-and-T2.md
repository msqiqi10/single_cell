# 白血病 CML / AML 与 “T2” 条件（Leukemia, CML, AML and T2）

> T2 = 肿瘤条件；本项目文件没有说明 T2 具体对应哪种白血病，只能说是 tumour。

## 定义

白血病（leukemia）是血细胞恶性增殖；慢性髓系白血病（chronic myeloid leukemia, CML）进展慢，急性髓系白血病（acute myeloid leukemia, AML）进展快。Rob 实验室用嵌合小鼠模型研究 CML（Yue 在转录稿约 00:34:00 的介绍）。

## CS 类比

Ctrl vs T2 ≈ 对照组 vs 实验组（A/B test）；T2 是“加了肿瘤压力”的那一臂。

## 常见误读与注意

- 不确定：项目文件未记录 T2 的具体模型、时间点和处理（`iNKT_by_date/2026-09-25/results/README.md` 第 7 节列为“尚需合作者提供”）。不要把 T2 直接写成 CML 或 AML。

## 在本项目中

- ZERU 在 9/24 会议上说 T2 是 “the tumors dataset”（转录稿 00:34:51）；早期文档把 Ctrl/T2 写成 control/tumor（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 4）。
- AML 患者与小鼠单细胞数据（全骨髓与 NK-only，约 40 位患者），Rob 愿意通过 Box 分享 FASTQ，尚未收到（转录稿 01:31:12–01:34:42；`iNKT_by_date/2026-09-25/README.md` 外部依赖）。
- Rob 提出可比较 CML 与 AML 的 NK 基线差异（转录稿 01:32:52）。

## 相关概念

- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]] — 每个样本标签由 3 只小鼠混合而成：每组只有 n=1 个 pool。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
- [[人物与各自关心的点（Who Is Who）|Who-is-who]] — Rob / Yue / Rajesh 分别关心什么。
