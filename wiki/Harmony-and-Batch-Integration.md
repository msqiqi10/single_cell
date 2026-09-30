# 批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）

> 批次效应是技术带来的系统偏差；Harmony 把不同批次的细胞在低维空间“对齐”。

## 定义

批次效应（batch effect）指不同实验批次带来的非生物学差异。Harmony 等整合（integration）方法在 PCA 空间中迭代调整，使不同批次细胞混合。

## CS 类比

≈ 分布式系统里的时钟漂移校正：能对齐，但也可能把真实的差异一并抹掉。

## 常见误读与注意

- 本项目没有独立的技术批次元数据，因此“去批次”无法证明被移除的都是技术偏差；因为每个 tissue×condition 只有一个文库，条件效应与文库效应无法拆开（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 5.4）。

## 在本项目中

- R01：已尝试 Harmony，使组织更混合，同时改变了分群和条件结构；纯技术 batch 解释未解决（`docs/inkt_history_and_supervisor_requirements_20260915.md` 第 3 节表格；37 页版第 30 页 “Batch sensitivity”）。
- 9/5 会议的 Yue 要求尝试 batch removal 并比较前后结构（会议时刻 10:14–10:18，同 H 文档）。
- Borra et al. 2026 用 scVI 做批次整合（`docs/Borra_et_al_2026_GAFA_CML_NK.pdf`）；本项目没有运行 scVI/SCORPION（`iNKT_by_date/2026-09-19/README.md`）。

## 相关概念

- [[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]] — 每个样本标签由 3 只小鼠混合而成：每组只有 n=1 个 pool。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
- [[UMAP 降维可视化（UMAP: why position is not physical）|UMAP]] — UMAP 把近邻图摊平成二维图，位置没有物理意义。
- [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] — 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。
