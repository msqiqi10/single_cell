# 人物与各自关心的点（Who Is Who）

## Rob（湿实验/细胞生物学合作者）

- 身份：Yue 介绍他是长期合作者，血液肿瘤专家，实验室用嵌合小鼠模型研究慢性髓系白血病（转录稿 00:33:18–00:34:00）。任务书称其在 UAB；参考论文 Borra et al. 2026 的作者中有 Robert S. Welner（UAB 血液肿瘤科，`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`），Rob 很可能就是他——这是**推断**，转录稿里没有出现全名。
- 他自称“不是计算生物学视角，是细胞生物学视角”，关心的是**能在台面上测的东西**：
  - **MAPK 信号、AP-1 转录因子（Fos/Jun/Fosb）、iNKT17/IL-17**：与他们的实验认识一致，值得优先（转录稿 01:22:19、01:23:19、01:33:23）。[[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]]
  - 核糖体蛋白（Rpl/Rps）的大量变化：可能是蛋白合成受损，**也可能只是细胞“angry”（处在炎症应激环境）**，不是新发现，也不是他优先的方向（01:21:52–01:22:13）。[[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]]
  - 组织比较很多只是“邻域/环境”差异，应聚焦“能反映细胞功能”的信号（01:22:42–01:23:10）。
  - 通路名字如 “measles” 只是免疫/炎症基因集，不是真的麻疹（01:26:35–01:28:20）。[[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]]
  - 点名基因：Dusp、Fos、Il1r1、Jun、Nr4a1、Fosb、Fosl2、Cd8a、Ciita、Il6ra（01:24:55–01:25:14）。
- 数据合作：AML 的全骨髓和 NK-only 单细胞数据（小鼠+人，提到约 40 位患者），通过 Box 分享 FASTQ，没有 BAM（01:31:12–01:34:42）；也提到刚被 JCI 接收的 NK 数据集（01:31:32–01:31:42）。
- 确认：每个样本标签混合 3 只小鼠（00:37:00）。[[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]]

## Yue（Zongliang Yue，PI / 计算方法负责人）

- 身份：PI（任务书）；PAGER-scFGA 与 Borra et al. 2026 的通讯作者（`docs/references/pager-scFGA.pdf`、`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`；Borra 论文署名 Auburn University）。
- 关心的点：
  - 完整候选基因列表，而不是只看 Top 基因；被 GO/通路覆盖了多少（01:39:30–01:42:06）。[[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]]
  - **GO/KEGG → 功能组（super-PAG）→ 关系网络**，GOLDEN / GOLDEN Fusion / GNN，且要能追溯到基因。[[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] [[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]]
  - 汇报要能解释“为什么这样算、能让实验做什么”（01:45:10–01:47:57）。
  - 验证：GNN 的增益如何衡量、独立数据/实验验证。[[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]]
  - 9/29：再做 GO MF 和 KEGG；交互式 HTML；基因层面展开；核糖体/应激“不新颖”。[[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]]
- 事实来源：`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md`；09-29 转录稿。

## Rajesh（Rajesh Amin，Auburn 药物发现与开发）— 独立项目

- Yue 介绍他来自 Auburn University 药物发现与开发系，专长新药设计、疗效与毒性（转录稿 00:34:10–00:34:43）。
- 他在会末提到自己的项目：NIH 提供的额颞叶痴呆（frontotemporal dementia）人脑样本，药物改变核内转录本，想做信息学（01:36:06–01:38:17）。**与 iNKT 项目无关。**
- 另：`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 第 3 节末提到 09-05 会议另有一项肝/肾毒性数据任务（约 12 篇研究清单），同样不属于 iNKT R01–R13，当前没有开始。

## ZERU Zhang（分析者）

- 把数据分析、报告和汇报串起来，会上负责讲结果。本 wiki 面向的就是想接手/审阅这些分析的计算机背景读者。
