# MAPK 级联（MAPK Cascade: ERK, p38, JNK）

> MAPK 是一条三层激酶级联，末端 ERK/p38/JNK 会激活 AP-1 等转录因子。

## 定义

丝裂原活化蛋白激酶（mitogen-activated protein kinase, MAPK）级联由 MAPKKK→MAPKK→MAPK 逐层磷酸化构成，三条主要臂为 ERK、p38、JNK，末端可激活 Fos/Jun（AP-1）。Dusp1 是使 MAPK 去磷酸化的反馈基因。（背景知识。）

## CS 类比

MAPK ≈ 三级中间件管道；磷酸化 ≈ 每层的“acknowledge”标志；Dusp ≈ 撤销处理器。

## 常见误读与注意

- AP-1 基因升高与 MAPK 通路显著是两件事：骨髓 C0 的 GO “MAPK cascade”（GO:0000165）家族 FDR=1，旧 KEGG MAPK 在 paper 规则下约 0.078、robust 规则下约 0.627——均不显著。
- 不能因 AP-1 显著就说 MAPK 被激活。

## 在本项目中

- Rob：“MAP kinase signaling, super important”，与台面实验一致（转录稿 01:22:19、01:28:20）。
- `iNKT_by_date/2026-09-25/results/README.md` 第 3 节；`iNKT_by_date/2026-09-25/results/BRIEF_EN.md`：MAPK activation remains unproven。
- 09-30 分析把 MAPK cascade、ERK1/ERK2 cascade、p38MAPK cascade、JNK cascade 和 KEGG MAPK signaling pathway 设为 anchor（`iNKT_by_date/2026-09-30/code/common.py`）。结果：这些 MAPK 锚点全部不显著（KEGG MAPK q=0.64，GO 的 MAPK/ERK/p38/JNK cascade q=1），KEGG MAPK 与 AP-1 之间在基因层面只共享 Fos、Jun、Jund，MAPK 骨架基因（Mapk1/3/14 等）不是 DEG（`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv`、`iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv`）；见 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] 与 [[锚点节点（Anchor Nodes）|Anchor-Nodes]]。
- Dusp1 在骨髓 C0 T2 更高（log2FC 0.804，FDR 5.35e-19）：`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节。
- 实验建议：蛋白磷酸化读出由 Rob 团队确定（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md` A 节）。

## 相关概念

- [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] — AP-1 是由 Fos、Jun 家族二聚化形成的转录因子；本项目最扎实的富集结果之一。
- [[信号通路（Signalling Pathway）|Signalling-Pathway]] — 信号通路是“受体→级联→转录因子→基因”的信息传递链。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[人物与各自关心的点（Who Is Who）|Who-is-who]] — Rob / Yue / Rajesh 分别关心什么。
- [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] — 基因层面：MAPK 与 AP-1 只共享 Fos/Jun/Jund，骨架基因不是 DEG。
- [[锚点节点（Anchor Nodes）|Anchor-Nodes]] — MAPK 锚点为何不显著也保留。
