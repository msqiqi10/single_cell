# AP-1 转录因子复合体（AP-1: Fos/Jun family）

> AP-1 是由 Fos、Jun 家族二聚化形成的转录因子；本项目最扎实的富集结果之一。

## 定义

AP-1（activator protein 1）是 Fos 家族（Fos、Fosb、Fosl1、Fosl2）与 Jun 家族（Jun、Junb、Jund）成员形成的二聚体转录因子，控制免疫激活、增殖和应激相关基因。它也是 MAPK 通路的下游“输出”之一。

## CS 类比

AP-1 ≈ 一个由多个类共同实现的接口；Fos 与 Jun 是两类实现，必须成对才能调用。

## 常见误读与注意

- GO:0035976 “transcription factor AP-1 complex” 与 “AP-1 adaptor complex” 是不同条目（后者是膜运输蛋白复合体），已用准确 GO ID 区分。
- 该条目反映 AP-1 成员基因的转录富集，不测复合体形成、DNA 结合或蛋白活性。

## 在本项目中

- `iNKT_by_date/2026-09-25/results/tables/AP1_complex_GO_exact_ID.csv`：骨髓 C0 命中 Fos/Fosl2/Jun/Junb（背景 M=9562，词条大小 K=6，输入 n=225，命中 k=4；p=4.3e-6，家族 q=0.000261，全局 q=0.00286）；脾脏 C3 命中 Fos/Jun/Junb/Jund（家族 q=0.000193，全局 q=0.00229，见 README）。
- 骨髓 C4（Fos、Junb）和脾脏 C5-1（Jun、Junb）只在家族内通过（全局 q=0.159、0.102）。
- 20 次固定种子 QC 匹配后方向全部保持；最大 QC 标准化差异约 0.057（骨髓 C0）/0.070（脾脏 C3）：`iNKT_by_date/2026-09-25/results/README.md` 第 3 节。
- 基因层面：骨髓 C0 的 Junb（log2FC 0.688，FDR 6.53e-31）、Fos（0.932，1.97e-19）、Jun（0.614，8.48e-16）、Fosb、Fosl2、Nr4a1（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。
- 09-30 混合 GO/KEGG 网络：AP-1 复合体是 12 个锚点中**唯一显著**的（q≈1.7e-4；Fos/Fosl2/Jun/Junb/Jund，骨髓 C0、脾 C3；`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv`）。
- 基因层面：AP-1 与 KEGG MAPK / IL-17 / TCR / Th17 / TNF 各共享恰好 3 个基因——Fos + Jun，再加 Jund（MAPK、IL-17）、Junb（TNF）或 Nfatc2（TCR、Th17，非 DEG）之一；Fos/Jun 在 BM C0、Spleen C3 上调（Fos +0.93 / +1.78，Jun +0.61 / +1.24）（`iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv`）。Fosb、Dusp1、Nr4a1、Egr1 差异明显，但各自至多在一个锚点里，不构成桥。见 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。
- 解读（推断）：T2 特征是即刻早期 AP-1（Fos–Jun）反应；KEGG 通路只是因为基因集含 Fos/Jun 才与之“相连”；可能是生态位内激活，也可能是解离应激，需要台面实验区分。
- Rob：Fos/Jun/Fosb 是“对功能和生物学都重要”的 AP-1 转录因子（转录稿 01:23:19）。

## 相关概念

- [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] — MAPK 是一条三层激酶级联，末端 ERK/p38/JNK 会激活 AP-1 等转录因子。
- [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] — Fos/Jun 等即刻早期基因几分钟内就会被诱导，组织解离/处理本身也能诱导它们——这是 AP-1 结果的主要替代解释。
- [[转录因子（Transcription Factor, TF）|Transcription-Factor]] — 转录因子是决定“哪些基因被开启”的调控蛋白，相当于配置开关。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[目前的主要发现与证据等级（Key Findings so Far）|Key-Findings]] — 结论与证据等级汇总。
- [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] — AP-1 与 KEGG 锚点之间唯一有差异表达支持的桥是 Fos + Jun。
- [[锚点节点（Anchor Nodes）|Anchor-Nodes]] — AP-1 是唯一显著的锚点。
