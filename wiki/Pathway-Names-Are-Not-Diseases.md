# 为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）

> 通路/GO 名字是历史命名，命中的其实是共享的免疫、炎症、核糖体或热休克基因。

## 定义

数据库条目按最初的实验或注释命名。“Measles” 来自“用麻疹病毒感染细胞后差异基因的整理”，命中的是干扰素/炎症基因；“translation at synapse” 是 GO 里对“突触处翻译”的条目，但其成员基因大多是核糖体蛋白。（Rob 的解释见转录稿。）

## CS 类比

≈ 变量名叫 `userLogin`，实际被别的模块复用来存别的东西：读实现（基因列表），不要只读名字。

## 常见误读与注意

- 不要因名字推断小鼠感染了麻疹或神经元形成了突触。
- 应始终展开到实际驱动基因（driver genes）。

## 在本项目中

- Rob：“Clearly these cells do not have measles”，它反映强烈的病毒/炎症反应；名称源于早期微阵列研究（转录稿 01:26:35–01:28:20）。
- 该行实际显示的 driver 是 Fos、Hspa1a、Hspa1b、Hspa2、Jun、Tlr7（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。
- “synapse translation”可由核糖体基因驱动，不能推断形成神经突触（`iNKT_by_date/2026-09-19/README.md`）；G03 组即由 “translation at synapse” 等条目构成（`iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`）。
- 旧 PPT 的 “Prion disease/Measles/Estrogen signaling” 同理（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 10.4）。

## 相关概念

- [[干扰素反应（Interferon Response）|Interferon-Response]] — 干扰素是抗病毒/炎症信号；Rob 用它解释为什么 “measles” 会出现在结果里。
- [[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] — 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。
- [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] — 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。
- [[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]] — 很多 GO/KEGG 术语共享同一批基因；Jaccard = 交集/并集，用来量化两个术语有多像。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
