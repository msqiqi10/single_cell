# 干扰素反应（Interferon Response）

> 干扰素是抗病毒/炎症信号；Rob 用它解释为什么 “measles” 会出现在结果里。

## 定义

干扰素（interferon, IFN）分 I 型（IFN-α/β）和 II 型（IFN-γ）；IFN 通过 STAT1 等转录因子诱导一批干扰素刺激基因（ISG），使细胞进入抗病毒/炎症状态。

## CS 类比

IFN 反应 ≈ 收到“全局告警”后，全网节点一起进入防御模式。

## 常见误读与注意

- 数据库里的 “Measles” 通路命中的是干扰素/炎症相关基因，不是麻疹病毒感染。
- 就项目实际数据而言，“measles” 行 driver 是 Fos、Hspa1a、Hspa1b、Hspa2、Jun、Tlr7；Rob 讲的 STAT1/interferon 是免疫含义的解释，不是该行读出的基因（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。

## 在本项目中

- Rob 的解释：Ifng 和 Tbx21 高表达使 “measles” 出现在某些簇（转录稿 01:26:35–01:28:20）。
- 细胞元数据中含 `NK_type_I_IFN` 参考状态分数：`iNKT_by_date/2026-09-19/results/tables/cell_metadata_scores.csv.gz`。
- 参考状态覆盖不足，不能强制五分类（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 7）。

## 相关概念

- [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]] — 通路/GO 名字是历史命名，命中的其实是共享的免疫、炎症、核糖体或热休克基因。
- [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] — 细胞因子是免疫细胞之间传递指令的分泌蛋白。
- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
