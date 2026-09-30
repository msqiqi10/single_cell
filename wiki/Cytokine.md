# 细胞因子（Cytokine: IFN-γ, IL-4, IL-17）

> 细胞因子是免疫细胞之间传递指令的分泌蛋白。

## 定义

细胞因子（cytokine）是细胞分泌的小蛋白，用来在细胞之间传递信号。IFN-γ（interferon gamma，基因 Ifng）促炎、激活巨噬细胞；IL-4（Il4）偏向抗寄生虫/体液反应；IL-17（Il17a/Il17f）参与黏膜防御和炎症。（背景知识。）

## CS 类比

细胞因子 ≈ 微服务之间的消息队列消息；受体（如 Il1r1、Il6ra）是订阅者。

## 常见误读与注意

- 本数据只测 mRNA。Il17a/Il17f 的 RNA 检出极低（0 和 2 个细胞）并不表示细胞不分泌 IL-17——刺激后才会大量分泌，需要蛋白/分泌实验。

## 在本项目中

- `iNKT_by_date/2026-09-25/results/tables/full_feature_marker_status.csv`：Il17a 保留细胞中 0 个检出，Il17f 2 个。
- Rob 点名 Il1r1、Il6ra 与 IL-17 相关信号；骨髓 C0 中它们的基因 FDR 分别约 0.175、0.124，未过 0.05（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。
- IL-17 signaling pathway、Th17 cell differentiation 是 09-30 分析的 anchor KEGG 通路（`iNKT_by_date/2026-09-30/code/common.py`）。
- 实验读出建议：身份、数量、分泌三者分开（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md` B 节）。

## 相关概念

- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[干扰素反应（Interferon Response）|Interferon-Response]] — 干扰素是抗病毒/炎症信号；Rob 用它解释为什么 “measles” 会出现在结果里。
- [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] — 细胞毒性指杀伤靶细胞的能力，转录层面常用 Prf1/Gzma/Gzmb 三个基因的均值近似。
- [[基因、mRNA 与蛋白质（Gene / mRNA / Protein, Central Dogma）|Gene-mRNA-Protein]] — 基因是 DNA 上的“函数定义”，mRNA 是一次“调用的输出缓存”，蛋白质是最终干活的“进程”；scRNA-seq 只数 mRNA。
