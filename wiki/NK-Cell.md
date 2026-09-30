# 自然杀伤细胞（Natural Killer Cell, NK cell）

> NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。

## 定义

自然杀伤细胞（NK cell）无需预先致敏就能杀伤病毒感染或肿瘤细胞，杀伤机制包括释放穿孔素（perforin）和颗粒酶（granzyme），并分泌 IFN-γ 等细胞因子。

## CS 类比

NK 细胞 ≈ 常驻的巡检守护进程；T 细胞 ≈ 按需部署的专用 worker。

## 常见误读与注意

- 本项目的细胞是 iNKT，不是 NK。用 NK 论文的图式做参照，不等于两个研究的 cluster 一一对应（ZERU 在 9/24 会议上明确这一点，转录稿 00:35:30）。

## 在本项目中

- 参照论文：PAGER-scFGA（小鼠 NK CITE-seq，`docs/references/pager-scFGA.pdf`）、Borra et al. 2026 GAFA（CML 缓解相关 NK，`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`）、Blood Advances（慢性炎症损害 NK 适应度与细胞毒性，`docs/references/blooda_adv-2024-014592-main.pdf`）。
- NK 参考 marker：CD56（Ncam1）在全部六个 raw feature 表存在但累计读数为 0；CD94、IL-4 用于 C5 对照（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 5）。
- 09-19 报告采用 PAGER 的细胞毒性模块思路（`iNKT_by_date/2026-09-19/notes/method_adaptation.md`）。

## 相关概念

- [[iNKT 细胞（invariant Natural Killer T cell）|iNKT-Cell]] — iNKT 是表达“不变” TCR、同时带 NK 样特征的 T 细胞亚群，能快速分泌大量细胞因子。
- [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] — 细胞毒性指杀伤靶细胞的能力，转录层面常用 Prf1/Gzma/Gzmb 三个基因的均值近似。
- [[白血病 CML / AML 与 “T2” 条件（Leukemia, CML, AML and T2）|Leukemia-CML-AML-and-T2]] — T2 = 肿瘤条件；本项目文件没有说明 T2 具体对应哪种白血病，只能说是 tumour。
- [[免疫系统基础（Immune System Basics）|Immune-System-Basics]] — 免疫系统分先天免疫与适应性免疫；iNKT 处在两者之间。
