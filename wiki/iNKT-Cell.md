# iNKT 细胞（invariant Natural Killer T cell）

> iNKT 是表达“不变” TCR、同时带 NK 样特征的 T 细胞亚群，能快速分泌大量细胞因子。

## 定义

不变自然杀伤 T 细胞（invariant natural killer T cell, iNKT cell）是一类 T 细胞，用近乎固定的 TCR 识别脂类抗原（由 CD1d 呈递），活化后迅速分泌 IFN-γ、IL-4、IL-17 等细胞因子，桥接先天与适应性免疫。（背景知识。）

## CS 类比

iNKT ≈ 一个“多态接口”对象：同一个不变签名，运行时可以表现为三种不同的子类（iNKT1/2/17）。

## 常见误读与注意

- iNKT 的亚型是连续程序的混合，本项目只把 C5-2 描述为“iNKT17 相关特征”，并非已证实的离散亚型（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 阶段 12）。

## 在本项目中

- 数据：6 个样本标签共 15,532 个 iNKT 细胞（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md`）。
- 12 套 signature（文献、in-house、Wang 2022 各含 iNKT1/2/17，加 Circulatory、Direct TCR activation、Residency）：同上 6.2 节。
- 研究问题：肿瘤 T2 条件是否改变骨髓、脾、胸腺中的 iNKT 细胞（`iNKT_by_date/2026-09-25/README.md`）。

## 相关概念

- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[T 细胞与 T 细胞受体（T Cell, TCR）|T-Cell]] — T 细胞用 TCR 识别抗原；iNKT 用一种近乎固定的“不变” TCR。
- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] — 细胞因子是免疫细胞之间传递指令的分泌蛋白。
- [[细胞类型与细胞状态（Cell Type vs Cell State）|Cell-Type-vs-Cell-State]] — “类型”是相对稳定的身份，“状态”是同一类细胞当前的活动模式；本项目的 C0…C7 是算法分出的表达状态簇，不是已验证的细胞类型。
