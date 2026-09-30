# 免疫系统基础（Immune System Basics）

> 免疫系统分先天免疫与适应性免疫；iNKT 处在两者之间。

## 定义

免疫系统（immune system）负责识别并清除病原体和异常细胞。先天免疫（innate）反应快、识别模式固定（如 NK 细胞）；适应性免疫（adaptive）反应慢、能针对特定抗原学习（如 T、B 细胞）。淋巴细胞（lymphocyte）在骨髓（bone marrow）和胸腺（thymus）产生或成熟，在脾脏（spleen）等外周器官驻留和活化。（以上为通用背景知识。）

## CS 类比

先天免疫 ≈ 内核里写死的规则引擎（防火墙规则）；适应性免疫 ≈ 会根据流量日志更新的机器学习模型。

## 常见误读与注意

- 免疫细胞“活化”常伴随大量应激/炎症基因上调，这不等于抗肿瘤功能增强。

## 在本项目中

- 本项目细胞是 iNKT 细胞，兼有 T 细胞受体和 NK 样标志（marker panel：iNKT/T lineage 与 NK/cytotoxic，见 `docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 6.1 节）。
- 项目参照的论文均涉及 NK 细胞与白血病微环境（`docs/references/blooda_adv-2024-014592-main.pdf`、`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`）。
- 与之相关的 KEGG/GO 疾病式名称（如 measles）指向的是免疫/炎症基因集，见 [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]]。

## 相关概念

- [[T 细胞与 T 细胞受体（T Cell, TCR）|T-Cell]] — T 细胞用 TCR 识别抗原；iNKT 用一种近乎固定的“不变” TCR。
- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[iNKT 细胞（invariant Natural Killer T cell）|iNKT-Cell]] — iNKT 是表达“不变” TCR、同时带 NK 样特征的 T 细胞亚群，能快速分泌大量细胞因子。
- [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] — 细胞因子是免疫细胞之间传递指令的分泌蛋白。
- [[组织：骨髓、脾脏、胸腺（Bone Marrow, Spleen, Thymus）|Tissues-BM-Spleen-Thymus]] — 骨髓造血、胸腺培育 T 细胞、脾脏是血液过滤与免疫应答场所；三者是不同的“微环境”。
