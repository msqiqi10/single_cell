# T 细胞与 T 细胞受体（T Cell, TCR）

> T 细胞用 TCR 识别抗原；iNKT 用一种近乎固定的“不变” TCR。

## 定义

T 细胞（T cell）是适应性免疫的主力，在胸腺发育成熟；其表面的 T 细胞受体（T cell receptor, TCR）决定识别什么抗原，TCR 信号会激活转录因子并引发细胞因子产生。

## CS 类比

TCR ≈ 每个 T 细胞独有的“API 签名”；TCR 信号 ≈ 收到匹配请求后触发的回调链。

## 常见误读与注意

- “Direct TCR activation”是一个基因集分数，不是对 TCR 信号的直接测量。

## 在本项目中

- marker panel 中 iNKT/T lineage 基因：Trac、Trav11、Traj18、Cd3d、Cd3e、Zbtb16（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 6.1）。
- T2−Ctrl 的 signature 变化中，Direct TCR activation 全局约 +0.121，在骨髓 +0.129、脾脏 +0.131、胸腺 +0.099（同上 6.5 节）。
- KEGG “T cell receptor signaling pathway” 被列为 09-30 分析的 anchor 通路之一（`iNKT_by_date/2026-09-30/code/common.py`）。
- Rob 点名的 Cd8a、Ciita（抗原处理相关）在骨髓 C0 未过 FDR 0.05（Cd8a FDR≈0.0886）：`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节。

## 相关概念

- [[iNKT 细胞（invariant Natural Killer T cell）|iNKT-Cell]] — iNKT 是表达“不变” TCR、同时带 NK 样特征的 T 细胞亚群，能快速分泌大量细胞因子。
- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] — 细胞因子是免疫细胞之间传递指令的分泌蛋白。
- [[信号通路（Signalling Pathway）|Signalling-Pathway]] — 信号通路是“受体→级联→转录因子→基因”的信息传递链。
