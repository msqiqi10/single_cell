# 转录因子（Transcription Factor, TF）

> 转录因子是决定“哪些基因被开启”的调控蛋白，相当于配置开关。

## 定义

转录因子（transcription factor, TF）是结合在 DNA 上、控制某些基因转录多少的蛋白质。一个 TF 可以同时调控几十到上千个下游基因。

## CS 类比

TF ≈ 配置文件里的 feature flag / 中间件路由：它本身是一个程序，但作用是决定别的模块是否运行。

## 常见误读与注意

- 一个 TF 基因的 mRNA 升高，不一定意味着 TF 蛋白已入核、已结合 DNA。
- 同一个下游基因可被多个 TF 调控，反推 TF 需要额外证据。

## 在本项目中

- 本项目最核心的 TF 是 AP-1 家族（Fos/Jun，见 [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]]）；iNKT 亚型 marker 中的 Tbx21（T-bet）、Gata3、Rorc、Zbtb16 也是 TF（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 第 6 节的 marker panel 含 Tbx21、Zbtb16）。
- Rob 解释：Tbx21 又叫 T-bet，是诱导 interferon gamma 产生的转录因子（转录稿约 01:27:12–01:27:25）。
- Borra et al. 2026（`docs/references/Borra_et_al_2026_GAFA_CML_NK.pdf`）的 GAFA 流程含转录因子调控网络（SCORPION），本项目未运行（`iNKT_by_date/2026-09-19/README.md` “解释边界”）。

## 相关概念

- [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] — AP-1 是由 Fos、Jun 家族二聚化形成的转录因子；本项目最扎实的富集结果之一。
- [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] — iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。
- [[信号通路（Signalling Pathway）|Signalling-Pathway]] — 信号通路是“受体→级联→转录因子→基因”的信息传递链。
- [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] — Fos/Jun 等即刻早期基因几分钟内就会被诱导，组织解离/处理本身也能诱导它们——这是 AP-1 结果的主要替代解释。
