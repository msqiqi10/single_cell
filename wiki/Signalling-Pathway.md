# 信号通路（Signalling Pathway）

> 信号通路是“受体→级联→转录因子→基因”的信息传递链。

## 定义

信号通路（signalling pathway）指细胞把外界信号（细胞因子、抗原）传到细胞核的分子链：受体接收信号，一串蛋白依次磷酸化（激酶级联），最终激活转录因子改变基因表达。

## CS 类比

通路 ≈ 事件驱动的调用链：event → handler1 → handler2 → … → 配置更新。

## 常见误读与注意

- 数据库里的“通路”是人工整理的基因集合（见 [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]]、[[KEGG 通路数据库（KEGG Pathway Database）|KEGG]]），并不是一次测量；基因集富集不等于通路被激活（[[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]]）。

## 在本项目中

- Rob 优先的通路：MAPK 信号、AP-1 转录因子、iNKT17（转录稿 01:22:19、01:33:23）。
- 本项目区分了基因层面（Fos、Jun 显著）与通路层面（MAPK cascade FDR=1）的证据（`iNKT_by_date/2026-09-25/results/README.md` 第 3 节）。
- 19 页/28 页报告按 PAGER、GAFA 的“功能富集→气泡矩阵”逻辑展示（`iNKT_by_date/2026-09-20/package/README.zh.md`）。

## 相关概念

- [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] — MAPK 是一条三层激酶级联，末端 ERK/p38/JNK 会激活 AP-1 等转录因子。
- [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] — AP-1 是由 Fos、Jun 家族二聚化形成的转录因子；本项目最扎实的富集结果之一。
- [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] — GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
