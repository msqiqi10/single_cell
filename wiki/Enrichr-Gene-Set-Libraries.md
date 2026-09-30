# Enrichr 基因集库（Enrichr Gene-set Libraries, KEGG_2019_Mouse）

> Enrichr 是收集大量基因集库的网站；本项目从中下载 KEGG_2019_Mouse 并冻结成 GMT 文件。

## 定义

Enrichr 是提供数百个基因集库（gene-set library）的富集分析平台，库可以下载为 GMT 格式：每行是“术语名 \t 描述 \t 基因1 \t 基因2 …”。KEGG_2019_Mouse 是其中一个小鼠 KEGG 库。

## CS 类比

≈ 一份版本固定的 lookup table（term → set of genes）；“冻结”就是把它 pin 到某个 commit。

## 常见误读与注意

- Enrichr 的库会更新；不冻结就无法复现。
- GMT 里的基因符号大小写规则可能和数据不同，需要映射。

## 在本项目中

- 库来源标记 `downloaded_from_enrichr_and_frozen`，含 GMT SHA256（`notebooks/scripts/iNKT/run_inkt_c5_paper_followup.py`）。
- 37 页版第 23–25 页脚注：“our ORA: nominal DEG P ≤ 0.05, linear FC ≥ 1.5; frozen KEGG_2019_Mouse”（`iNKT_by_date/2026-09-20/package/notes/presenter_guide.zh.md`）。
- 09-30：KEGG GMT 大写符号按不区分大小写映射到测得的小鼠符号，映射统计见 `iNKT_by_date/2026-09-30/results/tables/KEGG_gene_mapping_summary.csv`。
- 本项目没有原始 PAGER/Reactome/KEGG GMT（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 10.2）。

## 相关概念

- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
- [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] — GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] — PAG = 通路/注释基因列表/基因签名；PAGER 定义 PAG 之间的关系；Yue 用 super-PAG 指“把许多 PAG 聚成的功能大组”。
