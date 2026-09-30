# KEGG 通路数据库（KEGG Pathway Database）

> KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。

## 定义

KEGG（Kyoto Encyclopedia of Genes and Genomes）把基因组织成通路（如 MAPK signaling pathway）和疾病条目（如 Measles）。与 GO 的树/图结构不同，KEGG 每个条目是一份基因列表加一张通路图。

## CS 类比

KEGG ≈ 一份手工维护的“功能模块清单”，每项是若干基因的集合，而非层次结构。

## 常见误读与注意

- 疾病名条目通常由共享的免疫/炎症/应激基因驱动（[[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]]）。
- 库版本、物种、背景基因不同，条目结果无法与旧 PPT 精确复现。

## 在本项目中

- 本项目用冻结的 Enrichr `KEGG_2019_Mouse` 库，共 298 个实际检验的条目（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 3 节；库文件来源见 `notebooks/scripts/iNKT/run_inkt_c5_paper_followup.py`）。
- 旧 PPT 130 条通路记录逐项核查，部分旧库/条目不能精确复现；Blood S3 48 条中 44 条可检验（`docs/inkt_history_and_supervisor_requirements_20260915.md` 阶段 7）。
- 9/29 Yue：除 GO BP 外，把 KEGG 和 GO MF 都做成网络，显著性 q≤0.05 的条目都进图（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:19:09–16:19:45；ASR 里 “password/cag” 是 pathway/KEGG 的误识，推测）。
- 09-30 分析（已完成）：把 KEGG 与 GO BP/MF 一起做 ORA 与网络；KEGG 锚点 IL-17 q=0.071、MAPK 0.64、TNF 0.56，见 [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]]。

## 相关概念

- [[Enrichr 基因集库（Enrichr Gene-set Libraries, KEGG_2019_Mouse）|Enrichr-Gene-Set-Libraries]] — Enrichr 是收集大量基因集库的网站；本项目从中下载 KEGG_2019_Mouse 并冻结成 GMT 文件。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]] — 通路/GO 名字是历史命名，命中的其实是共享的免疫、炎症、核糖体或热休克基因。
- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 09-30 已完成：GO BP/MF + KEGG 混合网络、锚点与分组结果。
