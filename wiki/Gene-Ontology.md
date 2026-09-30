# 基因本体（Gene Ontology, GO: BP / MF / CC）

> GO 是一棵（准确说是有向无环图）术语字典，用来给基因标注它参与的过程、分子功能和所在位置。

## 定义

基因本体（Gene Ontology, GO）把基因功能分为三个分支：生物过程（biological process, BP）、分子功能（molecular function, MF）、细胞组分（cellular component, CC）。术语之间用有向无环图（directed acyclic graph, DAG）连接，关系包括 is_a（是一种）和 part_of（是……的一部分）。基因注释会沿关系向上“传播”到祖先术语。注释带证据代码，如 IEA（自动推断的电子注释）、NOT（明确“不具有”该功能）、ND（无数据）。

## CS 类比

GO ≈ 一个带继承（is_a）和组合（part_of）关系的类型系统；给基因贴标签 ≈ 接口实现声明，向上传播 ≈ 子类自动满足父类接口。

## 常见误读与注意

- 术语的名称不一定符合直觉，见 [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]]。
- 术语很多共享同一批基因，条目数量 ≠ 独立机制数量（[[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]]）。
- GO 的“regulates”关系本项目不传播。

## 在本项目中

- GO 官方本体 `releases/2026-07-26`，仅小鼠 MGI 注释；排除 NOT 和 ND，保留 IEA；只沿 is_a/part_of 向上继承，不沿 regulates；词条大小 5–500（`iNKT_by_date/2026-09-19/README.md` GO 方法；来源与哈希 `iNKT_by_date/2026-09-19/sources/GO_provenance.json`）。
- 骨髓 C0 独立复算：BP/MF/CC 共 17,912 项 GO 检验（`iNKT_by_date/2026-09-25/notes/local_inkt_README.md`）。
- 关键条目：GO:0035976 transcription factor AP-1 complex（CC）；GO:0000165 MAPK cascade（BP）；GO:0042267 NK 细胞介导的细胞毒性（29 个基因）、GO:0001906 cell killing（66）、GO:0140507（10）（`docs/audits/inkt_cytotoxicity_source_audit_20260922.md`）。
- 祖先路径表：`iNKT_by_date/2026-09-25/results/network/ontology_ancestry_paths.csv`。
- 9/29：Yue 要求用 GO MF 再做一遍（以前只用 BP）（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:21:14–16:21:31）。

## 相关概念

- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[术语冗余与 Jaccard 相似度（Term Redundancy, Jaccard Index）|Term-Redundancy-and-Jaccard]] — 很多 GO/KEGG 术语共享同一批基因；Jaccard = 交集/并集，用来量化两个术语有多像。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] — Yue 的问题：输入了多少基因，图上只展示了多少，还有多少没被任何显著条目覆盖？
