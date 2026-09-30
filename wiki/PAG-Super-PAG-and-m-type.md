# PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）

> PAG = 通路/注释基因列表/基因签名；PAGER 定义 PAG 之间的关系；Yue 用 super-PAG 指“把许多 PAG 聚成的功能大组”。

## 定义

PAGER 平台把 pathways、annotated gene lists、gene signatures 统称 PAG（PAGs）。PAGER 还提供 PAG 与 PAG 之间的关系表（m-type、r-type）。super-PAG（super-PAG, SP）是 Yue 在 GOLDEN 中的说法：先把许多 PAG 聚成少数功能组，再看组内功能和组间联系。

## CS 类比

PAG ≈ 一个“模块”；super-PAG ≈ 把模块聚合后的“包/命名空间”；m-type 关系 ≈ 模块之间的依赖边。

## 常见误读与注意

- 不同项目里 SP1、SP4 等标签指的是不同的东西；Yue 展示的是 TCGA_BRCA_GOLDEN 示例，不是 iNKT 结果。
- super-PAG 是功能组，不是新的细胞亚型（N25 第 2 节）。
- m-type 的精确定义在项目文件中没有给出；只可确认 PAGER 表格描述 m-type/r-type 为 PAG-PAG 关系（`docs/pager-scFGA.pdf`），本项目用 Jaccard 网络替代（不确定其数学等价性，`fusion_manifest.json` 明确写明“不完全相同”）。

## 在本项目中

- PAGER-scFGA（`docs/pager-scFGA.pdf`）：PAGER 分析部分给出 m-type/r-type PAG-PAG 关系表。
- 本项目的替代：Jaccard≥0.25 且共享≥3 的冻结小鼠 GO 网络（`iNKT_by_date/2026-09-25/results/fusion_manifest.json`）。
- 9/23 会议：Yue 的“大象”比喻——只看一个簇像只摸到象的一只角，super-PAG 把几百个条目分块并高亮（`/Users/zeruzhang/Downloads/Zoom会议/zoom-0923/GMT20260923-180025_Recording.transcript.vtt` 约 00:39:51–00:41:15；ASR 中文，字词有误）。
- 9/29：Yue 建议构建“大图 + 小图”两层：小图为 super-PAG，大图为详细网络（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:19:45–16:20:10）。

## 相关概念

- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[GOLDEN 与融合聚类（GOLDEN, Fusion）|GOLDEN]] — GOLDEN 是 Yue 团队的“功能条目语义 + 关系图”聚类方法；本项目做的是小鼠 GO 本地适配版。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[Enrichr 基因集库（Enrichr Gene-set Libraries, KEGG_2019_Mouse）|Enrichr-Gene-Set-Libraries]] — Enrichr 是收集大量基因集库的网站；本项目从中下载 KEGG_2019_Mouse 并冻结成 GMT 文件。
