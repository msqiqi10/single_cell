# 交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）

> Yue 想要的最终交付：可点击的网络，节点大小 = 显著性，颜色 = 功能组，点开能看到基因。

## 定义

知识图谱（knowledge graph, KG）在本项目里指把术语（GO/KEGG）、功能组、基因和 DE 结果连成可交互的网络（HTML 文件），让合作者自己点击探索。

## CS 类比

≈ 前端仪表盘：节点/边是数据，点击是查询。

## 常见误读与注意

- 必须能追溯：功能组 → 原始条目 → 贡献基因 → 组织/簇/条件效应。

## 在本项目中

- Yue 9/29 的要求（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt`，ASR 中文）：节点大小 = −log10 p，颜色 = 功能组；点击节点 highlight 其邻居和边（“现在错综复杂，太乱”）；可选边的类型/颜色；展开到基因层面，重叠基因加星或红色边框；同时给表达图与网络图；HTML 交互版；表格要含统计和边数、被注释的边数；可用 GO BP、GO MF、KEGG（约 16:13:27–16:15:30；16:27:52–16:29:02；16:40:56–16:41:29）。
- ZERU 已有本地网页版，但新版渲染不能拖拽、看不到边连到哪个节点，需要重新生成（同一转录稿约 16:12:17–16:14:40）；把所有 table 打包一起发给 Yue。
- **09-30 两条线都已交付**：
  1. 本项目的通路/基因网络浏览器：`iNKT_by_date/2026-09-30/drilldown/html/inkt_network_explorer.html`（`results/html/` 仍为空），见 [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]]。
  2. GOLDEN 论文（补充材料）的 14 个 KG 页面：点击节点后该节点及直接邻居保持原色、其余节点和边淡化，相连边加粗；点击空白处或按 Esc 恢复；右侧面板按证据通道（gene overlap / in the same study with published PMID / shared GO BP）分组列出邻居，可跳转；悬停节点显示名称、区域和度，悬停边显示证据通道与两端节点。颜色、曲率、图例、节点/边数量与 09-29 版一致。邻居/度按当前显示的边计算（GO BP 默认 top 100）。需联网加载 vis-network。新压缩包 `archive/2026-09-30/Supplement_S4_and_KGs_20260930.zip`，说明见 `/Users/zeruzhang/Documents/zhou2_latest_20260929/deliverables_clean_latest/README_for_Yue_20260930.md`，源码目录 `/Users/zeruzhang/Documents/zhou2_latest_20260929/deliverables_clean_latest/project/verify/supplement_top5_kg_20260930_interactive/`。
- 09-25 已有节点/边/模块证据表：`iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`、`iNKT_by_date/2026-09-25/results/tables/module_gene_evidence.csv.gz`。
- 时间：Yue 说最晚 11 月初要完成，下次与 UAB 开会约在下月底（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:17:29–16:18:05、16:31:23–16:31:37）；11 月初对应哪项截止，转写不清，不确定。

## 相关概念

- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 把每个显著 GO 术语当作节点，术语之间共享的基因当作边，得到一张功能关系图。
- [[图中的富集分数：−ln p 与 −log10 q（Enrichment Scores in Figures）|Enrichment-Scores-in-Figures]] — 不同图用了不同的“分数”：旧图 −ln(名义 p)，新 GO 气泡图 −log10(BH q)，GSEA 图用带符号 NES。
- [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] — Yue 的问题：输入了多少基因，图上只展示了多少，还有多少没被任何显著条目覆盖？
- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 09-30 已完成的混合 GO/KEGG 网络分析。
- [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] — 通路网络 + 基因网络的交互页面。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
