# 交互式网络浏览器（Pathway Network Explorer HTML）

> 09-30 交付的单文件 HTML：通路网络 + 基因网络，点击高亮邻居，可从通路跳到基因。

## 定义

`iNKT_by_date/2026-09-30/drilldown/html/inkt_network_explorer.html` 是回应 Yue 9/29 “点击 highlight 邻居、展开到基因”要求的交互页面，用 vis-network 绘制（库从 CDN `unpkg.com` 加载，**需要联网**）。

**打开方式**：用浏览器直接打开该 HTML 文件（双击即可，无需服务器）。断网时页面空白，因为图形库没加载。

**内容**：
- 标签页 (i) **通路网络**：版本 a / b 切换、按功能组过滤、点击节点高亮该节点及其邻居；节点大小 = 显著性，颜色 = 功能组。
- 标签页 (ii) **基因网络**：比较下拉（BM C0 等）、通路过滤、点击高亮；DEG 加标记，边为共成员（不是 PPI）。
- 在 (i) 点击节点后，按钮可跳转到 (ii) 并按该通路的基因过滤。

## CS 类比

≈ 一个自带前端的只读仪表盘：数据内嵌在 HTML 里，图形库靠 CDN，所以是“离线数据 + 在线渲染”。

## 常见误读与注意

- 需要联网（CDN）。
- 渲染依赖浏览器；已用 headless Chrome 截图检查（`iNKT_by_date/2026-09-30/drilldown/figures/html_screenshot_terms.png`、`iNKT_by_date/2026-09-30/drilldown/figures/html_screenshot_genes.png`；交互步骤截图在 `iNKT_by_date/2026-09-30/drilldown/html/qa/`）。
- 与 GOLDEN 论文的 14 个 KG 页面（[[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]]）是两套不同的东西。
- 页面中的边是注释共成员，节点大小反映显著性，都不代表激活。

## 在本项目中

- 构建代码：`iNKT_by_date/2026-09-30/drilldown/code/08_html.py`，模板 `iNKT_by_date/2026-09-30/drilldown/code/explorer_template.html`。
- 数据来源：[[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] 与 [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] 的表。

## 相关概念

- [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] — 通路网络的来源与结论。
- [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] — 基因标签页的数据。
- [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] — 另一线：论文补充的 KG。
- [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] — 通路网络的概念基础。
