# 2026-09-30｜GO/KEGG 功能网络与锚点分组（PI Yue 的问题）

接 2026-09-25。复用既有 DE（20 个合格比较）、冻结 GO（BP/MF）与 GOLDEN/GNN 适配代码；新增 KEGG_2019_Mouse。未联系任何人、未上传或下载数据、仅 CPU。

## 做了什么
- ORA：基因 BH FDR≤0.05 且 |log2FC|≥0.25，T2-up / Ctrl-up 分开；GO BP、MF（另算 CC 仅用于 AP-1 锚点）、KEGG；词条大小 5–500；BH 在 比较×方向×来源 内。GO 部分与 09-19 结果逐项一致（312,550 条，q 最大差 1e-16）。
- 两个版本：(a) 全部基因；(b) **主版本**：query 和背景同时剔除，正则 `^(Rp[ls]\d|Rplp\d|Rpsa$|Mrp[ls]\d|mt-|Hsp(?!g)|Dnaj)`（242/10670 基因；保留 Rps6k*激酶和 Hspg2）。
- 节点：q≤0.05 且命中≥3；加 12 个固定锚点；Jaccard≥0.75 合并（锚点不会被并入他人）。边：Jaccard≥0.25 且共享≥3。
- 分组：Louvain（仅重叠）、语义（MiniLM）、GOLDEN 融合（α=0.5）、GraphSAGE+语义（seed 42/43/44）；Ward，k 由轮廓系数在 2–15 选择。

## 关键数字
| | (a) 全部基因 | (b) 主版本 |
|---|---|---|
| 合并前/后节点 | 159 / 112 | 94 / 66 |
| 边（孤立点） | 231（32） | 157（21） |
| 组数 Louvain/语义/GOLDEN/GNN+语义 | 49/15/2/15 | 31/15/3/4 |

- 锚点：仅 AP-1 复合体（GO:0035976）显著（q≈1.7e-4，Fos/Fosl2/Jun/Junb/Jund，骨髓C0、脾C3）。KEGG IL-17 q=0.071，MAPK 0.64，TNF 0.56，其余 GO 锚点 q=1；“interleukin-17 production”在冻结注释中无实测基因，为孤立空节点。
- **锚点是否同组（主版本）**：GOLDEN 与 GNN+语义（按轮廓 k）把全部 AP-1/MAPK/IL-17/TCR 锚点放入同一个大组（G03，36 节点，泛信号/免疫杂项）；语义方法和 Louvain 则分散。这个结论依赖 k：GNN+语义 k≥7 时 MAPK 与 AP-1/IL-17/TCR 分开（见 `anchor_k_sensitivity_*.csv`）。版本 (a) 中 GNN+语义 k=15 时 AP-1、MAPK、IL-17/TCR 三处分离。
- 方法一致性（ARI，a）：Louvain–GNN 0.58、语义–GNN 0.53、GOLDEN–其他≈0.03–0.07；(b) GOLDEN–GNN 0.70，其余 0.10–0.35。GNN 跨种子 ARI：(a) 0.61–0.71，(b) 0.32–0.55（不稳定）。
- 留出边 AUC（3 种子均值）：(a) GNN 0.973 vs 语义 0.840 / 共同邻居 0.922 / Adamic-Adar 0.928；(b) 0.981 vs 0.723 / 0.922 / 0.922。GNN 比图基线高约 0.05，但测试边仅 23（a）/15（b）对，波动大；分组层面 GNN 相对 GOLDEN 融合的增益有限（(b) 中两者 ARI 0.70）。
- 热图显示氧化磷酸化/嘌呤核苷酸合成类组在 BM C0、Spleen C3 的 T2-up 驱动基因最多；BM/脾 C5-2 无显著词条。

## 局限
细胞层面探索性分析；每标签为 3 只小鼠 pooled（每组 n=1 pool）；富集不等于激活或分泌；k 由轮廓系数决定，极端值（k=2、15 边界）使分组偏粗或偏碎；GNN 为小图上的转导式重建，不是生物学验证；Louvain 将每个孤立点当作单独组。

## 偏离规格
锚点在合并时受保护（不与最小 q 规则竞争）；用 Louvain（networkx）代替 Leiden；额外加了 GNN-only 和 k 敏感性表；未做基因层面钻取与 HTML。

## 文件
`code/`（01_ora → 05_validate）、`results/tables/`（ORA 全表、node_attributes_*、anchor_by_method_*、group_summary_*、ARI_*、gnn_vs_baselines_*、partition_quality_*、group_driver_genes_*）、`results/network/{all,excl}/`（nodes/edges/embeddings/GNN）、`results/figures/`（全网络、super-PAG、热图；PNG+SVG）、`results/validation.json`。

## Drill-down（基因层面，2026-09-30 追加）
目录 `drilldown/`（code/、tables/、figures/、html/）；复用既有 ORA 与网络，未重算 ORA/GNN，未新增数据，仅 CPU。
- **背景**：版本 (b) 中锚点同在 G03 不是"锚点相关"的证据（G03 是 36 节点的杂项组，多为孤立点与 2 节点碎片），因此在基因层面回答 AP-1 / MAPK / IL-17·Th17 / TCR / TNF 是否功能相关。
- **钻取对象**：KEGG MAPK、IL-17、TCR、Th17、TNF，GO:0035976 AP-1 复合体，以及版本 (b) GNN+语义 G01（氧化磷酸化）作为"能量"对照（成员=在 G01 的 15 个词条中出现≥3 次的基因，共 181；并集 773 含疾病类词条，过宽）。成员为冻结注释的全部基因（KEGG 用冻结 gmt；GO 用冻结 GAF+is_a/part_of 传播），含未检测基因。
- **表**：`tables/pathway_gene_membership.csv`（每基因：是否检出、6 个比较的 log2FC/FDR/DEG、注释频次、shared/specific）、`pathway_drill_summary.csv`、`anchor_bridge_genes.csv`（锚点两两共享基因）、`bridge_candidate_genes.csv`、`gene_comembership_edges.csv`。
- **结论**：锚点间共享基因很少，且真正被差异表达支持的只有 Fos + Jun（MAPK、IL-17 另有 Jund，TNF 另有 Junb；TCR、Th17 与 AP-1 重叠中的 Nfatc2 仅为注释重叠、非 DEG）。AP-1 与各 KEGG 通路的重叠仅 3 个基因（Fos、Jun 加 Jund/Junb/Nfatc2）；5 个 KEGG 通路两两共享 19–43 个基因，其中 DEG 只有 3–6 个（Fos、Jun、Chuk 等，主要在 BM C0、脾 C3，Fos/Jun 上调而 Chuk/Traf2/Akt2/Map2k7/Grb2 在脾 C3 下调）。共享的多为 MAPK/NF-κB 骨架基因（Mapk1/3/14、Nfkb1、Rela、Ikbkb 等），未见差异表达。Fosb、Dusp1、Nr4a1、Egr1 差异明显，但 Fosb 仅在 IL-17，Dusp1/Nr4a1 仅在 MAPK，不构成桥。能量模块与 AP-1 共享 0 个基因、与各 KEGG 通路仅共享 0–3 个基因且均非 DEG。BM/Spleen C5-2 无任何 DEG。→ 以 Fos/Jun 为核心的"即刻早期基因"信号是唯一有 DE 支持的桥，其余是注释重叠。
- **静态图**：`figures/gene_network_*.png/svg`（每个通路一幅，BM C0 与 Spleen C3 双面板，DEG 红框加星，未检出灰色空心；MAPK 仅 DEG 及其一阶邻居不淡化）、`bridge_genes_heatmap`、`network_full_b_primary_excl_ribo_mito_hsp_v2`（标签避让）。
- **交互 HTML**：`html/inkt_network_explorer.html`（vis-network 从 CDN 加载，需联网；标签页 (i) 通路网络 a/b 切换、组过滤、点击高亮；(ii) 基因网络、比较下拉、通路过滤、点击高亮；在 (i) 点击节点后按钮跳转 (ii) 并按该通路基因过滤）。已用 headless Chrome 截图（`figures/html_screenshot_*.png`）。
- **限制**：未使用 STRING/PPI，边仅为注释共成员（≥2 个显著词条/锚点共有，权重=共享词条数），不是物理相互作用；细胞层面探索性分析，每个组织×条件为 3 只小鼠 pooled（n=1 pool）；富集/注释不等于通路激活。基因 "检出" = 出现在 DE 表（10,670 基因）。
