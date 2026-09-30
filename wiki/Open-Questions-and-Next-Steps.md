# 未决问题与下一步（Open Questions and Next Steps）

## 来自 09-24 会议（Rob / Yue）

| 项目 | 状态 | 说明 / 来源 |
|---|---|---|
| **AML 数据（NK-only 与全骨髓，小鼠+人）** | 未收到 | Rob 愿意通过 Box 分享 FASTQ，没有 BAM（转录稿 01:34:19–01:34:42）。`iNKT_by_date/2026-09-25/results/REQUIREMENTS_CHECKLIST.csv` R12：`external_pending` |
| **FASTQ → RNA velocity** | 未做 | 需要 spliced/unspliced 计数；FASTQ 不保证一定能得到可用结果（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 5 节）。[[拟时序（Diffusion Map, DPT, PAGA）与 RNA 速度（Pseudotime vs RNA Velocity）|Pseudotime-DPT-vs-RNA-Velocity]] |
| **每个 pool 的动物组成与处理记录** | 已知 3 只/标签，其余待提供 | 跨组织是否同一批动物、批次、T2 处理与取样时点、既有 IL-17 实验条件（`iNKT_by_date/2026-09-25/results/README.md` 第 7 节）。[[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]] |
| **台面实验验证** | 未做（R13 `external_pending`） | 优先：骨髓 C0/脾脏 C3 的 AP-1 相关响应；C5-2 的身份、比例、分泌三读出分开；AP-1 与效应功能是否脱钩（`iNKT_by_date/2026-09-25/results/EXPERIMENT_CANDIDATES.md`）。C0/C3 是计算编号，需要先映射到可分选表型。 |
| **T2 的具体模型/时间点** | 未在项目文件中找到 | [[白血病 CML / AML 与 “T2” 条件（Leukemia, CML, AML and T2）|Leukemia-CML-AML-and-T2]] |

## 来自 09-29 会议（Yue）

| 项目 | 说明（来源：09-29 转录稿，ASR 中文，术语有误处以推测标注） |
|---|---|
| **GO MF + KEGG 网络** | 之前只用 GO BP；再用 MF 和 KEGG（ASR 里的 “cag/password” 推测为 KEGG/pathway），q≤0.05 的条目都进图（约 16:19:09–16:19:45；16:21:14–16:21:31；16:24:36–16:25:36）。**已完成（09-30）**：GO BP + MF + KEGG，两个版本（a 全部基因 / b 主版本），[[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] |
| **节点大小 = 显著性（−log10 p）、颜色 = 功能组** | **已完成（09-30）**：网络图与 HTML 均按此着色/缩放。约 16:28:05–16:29:02。[[图中的富集分数：−ln p 与 −log10 q（Enrichment Scores in Figures）|Enrichment-Scores-in-Figures]] |
| **基因层面 drill-down** | **已完成（09-30）**：[[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]]、[[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]。每个 PAG 里有多少基因、哪些与 DE 重叠（星标/红色边框）、每个基因被注释多少次（共享基因 vs 特异功能基因）（约 16:30:27–16:35:52；16:42:27–16:43:16）。 |
| **交互 HTML（点击 highlight 邻居）** | **已完成（09-30）**：网络浏览器 [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]]；GOLDEN 论文的 14 个 KG 页面也已加点击高亮邻居与悬停提示。约 16:13:27–16:14:40；16:40:56–16:41:29。[[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] |
| **super-PAG“大图 + 小图”** | **已完成（09-30）**：`iNKT_by_date/2026-09-30/results/figures/superPAG_panels_b_primary_excl_ribo_mito_hsp.png`。约 16:19:45–16:20:10。[[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] |
| **把 GNN 融合进去并反复推敲、加入验证信息** | **部分完成（09-30）**：已加入 GNN+语义 与 GNN-only、3 种子、基线对比、k 敏感性；结论是增益有限且不稳定，验证信息（外部数据、更大图）仍待做。约 16:31:40–16:31:55。[[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] |
| **与 UAB 的下次会议** | 约下月底；最晚 11 月初需完成某项截止（对应事项转写不清，不确定）（约 16:17:29–16:18:05；16:31:23–16:31:37）。 |
| **重新生成表格并打包发给 Yue** | 表已生成（`iNKT_by_date/2026-09-30/results/tables/`、`iNKT_by_date/2026-09-30/drilldown/tables/`），是否已发送本 wiki 不知道。约 16:14:42–16:15:52；16:45:16–16:46:14。 |

## 其它待办与边界

- 每个 pool 缺独立重复：任何结论保持“探索性”（[[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]]）。
- 外部数据（仅检索，未下载）：GSE296020（小鼠 iNKT1 细胞毒性亚群 bulk）、GSE298293（胸腺 iNKT scRNA 参考）、GSE306154（结肠 iNKT）、Mouse MSigDB 2026.1.Mm（`docs/audits/inkt_cytotoxicity_source_audit_20260922.md`）。
- 肝/肾毒性数据任务不属于 iNKT（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 第 3 节末）。
- 建议的科学问题（不是已证实结论）：AP-1 相关响应与效应/折叠状态是否脱钩？（`iNKT_by_date/2026-09-25/results/README.md` 第 6 节）

## 09-30 之后新增的待办（推断）

- **台面验证 AP-1（Fos–Jun）反应是生态位内激活还是解离应激**：这是 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] 的核心解读，需要在骨髓 C0/脾 C3 对应的可分选表型上做，并记录处理时长/批次。
- 分组结果依赖方法和 k，需要 Yue 决定报告哪个 k、哪种分组（[[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]]）。
- GNN 增益要更有说服力，需要更大的图、共享基因之外的特征和外部验证（[[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]]）。
- GOLDEN 论文 KG 页面：包内 `README_KG.md` 仍是 09-29 版文字，需要更新；邻居/度按当前显示的边计算（GO BP 默认 top 100）（`/Users/zeruzhang/Documents/zhou2_latest_20260929/deliverables_clean_latest/README_for_Yue_20260930.md`）。
- 论文方法待确认：筛选在前驱输入（384/390 维）上进行，是否要写明；Leiden 分辨率网格 8 个还是 9 个值；DropEdge 新旧措辞并存（同上）。
