# 项目时间线（Project Timeline, 09-05 … 09-30）

每一行：日期 → 做了什么/改了什么 → 为什么 → 相关概念。“来源”给出项目文件；文件不在本地的，只能引用记录它的文档。

## 前期（背景）

| 日期 | 内容 | 来源 |
|---|---|---|
| 07-01/02 | 部署 Scanpy 环境，跑通 iNKT 基础流程：18,204 × 15,741，9 群（宽松 QC 基线，不同于后来的口径） | `docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 1 |
| 08-18/19 | 按旧 PPT 重建 QC：**15,532 × 10,670**，六个样本细胞数与旧 PPT 一致；旧 11 群 vs 新 8 群；旧 tissue DEG 恢复 胸腺 31/34、骨髓 35/39、脾脏 44/50 | 同上 阶段 3；`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` |
| 08-25 | 50/100/200 PCs 的 UMAP 比较（选 50）；参考 marker；按组织分层比例 | 同上 阶段 4 |
| 08-30 | C5 拆成 C5-1/C5-2（五种子平均 ARI=0.9868）；Blood Figure 3D/E/F 方法迁移（火山图、GSEA、Venn/UpSet） | 同上 阶段 5 |

## 09 月

| 日期 | 做了什么 / 改了什么 | 为什么 | 概念 |
|---|---|---|---|
| 09-05 | 会议要求落实：Harmony、单组织重分析、旧 PPT 通路核查、heat-shock 排除；第二轮“用当前 pipeline 找新信息”：骨髓 C4、脾脏 C3 CCT/TriC 等候选 | Yue 会议要求（R01–R13）；RNA velocity 因缺输入未完成 | [[批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）|Harmony-and-Batch-Integration]] [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] [[拟时序（Diffusion Map, DPT, PAGA）与 RNA 速度（Pseudotime vs RNA Velocity）|Pseudotime-DPT-vs-RNA-Velocity]] |
| 09-06 | R01–R13 逐项回应文档（16 页中文 PDF 等） | 每项写明做没做、没做为什么 | [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] |
| 09-09 | 20–30 分钟组会准备；讨论 R01 UMAP 含义 | 汇报 | [[UMAP 降维可视化（UMAP: why position is not physical）|UMAP]] |
| 09-15 | 参考论文格式的 37 页汇报（仓库 README 链接了 2026-09-15 目录，但该文件夹本地不存在，见 [[文档之间看起来不一致或需要核对的地方（Document Inconsistencies and Unsure Points）|Document-Inconsistencies]]） | 对齐 Blood 论文图式 | [[火山图（Volcano Plot）|Volcano-Plot]] [[基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）|GSEA-and-NES]] |
| 09-19 | **细胞毒性 + GO + 拟时序**（19 页）：Prf1/Gzma/Gzmb 核心模块，GO BP/MF/CC 全库 ORA，局部 DPT | 采用 PAGER-scFGA / Borra 2026 的展示逻辑 | [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] [[基因本体（Gene Ontology, GO: BP / MF / CC）|Gene-Ontology]] [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] |
| 09-20 | 可编辑资料包：28 页主报告、58 页整合版、37 页转换版、19 页转换版 | 方便 Yue 直接修改 PPT | 见 `iNKT_by_date/2026-09-20/README.md` |
| 09-22 | 细胞毒性签名选取与外部数据核查（GSE296020 等，仅检索、未下载） | 回应“签名怎么选” | [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] |
| 09-23 | 与 Yue 的会议（转录稿）：GNN 增益如何衡量、super-PAG “大象”比喻 | 方法讨论 | [[PAG、super-PAG 与 m-type 关系（PAGER Terms: PAG, super-PAG, m-type）|PAG-Super-PAG-and-m-type]] [[基线：共同邻居与 Adamic-Adar（Baselines: Common Neighbours, Adamic-Adar）|Baselines-Common-Neighbours-Adamic-Adar]] |
| **09-24** | **与 Rob、Yue（及 Rajesh）开会**：Rob 建议聚焦 MAPK、AP-1、iNKT17；核糖体/应激 = “angry cells”；提供 AML 数据（Box）；Yue 要完整候选、GO→功能组→GOLDEN | 见 [[人物与各自关心的点（Who Is Who）|Who-is-who]] | [[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] [[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] [[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] |
| **09-25** | 会后跟进（PPT 第 20–31 页）：AP-1 复合体 GO:0035976（骨髓 C0、脾脏 C3）；MAPK 不显著；Ccr6 从原始矩阵找回；97 节点/577 边功能网络、G01–G07、GNN | 对应 Rob 与 Yue 的意见（`iNKT_by_date/2026-09-25/results/REQUIREMENTS_CHECKLIST.csv` R01–R13） | [[术语网络（Term Network: nodes = terms, edges = shared genes）|Term-Network]] [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] [[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]] |
| 09-26 | 又一次会议（转录稿 `/Users/zeruzhang/Downloads/Zoom会议/zoom-0927/GMT20260926-205938_Recording.transcript.vtt`），我们只检索过关键词，主要是论文稿措辞（如 self-supervised）；未深读 | — | — |
| **09-29** | 与 Yue 会议（Zoom 转录稿，见文末证据）：要求再用 **GO MF + KEGG** 建网络；节点大小 = 显著性、颜色 = 功能组；基因层面 drill-down；交互 HTML；ribosome/stress 结果“不新颖” | Yue 认为 G01 的应激不是合作者关心的；G04–G06 提供额外信息 | [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] |
| **09-30** | 混合 GO/KEGG 网络分析**完成**：两个版本（a 全部基因 159 → 112 节点 / 231 边；b 主版本 94 → 66 节点 / 157 边）；锚点仅 AP-1 显著；基因层面钻取（AP-1 与 KEGG 的桥 = Fos + Jun 加一个）；交互 HTML；幻灯片第 32–43 页（`iNKT_by_date/2026-09-30/presentation/iNKT_pathway_network_drilldown_20260930.pptx`） | 落实 09-29 要求 | [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] [[锚点节点（Anchor Nodes）|Anchor-Nodes]] [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] [[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]] [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] |
| 09-30（GOLDEN 论文侧） | 补写方法部分“最终模型启用了哪些机制”（协方差正则 λ=1.0、DropEdge、PairNorm，input-skip 未启用；来自 160 组配置初筛）；14 个交互式 KG 页面加点击高亮邻居与悬停提示，新压缩包 `archive/2026-09-30/Supplement_S4_and_KGs_20260930.zip` | 回应 Yue 的 KG 交互要求；见 `/Users/zeruzhang/Documents/zhou2_latest_20260929/deliverables_clean_latest/README_for_Yue_20260930.md` | [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] |

证据：09-24 会议 `/Users/zeruzhang/Downloads/Zoom会议/zoom-9024/GMT20260924-133029_Recording.transcript.vtt`；09-29 会议 `/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt`；09-23 会议 `/Users/zeruzhang/Downloads/Zoom会议/zoom-0923/GMT20260923-180025_Recording.transcript.vtt`。08-05 另有一次会议录音（`/Users/zeruzhang/Downloads/Zoom会议/zoom-08-10/GMT20260805-173214_Recording.transcript.vtt`），本 wiki 仅检索，未细读。
