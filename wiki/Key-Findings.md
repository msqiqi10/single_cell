# 目前的主要发现与证据等级（Key Findings so Far）

证据等级（本 wiki 自定义，便于对话）：

- **A**：计数/QC/旧结果在数据内被精确重现或多层分析一致，但仍是细胞层面；
- **B**：单一分析或依赖所选签名/规则，做过敏感性检查；
- **C**：仅描述性或探索性；
- **N**：阴性/未测。

**所有结果均没有独立动物重复（n=1 pool）；所有 p/FDR 是细胞层面探索性统计。**

| # | 发现 | 数字与来源 | 证据 | 主要 caveat |
|---|---|---|---|---|
| 1 | **AP-1 复合体富集：骨髓 C0 与脾脏 C3**（T2 上调） | GO:0035976；C0 命中 Fos/Fosl2/Jun/Junb，家族 q=0.000261、全局 q=0.00286；C3 命中 Fos/Jun/Junb/Jund，0.000193 / 0.00229；20 次 QC 匹配方向全保留；剔除 151 个 Rpl/Rps/Hsp/Dnaj 基因后仍通过。`iNKT_by_date/2026-09-25/results/tables/AP1_complex_GO_exact_ID.csv`；`iNKT_by_date/2026-09-25/results/README.md` 第 3 节 | A/B | 转录富集≠蛋白活性；可能是解离应激（无处理时长/批次）。[[AP-1 转录因子复合体（AP-1: Fos/Jun family）|AP-1-Fos-Jun]] [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] |
| 2 | **MAPK 不显著** | 骨髓 C0：GO MAPK cascade（GO:0000165）家族 q=1；旧 KEGG MAPK paper 规则约 0.078，robust 规则约 0.627。同上 README；`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节 | N | 不能因 AP-1 显著就说 MAPK 被激活。[[MAPK 级联（MAPK Cascade: ERK, p38, JNK）|MAPK-Cascade]] |
| 3 | **iNKT17 身份特征在 C5-2；IL-17 未检出** | C5-2 的 Rorc/Il23r 特征突出；Ccr6 原始矩阵 92 个细胞（被 min_cells=100 滤掉），保留细胞 60 个、其中 54 个在 C5-2；Il17a 0 个，Il17f 2 个（均不在 C5-2）。`iNKT_by_date/2026-09-25/results/README.md` 第 4 节；`iNKT_by_date/2026-09-25/results/tables/full_feature_marker_status.csv` | B/C | 身份 ≠ 数量 ≠ 分泌；胸腺 C5-2 Ctrl 只有 11 个细胞。[[iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）|iNKT-Subsets]] [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] |
| 4 | **核糖体/应激信号** | Rpl/Rps 反复出现在上调基因，Hspa8/Dnaja1 下调（转录稿 00:46:59–00:47:28）；功能组 G02/G03；旧 PPT 的 IEG/HSP/线粒体驱动模块被高度保留（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 9.4）。 | A | Rob：可能是“angry cells”；Yue：不新颖。[[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] |
| 5 | **OXPHOS 在排除后仍在** | 骨髓 C0、脾脏 C3 的 T2 上调基因富集氧化磷酸化（`iNKT_by_date/2026-09-19/README.md`）；剔除 Rpl/Rps/Hsp/Dnaj 后仍通过（`iNKT_by_date/2026-09-25/results/README.md`）；G04–G06。 | B | 表达富集≠代谢通量；线粒体基因也是 QC 信号。[[线粒体、氧化磷酸化与 ATP（Mitochondria, Oxidative Phosphorylation, ATP: OXPHOS）|Mitochondria-and-OXPHOS]] |
| 6 | **细胞毒性核心模块：骨髓 C4 升高，脾脏 C3 降低** | 骨髓 C4 差 +0.08620（q=0.02513）；脾脏 C3 −0.02297（q=0.001538）；加 Ctla2a 后 q=0.93、0.53；杀伤 GO 最小 q=0.632。`docs/audits/inkt_cytotoxicity_source_audit_20260922.md`；`iNKT_by_date/2026-09-19/README.md` | B | 依赖签名；不是杀伤实验。[[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] |
| 7 | **C5-1 在 T2 脾脏中占比更高** | 脾脏 C5-1：Ctrl 239/3,422（7.0%）→ T2 345/3,371（10.2%）；我们对 `iNKT_by_date/2026-09-19/results/tables/cell_metadata_scores.csv.gz` 计数；ZERU 也在会上指出（转录稿 00:44:36） | C | 捕获比例≠组织内绝对数量；无动物层面误差。[[细胞组成与比例（Cell Composition, Proportions）|Cell-Composition-and-Proportions]] |
| 8 | 其它保留候选：骨髓 C4 的 Tcf7/Itga4/Emb 下调、Ly6a/S100a6/Cdca4 上调；脾脏 C3 的 CCT/TriC 下调；骨髓 C0 剪接（次级） | `docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 8 | C | 探索性；GSEA 不支持整个 spliceosome。[[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] [[DNA 修复与 RNA 剪接（DNA Repair, RNA Splicing: Shld1, Sfpq）|DNA-Repair-and-RNA-Splicing]] |
| 9 | 复现层面：QC 后 15,532 × 10,670 精确重建；旧 tissue DEG 恢复 胸腺 31/34、骨髓 35/39、脾脏 44/50 | `docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 3 | A | 下游簇不是逐位复刻。 |
| 10 | T2 更像 activation/residency 程序变化，而不是整体亚型转换 | Direct TCR activation 全局约 +0.121（骨髓 +0.129、脾脏 +0.131、胸腺 +0.099）；亚型分数全局 T2−Ctrl 很小（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 6.5） | B | signature 依赖。 |
| 11 | 网络：97 个 BP 节点、577 条边、7 个探索性功能组；GNN 内部 AUC 0.982–0.993 | `iNKT_by_date/2026-09-25/results/network_validation.json` | C | 小型内部图重建，不是验证。[[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] |
| 12 | **混合 GO/KEGG 网络：显著核心是能量代谢，信号通路锚点几乎不显著** | 主版本 (b)（剔除 ribo/mito/Hsp/Dnaj，242 基因）94 → 66 节点、157 边；版本 (a) 159 → 112 节点、231 边；能量代谢（OXPHOS / ATP / 呼吸链）在排除后仍是显著核心。锚点里只有 GO AP-1 复合体显著（q≈1.7e-4）；KEGG IL-17 q=0.071，MAPK 0.64，TNF 0.56，其余 1；“interleukin-17 production”无实测基因。`iNKT_by_date/2026-09-30/results/network/network_summary.json`；`iNKT_by_date/2026-09-30/results/tables/anchor_by_method_excl.csv` | B | 细胞层面探索性；富集≠激活。 [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]] [[锚点节点（Anchor Nodes）|Anchor-Nodes]] |
| 13 | **通路层面“锚点是否同组”取决于方法和 k；(b) 中“全在一组”是杂项组** | (b) 中 GOLDEN 与 GNN+语义（按轮廓 k）把锚点放进同一个 36 节点的 G03，那是孤立点与碎片组成的 catch-all 组；k≥7 时 MAPK 与 AP-1/IL-17/TCR 分开；锚点彼此及与能量核心几乎没有共享基因边。`iNKT_by_date/2026-09-30/results/tables/anchor_k_sensitivity_excl.csv` | C | 不能把同组当作功能联系。 |
| 14 | **基因层面：AP-1 与 KEGG 锚点之间唯一有 DE 支持的桥是 Fos + Jun（再加 Jund/Junb/Nfatc2 之一，视通路对而定）** | KEGG 五个通路两两共享 19–43 个基因，其中 DEG 仅 3–6 个；骨架基因 Mapk1/3/14、Nfkb1、Rela、Ikbkb 不是 DEG；Chuk/Traf2/Akt2/Map2k7/Grb2 在 Spleen C3 下调；Fosb/Dusp1/Nr4a1/Egr1 差异明显但各自至多在一个锚点里；能量模块与 AP-1 共享 0 个基因。`iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv` | B | 共成员边不是 PPI。解读（推断）：T2 特征是即刻早期 AP-1（Fos–Jun）反应，可能是生态位内激活或解离应激，需要台面实验。 [[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]] [[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]] |
| 15 | **GNN 在 09-30 图上增益有限且不稳定** | 留出边 AUC：GNN 0.973（a）/0.981（b），共同邻居 0.922 / 0.922，Adamic-Adar 0.928 / 0.922，语义 0.840 / 0.723；测试边仅 23 / 15 对；GNN 跨种子 ARI (a) 0.61–0.71，(b) 0.32–0.55。`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_all.csv`、`iNKT_by_date/2026-09-30/results/tables/gnn_vs_baselines_excl.csv` | C | 转导式重建，非验证。 [[为什么在 97 节点图上 GNN 增益有限（Why GNN Adds Little on a 97-node Graph）|Why-GNN-Adds-Little]] |

## 覆盖：骨髓 C0 与脾脏 C3

骨髓 C0 共 274 个 T2 上调候选，仅 66 个落入显著 BP 条目；脾脏 C3 为 96/275（`iNKT_by_date/2026-09-25/README.md`）。[[候选基因被 GO 覆盖了多少（Candidate Coverage by GO）|Candidate-Coverage-by-GO]]

## 09-30 网络与钻取

详见 [[混合 GO/KEGG 网络分析（Mixed GO/KEGG Network, 2026-09-30）|Mixed-GO-KEGG-Network-0930]]、[[锚点节点（Anchor Nodes）|Anchor-Nodes]]、[[桥接基因（Bridge Genes: Fos, Jun）|Bridge-Genes]]、[[基因层面钻取（Gene-level Drill-down）|Gene-Level-Drilldown]]、[[交互式网络浏览器（Pathway Network Explorer HTML）|Pathway-Network-Explorer-HTML]]。
