# iNKT：旧 PPT 与本地重跑结果逐页内容对照

> 用途：为下一版 PPT 准备内容；本文件只整理证据、共同点、差异和建议表达，不制作幻灯片。
>
> 对照对象：`input/iNKT/iNKT.pptx`（33 页）与本地运行 `output/iNKT_scanpy_tutorial_run`、`output/iNKT_extended_runs/20260818_105349`。
>
> **术语限定：** PPT 中 6 个样本的名称、原始细胞数和 feature 维度与本地输入逐项吻合，强烈支持二者来自同一 intended dataset；但 PPT 没有矩阵哈希，不能证明 byte-for-byte identity。本文的“复现/重现”只指对该数据做本地重分析后恢复旧信号（computational recovery），不是独立动物、独立队列或独立实验中的 biological replication。

## 一、先定下来的总叙事

1. **可见输入元数据高度一致。** 两边都显示 6 个样本、18,458 个细胞、32,285 个 Gene Expression features；PPT 缺少矩阵哈希，所以只能写 “strongly consistent with the same intended input”。
2. **本地重跑保留了更多细胞和基因。** 旧 PPT 为 15,532 cells × 10,670 genes；本地为 18,204 cells × 15,741 genes。
3. **组织主导的结构可跨分析流程稳定恢复。** 两套流程都得到 BM、spleen、thymus 主导的细胞群，以及若干混合群。
4. **旧、新 cluster 编号不能直接对应。** 旧稿使用 `legacy_ppt_c1_c11`，本地使用 `current_leiden_res_0_5_c0_c8`。只能比较组织组成、signature 和 DEG 响应程序，不能写“旧 c1 = 新 c1”或“旧 c1 = 新 c5”。
5. **旧稿的主要 T2 响应程序可在本地流程中大量恢复。** 尤其是 BM 主导群和 spleen 主导群；若按 PPT 标题将旧 logFC 解释为 T2−Ctrl，已恢复基因的符号也高度一致，但旧 workbook 的 group order 仍需原代码确认。
6. **工作假设更偏向 activation/proteostasis/mitochondrial remodeling，而非稳定的 subtype conversion。** 这一判断来自 IEG/HSP/mitochondrial DE、curated ORA 和跨来源 subtype assignment 的综合一致性；不同 gene sets 的原始 `score_genes` delta 量纲并不等价，不能仅凭数值大小排序生物效应。
7. **旧 pathway 中的疾病名称不宜当成疾病机制。** “Prion disease / Measles / Estrogen / Atherosclerosis”等结果主要由共享的 HSP、即时早期和线粒体基因驱动，应该归并为机制模块。
8. **轨迹结论必须降级。** 本地 DPT 主要分开了小型 root c8 与其余细胞，不能作为精细分支发育或 lineage tracing 的证据。
9. **所有 T2 vs Ctrl 推断均为细胞层面的探索性结果。** 每个 tissue × condition 只有一个样本，不能把细胞级 p 值当作独立动物重复。

## 二、比较口径

| 项目 | 旧 PPT | 本地重跑 | 新 PPT 应采用的口径 |
|---|---|---|---|
| 样本 | 6 个：Ctrl/T2 × BM/Spleen/Thymus | 样本名和逐样本维度均吻合 | 写“高度支持同一 intended input”；PPT 无矩阵哈希，不能证明 byte-level identity |
| 原始规模 | 18,458 cells × 32,285 genes | 18,458 cells × 32,285 Gene Expression features；另检测到 3 个 Multiplexing Capture features | 把 3 个 multiplex features 作为技术注释，不混入基因数 |
| 过滤后 | 15,532 cells × 10,670 genes | 18,204 cells × 15,741 genes | 两套流程并列，不把其中一套称为“绝对正确” |
| 细胞过滤 | 页面只写下限；用原矩阵可精确重构为 `200≤n_genes<2500`、`pct_counts_mt<5%` | `n_genes≥200`、`pct_counts_mt<20%`；无 `max_genes` | 旧流程同时剔除高复杂度/潜在 doublet 和 5%–20% mt 细胞，因此明显更严格 |
| 基因过滤 | 全体 18,458 cells 上 `min_cells=100`，恰得 10,670 genes；随后在该矩阵上重算 QC | 细胞 QC 后 `min_cells=3`，得 15,741 genes | 基因阈值和操作顺序共同改变后续 cell-level QC universe |
| 归一化/HVG | 参数不完整 | 10,000 counts/cell、log1p、按 sample 选 3,000 HVGs | 当前参数完整写入方法页 |
| 降维 | UMAP+t-SNE；50/100/200 PCs；`n_neighbors=20` | 50 PCs；邻居图用前 30 PCs；`n_neighbors=15`；UMAP | 不声称本地复刻了旧参数 sweep |
| 聚类 | 11 群，c1–c11 | Leiden resolution 0.5，9 群，c0–c8 | 始终带 `legacy` / `current` 前缀 |
| subtype score | literature 与 in-house，量纲不统一 | 12 套 signatures，用统一 `scanpy.tl.score_genes` | 比较各 signature 内部的空间与条件模式；跨 signature 的效应量须先标准化，不能直接比较原始 score、winner margin 或 delta 大小 |
| DE | Wilcoxon；旧页约 `abs(logFC)≥1.2`、padj≤0.05 | Wilcoxon；摘要采用 FDR≤0.05、`abs(logFC)≥0.25` | 重点比较方向和旧集合 recovery，不比较 DEG 总数或 logFC 绝对大小 |
| pathway | PAGER/Reactome/KEGG 等旧结果；可见 `KEGG_2021_HUMAN` 与 human-labeled Reactome terms | 无离线 Reactome/KEGG GMT；仅对 12 套 curated signatures 做 ORA | 数据为 mouse，旧 human-library orthology/mapping provenance 不明；旧 pathway 只能恢复表格并重新解释，不能说已同库复算 |
| trajectory | t-SNE/cluster graph 加手工组织命名 | PAGA + diffusion map + DPT | 只称 exploratory expression-state ordering |

### 本地当前 cluster 基本信息

| Current cluster | 细胞数 | 主要组织 | 稳定解释 |
|---|---:|---|---|
| c0 | 5,051 | BM 95.1% | BM-dominant response anchor |
| c1 | 175 | BM 53.7%、Spleen 44.0% | 小型 BM/Spleen mixed cluster |
| c2 | 1,563 | Thymus 47.2%、BM 34.0% | mixed，Ctrl 富集 |
| c3 | 1,680 | Spleen 50.7%，其余混合 | subtype source-dependent；literature/Wang 偏 iNKT2，in-house 偏 iNKT17 |
| c4 | 1,470 | 三组织近似混合 | mixed cluster |
| c5 | 1,997 | Thymus 98.3% | Thymus-dominant response anchor |
| c6 | 5,977 | Spleen 92.9% | Spleen-dominant response anchor；三来源均偏 iNKT1 |
| c7 | 174 | 混合 | 小群，T2 富集，部分比较 0 DEG |
| c8 | 117 | Thymus 97.4% | 三来源一致 iNKT2-like；DPT root，但规模很小 |

## 三、逐页对照

### 第 1 页｜空白页

**旧 PPT**

- 无可提取文字、图或数据。

**本地结果**

- 无对应分析页。

**共同点与区别**

- 没有需要复现的分析内容；这是可以直接重构的一页。

**新 PPT 建议**

- 改成封面。
- 建议标题：`iNKT single-cell reanalysis: legacy audit and reproducible local analysis`。
- 副标题：6 samples；18,204 cells；15,741 genes；run `20260818_105349`。

---

### 第 2 页｜材料、样本与 QC

**旧 PPT**

- 6 个样本的原始细胞数：Ctrl_Thymus 2,101、Ctrl_BM 3,998、Ctrl_Spleen 3,829、T2_Thymus 1,639、T2_BM 3,131、T2_Spleen 3,760；总计 18,458。
- 每个样本 32,285 genes。
- 页面写明 `min_genes=200`、`min_cells=100`。
- 过滤后为 15,532 cells × 10,670 genes，细胞保留率 84.1%。
- 旧图例中的过滤后样本数合计也为 15,532。
- 虽然 PPT 未写全，但对原矩阵做以下流程可把 6 个样本的旧细胞数逐个精确重现：先在全体 18,458 cells 上保留 `min_cells=100` 的 genes（正好 10,670），重算 QC，再保留 `200≤n_genes<2500` 且 `pct_counts_mt<5%` 的 cells。

**本地结果**

- 原始样本名、逐样本细胞数和 Gene Expression feature 数逐项吻合；另识别到 3 个 Multiplexing Capture features。
- 当前规则：每个细胞检测基因数 ≥200、线粒体比例 <20%；基因至少在 3 个细胞中检测到。
- 过滤后为 18,204 cells × 15,741 genes，细胞保留率 98.6%。
- 样本细胞数：Ctrl_BM 3,978、Ctrl_Spleen 3,808、Ctrl_Thymus 2,033、T2_BM 3,095、T2_Spleen 3,737、T2_Thymus 1,553。

**共同点**

- 可见元数据强烈支持二者来自同一 intended dataset；但没有旧矩阵哈希，不能彻底排除内容级差异。

**关键区别**

- 本地多保留 2,672 个细胞和 5,071 个基因。
- 基因差异来自：legacy 在所有原始 cells 上先做 `min_cells=100`，直接得到 10,670 genes；current 在细胞 QC 后只要求 `min_cells=3`，保留 15,741 genes。
- 细胞差异来自：legacy 还有 `n_genes<2500` 上限和 `mt<5%` 严格阈值；current 没有 gene-count 上限且使用 `mt<20%`。Legacy 删除 2,926 cells，current 只删除 254 cells。
- 该 legacy 流程重建得到的逐样本细胞数为 Ctrl_BM 3,379、Ctrl_Spleen 3,422、Ctrl_Thymus 1,372、T2_BM 2,731、T2_Spleen 3,371、T2_Thymus 1,257，与 PPT 6/6 完全吻合。由于 PPT 没有原代码，正式措辞应为 `legacy QC inferred by exact count reconstruction`。
- 旧稿中同时出现 “CML” 和文件名/本地元数据 “T2”。新稿必须在研究背景页定义 T2 的确切生物学含义。

**新 PPT 建议**

- 左边显示“高度一致的原始输入元数据”，右边显示“两套过滤输出”。
- 主句：`Input metadata strongly concordant; retained cell/gene universe differs`。
- 不用两张大 violin 占满页面；用 98.6% vs 84.1% retention、两条过滤流程和 6/6 exact reconstruction 更清楚。

---

### 第 3 页｜高变基因与 PCA

**旧 PPT**

- 一幅 top-expression 图，主要为 Malat1、Tmsb4x、Tpt1、Rps24、Rpl13、Rpl19、Actb 等高丰度/housekeeping/ribosomal genes。
- 一幅 HVG mean–dispersion 图。
- 四幅按样本着色的 PCA：PC1/2、PC3/4、PC5/6、PC7/8。
- 标题 `Principle components` 应更正为 `Principal components`。

**本地结果**

- 每细胞归一化到 10,000 counts 后 log1p。
- 以 sample 为 batch 选择 3,000 HVGs。
- 在 HVGs 上计算 50 PCs，并保存 HVG 图与 50-PC explained-variance 图。

**共同点**

- 都使用 HVG → PCA 的标准分析路径。
- 低维空间中均能看到强组织/样本结构。

**关键区别**

- 旧稿未记录 HVG 数量、选择方法和 PCA 精确参数。
- 本地用可复现的 explained-variance 曲线代替旧 PC1–PC8 散点矩阵。
- top-expression 图主要说明测序丰度，不足以支持 subtype 或 condition 结论。

**新 PPT 建议**

- 主图只留 `HVG selection + PCA variance ratio`。
- 参数条写明：`3,000 sample-aware HVGs | 50 PCs`。
- top-expression 图移到附录。

---

### 第 4 页｜UMAP/t-SNE 参数探索

**旧 PPT**

- `n_neighbors=20`。
- 上排 UMAP、下排 t-SNE；比较 50、100、200 PCs。
- 备注认为 “200 PCs is the best”。

**本地结果**

- 计算 50 PCs；邻居图使用前 30 PCs。
- `n_neighbors=15`、UMAP `random_state=0`。
- Leiden 分辨率测试 0.2、0.5、1.0，主分析用 0.5。
- 未重跑 t-SNE，也未做 50/100/200-PC sweep；scVI integration 暂缓。

**共同点**

- 两边都恢复出 thymus 独立区、BM 与 spleen 主区等强组织结构。

**关键区别**

- PC 数、邻居数、过滤后的细胞集合和嵌入算法均不匹配。
- 旧稿没有定量稳定性指标，所以“200 PCs 最佳”不能由本地结果确认。

**新 PPT 建议**

- 用本地 `sample / condition / tissue / cluster` UMAP 四联图。
- 页脚明确当前参数；旧 sweep 若保留，只标为 `legacy, not parameter-matched`。

---

### 第 5 页｜聚类与“轨迹”概览

**旧 PPT**

- 旧 t-SNE 上有 11 个群，并使用 T、B、S、T_u1、TBS、BS_u1、BS_u2、BS_l1、BS_l2、BS_b1、BS_b2 等组织型名称。
- 同页有样本 cluster count、cluster graph 和 cluster-correlation heatmap。

**本地结果**

- Leiden resolution 0.5 得到 9 群：c0–c8，细胞数见前表。
- PAGA 描述 cluster 连通性。
- Root cluster 按 cluster mean `Cd27 − Itgam` 选择为 c8；Il2rb/Klrb1c/Zbtb16 虽被导出检查，但不参与这一步选择。
- 实际 `iroot` 不是 c8 的 medoid 或最高 marker cell，而是代码取 `root_cells[0]`；该细胞来自 Ctrl_Spleen，尽管 c8 的 117 个细胞中有 114 个来自 thymus。
- c8 的 DPT 中位数约 0.072；其余群约 0.945–0.965，显示 pseudotime 被压缩成近似“root vs 其余”的两段结构。

**共同点**

- 两套分析都得到组织主导群和混合群，并显示它们之间的连通关系。

**关键区别**

- 11 群变成 9 群；旧名和当前编号不可一一映射。
- 旧页的“trajectory”主要是 graph/correlation；当前才真正计算 PAGA、diffusion map 和 DPT。
- 即便如此，DPT 仍只是表达相似性排序，不是 lineage tracing。
- 所谓 alternative biomarker-root 分析再次选中 c8 并取同一个 first cell，因此 rho=1.0 不能作为有效的多-root 敏感性验证。

**新 PPT 建议**

- 使用 `current UMAP + PAGA + cluster size/composition`。
- 把轨迹称为 `exploratory maturation/state ordering`，并直接标出 c8 只有 117 个细胞。

---

### 第 6 页｜iNKT1 markers

**旧 PPT**

- iNKT1 score map、按旧 cluster 的 violin、feature maps 和 dotplot。
- 核心基因包括 Xcl1、Nkg7、Klrd1、Klrb1c、Klrk1、Fasl、Gzmb、Cxcr3、Il2rb、Tbx21、Ifng、Gzma 等。
- 页面备注：`Xcl1 is uniquely expressed in case (CML)`。

**本地结果**

- Literature 25/26、in-house 35/35、Wang 2022 46/47 个基因可测。
- 三套来源都把 c0、c1、c2、c4、c5、c6、c7 判为 iNKT1-dominant。
- Literature–in-house Jaccard 0.20；literature–Wang 约 0.42；说明结论有共同核心，但 gene-list 本身并不高度重叠。
- Xcl1 并非 T2 特异：Ctrl 表达比例约 61.0%，T2 约 64.4%；全局描述性 log2FC 仅 +0.155，而且 thymus 方向略为负。

**共同点**

- 广泛而强的 iNKT1/NK-cytotoxic program 稳定存在。

**关键区别**

- 旧稿“Xcl1 仅在 case 表达”被本地数据否定。
- 正确表述应为：Xcl1 在两组均广泛表达，T2 全局仅轻度增加且具有组织依赖性。

**新 PPT 建议**

- 用三来源 cluster-level evidence heatmap 代替几十张 feature plot。
- Xcl1 单独做一个 claim-correction 小框：`widely expressed; modest and tissue-dependent T2 change`。

---

### 第 7 页｜iNKT2 markers

**旧 PPT**

- iNKT2 score map、violin、大量 feature maps 和约 80-gene dotplot。
- 基因同时包含 Il4/Gata3/Zbtb16/Rora/Il17rb 等 iNKT2 markers 与 Birc5/Rrm2/Cdca8/Tyms 等增殖基因。

**本地结果**

- Literature 80/86、in-house 18/18、Wang 37/39 个基因可测。
- 同 subtype 的 gene-list overlap 偏低：literature–in-house Jaccard 0.054；literature–Wang 0.182。
- c8 被三套来源一致判为 iNKT2-like。
- c3 被 literature/Wang 判为 iNKT2，但被 in-house 判为 iNKT17。
- 三套 iNKT2 scores 的全局 T2–Ctrl mean delta 很小：约 +0.014、+0.001、+0.001。

**共同点**

- iNKT2-related program 可检测，c8 是最稳定的对应群。

**关键区别**

- c3 的标签依赖 signature source，不能硬命名为 iNKT2。
- T2 没有显示大幅全局 iNKT2 score 上升。

**新 PPT 建议**

- 结论写成：`c8 robust iNKT2-like; c3 source-dependent/ambiguous`。
- 并列显示 gene-set overlap，提醒读者 subtype 定义不是唯一标准。

---

### 第 8 页｜iNKT17 markers

**旧 PPT**

- iNKT17 score map、violin、feature maps 和 dotplot。
- 包含 Rorc、Tmem176a/b、Il17a、Il17re、Ccr2/Ccr6、Il23r、Il1r1、Selenop、Serpinb1a 等。

**本地结果**

- Literature 34/36、in-house 21/22、Wang 43/46 个基因可测。
- Literature–in-house Jaccard 0.222；literature–Wang 0.351。
- Literature/Wang 没有把任何 current cluster 判为 iNKT17 winner；只有 in-house 把 c3 判为 iNKT17。
- c3 的 literature iNKT17 score 可检测，但低于其 iNKT2 score。
- marker 效应不完全一致：Rorc +0.404，而 Il17a −0.160、Il17f −0.924；Il17a/f 表达细胞很少。

**共同点**

- Rorc/IL-17-related signal 仍可检测，且集中在有限区域。

**关键区别**

- 本地不支持一个跨 signature 来源稳定、独立的 iNKT17 cluster。
- 转录因子 Rorc 与低丰度效应因子 Il17a/f 的方向不一致，不宜用单一 marker 定义整群。

**新 PPT 建议**

- 使用 `iNKT17-like program/gradient`，不使用确定性硬分群。

---

### 第 9 页｜Literature 与 in-house signatures

**旧 PPT**

- 上排为 literature iNKT1/2/17，下排为 in-house iNKT1/2/17，并配 6 幅 violin。
- 不同 gene set 的 score 尺度明显不同，不能直接横比绝对值。

**本地结果**

- 在相同 log-normalized raw expression 上，用 `scanpy.tl.score_genes` 统一计算 12 套 signatures。
- 三来源 subtype list 的 Jaccard 已定量计算。
- 全局 T2–Ctrl mean-score delta：literature iNKT1 −0.013、iNKT2 +0.014、iNKT17 +0.015；Direct TCR activation +0.109；residency +0.029；circulatory −0.002。
- Direct TCR 在 BM、spleen、thymus 都增加约 +0.118 至 +0.133。

**共同点**

- 不同来源在大体空间模式上存在一致性。

**新发现/区别**

- 多来源 subtype winners 没有显示稳定的全局 iNKT1 → iNKT2/17 转换；IEG DE 与 curated ORA 则反复指向 activation/residency，因而后者是更一致的工作假设。
- Direct TCR 的 raw score delta 数值最大，但 `score_genes` 会受 gene-set size、control genes 和各集合分布影响；不同 signatures 的 raw delta 不能直接比较幅度。正式比较需先做 per-signature standardization 或统一效应量。
- Measured Direct TCR 与 Residency sets 共享 7 个 IEG genes，两者的 score/ORA 证据部分重复，不能当作两条完全独立验证。

**新 PPT 建议**

- 一页组合：signature Jaccard、cluster subtype evidence、T2 effect heatmap。
- 推荐结论标题：`Activation-related signals recur across complementary analyses`。

---

### 第 10 页｜全局 T2 vs Ctrl DEGs

**旧 PPT**

- Wilcoxon；页面写 `logFC >1.2`、padj≤0.05，但表中有负值，实际应理解为约 `|logFC|≥1.2`。
- 展示 39 genes。
- 上调代表：Tmsb10、Cd52、Junb、Jun、Fos、Hspa1b、Iglc2。
- 下调代表：Slc15a2、Mxd4、Dynll1、Sfpq、Dnaja1、Hspa8。

**本地结果**

- Base global DE 在 3,000 HVGs 上做 Wilcoxon；Ctrl 9,819 cells，T2 8,385 cells。
- 252 个 HVGs 为 FDR≤0.05；其中 91 个同时满足 `|logFC|≥0.25`，3 个满足 `|logFC|≥1.2`。
- 仅按 FDR≤0.05，与旧 39-gene list 共有 15 个：Actb、Eps8l1、Fos、Fosb、Hspa1b、Iglc2、Jun、Junb、Mxd4、Slc15a2、mt-Co3、mt-Cytb、mt-Nd1、mt-Nd2、mt-Nd3。
- 同时要求当前 `|logFC|≥0.25` 时，共有 11/39；要求 `|logFC|≥1.2` 时，共有 3/39：Iglc2、Eps8l1、Slc15a2。

**共同点**

- Immediate-early/stress、线粒体基因以及 Iglc2/Eps8l1/Slc15a2 信号重复出现。

**关键区别**

- 当前 global DE 受 3,000-HVG 背景限制；旧表看起来来自更广基因背景。
- 过滤细胞、基因背景和 logFC 尺度均不同，不能直接比较效应绝对值或显著基因数量。

**新 PPT 建议**

- 用 `15/39 shared` 的 overlap 卡片和 6–10 个代表基因替代 39 行表。
- 主句：`Qualitative program concordance, not numeric replication`。

---

### 第 11 页｜按组织的 T2 vs Ctrl DEGs

**旧 PPT**

- 可见列表：Thymus 34、BM 39、Spleen 50 genes。
- 三组织共享核心包含 Tmsb10、Cd52、Rplp2、Gm10076、Atp5e、Cirbp、Junb、Fos、Hspa1b、Hspa8、Dnaja1、Sfpq、Dynll1 等。

**本地结果**

- 相同阈值口径 FDR≤0.05、`|logFC|≥0.25`：BM 711（509 up/202 down）；Spleen 777（431/346）；Thymus 1,638（1,596/42）。
- 对旧 PPT 可见列表的方向无关 recovery：BM 36/39（92.3%）；Spleen 45/50（90.0%）；Thymus 29/34（85.3%）。
- 现有 `legacy_vs_current_tissue_gene_overlap.csv` 给出的 BM 65、Spleen 54 并不是更宽的 tissue-only 表；生成逻辑按 tissue 过滤但未限制 `legacy_scope`，把同组织的 legacy cluster×tissue rows 也并入了 union。逐页比较必须使用严格 `legacy_scope=legacy_tissue` 的 36/39、45/50、29/34，不得把 65/54 当成 PPT tissue 分母。
- 若在当前结果上也使用 `|logFC|≥1.2`，共享数急剧降低；这说明两套 logFC 的尺度/预处理不能直接等价。

**共同点**

- 旧组织 DE 的核心基因在本地大多仍显著，是最强的跨流程复现证据之一。

**关键区别**

- 当前阈值更宽，显著集合膨胀到 711–1,638，导致 Jaccard 低但 legacy coverage 高。
- Thymus 的 1,596 个 T2-up genes 尤其需要谨慎；可能同时受阈值、归一化、cluster composition 和 cell-level pseudoreplication 影响。
- Current global DE 只测试 3,000 HVGs，而 extended tissue/cluster DE 测试 15,741 genes。两者的 DEG 总数不能直接并列；在补做 full-gene global DE 前，也不能把 “tissue DEG but absent from global” 全部解释为真正 tissue-specific。

**新 PPT 建议**

- 使用 recovery bar 或 UpSet，不再堆三张长表和 Venn。
- 标题：`Legacy tissue DEG core is largely recovered; current list size is threshold-dependent`。

---

### 第 12 页｜Cluster 与组织组成

**旧 PPT**

- c1–c11 t-SNE、cluster graph、相关热图和 11 个样本组成饼图。
- 明显组织优势群：legacy c1/c11 偏 thymus，c5 偏 BM，c6/c8 偏 spleen；其余为混合群。

**本地结果**

- c0 95.1% BM；c5 98.3% thymus；c6 92.9% spleen；c8 97.4% thymus。
- c1/c2/c3/c4/c7 为不同程度混合群。
- 描述性 condition shift：c2 为 Ctrl 72.5%；c7 为 T2 69.0%。

**共同点**

- BM、spleen、thymus 主导的大群在两套流程中稳定出现，是最可靠的结构性共同点。

**关键区别**

- 群数、边界和编号改变，旧 c1–c11 与 current c0–c8 不能按数字匹配。

**新 PPT 建议**

- 用 `cluster × sample/tissue` 100% stacked bar 或 heatmap 代替 11 个小饼图。
- 每群标 N，并圈出 current c0/c5/c6 三个组织主导 response anchors。

---

### 第 13 页｜Legacy c1 thymus DE

**旧 PPT**

- Legacy c1 thymus：30 DEGs，23 up、7 down；12 个标记为重点。
- 核心模式：Tmsb10/Cd52/Hspa1b 上调；Ccl5/Cxcr6/Dnaja1/Hspa8/Sfpq/Dynll1 下调。

**本地结果**

- 按“相同组织主导 + 最大旧集合 recovery”比较，最佳响应匹配为 current c5 thymus：999 T2、964 Ctrl，190 DEGs。
- 旧 30 genes 中 23 个在 current c5 重现；按 PPT 的 `T2 vs Ctrl` 标题解释旧符号时，23/23 方向一致。
- 代表性旧→新 logFC：Tmsb10 `+2.5→+0.45`；Hspa1b `+3.0→+0.70`；Hspa8 `−1.5→−0.37`；Dnaja1 `−1.6→−0.47`；Ccl5 `−1.8→−0.52`。

**共同点**

- Thymus 主导群的 T2 response program 高度保留；在上述待确认的 legacy sign interpretation 下，已恢复基因方向一致。

**关键区别**

- 这是“响应程序匹配”，不是 cluster identity 映射。
- Current c1 thymus 只有 1 个 T2 与 3 个 Ctrl 细胞，不能进行同编号验证。
- 旧表中的 Fos/Dusp1 在整个 thymus 中仍显著，但在 current c5 内被分层稀释；mt-Co3 在整个 thymus 中方向相反，提示少数基因不稳定。

**新 PPT 建议**

- 标题：`Legacy c1 thymus response is retained in current thymus-dominant c5 (23/30 recovered; signs agree under the PPT T2-vs-Ctrl interpretation)`。
- 图中明确画虚线“program match”，不要画等号。

---

### 第 14 页｜Legacy c1 thymus pathways

**旧 PPT**

- 38 个显著 terms，但只由 8 个独特 driver genes 支撑。
- 主要标签：HSP70/HSP40/HSF1、MAPK/NFκB、ER protein processing，以及 prion/measles/atherosclerosis 等疾病名。

**本地结果**

- Current c5 保留 6/8 drivers：Dnaja1、Hspa1a、Hspa1b、Hspa8、Cox8a、Ccl5；Fos/Dusp1 在整个 thymus 中显著。
- 当前没有同版本 Reactome/KEGG GMT，所以没有正式复算旧 term。

**共同点**

- HSP/chaperone、即时早期和线粒体相关 drivers 大量保留。

**关键区别**

- 诱导型 Hspa1a/b 上调，而 constitutive Hspa8/Dnaja1 下调；这不是简单的“heat-shock pathway 全面激活”。
- 疾病名称来自共享 genes，不能作为 prion、measles 等疾病机制结论。
- Legacy 页面使用 human-labeled KEGG/Reactome libraries 处理 mouse genes，但未记录 orthology/gene-symbol mapping；这进一步限制了 term-level 解释。

**新 PPT 建议**

- 压缩为 `proteostasis/chaperone remodeling`，再分 up/down 两支。
- 旧 disease terms 只放附录，并写 `database label driven by shared HSP/mitochondrial genes`。

---

### 第 15 页｜Legacy c2：BM 与 spleen DE

**旧 PPT**

- BM：18 DEGs（15 up/3 down）。
- Spleen：27 DEGs（23 up/4 down）。
- 两组织都有 Tmsb10/Cd52/线粒体基因上调、Hspa8/Dnaja1 下调；spleen 另有 Jun 上调、Hsp90aa1 下调。

**本地结果**

- BM 最佳响应匹配为 current c0：旧 18 genes 中 15 个重现；按 PPT 标题解释旧符号时方向一致。
- Spleen 最佳响应匹配为 current c4：旧 27 genes 中 23 个重现；按 PPT 标题解释旧符号时方向一致。
- Spleen 代表性旧→新：Jun `+2.9→+1.60`；Hspa8 `−3.1→−0.68`；Dnaja1 `−2.9→−0.85`；Hsp90aa1 `−2.0→−0.73`。

**共同点**

- 两个组织中的 activation/chaperone response 都高度保留。

**关键区别**

- 同编号 current c2 在 BM 仅 2 个显著 DEG，在 spleen 仅 13 个，且与旧表几乎无重叠；这是“不允许按编号匹配”的直接证据。

**新 PPT 建议**

- BM/Spleen 两栏，各显示 old-set recovery、cell N 和 4–6 个核心基因；不再显示整张长表。

---

### 第 16 页｜Legacy c2：BM 与 spleen pathways

**旧 PPT**

- BM：4 terms、4 个独特 drivers，集中于 HSP70/HSP40 与“prion disease”。
- Spleen：13 terms、7 drivers，集中于 HSP90/HSP70、线粒体和 Jun。
- 原页面板标签不够清楚。

**本地结果**

- BM/current c0 保留 3/4 drivers；Ndufa4 logFC +0.239，仅因略低于 +0.25 阈值未入选。
- Spleen/current c4 保留 6/7 drivers。
- 当前 curated ORA：c0 BM 的 Direct TCR activation FDR `1.45×10⁻7`、Residency `1.24×10⁻6`；c4 spleen 的 Residency `0.00144`、Direct TCR activation `0.00215`。
- Measured Direct TCR set（11 genes）与 Residency set（37 genes）共享 7 个 genes：Fos、Icos、Jun、Junb、Ppp1r15a、Prdx6、Ptp4a1；两者共显著不是两条独立机制证据。

**共同点**

- HSP/线粒体 drivers 与 activation/residency 方向均支持组织依赖的 T2 response。

**关键区别**

- 旧 pathway 和当前 ORA 使用不同 gene-set library，不能称为同一 pathway 的直接复算。
- Legacy library 带 human 标签而输入为 mouse；orthology/mapping provenance 未记录。

**新 PPT 建议**

- 明确分开 `legacy driver retention` 与 `current curated ORA` 两个证据层。

---

### 第 17 页｜Legacy c5 BM DE

**旧 PPT**

- 61 DEGs（50 up/11 down）。
- T2 中线粒体/核糖体、Jun/Junb/Fos/Dusp1/Hspa1a/b 上调；Hspa8/Dnaja1/Hsph1/Mxd4/Sfpq/Dynll1 下调。

**本地结果**

- 最佳响应匹配为 current c0 BM：53/61 重现；按 PPT 标题解释旧符号时 53/53 方向一致。
- 代表性旧→新：Hspa1b `+2.1→+1.42`；Dusp1 `+1.8→+0.78`；Fos `+1.7→+0.88`；Hspa8 `−1.9→−0.48`。
- 同编号 current c5 在 BM 只有 21 T2/10 Ctrl，因细胞不足跳过。

**共同点**

- BM activation/proteostasis program 是全稿最强复现信号之一。

**关键区别**

- Current c0 汇入了多个旧 BM response programs，提示聚类合并/边界变化，而不是旧 c5 身份原样保留。

**新 PPT 建议**

- 主标题：`BM activation/proteostasis program strongly recovered across analyses (53/61)`。

---

### 第 18 页｜Legacy c5 BM pathways

**旧 PPT**

- 32 terms，只由 15 个独特 drivers 支撑。
- 首项“Prion disease”主要由 COX/NDUFA/HSP genes 驱动；其余为 HSP70、MAPK、ER processing。

**本地结果**

- Current c0 中 15/15 pathway drivers 全部显著。
- 当前 curated ORA 同时命中 Direct TCR activation 与 Residency；两套 measured genes 共享 7 个 IEG drivers，因此不是两条完全独立证据。

**共同点**

- 这是 BM pathway-driver 层面的最强复现证据。

**关键区别**

- “Prion disease”不是 prion-specific signal，而是 respiratory-chain/proteostasis genes 的数据库共享标签。
- Legacy human-library 到 mouse genes 的映射未记录，term-level 结论需进一步降级。

**新 PPT 建议**

- 改写成三个机制模块：
  1. Mitochondrial respiratory-chain remodeling
  2. HSP/proteostasis remodeling
  3. Immediate-early/TCR activation

---

### 第 19 页｜Legacy c6 spleen（旧称 iNKT2）DE

**旧 PPT**

- 29 DEGs（24 up/5 down）。
- Jun/Junb/Ndufv3/Tmsb10 上调；Hspa8/Dnaja1/Hsp90ab1/Dynll1 下调。

**本地结果**

- Current c3 spleen：506 T2/346 Ctrl，97 DEGs；旧 29 genes 中 23 个重现，按 PPT 标题解释旧符号时方向一致。
- Current c6 也重现 22/29，但三套 subtype signatures 均把 current c6 判为 iNKT1。
- Current c3 的 literature/Wang 偏 iNKT2，in-house 偏 iNKT17。
- 代表性旧→新：Jun `+2.9→+1.27`；Ndufv3 `+1.8→+0.63`；Hspa8 `−2.1→−0.49`；Dnaja1 `−2.3→−0.71`。

**共同点**

- Spleen 中的即时早期上调与 constitutive chaperone 下调模式高度保留。

**关键区别**

- 旧的“c6 (iNKT2)”不能转写成 current c6；最佳 response match 是 current c3，且 subtype evidence 仍有来源依赖。

**新 PPT 建议**

- 标题：`Legacy c6 response aligns best with current c3; subtype assignment remains provisional`。

---

### 第 20 页｜Legacy c6 spleen pathways

**旧 PPT**

- 13 terms、6 个独特 drivers：DNAJA1/HSPA8/HSP90AB1 为核心，另含 JUN/COX8A/NDUFV3。

**本地结果**

- Current c3 保留 5/6 drivers；Cox8a +0.218，仅略低于当前效应阈值。
- Current c3 spleen curated ORA：Direct TCR activation up，FDR `6.91×10⁻7`；iNKT1 down，`4.23×10⁻6`；Residency up，`6.43×10⁻5`。Direct TCR/Residency 的 shared IEG drivers 使这两个结果部分重复。

**共同点**

- 旧 driver-level response 与当前 activation/residency enrichment 相互支持。

**关键区别**

- 方向是 `chaperone down + immediate-early activation up` 的双向重塑，不应概括成单向 heat-shock activation。
- Legacy pathway 使用 human-labeled library，且未记录 mouse orthology mapping；当前只验证 driver-level recovery 与独立 curated ORA。

**新 PPT 建议**

- 用双向箭头模块图取代冗长 pathway table。

---

### 第 21 页｜Legacy c7 BM（旧称 iNKT17）DE

**旧 PPT**

- 只有 4 个 DEGs：Cox7c +2.4、Gm10076 +2.7、Rplp2 +1.6、Sfpq −2.6。

**本地结果**

- Current c3 BM 完整重现 4/4：约 +0.66、+0.69、+0.35、−0.63。
- Current c3 BM 共有 84 DEGs；current c7 BM 为 0 DEG。
- Subtype evidence 对 current c3 不一致；只有 in-house 支持 iNKT17，literature/Wang 更偏 iNKT2。

**共同点**

- 旧四基因均被恢复；按 PPT 标题解释旧符号时四者方向一致。

**关键区别**

- 四个通用线粒体/核糖体/RNA-processing genes 不足以证明 iNKT17-specific response。

**新 PPT 建议**

- 表述：`Sparse legacy response recovered across analyses, but iNKT17 identity is not supported consistently across signatures`。

---

### 第 22 页｜Legacy c8 spleen DE

**旧 PPT**

- 84 DEGs（72 up/12 down），是旧稿最大的 cluster × tissue responses 之一。

**本地结果**

- Current c6 spleen：2,687 T2/2,865 Ctrl，657 DEGs；旧 84 genes 中 73 个重现，按 PPT 标题解释旧符号时 73/73 方向一致。
- 同编号 current c8 在 spleen 只有 0 T2/3 Ctrl，无法计算。
- Current c6 的三套 subtype signatures 均指向 iNKT1。

**共同点**

- 这是 spleen 最强的跨分析 response-program recovery。

**关键区别**

- 这是 legacy c8 → current c6 的功能响应类比，不代表 cluster identity 保留。

**新 PPT 建议**

- 主标题：`Dominant spleen response program recovered in current c6 (73/84; signs agree under the PPT T2-vs-Ctrl interpretation)`。

---

### 第 23 页｜Legacy c8 spleen pathways

**旧 PPT**

- 15 terms、20 个独特 drivers；包含“Prion disease”、ER processing、MAPK、HSP70/HSP90。

**本地结果**

- Current c6 中 20/20 drivers 全部显著。
- 典型模式：Dnaja1/Hspa8/Hsp90ab1/Hsph1 下调；Fos/Jun/Dusp1/Ppp1r15a 上调；多个 COX/NDUFA/SEC61 genes 改变。
- Current curated ORA：Direct TCR activation 与 Residency 均 FDR `5.19×10⁻8`；iNKT1 up FDR `0.032`。
- Direct TCR 与 Residency 的共显著部分由 7 个共享 IEG drivers 推动，不能把它们当作两次独立复制。

**共同点**

- Spleen pathway-driver 复现极强。

**关键区别**

- Current c6 是 iNKT1-like；旧 disease labels 仍需降级为共享基因标签。
- 旧 human-labeled enrichment 未记录 mouse orthology mapping，不能视为同物种、同库复算。

**新 PPT 建议**

- 用 activation、proteostasis、mitochondrial 三模块替代 15-row pathway table。

---

### 第 24 页｜Legacy c9：BM 与 spleen DE

**旧 PPT**

- BM 表：11 DEGs（9 up/2 down）。
- Spleen 表：28 DEGs（23 up/5 down）。

**本地结果**

- BM/current c0：9/11 重现；按 PPT 标题解释旧符号时 9/9 方向一致。
- Spleen/current c6：19/28 重现；按 PPT 标题解释旧符号时 19/19 方向一致。
- Current 已无 c9。

**共同点**

- 两组织都保留 Tmsb10/Cd52 上调和 Dnaja1/Hspa8 下调；spleen 还保留 Jun/Junb/Hsp90ab1/Sfpq。

**关键区别**

- Legacy c9 的响应被分别吸收到 current BM c0 和 spleen c6，支持“聚类合并/边界重划”的解释。

**新 PPT 建议**

- 不再保留 c9 作为当前实体；并入 `legacy clusters consolidated into current tissue anchors` 总结页。

---

### 第 25 页｜Legacy c9 pathways

**旧 PPT**

- 页面标题只写 BM pathways，但实际上上方 3-term 表对应 BM，下方 12-term 表按第 24 页布局应对应 spleen；下方缺少标签。
- BM 仅由 DNAJA1/HSPA8 驱动；spleen 由 DNAJA1/HSPA8/HSP90AB1/JUN 驱动。

**本地结果**

- BM/current c0：2/2 drivers 显著。
- Spleen/current c6：4/4 drivers 显著。

**共同点**

- 两个组织的核心 driver axis 复现。

**关键区别**

- 旧页标签错误；BM 只有两个 drivers，不足以被展开成多个独立机制。
- 旧 pathway 仍受 human-library/mouse-input mapping provenance 缺失限制。

**新 PPT 建议**

- 如保留，明确标 BM/Spleen；更推荐并入机制总览，并命名为 `DNAJA1–HSPA8 chaperone axis perturbed`。

---

### 第 26 页｜Legacy c10/c11：无 DEG

**旧 PPT**

- 标题：`T2 vs. Ctrl DEGs in c10 (No DEGs) and c11 (No DEGs)`。
- 同时展示 c3、c4、c10、c11 的 condition/tissue t-SNE 子图。
- 没有给细胞数、阈值或统计表。

**本地结果**

- 当前也存在“完成检验但 0 个显著 DEG”的组合：BM c1（40 T2/54 Ctrl）、BM c7（45/28）、Spleen c7（54/20）、Thymus c8（53/61）。
- 另有若干 `<20 cells/group` 的组合被跳过；“未检验”不能写成“无 DEG”。

**共同点**

- T2 response 明显依赖 cluster；部分群即使两组均达到最低细胞数，也没有可检测的显著变化。

**关键区别**

- 旧 c10/c11 不存在于当前命名空间，不能指定某个 current cluster 为其直接对应物。
- 0 DEG 既可能是低响应，也可能来自 power、dropout 或多重校正。

**新 PPT 建议**

- 标题：`Not all iNKT subsets show a detectable T2 transcriptional response`。
- 用状态表区分：`tested + significant`、`tested + 0 DEG`、`not tested: insufficient cells`。

---

### 第 27 页｜Tissue DEG 与 cluster-specific DEG overlap

**旧 PPT**

- 左侧 tissue-level Venn：Thymus 34、BM 39、Spleen 50；三组织共有 19。
- 左侧精确分区：Thymus-only 10、BM-only 4、Spleen-only 16、Thymus∩BM only 3、Thymus∩Spleen only 2、BM∩Spleen only 13、三者共有 19。
- 右侧 cluster Venn：legacy c1_Thymus 30、c5_BM 41、c8_Spleen 50；三者共有 17。
- 右侧精确分区：c1-only 8、c5-only 6、c8-only 16、c1∩c5 only 3、c1∩c8 only 2、c5∩c8 only 15、三者共有 17。
- 星号定义：tissue-specific = 不在 global DEG；cluster-specific = 不在对应 tissue DEG。
- 注意：本页 c5/c8 的 41/50 与第 17/22 页完整表的 61/84 不同。Legacy c5-BM 的 61 行中，按精确 `abs(logFC)>1.2` 正好得到 41（已确认）；c1 的 30 行均超过该阈值。c8 的 50 推测为同类严格 effect subset，但 PPT 显示值经过四舍五入，缺少原始精度，无法精确重建其 50-gene rule（待确认）。新稿应把 full table 与 strict Venn subset 分开标注。

**本地结果**

- Tissue DE：BM 711、Spleen 777、Thymus 1,638。
- 旧 tissue list recovery：BM 36/39（92.3%）、Spleen 45/50（90.0%）、Thymus 29/34（85.3%）。
- Tissue vs current tissue-dominant cluster：
  - BM c0：up Jaccard 0.489，down 0.458。
  - Spleen c6：up 0.526，down 0.544。
  - Thymus c5：up 0.064，down 0.147；但 c5 up 的 overlap coefficient 0.963，低 Jaccard 主要因为 thymus tissue list 极大。

**共同点**

- 多数 cluster response 同时存在于对应组织整体中，但仍有 cluster-only genes。
- 旧 tissue DEG 的 gene-level recovery 很高。

**关键区别**

- 当前集合远大于旧集合，Venn 不再易读；Jaccard 低不等于旧核心消失。
- Thymus 1,596 个 T2-up genes 使 union 极大，必须和 overlap coefficient/legacy coverage 一起看。
- 旧页“tissue-specific = absent from global”的定义不能直接套到当前结果：current global 仅测试 3,000 HVGs，而 tissue/cluster DE 测试 15,741 genes；需要先补做同一 gene universe 的 full-gene global DE。

**新 PPT 建议**

- 用当前 UpSet/Jaccard 图替代 Venn；旁边放 92.3%/90.0%/85.3% recovery 卡片。
- 标题：`Most legacy tissue signals reproduce, while cluster context adds distinct DEGs`。

---

### 第 28 页｜Legacy c1/c5/c8 cluster-associated DEG 验证

**旧 PPT**

- 三组 violin：c1_Thymus、c5_BM、c8_Spleen。
- 标出 cluster-only genes，例如 c1 的 Klf6/Cxcr6/Rhob，c5 的 Tomm6/Pink1/Fus，c8 的 mt-Co2/Atp5j2/Actg1/Xist/Xcl1/Snrnf/Ccnd2。

**本地结果**

- 按组织主导群做功能类比，而非编号对应：
  - Legacy c1-Thymus → current c5：23/30（76.7%）。
  - Legacy c5-BM 的本页 41-gene subset → current c0：37/41（90.2%）。
  - Legacy c8-Spleen 的本页 50-gene subset → current c6：43/50（86.0%）。
- 当前组织主导群：c0 95.1% BM、c5 98.3% thymus、c6 92.9% spleen。

**共同点**

- 旧 cluster-associated gene sets 大部分能在当前相同组织主导群中找到。

**关键区别**

- 不少复现项为 housekeeping、stress、ribosomal 或 mitochondrial genes，统计稳定不等于亚型特异。
- 当前 cluster-only gene 数更多：BM c0 142、Spleen c6 157、Thymus c5 71（up+down）。

**新 PPT 建议**

- 左侧用 old-set recovery bars；右侧显示 current cluster-only gene counts。
- 标题：`Legacy cluster-associated signatures are recoverable, but not cluster-ID equivalent`。

---

### 第 29 页｜DEG interaction network

**旧 PPT**

- Legacy c1-Thymus、c5-BM、c8-Spleen 的 DEG interaction network。
- 中心常见 Jun、Actb、mt-Cytb、Cox6b1、Ndufa4、Ndufa13、Cox6c、Cox7c、Atp5e、Hspa1b。
- 页面未写 interaction database、物种版本、confidence cutoff 或背景集。

**本地结果**

- 当前没有 DEG/PPI interaction network 输出。
- 现有的是 UpSet、Jaccard、signature overlap 和 PAGA；PAGA 是细胞状态图，不是 gene interaction network，不能替代本页。

**共同点**

- 旧网络中的多数 DEG 在当前组织主导群仍能找到；stress/mitochondrial genes 仍是高频核心。

**关键区别**

- 网络本身没有复现，且旧网络 provenance 不完整，degree 和 edge score 无法审计。
- Housekeeping/mitochondrial genes 可能因广泛连接而人为占据网络中心。

**新 PPT 建议**

- 正式主稿删除旧 PPI 网络，改用可审计的 UpSet + module summary。
- 若后续重建，必须写明 Mus musculus 数据库及版本、score cutoff、背景集，并单独评估移除 mt/ribosomal/泛 stress genes 后的稳定性。

---

### 第 30 页｜Clusters and underlying trajectories

**旧 PPT**

- 按 tissue composition 把 t-SNE 区域命名为 Thymus、BM、Spleen、BM+spleen、upr-mixture、btm-mixture。
- 代表性组成：Thymus 区约 Ctrl 51.2%/T2 46.7%；BM 区 51.8%/43.9%；Spleen 区 43.2%/40.3%。
- BM+spleen 区约 51.1% BM、45.6% spleen；upr-mixture 的 6 个样本各约 9.5%–23.8%；btm-mixture 约 61.5% spleen、31.6% BM。
- 右侧网络和手工箭头被解释为 trajectory，但未报告 pseudotime、root 或 lineage 方法。

**本地结果**

- Tissue topology 稳定复现：c0 BM、c5 thymus、c6 spleen；c1/c2/c4 为不同类型 mixed clusters。
- PAGA 最强 edges：c4–c5 0.498、c0–c1 0.275、c3–c8 0.260、c2–c8 0.224。
- Root cluster 按 `mean(Cd27) − mean(Itgam)` 选为 c8：117 cells，其中 114 thymus、3 spleen。
- 实际 `iroot` 是 c8 列表中的第一个细胞 `Ctrl_Spleen_ACGAT...`，不是 thymus medoid/最高 root-score cell；alternative run 又使用同一 cell，所以 rho=1.0 只是同 root 重算。
- c8 median DPT 0.072；其余 clusters 0.945–0.965，说明全局排序高度压缩。
- Program–DPT correlations：iNKT1 Wang 0.503、iNKT1 literature 0.496、stage3 cytotoxic mature 0.461、residency 0.301、Direct TCR activation 0.299；stage0 immature −0.075、in-house iNKT2 −0.112。
- Tissue 内 Ctrl/T2 的平均 DPT 差异接近 0。

**共同点**

- 组织来源是低维结构的主轴；c8 与成熟/cytotoxic programs 的相对分离可见，但当前 root-cell 选取不足以把它命名为确定的 thymus developmental root。

**关键区别**

- 旧“trajectory”是组织组成 adjacency；当前虽然计算了 DPT，但主要是 c8 与其余细胞的二段式排序。
- 没有 branch labels，且 DPT 不是 lineage tracing。
- Root cell 由 first-in-cluster 规则确定，现有 alternative run 没有改变 root；正式汇报前应以 c8 medoid、最高 `Cd27−Itgam` thymus cell 和多个候选 roots 做敏感性分析。

**新 PPT 建议**

- 标题：`Tissue-associated topology is robust; the DPT maturation axis remains provisional`。
- 三联图：UMAP/PAGA、root 与 top edges、program-vs-DPT heatmap；不画确定性多分支发育箭头。

---

### 第 31 页｜Clusters and tissue-specific/condition markers

**旧 PPT**

- 全局 t-SNE 标出组织及 mixed 区域。
- 分组织叠加 Ctrl/T2。
- 三个 marker heatmaps：Thymus 10、BM 22、Spleen 27 genes。
- Notes 表明实际对比是在每个 tissue subset 内做 T2 vs Ctrl，而非简单 tissue-vs-rest。

**本地结果**

- 对旧 marker panel 的当前 tissue-level recovery：Thymus 8/10（80%）、BM 13/22（59%）、Spleen 21/27（78%）。
- Curated iNKT markers 的描述性变化：全局 Rorc +0.404、Xcl1 +0.155、Il4 +0.100；Gzma −0.699、Fasl −0.264。
- 组织例子：BM Rorc +0.406、Gzmb +0.455、Xcl1 +0.183；Spleen Xcl1 +0.229、Rorc +0.229；Thymus Rorc +0.700、Klf2 +0.231。

**共同点**

- 三组织具有不同表达背景；多数旧 marker panel 仍显示 condition 差异。

**关键区别**

- 旧 panel 混有 Iglc2、S100a8/a9、线粒体/应激与低检测率 genes；当前 curated panel 更贴近 subtype、迁移、细胞毒和激活生物学。
- 旧 heatmap 未同时显示 effect、FDR 和 expressing fraction。

**新 PPT 建议**

- 主图用 curated marker/signature effect heatmap；小图给 legacy panel recovery 80%/59%/78%。
- 每个基因同时显示 logFC 与 fraction expressing。

---

### 第 32 页｜Tissue marker Venn 与热图

**旧 PPT**

- Venn 共列 49 genes。
- Thymus-only 5：Mtus2、Mss51、P2rx7、Mttp、Lig4。
- BM-only 17：Rasgef1a、Hist1h1a、Cenps、Aqp11、Cryzl2、Smim40、Samd13、Gm17354、Kif3c、Zfp300、Gm26520、Gm19325、Gm43581、St3gal3、Slc6a20b、Myo19、Pclaf。
- Spleen-only 21：Usp49、Gm37469、1810006J02Rik、E230034D01Rik、Cd200、Ccdc146、Fos、9930022D16Rik、Gm16675、Jchain、Hsph1、Sh3rf1、Dtl、Tctex1d4、S100a8、Gm29336、Ccna2、Gm11579、Prdm16、G0s2、Nek3。
- 四个跨组织共有 genes：Iglc2、Slc15a2、Xlr、3830403N18Rik。
- 组织特异候选包括 thymus P2rx7、BM St3gal3、spleen S100a8；Thymus∩Spleen 有 S100a9。
- BM∩Spleen 另有 Vill。

**本地结果**

- 四个共享 genes 在当前三组织中都显著，且 current T2−Ctrl 方向跨组织一致：Iglc2 上调；Slc15a2、Xlr、3830403N18Rik 下调。
- 代表性 logFC：Iglc2 在 BM/Spleen/Thymus 为 +4.93/+5.76/+26.05；Slc15a2 为 −2.68/−2.64/−2.25。
- 旧 pathway drivers 的稳定性较差：Thymus P2rx7 +0.50、FDR 0.463；BM St3gal3 −0.85、FDR 0.163，均未重现。
- Spleen S100a8/S100a9 强烈下降，但极稀有：Ctrl 表达比例约 0.97%/0.89%，T2 约 0.027%/0%。

**共同点**

- 四个 genes 在三组织中呈方向一致的低频 gene-level signal；这是 analysis-level recovery，不是广泛细胞群 program。

**关键区别**

- P2rx7/St3gal3 不稳定；S100a8/a9 由极少数细胞驱动，不能直接解释为主体 iNKT program。
- Iglc2 的 T2 expressing fraction 仅约 4.64%–4.74%，Ctrl 为 0%–0.20%；Slc15a2、Xlr、3830403N18Rik 的 T2 expressing fraction 均 ≤1.29%，Ctrl 均 ≤4.97%。
- Thymus Iglc2 的极大 logFC 来自 Ctrl mean 接近 0，不代表在主体细胞中广泛增强；图中必须同时显示 fraction expressing。
- S100a8/a9 的 rare-cell、ambient RNA 或细胞组成来源均不能排除，不能直接判定为污染或主体 iNKT regulation。

**新 PPT 建议**

- 标题：`A directionally concordant, low-frequency four-gene signal is recovered across tissues`。
- 主表列 Iglc2 ↑、Slc15a2 ↓、Xlr ↓、3830403N18Rik ↓；旁栏列 `not reproduced` 和 `rare-cell signal`。

---

### 第 33 页｜Tissue marker pathways

**旧 PPT**

- Thymus-like：多个 P2X7/P2RX4/7 terms，绝大多数 overlap=1。
- BM：ST3GAL3 与 peptide-transport terms，overlap=1。
- Spleen：S100A8:S100A9 结合 Zn²⁺/Mn²⁺/TLR4:LY96，overlap=2。
- 页面报告的是 PVAL，不是 FDR；很多 gene set size 只有 1–4。

**本地结果**

- 旧 pathway 结果已从 workbook 恢复，但只属于 `recovered_result_only_not_recomputed`；scope 未完整解析。
- 当前没有同版本 Reactome GMT，不能重算旧 terms。
- 当前 curated signature ORA：
  - Thymus T2-up：iNKT17，22/34，FDR `2.34×10⁻13`。
  - BM T2-up：Direct TCR 7/11、Residency 11/37，均 FDR `8.69×10⁻8`。
  - Spleen T2-up：Direct TCR 6/11，FDR `2.01×10⁻6`；Residency 9/37，`3.01×10⁻6`。
  - Spleen-dominant c6：Direct TCR 与 Residency 均 FDR `5.19×10⁻8`。

**共同点**

- 两边均支持 tissue-dependent activation/stress differences。
- S100a8/a9 gene-level 方向仍有信号；Slc15a2 跨组织下降也稳定，但对应 Reactome term 没有复算。

**关键区别**

- P2rx7 与 St3gal3 作为旧单基因 pathway drivers 未通过当前 FDR。
- S100 terms 由两个极稀有 genes 驱动。
- 旧库与当前 curated library 不同，不能比较 pathway 名次或说“原通路已复现”。
- Legacy enrichment 使用 human-labeled libraries，而输入是 mouse；orthology/mapping 方法和版本未记录。
- Current Direct TCR 与 Residency sets 有 7/11 measured Direct TCR genes 重叠，因此两者同时显著不代表两条独立证据链。
- Thymus T2-up 的 iNKT17 ORA 表示 DEG list 与该 signature 重叠，不等于 iNKT17 细胞比例增加；cluster abundance 与多来源 subtype scores 并未给出同等强度的证据。

**新 PPT 建议**

- 正式主页面改成 `T2-associated programs are tissue dependent`：Thymus iNKT17 enrichment；BM 与 Spleen/c6 为 Direct TCR + Residency。
- 旧 P2X7/ST3GAL3/S100 tables 移附录，并标：`legacy recovered; scope unresolved; not recomputed; single/dual-gene terms`。

## 四、把 33 页压缩成可讲清楚的新故事

### 4.1 Legacy cluster response 的事后最佳功能匹配

下表只表示 **response-program similarity**，不表示 cluster identity。

方向一致率按旧页面标题 `T2 vs Ctrl` 将正 logFC 解释为 T2-up；但 legacy workbook 未保存完整的 group-order 元数据，正式发表前仍应从原分析代码确认该符号约定。第 11/27 页的宽口径 tissue overlap 因此继续使用 direction-agnostic recovery。

候选 current cluster×tissue 单元是在看到结果后，以旧集合 recovery 最大者选出；current recovery 使用 FDR≤0.05 且 `abs(logFC)≥0.25`。因此这些是 **post-hoc descriptive matches**，比例会受选择偏倚影响，不能称为 validation，也不能替代独立数据验证。

| Legacy 单元 | Current 最佳响应单元 | 旧基因重现 | 表观方向一致性（按 PPT 标题解释） | 解释 |
|---|---|---:|---:|---|
| c1 thymus | c5 thymus | 23/30 | 23/23 | Thymus-dominant response |
| c2 BM | c0 BM | 15/18 | 15/15 | 汇入 BM anchor |
| c2 spleen | c4 spleen | 23/27 | 23/23 | Mixed-cluster spleen response |
| c5 BM | c0 BM | 53/61 | 53/53 | 最强 BM recovery |
| c6 spleen，旧称 iNKT2 | c3 spleen | 23/29 | 23/23 | Subtype 仍 source-dependent |
| c7 BM，旧称 iNKT17 | c3 BM | 4/4 | 4/4 | 基因过少，不足以确认 iNKT17 |
| c8 spleen | c6 spleen | 73/84 | 73/73 | 最强 spleen recovery；current c6 偏 iNKT1 |
| c9 BM | c0 BM | 9/11 | 9/9 | 汇入 BM anchor |
| c9 spleen | c6 spleen | 19/28 | 19/19 | 汇入 spleen anchor |

从这张表可得到比“强行映射 cluster 编号”更稳妥的结论：

- 多个 legacy BM responses 在 current c0 中合并。
- 多个 legacy spleen responses 在 current c6 中合并。
- Thymus 的主要 legacy response 在 current c5 中保留。
- 聚类边界变化很大，但主要 T2 transcriptional programs 并未消失。

### 4.2 三个可以作为主线的机制模块

#### A. Immediate-early / TCR activation

- 代表 genes：Fos、Fosb、Jun、Junb、Dusp1、Ppp1r15a、Cd69。
- 证据：旧 global/tissue/cluster tables 与当前 DE 重复出现；Direct TCR activation score 的全局及三组织 T2–Ctrl delta 均为正；BM/Spleen/c6 ORA 显著。
- Direct TCR 与 Residency 的 measured sets 共享 7 个 IEG genes，所以两者共显著应视为部分重复的证据，而非两条独立机制。
- 可讲结论：An immediate-early/TCR-response program is consistently recovered across the legacy and local analyses of this dataset。
- 不可讲结论：没有时间序列与独立重复，不能证明反应持续时间或因果通路。

#### B. Proteostasis / chaperone remodeling

- 诱导型 Hspa1a/Hspa1b 往往上调。
- Constitutive Hspa8、Dnaja1、Hsp90ab1、Hsph1 往往下调。
- 可讲结论：Chaperone/proteostasis program is remodeled bidirectionally。
- 不可讲结论：不能概括成 “HSP70 pathway activated”，因为同一家族成员方向相反。

#### C. Mitochondrial / respiratory-chain remodeling

- 代表 genes：COX、NDUFA/NDUFV、ATP5 相关 genes，以及多个 mt genes。
- 证据：旧、新 DE 和旧 pathway drivers 中高度重复。
- 限定：线粒体 genes 也受 QC、细胞应激和文库质量影响，必须与 mt fraction/QC 并列解释。

### 4.3 新结果中值得单列的描述性变化

- 为追踪性可列 raw `score_genes` deltas，例如 literature iNKT1 −0.013、iNKT2 +0.014、iNKT17 +0.015、Direct TCR +0.109、Residency +0.029；这些值因 gene-set size、control genes 和分布不同而**不可横向比较幅度**。主叙事应依赖跨来源 subtype winners、gene-level DE 与 ORA 的一致性，而不是对 raw deltas 排名。
- Current cluster 比例的 T2–Ctrl 描述性变化：
  - BM：c0 +4.84 percentage points；c2 −4.27；c6 −3.18。
  - Spleen：c3 +4.45；c6 −3.33；c2 −2.44。
  - Thymus：c5 +16.91；c2 −16.72；c4 −6.42。
- 这些只是一个样本/condition/tissue 的细胞组成差异，不能称为可推广的 abundance change。

## 五、建议的新 PPT 主稿结构（18 页）

| 新页 | 建议标题 | 一句话结论 | 主要内容/图 |
|---:|---|---|---|
| 1 | iNKT single-cell reanalysis | Legacy results were audited against a reproducible local Scanpy analysis | 封面；6 samples、18,204 cells、run ID |
| 2 | Executive summary | Input metadata strongly concordant; tissue structure and response programs recover across analyses; cluster IDs do not | 4 个结论卡：input、structure、response、limitations |
| 3 | Concordant input metadata, different retained universe | 样本名和逐样本维度吻合；过滤差异改变细胞/基因宇宙 | 18,458 起点；15,532×10,670 vs 18,204×15,741；注明无旧矩阵哈希 |
| 4 | Reproducible current workflow | 当前 pipeline 参数完整可追踪 | HVG 图、PCA variance、参数条；旧 200-PC sweep 放附录 |
| 5 | Tissue-associated structure is stable across analyses | BM/spleen/thymus 主导群可跨流程恢复 | Current UMAP 四联图 + cluster composition heatmap |
| 6 | Cluster labels are analysis-specific | 旧 c1–c11 与 current c0–c8 只能做功能比较 | namespace 示意 + response-program matching 方法 |
| 7 | Subtype calls depend on signature source | c8 稳定 iNKT2-like；c3 标签不稳定；多数大群偏 iNKT1 | Signature Jaccard + cluster subtype evidence |
| 8 | Activation-related signals recur across analyses | IEG DE 与 curated ORA 均指向 activation/residency；多来源 subtype winners 不支持稳定的全局 conversion | Signature effect heatmap + ORA；标注不同 raw score deltas 不可直接比较 |
| 9 | Legacy global and tissue DEG cores are computationally recovered | 旧 tissue lists 的 85–92% 在同数据本地重分析的显著集合中重现 | Global 15/39（FDR-only）、11/39（另加 abs(logFC)≥0.25）；tissue recovery 92.3/90.0/85.3；显式标 3,000-HVG vs 15,741-gene universe |
| 10 | Post-hoc response-program matches across clusterings | 多个旧 BM/spleen responses 的最大-recovery 候选分别落在 current c0/c6 | 9-row match table；标 `post-hoc, FDR≤0.05, abs(logFC)≥0.25`；不要画 cluster identity 等号 |
| 11 | BM response is strongly recovered across analyses | Legacy c5 BM 53/61 在 current c0 中恢复；符号一致性依赖待确认的旧 group order | BM overlap + 代表基因 slope/dot plot |
| 12 | Spleen response is strongly recovered across analyses | Legacy c8 spleen 73/84 在 current c6 中恢复；符号一致性依赖待确认的旧 group order | Spleen overlap + 代表基因；标 current c6 iNKT1-like |
| 13 | Thymus response is recovered but broader | Legacy c1 thymus 23/30 在 current c5 中保留；当前 thymus list 极大 | Thymus overlap + c5/tissue distinction + caution |
| 14 | Three response modules summarize legacy pathways | Activation、proteostasis、mitochondrial remodeling 比疾病名更可解释 | 三模块图；up/down drivers；旧 disease labels 灰色附注 |
| 15 | Tissue and cluster contexts share—and add—signals | BM/Spleen anchors 与 tissue DE 高重叠，同时存在 cluster-only genes | Current Jaccard/UpSet + legacy recovery cards |
| 16 | A low-frequency four-gene signal is directionally concordant | Iglc2 ↑；Slc15a2/Xlr/3830403N18Rik ↓，但表达细胞比例低 | Effect + fraction-expressing table；P2rx7/St3gal3 fragile；S100 rare/ambient/composition caution |
| 17 | Tissue topology is stable; DPT is provisional | PAGA 支持组织拓扑，但 DPT 主要是 c8 vs rest，且当前 iroot 是 c8 的 first cell | PAGA + DPT diagnostic；把多-root/medoid sensitivity 列为未完成，不画 lineage 箭头 |
| 18 | Conclusions and evidence limits | Major programs recover across analyses of this dataset; biological inference still needs independent animals | Keep/rewrite/drop 总结；下一步实验与分析 |

### 建议附录

1. 旧 33 页逐页审计索引。
2. 样本级 raw/filtered counts 与 QC plots。
3. PCA/UMAP 参数细节和 Leiden resolution 0.2/0.5/1.0。
4. 12 套 signature gene coverage 和 pairwise Jaccard。
5. 所有 legacy cluster response match 的完整表。
6. `tested + 0 DEG` 与 `insufficient cells` 状态表。
7. 旧 pathway table provenance；明确 `recovered, not recomputed`。
8. DPT root、occupancy 和压缩诊断。

## 六、旧稿结论的保留、改写与删除清单

### 可以保留

- 两边的 6 个样本名及逐样本原始维度吻合，高度支持同一 intended dataset。
- 组织来源是嵌入和 cluster composition 的主轴。
- Immediate-early、HSP/proteostasis 和 mitochondrial response 在多个组织/cluster 中重复出现。
- 旧 tissue DEG core 大部分可在本地重现。
- BM 与 spleen 的 cluster-associated T2 response programs 可在本地分析流程中强力恢复。

### 必须改写

| 旧表达 | 新表达 |
|---|---|
| `Xcl1 uniquely expressed in case` | Xcl1 在 Ctrl/T2 均广泛表达；T2 全局仅轻度增加且组织依赖 |
| `c6 is iNKT2` | Legacy c6 response 最接近 current c3；current c3 subtype evidence source-dependent |
| `c7 is iNKT17` | Cox7c/Gm10076/Rplp2/Sfpq 的四基因 response 可跨流程恢复，但不足以建立稳定 iNKT17 cluster identity |
| `iNKT17 cluster` | iNKT17-like program/gradient；跨来源 winner 不稳定 |
| `Prion/Measles/... pathway activated` | 由共享 HSP、IEG、mitochondrial genes 驱动的机制模块 |
| `underlying trajectory` | Tissue-associated topology + provisional expression-state ordering |
| `No DEGs` | Tested with adequate cell count and 0 significant DEG；与 insufficient cells 分开 |

### 主稿应删除或暂缓

- 旧 cluster 编号与 current cluster 编号的一对一等号。
- “200 PCs is the best”——当前没有参数匹配的稳定性验证。
- 旧 PPI interaction network——数据库、版本、阈值和背景不明，当前也未重建。
- P2RX7/ST3GAL3 单基因 pathway 作为主结论——当前 driver 不显著。
- 由 S100a8/a9 得出的主体 iNKT pathway 结论——信号由极少数细胞驱动。
- 把极小 cell-level p 值解释成动物层面的可重复生物效应。
- 把 DPT 解释成真实 lineage 或确定的分支发育方向。

## 七、统一页脚与术语

### 所有 DE、signature effect、cluster proportion 页统一页脚

`Exploratory cell-level comparison; one sample per condition within each tissue; p values do not represent biological-replicate inference.`

中文版本：

`细胞层面探索性比较；每个组织内每个条件仅 1 个样本，p 值不代表独立生物重复推断。`

### 术语

- `Legacy cluster`：旧 PPT 的 c1–c11。
- `Current cluster`：本地 Leiden resolution 0.5 的 c0–c8。
- `Response-program match`：在看到 current DE 后按组织背景、旧基因 recovery 和表观符号一致性挑选的 post-hoc 功能类比；**不是 cluster identity mapping，也不是独立 validation**。
- `Recovery`：旧集合中在当前显著集合内找到的比例。
- `Apparent directional concordance`：按 PPT 的 `T2 vs Ctrl` 标题解释旧 logFC 时的符号一致比例；legacy group order 未由原代码确认。
- `Current curated ORA`：只使用本地 12 套 iNKT-related gene sets，不等于旧 Reactome/KEGG enrichment。

## 八、制作 PPT 时可直接使用的本地图表

### QC、降维与组成

- `output/iNKT_scanpy_tutorial_run/figures/qc_violin_by_sample.png`
- `output/iNKT_scanpy_tutorial_run/figures/qc_scatter_counts_genes_sample.png`
- `output/iNKT_scanpy_tutorial_run/figures/highly_variable_genes.png`
- `output/iNKT_scanpy_tutorial_run/figures/pca_variance_ratio.png`
- `output/iNKT_scanpy_tutorial_run/figures/umap_sample_condition_tissue_cluster.png`
- `output/iNKT_scanpy_tutorial_run/figures/cluster_percent_by_sample.png`
- `output/iNKT_scanpy_tutorial_run/figures/cluster_percent_by_condition.png`

### Signatures 与 markers

- `output/iNKT_extended_runs/20260818_105349/signatures/gene_sets/signature_jaccard_heatmap.png`
- `output/iNKT_extended_runs/20260818_105349/signatures/gene_sets/subtype_signature_membership.png`
- `output/iNKT_extended_runs/20260818_105349/signatures/signature_scores/cluster_subtype_evidence.png`
- `output/iNKT_extended_runs/20260818_105349/signatures/signature_scores/signature_t2_ctrl_effect_heatmap.png`
- `output/iNKT_extended_runs/20260818_105349/signatures/marker_validation/marker_t2_ctrl_effect_heatmap.png`
- `output/iNKT_extended_runs/20260818_105349/signatures/marker_validation/marker_expression_condition_tissue_dotplot.png`

### DEG overlap 与 pathway

- `output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/current_tissue_vs_cluster_overlap_jaccard.png`
- `output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/current_upset_top_patterns__bone_marrow.png`
- `output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/current_upset_top_patterns__spleen.png`
- `output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/current_upset_top_patterns__thymus.png`

### Trajectory

- `output/iNKT_scanpy_tutorial_run/figures/trajectory_paga_cluster_graph.png`
- `output/iNKT_scanpy_tutorial_run/figures/trajectory_paga_compare_umap_pseudotime.png`
- `output/iNKT_extended_runs/20260818_105349/trajectory/figures/program_scores_by_cluster.png`
- `output/iNKT_extended_runs/20260818_105349/trajectory/figures/top_curated_dynamic_genes_along_dpt.png`

## 九、主要数值证据文件

- 总体运行：`output/iNKT_scanpy_tutorial_run/summary.json`
- 输入与过滤：`output/iNKT_scanpy_tutorial_run/tables/input_summary.json`、`filter_summary.json`
- Legacy QC 精确重构：`docs/audits/inkt_qc_exact_reconstruction.csv`；重构规则为全体 cells 上 `min_cells=100` → 重算 QC → `200≤n_genes<2500` 且 `pct_counts_mt<5%`
- Cluster 组成：`output/iNKT_scanpy_tutorial_run/tables/cluster_percent_by_sample.csv`、`cluster_percent_by_condition.csv`
- Global DE：`output/iNKT_scanpy_tutorial_run/tables/rank_genes_condition_t2_vs_ctrl_top500.csv` 只含 top 500；完整 3,000-HVG ranking 位于 `output/iNKT_scanpy_tutorial_run/inkt_scanpy_tutorial_processed.h5ad` 的 `uns/rank_genes_condition_t2_vs_ctrl`
- Current tissue/cluster DE：`output/iNKT_extended_runs/20260818_105349/de_pathway/current_de/`
- Legacy 表恢复：`output/iNKT_extended_runs/20260818_105349/de_pathway/legacy/`
- Tissue/cluster overlap：`output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/`
- PPT 严格 tissue-only overlap：`docs/audits/inkt_legacy_ppt_tissue_overlap_strict.csv`；现有 output 中的 `legacy_vs_current_tissue_gene_overlap.csv` 混入同组织 legacy cluster×tissue rows，只能作为宽 union 审计，不能支持 39/50/34 分母
- Signature QC：`output/iNKT_extended_runs/20260818_105349/signatures/gene_sets/signature_qc.csv`
- Subtype evidence：`output/iNKT_extended_runs/20260818_105349/signatures/signature_scores/cluster_subtype_evidence.csv`
- Signature effects：`output/iNKT_extended_runs/20260818_105349/signatures/signature_scores/signature_t2_vs_ctrl_effects.csv`
- Marker effects：`output/iNKT_extended_runs/20260818_105349/signatures/marker_validation/marker_t2_vs_ctrl_effects.csv`
- Trajectory audit：`output/iNKT_extended_runs/20260818_105349/trajectory/summary.json` 与 `trajectory/tables/`
- Pathway/overlap limitations：`output/iNKT_extended_runs/20260818_105349/de_pathway/overlap_pathway/overlap_pathway_manifest.json`

## 十、制作正式 PPT 前仍需确认的六件事

1. **T2 与旧稿 “CML/case” 的准确生物学定义。** 这会决定标题、图例和结论措辞。
2. **是否存在独立动物 ID 或更多生物重复。** 如果存在，应改为按 animal 做 pseudobulk/replicate-level inference；如果不存在，所有显著性必须保留探索性限定。
3. **是否要求完全复算旧 pathway/PPI。** 若要求，需要明确 Reactome/KEGG/PAGER 或 STRING 的数据库版本、物种、背景集和阈值；旧稿使用 human-labeled libraries 分析 mouse genes，orthology/mapping provenance 也必须补齐。
4. **补做 full-gene global DE。** 当前 global 只测试 3,000 HVGs，而 tissue/cluster 测试 15,741 genes；统一 gene universe 后才能严格定义 current tissue-specific/global overlap。
5. **Signature 跨集合比较先标准化。** 不同 `score_genes` raw deltas 不可直接比较大小；若要画 effect-ranking，应预先定义 per-signature z-score、standardized mean difference 或其他统一效应量。
6. **重做 DPT root sensitivity。** 以 c8 medoid、最高 `Cd27−Itgam` 的 thymus cell 和多个候选 roots 重跑；当前 first-in-cluster root 与“alternative”是同一细胞，不能视为敏感性验证。
