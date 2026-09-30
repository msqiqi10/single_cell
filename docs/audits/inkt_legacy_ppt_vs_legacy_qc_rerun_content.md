# iNKT：原 PPT 与 Legacy-QC 全流程重跑结果逐页内容对照

> 用途：为下一版 PPT 准备内容。本文逐页整理原 PPT、新重跑结果、共同点、差异、建议标题和建议图表；暂不制作幻灯片。
>
> 对照对象：`input/iNKT/iNKT.pptx`（33 页）与 `output/iNKT_legacy_ppt_qc_runs/20260818_125231`。
>
> 本文是新版本，不覆盖此前基于 18,204 × 15,741 宽松 QC 运行撰写的 `docs/inkt_legacy_vs_local_ppt_content.md`。

## 一、这次更新后可以确定的总叙事

1. **旧版 QC 已在计数层面精确复建。** 原始 18,458 cells × 32,285 Gene Expression features、过滤后 15,532 cells × 10,670 genes，以及 6 个样本的过滤后细胞数均与原 PPT 逐项一致。
2. **此前最大的 cell/gene universe 差异已被移除。** 因此本轮比较比 18,204 × 15,741 的宽松 QC 运行更接近原 PPT。
3. **精确复建的是 QC cohort，不是整套旧分析。** 原 PPT 的完整代码、软件版本、100/200-PC 设置细节、聚类算法参数和 pathway library provenance 不齐；后续降维、聚类、DE、signature 与 trajectory 是在同一 QC cohort 上重新分析。
4. **组织结构是最稳定的结果。** 新的 8 群中，c0 为 95.8% BM、c3 为 93.7% spleen、c6 为 97.3% thymus、c7 为 96.7% thymus；原 PPT 同样存在 BM、spleen、thymus 主导群和混合群。
5. **旧、新 cluster 编号仍不能直接对应。** 原 PPT 为 `legacy c1–c11`，新分析为 `current c0–c7`。本文只比较相同组织中的 response-program recovery，不写“旧 c1 = 新 c6”。
6. **旧 tissue DE 核心高度保留。** Thymus 31/34、BM 35/39、spleen 44/50 达到新阈值；在把旧页标题解释为 T2−Ctrl 的前提下，所有恢复基因方向一致。
7. **QC 选择显著影响 thymus response 的范围与方向平衡。** 宽松 QC 时 thymus 为 1,638 DEGs（1,596 up/42 down）；Legacy-QC 重跑为 352（195 up/157 down），方向更平衡，且 thymus-dominant c6 与整个 thymus 的 up/down Jaccard 均约 0.66。
8. **主要旧 cluster-response programs 仍大量恢复。** 代表性结果包括旧 c1-thymus 26/30、旧 c5-BM 53/61、旧 c8-spleen 73/84；已恢复基因的表观方向均一致。
9. **driver-level pathway recovery 很强，但旧 pathway term 没有同库复算。** 多页的 HSP、IEG 与 mitochondrial drivers 可达到 8/8、15/15、20/20；“Prion disease / Measles / Estrogen”等旧数据库名称不应直接解释为疾病机制。
10. **更一致的工作模型仍是 activation/proteostasis/mitochondrial remodeling。** BM 和 spleen 的 Direct TCR activation 与 Residency ORA 反复显著；同时 Hspa1a/b、Fos/Jun/Dusp1 与 Hspa8/Dnaja1 的双向变化提示应激/蛋白稳态重塑，而不是简单的“某条通路全面上调”。
11. **统一口径识别到 c5 的相对 iNKT17-list enrichment，但仍不支持把它升级为稳定、独立的 iNKT17 subtype。** c5 的 literature-list 标准化均值为 iNKT17 `+1.54 SD`、iNKT2 `+1.15 SD`，iNKT17 相对最大竞争项的 margin 为 `+0.39 SD`（cell-bootstrap 95% CI `+0.24` 至 `+0.54`）。然而严格 `min_cells=100` 移除了 Il17a/Il17f，Ccr6 也不在最终 gene universe，且 raw 多来源结果不稳定。因此正文将 c5 保守解释为 mixed iNKT2/iNKT17 program，而不是已确认的离散 subtype。
12. **所有 T2 vs Ctrl 统计仍是探索性的 cell-level comparison。** 每个 tissue × condition 只有一个样本；细胞级 p 值不能替代独立动物重复。

## 二、比较口径与命名

| 项目 | 原 PPT | Legacy-QC 重跑 | 新 PPT 应采用的表达 |
|---|---|---|---|
| 原始输入 | 6 samples；18,458 cells；32,285 genes | 6 samples；18,458 cells；32,285 Gene Expression features，另有 3 个 Multiplexing Capture features | `Input dimensions and per-sample counts match the legacy presentation exactly` |
| 过滤后 | 15,532 × 10,670 | 15,532 × 10,670；6/6 样本计数吻合 | `Legacy QC reconstructed exactly at the cell/gene-count level` |
| QC | 页面只明确写 `min_genes=200`、`min_cells=100` | 计数反推并复现：先在全部 cells 上 `min_cells=100`，再保留 `200≤n_genes<2500`、`pct_mt<5%` | 新方法页补齐旧页未完整记录的上限与顺序 |
| 降维 | UMAP+t-SNE；50/100/200 PCs；neighbors=20 | 3,000 HVGs；50 PCs；前 30 PCs；neighbors=15；UMAP | 不能写成旧参数 sweep 的严格复现 |
| 聚类 | 11 群，legacy c1–c11 | Leiden resolution 0.5，current c0–c7，共 8 群 | 所有图例显式带 `legacy` 或 `current` 前缀 |
| Global DE | 看起来是较宽 gene universe；约 `|logFC|>1.2` | 只在 3,000 HVGs 中排名；默认摘要 FDR≤0.05、`|logFC|≥0.25` | 先区分 `tested` 与 `not tested`，不要把未进入 HVG 的旧基因写成不显著 |
| Tissue/cluster DE | Wilcoxon；阈值和 logFC 实现记录不完整 | 全 10,670 genes；FDR≤0.05、`|logFC|≥0.25` | 比较方向、旧集合 recovery 和新集合范围，不直接比较 logFC 绝对值 |
| Cluster matching | 无跨流程映射 | 同组织内比较旧 DEG 指纹与所有当前 cluster×tissue 单元 | 报告 `maximum-coverage match`；有更紧凑候选时同时报告，始终称 post-hoc response match |
| Pathway | PAGER/Reactome/KEGG 等；存在 human-labeled libraries | 旧 term 仅恢复表格；新 ORA 只用本地 12 套 curated signatures | 只能说 driver/module recovery，不能说旧 pathway 已同库复算 |
| Trajectory | cluster graph + 手工组织箭头 | PAGA + diffusion map + DPT | 只称 exploratory expression-state ordering |

### 两种 response match 的读法

- **最大覆盖匹配**：在相同 tissue 内，选择恢复旧 DEG 数最多的当前单元。这是自动生成的主比较口径。
- **紧凑/特异性候选**：若一个较小的新 DEG 集合获得更高 Jaccard 或 enrichment odds ratio，则作为次级候选列出；其 coverage 可能明显较低，必须同时报告。
- 两种口径不一致时，结论不是“选一个正确 cluster”，而是旧响应被当前聚类拆分或合并；这正说明 cluster identity 不能一一映射。

### 当前 8 个 cluster 的组成

| Current cluster | 细胞数 | 组织组成 | T2 占比 | 保守解释 |
|---|---:|---|---:|---|
| c0 | 4,828 | BM 95.8%，spleen 3.8%，thymus 0.4% | 45.5% | BM-dominant response anchor |
| c1 | 162 | spleen 50.0%，BM 48.1% | 48.1% | 小型 BM/spleen mixed |
| c2 | 378 | BM 38.4%，spleen 31.5%，thymus 30.2% | 50.5% | 三组织 mixed |
| c3 | 5,677 | spleen 93.7% | 47.2% | Spleen-dominant response anchor |
| c4 | 887 | BM 52.6%，spleen 45.5% | 48.8% | BM/spleen mixed；多个旧响应的紧凑候选 |
| c5 | 1,129 | spleen 57.9%，BM 36.2%，thymus 5.8% | 54.4% | BM/spleen mixed；统一 literature-list 口径为 iNKT17-enriched，但 iNKT2 也高，核心 markers 不完整 |
| c6 | 2,379 | thymus 97.3% | 47.2% | Thymus-dominant response anchor |
| c7 | 92 | thymus 96.7% | 47.8% | 小型 thymus 群；统一 literature-list 三项均低于全局均值，保留为 unclassified；DPT root |

## 三、原 PPT 逐页对照与新版内容

### 第 1 页｜封面

**原 PPT**

- 空白页，无分析内容。

**新结果与比较**

- 已获得完整且审计通过的 Legacy-QC 全流程运行。
- 本页可以直接重建，不存在需要保留的旧图。

**建议标题**

`iNKT single-cell reanalysis: exact legacy QC reconstruction`

**建议副标题**

`Ctrl vs T2/CML across bone marrow, spleen and thymus`

**页脚**

`15,532 cells · 10,670 genes · 6 samples · run 20260818_125231`

**必须限定**

- 写“exact legacy QC reconstruction”，不要写“complete replication of the legacy analysis”。

---

### 第 2 页｜样本、输入与 QC

**原 PPT**

| 样本 | 原始细胞 | 过滤后细胞 | 保留率 |
|---|---:|---:|---:|
| Ctrl_BM | 3,998 | 3,379 | 84.52% |
| Ctrl_Spleen | 3,829 | 3,422 | 89.37% |
| Ctrl_Thymus | 2,101 | 1,372 | 65.30% |
| T2_BM | 3,131 | 2,731 | 87.22% |
| T2_Spleen | 3,760 | 3,371 | 89.65% |
| T2_Thymus | 1,639 | 1,257 | 76.69% |
| 总计 | 18,458 | 15,532 | 84.15% |

- 原页明确写 `min_genes=200`、`min_cells=100`，图中还能看到 2,500-gene 和 5% mitochondrial 上限。

**Legacy-QC 重跑**

- 先在全部 18,458 cells 上执行 gene `min_cells=100`，恰得 10,670 genes。
- 在该 gene universe 上重算 QC，再保留 `200≤n_genes<2500` 且 `pct_counts_mt<5%` 的 cells。
- 总规模和 6 个样本过滤后细胞数与原 PPT 6/6 完全一致。
- 原始数据另含 3 个 Multiplexing Capture features；不计入 32,285 个 Gene Expression features。

**共同点**

- 原始与过滤后计数均完全相同，可称为旧 QC cohort 的精确计算复建。

**区别与解释**

- 完整过滤顺序是通过精确计数重构得到，而非来自旧代码；新 PPT 应明确写 `inferred by exact count reconstruction`。
- 旧稿把 condition 写作 CML，而文件和本地元数据使用 T2；正式稿应先定义 T2 与 CML 的关系，再统一为 `T2/CML` 或一个名称。

**建议标题**

`Exact reconstruction of the legacy QC cohort`

**建议图表**

- 左侧样本漏斗表；右侧 QC violin/scatter。
- 页面主数字：`18,458 → 15,532 cells`、`32,285 → 10,670 genes`、`6/6 sample counts matched`。

---

### 第 3 页｜HVG 与 PCA

**原 PPT**

- Top-expressed genes 图，以 Malat1、Tmsb4x、Tpt1、Rps/Rpl、Actb 等高丰度基因为主。
- HVG mean–dispersion 图。
- PC1/2、PC3/4、PC5/6、PC7/8 按样本着色。
- `Principle components` 应更正为 `Principal components`。

**Legacy-QC 重跑**

- 归一化到 10,000 counts/cell 后 log1p。
- 选择 3,000 sample-aware HVGs。
- 计算 50 PCs；邻接图使用前 30 PCs。
- PC1/PC2 分别解释约 2.169%/1.339%；前 30 PCs 累计约 12.447%，前 50 PCs 累计约 15.105%。

**共同点**

- 都采用 HVG → PCA → neighbor graph 的标准路径。

**区别与解释**

- 原 PPT 未记录最终 HVG 数、HVG flavor 或用于邻接图的 PC 数。
- Top-expression 图主要反映 abundance/housekeeping，不支持 subtype 或 condition 机制。

**建议标题**

`HVG selection and PCA representation`

**建议图表**

- 主图用 HVG 和 PCA explained-variance；top-expression 图移附录。
- 参数条：`3,000 HVGs · 50 PCs computed · 30 PCs used downstream`。

---

### 第 4 页｜UMAP/t-SNE 参数探索

**原 PPT**

- `n_neighbors=20`。
- UMAP 与 t-SNE 比较 50、100、200 PCs；备注认为 200 PCs 最佳。

**Legacy-QC 重跑**

- `n_neighbors=15`、`n_pcs=30`、Euclidean metric、`random_state=0`。
- 正式输出 UMAP；没有重跑旧 t-SNE 或 50/100/200-PC sweep。
- 当前只计算了 50 PCs，不能直接复刻 100/200-PC 图。

**共同点**

- 都从 PCA 建立 kNN graph；当前 UMAP 上可见明显组织分层。

**区别与解释**

- 邻居数、PC 数、嵌入算法和软件版本不同，所以二维形状与 cluster 边界不同并不意味着数据/QC 不一致。
- 原页没有定量稳定性指标，不能验证“200 PCs is the best”。

**建议标题**

`Current embedding settings and reproducibility`

**建议图表**

- 使用 sample/condition/tissue/current-cluster UMAP 四联图。
- 旧参数 sweep 如保留，明确标 `legacy, not parameter-matched`。

---

### 第 5 页｜聚类与“轨迹”概览

**原 PPT**

- 11 个旧群，命名为 T、B、S、T_u1、TBS、BS_u1、BS_u2、BS_l1、BS_l2、BS_b1、BS_b2 等。
- 同页有 t-SNE、cluster counts、cluster graph 和相关性热图。
- 名称主要由组织组成推断。

**Legacy-QC 重跑**

- Leiden resolution 0.5 得到 8 群，current c0–c7。
- 新增 PAGA、diffusion map 与 DPT。
- DPT root 为 c7；c7 有 92 cells，其中 89 thymus、3 spleen。
- 实际 `iroot` 是 c7 中的首个 cell `Ctrl_Spleen_ACGATGTAGCTGACCC-1`，不是 thymus medoid 或最高 biomarker-score cell。
- `Itgam` 被 `min_cells=100` 过滤，原 Cd27−Itgam 规则退化为主要由 Cd27 决定。
- c7 median DPT 约 0.103；其余群约 0.977–0.992，几乎是 root-vs-rest 的二段排序。

**共同点**

- 都识别出组织主导状态和 mixed states，并用图结构描述状态邻接。

**区别与解释**

- 旧 11 群变成当前 8 群，反映下游图构建与聚类参数差异；不是 QC cohort 差异。
- 旧“trajectory”主要是 adjacency/相关图；新 DPT 虽更正式，但 root 选择和饱和分布限制了发育解释。

**建议标题**

`Eight current Leiden states and exploratory graph topology`

**建议图表**

- 正文：current UMAP + PAGA + cluster composition。
- DPT 放附录或明确标 `exploratory expression-state ordering; not lineage tracing`。

---

### 第 6 页｜iNKT1 markers

**原 PPT**

- iNKT1 score map、violin、feature maps 和 dotplot。
- 核心基因包括 Xcl1、Nkg7、Klrd1、Klrb1c、Klrk1、Fasl、Gzmb、Cxcr3、Il2rb、Tbx21、Ifng、Gzma。
- 备注：`Xcl1 is uniquely expressed in case (CML)`。

**Legacy-QC 重跑**

- 正文采用统一标准化的 literature lists：先对每个 measured gene 跨全部 15,532 cells 做 z-score，再对 subtype 内基因等权平均，最后把每套 subtype score 重新标准化为全局均值 0、标准差 1。
- 在这把共同尺子下，c4（`+0.72 SD`）和 c6（`+0.81 SD`）具有最强 iNKT1 enrichment；c0 为较弱正向 enrichment（`+0.20 SD`），c1 因竞争 program 未充分分离而标为 mixed。
- c2、c3 的 iNKT1 均低于全局均值，不能因为 raw `score_genes` 的相对 argmax 而称为 iNKT1-enriched。
- Xcl1 并非 T2 特异：Ctrl/T2 expressing fraction 为 67.15%/68.11%，描述性 log2FC +0.166。
- Xcl1 在 BM 和 spleen 略上升，在 thymus 略下降，具有 tissue dependence。
- 其他全局描述性效应：Ifng −0.053、Gzma −0.621、Fasl −0.181、Gzmb +0.058、Prf1 −0.099。

**共同点**

- iNKT1/NK-cytotoxic program 是数据中的主要表达状态，两版均稳定看到。

**关键区别**

- 原页的 `Xcl1 case-specific` 被否定；正确表达是两组广泛表达、T2 仅轻度增加且组织依赖。

**建议标题**

`Unified iNKT1 enrichment is strongest in c4 and c6`

**建议图表**

- 统一 literature-list iNKT1 UMAP 与 cluster violin；三个 subtype 共用 `−3.5` 至 `+3.5 SD` 色标；旁列 Xcl1 Ctrl/T2 expressing fraction。

---

### 第 7 页｜iNKT2 markers

**原 PPT**

- iNKT2 score、violin、feature maps 与大型 dotplot。
- 同时包含 Il4/Gata3/Zbtb16/Rora/Il17rb 等 subtype markers 和 Birc5/Rrm2/Cdca8/Tyms 等增殖基因。

**Legacy-QC 重跑**

- 可测基因：literature 76/86、in-house 16/18、Wang 36/39。
- 统一 literature-list 口径中，c5 的 iNKT2 score 最高（`+1.15 SD`），但 iNKT17 更高（`+1.54 SD`），因此生物学叙事保留为 mixed iNKT2/iNKT17。
- c7 的 iNKT2 score 为 `−1.62 SD`，而且三套 literature-list subtype scores 全部低于全局均值；因此撤回旧的 iNKT2-like 命名，标为 unclassified。
- 原始三来源 `score_genes` 的 argmax 结果只保留在附录作为方法敏感性分析。它不能在不同 gene-set 大小和 control-gene 背景之间形成共同量尺。
- 全局描述性效应：Il4 +0.060、Zbtb16 +0.043；三套 iNKT2 score 的 T2−Ctrl mean delta 均很小（约 +0.001–+0.008）。

**共同点**

- iNKT2-related program 可检测，且 PLZF/IL-4 相关信号仍存在。

**关键区别**

- 旧页把 subtype 与 proliferation 混在一起；新稿应分开 cell-cycle 和 subtype evidence。
- c5 同时承载 iNKT2/iNKT17-related signal；c7 没有在统一主分析中获得正向 subtype enrichment。

**建议标题**

`iNKT2 signal is strongest in c5; c7 is unclassified`

**建议图表**

- 统一 literature-list iNKT2 UMAP 与 cluster violin；使用和 iNKT1/iNKT17 完全相同的色标与统计单位。

---

### 第 8 页｜iNKT17 markers

**原 PPT**

- iNKT17 score、violin、feature maps 和 dotplot。
- 包含 Rorc、Tmem176a/b、Il17a、Il17re、Ccr2/Ccr6、Il23r、Il1r1、Selenop、Serpinb1a。

**Legacy-QC 重跑**

- 可测基因：literature 30/36、in-house 15/22、Wang 31/46。
- Il17a、Il17f 以及若干低频核心 marker 因 legacy `min_cells=100` 不在最终 gene universe。
- 统一 literature-list 口径中，c5 的 iNKT17 score 为 `+1.54 SD`，高于 iNKT2 的 `+1.15 SD`；其动态竞争 margin 为 `+0.39 SD`，10,000 次 sample-stratified cell bootstrap 的 95% CI 为 `+0.24` 至 `+0.54`。
- 该 bootstrap 只反映本数据集内的细胞重采样稳定性，不是动物层面的生物学重复推断。
- Rorc expressing fraction 很低：Ctrl 1.03%、T2 1.49%；描述性 log2FC +0.350。

**共同点**

- Rorc/Tmem176/Il1r1/Il23r 等部分 iNKT17-related signal 仍可检测。

**关键区别**

- c5 的 literature-list iNKT17-related enrichment 是可重复计算的连续 program 结果，但 Il17a/Il17f/Ccr6 不可测、其他来源结果不稳定，因此不足以确认一个离散 iNKT17 subtype。
- 这也不等于 iNKT17 生物学不存在；严格 gene filtering 和低丰度降低了身份判定能力。

**建议标题**

`c5 is legacy-list iNKT17-enriched, but identity is incomplete`

**建议图表**

- 统一 literature-list iNKT17 UMAP 与 cluster violin；旁边列 measured vs missing genes 和 c5 的 bootstrap margin。

---

### 第 9 页｜Literature、in-house 与 Wang signatures

**原 PPT**

- Literature 与 in-house 两套 iNKT1/2/17，共 6 套 score maps 和 violins。
- 不同 gene set 的原始 score 范围直接并列。

**Legacy-QC 重跑**

- 正文主比较只使用原 PPT 对应的 literature iNKT1/2/17 lists，并把三套 score 统一到同一标准差单位；增加 Wang 2022、in-house、circulatory、Direct TCR activation 与 Residency 作为补充/敏感性分析。
- 同 subtype 不同来源 measured-gene Jaccard 仅低到中等：iNKT1 约 0.20–0.44、iNKT2 约 0.06–0.19、iNKT17 约 0.22–0.39。
- 全局 T2−Ctrl mean-score delta：
  - literature iNKT1 −0.0112、iNKT2 +0.0018、iNKT17 +0.0079；
  - in-house iNKT1 −0.0221、iNKT2 +0.0076、iNKT17 +0.0083；
  - Wang iNKT1 −0.0049、iNKT2 +0.0010、iNKT17 +0.0019；
  - circulatory −0.0042、Direct TCR +0.1212、Residency +0.0339。
- Direct TCR 和 Residency 在测得基因中共享 Fos、Icos、Jun、Junb、Ppp1r15a、Prdx6、Ptp4a1；不是两条完全独立的证据。

**共同点**

- 两版都能看到非均匀的 subtype-related expression programs；统一主分析在 c0/c4/c6 看到相对 iNKT1 enrichment，在 c5 看到 iNKT2/iNKT17-related signal。

**关键区别**

- Gene lists 并不等价；raw `score_genes` 值受集合大小和 control genes 影响，不能跨 signature 直接比较绝对幅度。
- 因此 raw `score_genes` panels 移到附录，各自保留独立 colorbar；正文不再用 raw argmax 强制每个 cluster 获得 subtype winner。
- 数据更支持 activation-related shift，而不是稳定的 subtype conversion。

**建议标题**

`One common scale replaces forced subtype winners`

**建议图表**

- 统一 literature-list 3×8 cluster heatmap；三行共用 `−3.5` 至 `+3.5 SD` 色标并显示每格数值；原 PPT raw score UMAP 仅作为定性空间参考。

---

### 第 10 页｜全局 T2 vs Ctrl DEGs

**原 PPT**

- Wilcoxon；页面约写 `|logFC|>1.2`、padj≤0.05。
- 展示 39 genes：33 个按标题解释为 T2-up、6 个 T2-down。
- 代表性上调：Cd52、Jun/Junb、Fos/Fosb、Hspa1b、Iglc2、Eps8l1；下调：Slc15a2、Mxd4、Dynll1、Sfpq、Dnaja1、Hspa8。

**Legacy-QC 重跑**

- Global DE 只在 3,000 HVGs 上排名。
- 345 genes 为 FDR≤0.05；132 genes 同时满足 FDR≤0.05、`|logFC|≥0.25`。
- 旧 39 genes 中只有 20 个进入当前 HVG 检验；20/20 均 FDR≤0.05 且方向一致。
- 其中 16/20 还通过当前 `|logFC|≥0.25`：Cd52、Dnaja1、Eps8l1、Fos、Fosb、H3f3b、Hcst、Hspa1b、Iglc2、Jun、Junb、mt-Nd1、mt-Nd2、Mxd4、Sfpq、Slc15a2。
- 其余 19 个旧基因未进入 global HVG 检验，不能写成“不显著”或“未复现”。
- 机械套用 `|logFC|≥1.2` 后当前仅 Iglc2 与 Slc15a2 通过，说明两套 logFC 标度不应直接对齐。

**共同点**

- 在可检验的旧基因中，IEG/stress 与 Slc15a2/Mxd4/Sfpq/Dnaja1 的方向完全一致。

**关键区别**

- Global gene universe 不同；正确分母是“20 个进入检验的旧基因”，而不是把 19 个未测试基因算作失败。
- 代表性当前 logFC：Iglc2 +5.89、Slc15a2 −2.40、Eps8l1 +1.15、Hspa1b +0.90、Fos +0.77、Jun +0.73、Sfpq −0.70、Dnaja1 −0.66。

**建议标题**

`Legacy global DE directions are retained within the tested HVG universe`

**建议图表**

- 旧 39 genes 分类卡：`16 pass default / 4 FDR-only / 19 not tested`。
- 可配全局 volcano，但不要用“16/39 recovery”作为唯一主数字。

---

### 第 11 页｜按组织的 T2 vs Ctrl DEGs

**原 PPT**

- Thymus 34、BM 39、spleen 50 genes。
- 星号把“不在 global list”解释为 tissue-specific。

**Legacy-QC 重跑**

| 组织 | 当前显著 DEG | T2-up / T2-down | 旧列表恢复 | 表观方向一致 |
|---|---:|---:|---:|---:|
| Thymus | 352 | 195 / 157 | 31/34（91.2%） | 31/31 |
| Bone marrow | 718 | 350 / 368 | 35/39（89.7%） | 35/35 |
| Spleen | 864 | 356 / 508 | 44/50（88.0%） | 44/44 |

- 未达到当前阈值的旧基因：thymus 为 Fau、Rplp1、Ubb；BM 为 Actg1、mt-Cytb、Ptpn18、Rplp2；spleen 为 Actb、Fau、mt-Co3、Myl6、Rplp1、Rplp2。

**共同点**

- 原 PPT 的 tissue DE 核心被高度、同向恢复，是整套比较中最强的结果之一。

**关键区别**

- 当前列表较大，主要因为使用 `|logFC|≥0.25`、完整 10,670-gene universe 以及不同 logFC 实现；不能说发现了十倍更多的独立机制。
- 旧星号不能直接复核：当前 global 只测试 HVGs，tissue DE 测试全部 genes；`global not tested` 不等于 tissue-specific。
- 与此前宽松 QC 比较，thymus 从 1,638（1,596/42）降到 352（195/157），说明旧的“thymus 广泛上调”高度 QC-sensitive。

**建议标题**

`Tissue-stratified DE shows strong directional concordance with the legacy result`

**建议图表**

- 三组织 recovery bars + up/down counts；长基因表放附录。

---

### 第 12 页｜Cluster 与组织组成

**原 PPT**

- 11 个 legacy clusters 的 t-SNE、graph、相关热图和组成饼图。
- 旧 c1/c11 偏 thymus，c5 偏 BM，c6/c8 偏 spleen，其余为 mixed groups。

**Legacy-QC 重跑**

- 8 个 current clusters；组成见本文件前表。
- c0 95.8% BM、c3 93.7% spleen、c6 97.3% thymus、c7 96.7% thymus。
- c1、c2、c4、c5 为不同程度 mixed clusters。
- 各群 T2 占比为 45.5%–54.4%，没有 condition-exclusive cluster。

**共同点**

- 两套分析都呈现 `tissue-dominant anchors + mixed states`，组织来源是最稳定的结构信号。

**关键区别**

- 群数由 11 变 8；旧、新编号无一一关系。
- 新 c0/c3/c6 是后续组织 response 比较的 anchors，但不能分别称为旧 c5/c8/c1 的“同一个群”。

**建议标题**

`Tissue identity dominates cluster composition; no cluster is condition-exclusive`

**建议图表**

- cluster×tissue percentage heatmap 或 100% stacked bar；格内标 N。
- 全文固定写 `legacy c1–c11` 与 `current c0–c7`。

---

### 第 13 页｜Legacy c1 thymus DE

**原 PPT**

- Legacy c1 thymus：30 DEGs，23 up、7 down。
- 核心模式：Tmsb10/Cd52/Hspa1b 上调；Ccl5/Cxcr6/Dnaja1/Hspa8/Sfpq/Dynll1 下调。

**Legacy-QC 重跑**

- 最大覆盖且组织合理的匹配为 current c6 thymus：1,087 T2、1,227 Ctrl；304 DEGs（140 up/164 down）。
- 旧 30 genes 中恢复 26/30（86.7%）；在旧标题代表 T2−Ctrl 的假设下，26/26 方向一致。
- 恢复的关键基因包括 Hspa1a/b、Fos、Dusp1、Tmsb10、Cd52、Hspa8、Dnaja1、Ccl5、Cxcr6。
- 未达到当前阈值：Fau、Rplp1、mt-Co3、mt-Cytb。

**共同点**

- Thymus 主导群的 activation/proteostasis response 高度保留，而且比宽松 QC 运行的 23/30 更高。

**关键区别**

- Current c6 有 304 DEGs，响应范围比旧 30-gene 表更宽。
- 这是 response-program match，不是旧 c1 的身份标签转移。

**建议标题**

`Legacy c1 thymus response is retained in current thymus-dominant c6`

**建议主数字**

`26/30 recovered · 26/26 directionally concordant`

---

### 第 14 页｜Legacy c1 thymus pathways

**原 PPT**

- 38 个显著 terms，但只有 8 个独特 driver genes。
- 标签涉及 HSP70/HSP40/HSF1、MAPK/NFκB、ER protein processing，以及 prion/measles/atherosclerosis 等疾病名称。

**Legacy-QC 重跑**

- Current c6 thymus 中 8/8 drivers 全部达到当前 DEG 阈值并保持表观方向：Ccl5、Cox8a、Dnaja1、Dusp1、Fos、Hspa1a、Hspa1b、Hspa8。
- 相比此前宽松 QC 的 6/8，QC 对齐后 driver recovery 更完整。
- 当前没有相同版本的 Reactome/KEGG GMT，38 个原 terms 未同库重算。

**共同点**

- HSP/chaperone、即时早期和线粒体 driver module 被完整恢复。

**关键区别**

- Hspa1a/b、Fos/Dusp1 上调，而 Hspa8/Dnaja1/Ccl5 下调，应称双向 remodeling，不应概括成“heat-shock pathway 全面激活”。
- 旧疾病 term 主要由共享 HSP/mitochondrial genes 驱动，且 human-library 到 mouse 的映射方法未记录。

**建议标题**

`All eight legacy c1 pathway drivers are retained in current c6`

**建议图表**

- 8-gene direction heatmap；旧 38-term 表移附录并标 `legacy term not recomputed`。

---

### 第 15 页｜Legacy c2：BM 与 spleen DE

**原 PPT**

- BM：18 DEGs（15 up/3 down）。
- Spleen：27 DEGs（23 up/4 down）。
- 两组织共享 Tmsb10/Cd52/mitochondrial 上调与 Hspa8/Dnaja1 下调；spleen 另见 Jun 上调、Hsp90aa1 下调。

**Legacy-QC 重跑：BM**

- Current c4 BM：223 T2、244 Ctrl；83 DEGs（56 up/27 down）。
- 恢复 17/18（94.4%），方向 17/17 一致；仅 Ppp1r12a 未达到阈值。
- 这是覆盖和特异性均清晰的匹配。

**Legacy-QC 重跑：spleen**

- 最大覆盖为 current c3 spleen：2,557 T2、2,762 Ctrl；735 DEGs；恢复 19/27，方向 19/19。
- 更紧凑的候选为 current c4 spleen：208 T2、196 Ctrl；45 DEGs；恢复 18/27，方向 18/18，Jaccard 更高。
- 两个候选的 overlap 数量只差 1，但基因身份并不相同：两者共同恢复 14 个，另有 c3-only 5 个、c4-only 4 个。因此不应宣称唯一对应；保守叙述可把 c4 作为跨 BM/spleen 的紧凑匹配，同时注明 c3 的最大覆盖。

**共同点**

- 两组织共享的 Tmsb10/Cd52 与 Hspa8/Dnaja1 response 保留；Jun/Hsp90aa1 是 spleen 侧特征。

**关键区别**

- 同一个 legacy c2 response 在当前空间中可能被 c3/c4 分担，说明聚类边界已改变。
- Spleen recovery 由此前宽松 QC 的 23/27 降至最大覆盖 19/27；QC 对齐并非所有单元都提高 overlap。

**建议标题**

`Legacy c2 response is compactly represented by current c4, with broader spleen coverage in c3`

---

### 第 16 页｜Legacy c2：BM 与 spleen pathways

**原 PPT**

- BM：4 terms、4 个独特 drivers，集中于 HSP70/HSP40 与“prion disease”。
- Spleen：13 terms、7 个独特 drivers，集中于 HSP90/HSP70、线粒体与 Jun。

**Legacy-QC 重跑**

- BM/current c4：4/4 drivers 全部保留——Cox8a、Dnaja1、Hspa8、Ndufa4。
- Spleen 最大覆盖 current c3：7/7 drivers 全部保留。
- Spleen 紧凑 current c4：4/7，保留 Dnaja1、Hsp90aa1、Hspa8、Jun；Cox8a、Ndufa4、Ndufa7 未达阈值。
- Current c4 BM 的 curated ORA 支持 Direct TCR activation 与 Residency；c4 spleen 没有 FDR<0.05 的 curated ORA，而广谱 c3 spleen 强烈支持两者。

**共同点**

- BM driver recovery 由此前 3/4 提高至 4/4；spleen 在最大覆盖候选中由 6/7 提高至 7/7。

**关键区别**

- Driver completeness 与 DEG-set specificity 给出不同的最佳候选；应把两者并列，而不是伪造唯一 cluster identity。
- 原 pathway library 未同库复算。

**建议标题**

`Legacy c2 pathway drivers are fully recoverable, but distributed across current response states`

---

### 第 17 页｜Legacy c5 BM DE

**原 PPT**

- 61 DEGs（50 up/11 down）。
- T2 中 mitochondrial/ribosomal、Jun/Junb/Fos/Dusp1/Hspa1a/b 上调；Hspa8/Dnaja1/Hsph1/Mxd4/Sfpq/Dynll1 下调。

**Legacy-QC 重跑**

- 最大覆盖且 BM-dominant 的匹配为 current c0 BM：2,090 T2、2,533 Ctrl；592 DEGs（274 up/318 down）。
- 恢复 53/61（86.9%），方向 53/53 一致。
- 未达到阈值：Actg1、Fau、Ptpn18、Rplp1、Rplp2、mt-Co2、mt-Co3、mt-Cytb。
- Current c4 BM 是 lower-coverage、higher-specificity 的次级候选：34/61；其 Jaccard 与 enrichment odds ratio 更高，但 Fisher P 与 c0 近似且并不更小。

**共同点**

- BM activation/proteostasis/mitochondrial program 是全稿最强的跨分析信号之一。

**关键区别**

- Current c0 的 response 很宽，可能汇入多个旧 BM clusters；不能把它直接改名为 legacy c5。

**建议标题**

`Legacy c5 BM response is strongly recovered in the current BM-dominant anchor`

**建议主数字**

`53/61 recovered · 53/53 directionally concordant`

---

### 第 18 页｜Legacy c5 BM pathways

**原 PPT**

- 32 terms，只由 15 个独特 drivers 支撑。
- “Prion disease”等首要标签由 COX/NDUFA/HSP genes 驱动；其余涉及 HSP70、MAPK、ER processing。

**Legacy-QC 重跑**

- Current c0 BM 中 15/15 drivers 全部达到当前阈值并保持表观方向。
- Drivers 包括 Hspa1a/b、Hspa8、Hsph1、Dnaja1、Fos、Jun、Dusp1 与多个 Cox/Ndufa genes。
- P17 的 higher-specificity alternative c4 只保留 8/15 drivers；缺 Cox7a2、Cox7b、Hsph1、Jun、Ndufa13、Ndufa2、Ndufa7，显示紧凑性是以 driver coverage 为代价。
- Current curated ORA：Direct TCR activation 与 Residency T2-up，FDR 均约 10⁻6。

**共同点**

- BM 的 pathway-driver recovery 完整，是最强的 module-level concordance。

**关键区别**

- “Prion disease”应改写为 mitochondrial respiratory-chain/proteostasis shared module；不能推断神经退行性疾病机制。
- Direct TCR 与 Residency sets 共享多个 IEG drivers，并非两次独立验证。

**建议标题**

`A conserved BM activation–proteostasis–mitochondrial response module`

---

### 第 19 页｜Legacy c6 spleen（旧称 iNKT2）DE

**原 PPT**

- 29 DEGs（24 up/5 down）。
- Jun/Junb/Ndufv3/Tmsb10 上调；Hspa8/Dnaja1/Hsp90ab1/Dynll1 下调。
- 页面把 legacy c6 命名为 iNKT2。

**Legacy-QC 重跑**

- Current c5 spleen 与 c3 spleen 均恢复 22/29；current c5 只有 73 DEGs，因此是更紧凑、Jaccard 更高的候选。
- Current c5 spleen：387 T2、267 Ctrl；73 DEGs（44 up/29 down）；22/22 方向一致。
- 未达到阈值：Cox8a、H3f3b、mt-Co3、mt-Nd4、Rplp2、Uba52、Ubb。
- Current c5 的 literature/in-house signatures 偏 iNKT2，但 Wang 偏 iNKT1；亚型仍有来源依赖。

**共同点**

- Spleen 中 Jun/Junb activation 与 Hspa8/Dnaja1/Hsp90ab1 down-regulation 高度保留。

**关键区别**

- 不能把旧名称自动继承到 current c5；response similarity 和 subtype identity 是两套证据。

**建议标题**

`Legacy c6 spleen response has a compact match in current c5; subtype remains provisional`

---

### 第 20 页｜Legacy c6 spleen pathways

**原 PPT**

- 13 terms、6 个独特 drivers：Dnaja1、Hspa8、Hsp90ab1、Jun、Cox8a、Ndufv3。

**Legacy-QC 重跑**

- 紧凑 current c5 spleen 保留 5/6：Dnaja1、Hsp90ab1、Hspa8、Jun、Ndufv3；Cox8a 未达阈值。
- 广谱 current c3 spleen 可保留 6/6，但总 DEGs 为 735，特异性较低。
- Current c5 curated ORA：Residency 与 Direct TCR activation 均 T2-up，FDR 约 0.005。

**共同点**

- 旧 driver module 与当前 activation/residency enrichment 相互支持。

**关键区别**

- 模式是 `IEG activation up + constitutive chaperone down` 的双向重塑。
- 旧 13 个 pathway terms 未同库重算。

**建议标题**

`The legacy c6 HSP/IEG driver module is largely retained in current c5 spleen`

---

### 第 21 页｜Legacy c7 BM（旧称 iNKT17）DE

**原 PPT**

- 只有 4 个 DEGs：Cox7c、Gm10076、Rplp2 上调，Sfpq 下调。

**Legacy-QC 重跑**

- Current c0 BM 与 c5 BM 都恢复 3/4，且方向 3/3 一致；缺失项均为 Rplp2。
- c0 更 BM-dominant，但有 592 DEGs；c5 只有 39 DEGs，因此更紧凑。
- 没有唯一、稳定的新 cluster match。
- 统一 literature-list score 把 c5 判为相对 iNKT17-enriched，但 c5 同时有较高 iNKT2 score，且关键 iNKT17 markers 缺失；因此没有足够证据把它命名为稳定、独立的 iNKT17 subtype。

**共同点**

- Cox7c/Gm10076/Sfpq 的稀疏 response 可恢复。

**关键区别**

- 4 个通用 mitochondrial/ribosomal/RNA-processing genes 不足以支持 iNKT17-specific identity。

**建议标题**

`The sparse legacy c7 BM response is recoverable but not uniquely mappable`

**建议处理**

- 改成“不确定性/候选比较”页，或并入 response-summary；不要继续用确定性 iNKT17 标题。

---

### 第 22 页｜Legacy c8 spleen DE

**原 PPT**

- 84 DEGs（72 up/12 down），为旧稿最大的 cluster×tissue responses 之一。

**Legacy-QC 重跑**

- 最大覆盖且 spleen-dominant 的匹配为 current c3 spleen：2,557 T2、2,762 Ctrl；735 DEGs（275 up/460 down）。
- 旧 84 genes 中恢复 73/84（86.9%），方向 73/73 一致。
- 未达到阈值：Actb、Fau、mt-Co2、mt-Co3、Myl6、Ptpn18、Rplp1、Rplp2、Tmsb4x、Trac、Xist。
- Current c3 的统一 literature-list 三项 scores 均略低于全局均值，因此标为 unclassified；原 raw `score_genes` 的 iNKT1 argmax 不再作为主结论。

**共同点**

- Spleen 主响应的 activation/HSP/mitochondrial/TCR gene program 高度恢复。

**关键区别**

- 高覆盖来自一个 735-DEG 的广谱响应；应强调 coverage，不应声称高特异性或 cluster identity 保留。

**建议标题**

`The broad legacy c8 spleen response is highly covered by current spleen-dominant c3`

---

### 第 23 页｜Legacy c8 spleen pathways

**原 PPT**

- 15 terms、20 个独特 drivers；包含“Prion disease”、ER processing、MAPK、HSP70/HSP90。

**Legacy-QC 重跑**

- Current c3 spleen 中 20/20 drivers 全部达到阈值并保持表观方向。
- 包括 Cox6b1/Cox7a2/Cox7b/Cox8a、Ndufa2/4/7/13、Hspa8/Hsph1/Hsp90ab1/Dnaja1、Fos/Jun/Jund/Dusp1/Ppp1r15a 等。
- Current curated ORA 强烈支持 Direct TCR activation 与 Residency T2-up，FDR 约 10⁻6。

**共同点**

- Spleen driver-level recovery 完整。

**关键区别**

- 旧 disease labels 仍只能解释为共享 HSP/mitochondrial genes 的数据库标签。
- Current c3 的统一 subtype scores 接近全局基线并标为 unclassified；旧 c8 的身份名称不能继承。

**建议标题**

`All legacy c8 pathway drivers are retained in the current spleen anchor`

---

### 第 24 页｜Legacy c9：BM 与 spleen DE

**原 PPT**

- BM：11 DEGs（9 up/2 down）。
- Spleen：28 DEGs（23 up/5 down）。

**Legacy-QC 重跑：BM**

- 最大覆盖为 current c0 BM：9/11，方向 9/9；未达到阈值为 Actb、Rplp2。
- 更紧凑的 current c4 BM：8/11，方向 8/8；未达到阈值为 H2afj、Tma7、Rplp2。

**Legacy-QC 重跑：spleen**

- 最大覆盖为 current c3 spleen：19/28，方向 19/19。
- 更紧凑的候选为 current c4（12/28）与 c5（13/28）；c4 的集合富集更强，c5 的 pathway-driver coverage 更完整。
- 因此 spleen 没有唯一 cluster match。

**共同点**

- 两组织仍保留 Cd52/Tmsb10 上调与 Dnaja1/Hspa8 下调；spleen 还保留 Jun/Junb/Hsp90ab1/Sfpq。

**关键区别**

- Legacy c9 response 被多个 current states 吸收，支持 cluster consolidation/boundary redraw，而不是一一重命名。

**建议标题**

`Legacy c9 cross-tissue response is distributed across current BM/spleen states`

---

### 第 25 页｜Legacy c9 pathways

**原 PPT**

- 页面标题只写 BM pathways，但上方 3-term 表对应 BM、下方 12-term 表按上下文对应 spleen；下方标签缺失。
- BM 由 Dnaja1/Hspa8 驱动；spleen 由 Dnaja1/Hspa8/Hsp90ab1/Jun 驱动。

**Legacy-QC 重跑**

- 沿用第 24 页的 DEG maximum-coverage mapping：BM/current c0 为 2/2 drivers；spleen/current c3 为 4/4。
- 沿用第 24 页的紧凑候选：BM/current c4 为 2/2；spleen/current c4 为 3/4，current c5 为 4/4。
- 若只按 driver coverage，本页没有唯一最佳群：BM 的 c0/c4 都是 2/2，spleen 的 c3/c5 都是 4/4。

**共同点**

- Dnaja1–Hspa8 chaperone axis 在 BM 和 spleen 均稳定恢复；spleen 的 Hsp90ab1/Jun 也能在合适候选中恢复。

**关键区别**

- 原页标签需要修正。
- 多个 cluster 均可承载相同 2–4 gene driver set，driver recovery 不能建立 cluster identity。
- 原 pathway term 未同库复算。

**建议标题**

`The legacy c9 chaperone driver axis is retained across current response states`

---

### 第 26 页｜Legacy c10/c11：无 DEG

**原 PPT**

- 标题写 `T2 vs. Ctrl DEGs in c10 (No DEGs) and c11 (No DEGs)`。
- 没有报告各组细胞数、阈值或检验状态。

**Legacy-QC 重跑**

- 完成检验且 0 个显著 DEG 的 6 个单元：
  - c1-BM：32 T2 / 46 Ctrl；
  - c1-spleen：45/36；
  - c2-spleen：74/45；
  - c2-thymus：50/64；
  - c5-thymus：46/20；
  - c7-thymus：44/45。
- 另有 7 个单元因任一组少于 20 cells 而跳过：c0-thymus、c1-thymus、c3-thymus、c4-thymus、c6-BM、c6-spleen、c7-spleen。

**共同点**

- T2 response 具有 cluster context；确实存在完成检验后未检测到显著变化的群。

**关键区别**

- 新分析能严格区分 `tested + 0 DEG` 与 `not tested: insufficient cells`；原页没有这个区分。
- 0 DEG 不等于生物学上绝对无变化，仍受 power、dropout 和多重校正影响。

**建议标题**

`Not all rerun tissue × cluster strata show a detectable T2 response`

**建议图表**

- 状态矩阵：格内标 T2/Ctrl N 和 DEG 数；颜色区分 `DEG`、`tested-zero`、`skipped`。

---

### 第 27 页｜Tissue DEG 与 cluster-specific DEG overlap

**原 PPT**

- Tissue Venn：thymus 34、BM 39、spleen 50，三组织共有 19。
- Cluster Venn：legacy c1-thymus 30、c5-BM 41、c8-spleen 50，三者共有 17。
- Legacy c5 的 41 是完整 61 genes 经过严格 `|logFC|>1.2` 后的子集；c8 的 50 因旧图只给四舍五入数值，无法精确恢复筛选规则。

**Legacy-QC 重跑：tissue lists**

- BM 718、spleen 864、thymus 352；旧列表分别恢复 35/39、44/50、31/34。
- 当前三组织 union 为 1,426 genes：
  - BM-only 381；spleen-only 503；thymus-only 163；
  - BM∩spleen only 190；BM∩thymus only 18；spleen∩thymus only 42；
  - 三组织共有 129。

**Legacy-QC 重跑：组织 anchor 与整个 tissue**

| Anchor | Anchor DEG | 与 tissue 共享 | Up Jaccard | Down Jaccard | Anchor-only |
|---|---:|---:|---:|---:|---:|
| c0-BM | 592 | 489 | 0.629 | 0.566 | 103 |
| c3-spleen | 735 | 600 | 0.656 | 0.566 | 135 |
| c6-thymus | 304 | 261 | 0.667 | 0.655 | 43 |

**共同点**

- 旧 tissue 核心高度保留；current tissue-dominant anchors 的大部分响应也出现在对应组织总体。

**关键区别**

- 当前集合因阈值和 gene universe 不同而更大；旧 19 与新 129 不能被解释为“共享反应增强”。
- QC 对齐后 thymus anchor–tissue 的上下调 Jaccard 均约 0.66，明显比宽松 QC 时的 0.064/0.147 更一致。
- `legacy_vs_current_tissue_gene_overlap.csv` 混入了 legacy cluster-tissue rows；严格的 tissue 分母必须使用原页提取的 34/39/50 或 `legacy_ppt_de_best_matches.csv`。

**建议标题**

`Legacy tissue DEG cores are recovered; tissue anchors show balanced overlap`

**建议图表**

- 左侧 3 个 legacy recovery bars；右侧 anchor–tissue up/down Jaccard。
- 完整三组织 UpSet 放附录，避免大集合 Venn。

---

### 第 28 页｜Legacy c1/c5/c8 cluster-associated DEG 验证

**原 PPT**

- 三组 violin：c1-thymus、c5-BM、c8-spleen。
- 代表性 cluster-associated genes：c1 的 Klf6/Cxcr6/Rhob；c5 的 Tomm6/Pink1/Fus；c8 的 Atp5j2/Actg1/Xcl1/Ccnd2 等。

**Legacy-QC 重跑**

| Legacy 单元 | Current response anchor | 恢复 | 表观方向一致 | 代表性恢复基因 |
|---|---|---:|---:|---|
| c1-thymus，30 genes | c6-thymus | 26/30（86.7%） | 26/26 | Klf6、Cxcr6、Rhob |
| c5-BM，完整 61 genes | c0-BM | 53/61（86.9%） | 53/53 | Tomm6、Pink1、Fus |
| c5-BM，旧 Venn 严格 41 genes | c0-BM | 37/41（90.2%） | 37/37 | 同上 |
| c8-spleen，完整 84 genes | c3-spleen | 73/84（86.9%） | 73/73 | Atp5j2、Actg1、Xcl1、Ccnd2 |

**共同点**

- 三个主要组织 response programs 都被大量恢复。

**关键区别**

- 这些是 post-hoc response matches，不是旧 cluster identity 的独立验证。
- c8 的严格 50-gene membership 不能从旧页四舍五入数值精确恢复；正式数字应优先使用可审计的完整 84-gene 表。
- 多个恢复基因为 housekeeping/stress/ribosomal/mitochondrial genes；统计复现不等于 subtype specificity。

**建议标题**

`Legacy response programs map to current tissue anchors—not to matching cluster IDs`

---

### 第 29 页｜DEG interaction network

**原 PPT**

- 给出 legacy c1-thymus、c5-BM、c8-spleen DEG interaction networks。
- 常见 hub：Jun、Actb、mt-Cytb、Cox6b1、Ndufa4、Ndufa13、Cox6c、Cox7c、Atp5e、Hspa1b。
- 未记录 interaction database、物种/版本、confidence cutoff、背景集或 edge 定义。

**Legacy-QC 重跑**

- 本次没有重新生成 PPI/network；PAGA 是 cell-state graph，不能代替 gene interaction network。
- 对上述 10 个 hub 做 gene-level 核对：c6-thymus 5/10、c0-BM 7/10、c3-spleen 9/10 达到当前 DEG 阈值。
- Ndufa13、Cox6c、Atp5e、Hspa1b 在三个 anchors 中均达到阈值。

**共同点**

- 多数旧 hub genes 的 response 再次出现，尤其在 spleen anchor。

**关键区别**

- 只能说 gene-level hub response 恢复；旧 edge、degree 与 topology 没有复现。
- 泛连接的 stress/mitochondrial/housekeeping genes 容易成为网络中心，需单独做稳健性分析。

**建议标题**

`Gene-level response hubs recur, but the legacy interaction network was not reproduced`

**建议处理**

- 正文改为 10-gene effect/FDR dotplot；旧网络放附录。
- 若后续重建，必须注明 Mus musculus 数据库版本、score cutoff、background，并测试移除 mt/ribosomal/general-stress genes 后的稳定性。

---

### 第 30 页｜Clusters and underlying trajectories

**原 PPT**

- 按 tissue composition 手工命名 Thymus、BM、Spleen、BM+spleen、upr-mixture、btm-mixture 等区域。
- 右侧 graph 和手工箭头被解释为 trajectory，但没有 root、pseudotime 或 lineage 方法。

**Legacy-QC 重跑**

- 组织主轴稳定：c0 BM、c3 spleen、c6/c7 thymus。
- PAGA 较强连接包括 c4–c6 0.393、c0–c1 0.223、c4–c5 0.208、c1–c3 0.175、c2–c7 0.168。
- DPT root 为 c7；实际 root cell 是首个 c7 cell，而非 medoid/最高 score cell。
- c7 median DPT 0.103；其余群 0.977–0.992，近似 root-vs-rest。
- Tissue 内 T2−Ctrl mean DPT 很小：BM +0.00029、spleen +0.00041、thymus −0.00289。
- 与 DPT 正相关的 programs：iNKT1 Wang ρ≈0.466、iNKT1 literature 0.436、in-house iNKT1 0.374、Residency 0.364、stage3 cytotoxic 0.362、Direct TCR 0.321。
- Alternative root 又选择同一个 c7 和同一首细胞，ρ=1.0 不是独立 sensitivity test。

**共同点**

- 组织相关 topology 稳定；mixed states 连接组织主导群。

**关键区别**

- 新分析有正式 PAGA/DPT，但 DPT 严重饱和。两种规则都选择 c7 作为 root-cluster candidate；真正未验证的是 root cell 被机械设为 c7 首个细胞，且 alternative run 又使用同一细胞。因此 root-cell selection 仍近乎任意，敏感性尚未检验，不能支持确定性多分支发育模型。

**建议标题**

`Tissue topology is robust; DPT remains a saturated root-versus-rest ordering`

**建议图表**

- UMAP/PAGA + cluster DPT violin；不画确定性 lineage arrows。

---

### 第 31 页｜Clusters and tissue-specific/condition markers

**原 PPT**

- 标题称 `tissue-specific marker`，但 heatmaps 实际是在每个 tissue 内比较 T2 vs Ctrl。
- Panel 大小：thymus 10、BM 22、spleen 27。

**Legacy-QC 重跑：旧 panel recovery**

- Thymus 8/10（80.0%）。未恢复：Mttp、P2rx7。
- BM 16/22（72.7%）。未恢复：Aqp11、Gm19325、Zfp300、Kif3c、St3gal3、Cenps。
- Spleen 20/27（74.1%）。未恢复：Gm37469、Usp49、Ccna2、Cd200、G0s2、Sh3rf1、Nek3。

**Legacy-QC 重跑：curated markers**

- 全局描述性 log2FC：Rorc +0.350、Xcl1 +0.166、Klf2 +0.109、Il4 +0.060；Gzma −0.621、Fasl −0.181、Prf1 −0.099、Ifng −0.053。
- 组织例子：BM Gzmb +0.426、Gzma −1.298；spleen Xcl1 +0.241、Rorc +0.294；thymus Rorc +0.863、Gzma −0.679。
- Il17a/Il17f 因 `min_cells=100` 不在当前 gene universe。

**共同点**

- 旧 marker panels 中多数基因仍达到当前 within-tissue DE 阈值。

**关键区别**

- 原页名称不准确：这些不是 tissue-vs-rest markers，而是 within-tissue T2-vs-Ctrl genes。
- 旧 panel 混有 low-frequency、stress 与非 subtype-specific genes；平均表达图应补 fraction expressing。

**建议标题**

`Within-tissue T2 responses are more informative than “tissue-specific markers”`

**建议图表**

- Curated marker effect heatmap + fraction-expressing dotplot；角落保留旧 panel recovery 80.0%/72.7%/74.1%。

---

### 第 32 页｜Tissue marker Venn 与重点候选

**原 PPT**

- Venn 共列 49 genes。
- 四个跨组织共有 genes：Iglc2、Slc15a2、Xlr、3830403N18Rik。
- 组织候选包括 thymus P2rx7/Mttp、BM St3gal3、spleen S100a8；thymus∩spleen 有 S100a9。

**Legacy-QC 重跑：四个跨组织共同基因**

| Gene | BM logFC（FDR） | Spleen logFC（FDR） | Thymus logFC（FDR） |
|---|---:|---:|---:|
| Iglc2 | +5.851（7.3×10⁻32） | +5.456（1.2×10⁻32） | +25.767（4.2×10⁻13） |
| Slc15a2 | −2.734（1.6×10⁻10） | −2.613（4.5×10⁻17） | −1.784（4.3×10⁻5） |
| Xlr | −3.332（5.3×10⁻15） | −3.410（1.6×10⁻13） | −2.476（0.0042） |
| 3830403N18Rik | −2.074（2.6×10⁻4） | −2.538（2.3×10⁻6） | −3.827（0.017） |

**Legacy-QC 重跑：不稳定/稀有候选**

- Thymus P2rx7：+0.569，FDR 0.746；不显著。
- Thymus Mttp：+0.937，FDR 0.085；不显著。
- BM St3gal3：−1.062，FDR 0.084；不显著。
- Spleen S100a8：−5.682，FDR 6.3×10⁻7，但 T2/Ctrl expressing fraction 只有 0.03%/1.02%。
- Spleen S100a9：−23.454，FDR 2.8×10⁻6，但 T2/Ctrl 为 0%/0.85%。

**共同点**

- 四基因跨组织信号在新结果中方向完全一致并通过阈值。

**关键区别**

- 四基因信号是低频 gene-level recovery：Iglc2 在 T2 中也只有约 4.5% cells 表达；thymus 巨大 logFC 来自 Ctrl mean 接近 0。
- P2rx7、Mttp、St3gal3 未达到当前 tissue-DE 阈值；S100a8/a9 由极少数 cells 驱动，可能涉及 rare cells、ambient RNA 或组成变化。

**建议标题**

`A low-frequency four-gene cross-tissue signal is retained; several tissue candidates are unstable`

**建议图表**

- 用 `logFC + FDR + fraction expressing` bubble/table 替代 Venn。
- 标记 `not reproduced at the current tissue-DE threshold`：P2rx7、Mttp、St3gal3；标记 `rare-cell signal`：S100a8/a9。

---

### 第 33 页｜Tissue marker pathways

**原 PPT**

- Thymus：多个 P2RX4/7 terms；BM：ST3GAL3/peptide-transport terms；spleen：S100A8:S100A9 metal/TLR4 terms。
- 多数 term 只有 1–2 个 overlap genes，报告 PVAL 而非 FDR。
- 数据库版本、mouse→human mapping 与背景集不清楚。

**Legacy-QC 重跑：current curated ORA**

- BM T2-up：Direct TCR 6/11，FDR 5.76×10⁻6；Residency 9/36，1.02×10⁻5；Wang iNKT1 7/44，0.00207；literature iNKT1 4/25，0.0250。
- Spleen T2-up：Direct TCR 6/11，6.37×10⁻6；Residency 9/36，1.18×10⁻5。
- Thymus T2-up：literature iNKT17 6/30，1.71×10⁻4；Wang iNKT17 5/31，0.00134；Direct TCR 3/11，0.00356；Residency 4/36，0.0121。
- Thymus T2-down：Wang iNKT1 7/44，3.80×10⁻5；in-house iNKT1 6/35，6.31×10⁻5；literature iNKT1 4/25，0.00179；Wang iNKT17 3/31，0.0312（Cd7、Cxcr6、Itgae）。
- Dominant anchors：c0-BM 与 c3-spleen 的 Direct TCR/Residency T2-up FDR 约 10⁻6；c6-thymus 的 Wang/in-house iNKT1 T2-down FDR 约 0.0065/0.0117。

**共同点**

- 两版都显示 tissue-dependent marker/program differences；activation/stress 模型主要由跨页 driver recovery 与 current curated ORA 综合支持。

**关键区别**

- P2rx7、St3gal3 未通过当前 tissue-DE 阈值；S100a8/a9 通过 DEG 阈值但极低频。对应旧 Reactome terms 均未同库重算，因此 pathway 本身没有被重新检验。
- 没有可用于同库复算的外部 Reactome/KEGG GMT；current ORA 仅使用 12 套本地 curated signatures。
- Direct TCR 与 Residency gene sets 有较多共享 IEG genes，同时显著不是两条独立证据链。
- Thymus iNKT17 ORA 只表示 DEG-list overlap，不代表 iNKT17 cell fraction 增加。c5 的统一 literature-list iNKT17 enrichment 也是连续 program 证据；两者都不能单独证明离散 subtype 或 population expansion。

**建议标题**

`Curated activation and subtype programs are recoverable; legacy Reactome terms remain unrecomputed`

**建议图表**

- 正文放 current curated ORA dotplot。
- 旧 P2RX7/ST3GAL3/S100 表移附录，并标：`legacy result recovered; scope partly unresolved; not recomputed; one-/two-gene driven`。

## 四、最重要的“原 PPT vs 新结果”结论矩阵

| 旧结论/内容 | 新结果判断 | 新 PPT 处理 |
|---|---|---|
| 过滤后 15,532 × 10,670 | **在细胞数、基因数与 6 个样本计数层面精确复建** | 保留并升级为 QC audit 主结果；旧 PPT 无 cell/gene ID 或矩阵 hash |
| 6 个样本逐样本细胞数 | **6/6 完全一致** | 保留 |
| 组织主导结构 | **稳定恢复** | 作为结构主线 |
| 旧 11 个 cluster 编号 | **未复现为相同编号/边界** | 禁止按数字映射；用 current c0–c7 |
| iNKT1 为主导 program | **c0/c4/c6 在统一 literature-list 口径中相对富集；标准化分数不能用于断言全数据总体 subtype abundance** | 改写为 cluster-relative program enrichment |
| Xcl1 仅在 case 表达 | **不支持** | 改为广泛表达、T2 轻度且组织依赖变化 |
| 稳定 iNKT2 群 | **不支持 c7：统一 literature-list 三项均低于全局均值；c5 同时具有 iNKT2/iNKT17 signal** | c7 标 unclassified；c5 标 mixed program |
| 稳定 iNKT17 群 | **c5 有显著相对 enrichment，但核心 markers 与跨来源一致性不足** | 保留 continuous program，不能升级为离散 subtype identity |
| Global legacy DEG core | **20/20 tested genes 同向；16 过默认阈值** | 区分 tested/not-tested |
| Tissue legacy DEG core | **88.0%–91.2% 恢复且全部同向** | 强保留 |
| Legacy major cluster responses | **约 86.7%–86.9% 高覆盖** | 作为 response-program recovery，不能称 cluster validation |
| Legacy pathway drivers | **多数完整恢复** | 保留 module/drivers；删除确定性 disease interpretation |
| Legacy Reactome/KEGG term significance | **未同库复算** | 仅附录显示 recovered legacy result |
| P2RX7/ST3GAL3 pathways | **P2rx7/St3gal3 未通过当前 tissue-DE 阈值；旧 pathway 未同库检验** | 主稿删除；附录分开标 driver not retained 与 pathway not recomputed |
| S100A8/A9 pathway | **S100a8/a9 通过 DEG 阈值但极低频；旧 pathway 未同库检验** | 标 rare-cell gene signal，不作主体机制或 pathway replication |
| 多分支 developmental trajectory | **不支持** | 降级为 exploratory topology/DPT |
| T2 vs Ctrl 的生物学重复推断 | **不成立** | 全页统一写 cell-level exploratory |

> 表中所有“同向”结论都假设原 PPT 标题 `T2 vs Ctrl` 的正 logFC 表示 T2−Ctrl；旧代码的 group order 不可得，正式发表前仍需确认。

## 五、建议的新 PPT 主稿结构（18 页）

1. **Cover** — exact legacy QC reconstruction；15,532 cells、10,670 genes、6 samples。
2. **Question and dataset** — T2/CML 的定义、三组织、每个 tissue×condition 一个样本。
3. **Exact QC reconstruction** — 18,458→15,532 cells、32,285→10,670 genes、6/6 count match。
4. **What is and is not reproduced** — QC cohort exact；downstream parameters not identical。
5. **Current embedding and eight clusters** — UMAP、Leiden parameters、current namespace。
6. **Tissue composition dominates structure** — c0 BM、c3 spleen、c6/c7 thymus，且无 condition-exclusive cluster。
7. **Unified subtype-program evidence** — c0/c4/c6 相对 iNKT1 enrichment、c5 mixed iNKT2/iNKT17、c2/c3/c7 unclassified；raw multi-source scores 移入敏感性附录。
8. **Marker claim corrections** — Xcl1 not case-specific；Rorc sparse；Il17a/f filtered。
9. **Global DE audit** — 20/20 tested legacy genes directionally concordant；区分未测试 HVGs。
10. **Tissue DE concordance** — thymus 31/34、BM 35/39、spleen 44/50；方向全一致。
11. **Legacy response recovery map** — 所有旧 cluster×tissue 单元的 coverage/compact candidates。
12. **Thymus response module** — legacy c1→current c6：26/30，drivers 8/8。
13. **BM response module** — legacy c5→current c0：53/61，drivers 15/15。
14. **Spleen response module** — legacy c8→current c3：73/84，drivers 20/20。
15. **Three mechanistic modules** — immediate-early/TCR activation、proteostasis/chaperone remodeling、mitochondrial remodeling。
16. **Tissue-dependent curated ORA** — BM/spleen TCR+residency；thymus iNKT1-down/iNKT17-overlap，附严格解释。
17. **Tissue anchors and exploratory topology** — anchor–tissue overlap；PAGA；DPT 限制。
18. **Conclusions and limitations** — 哪些稳定、哪些不稳定、无 biological replicates。

### 建议附录

- 原 PPT 33 页逐页缩略图与本文件的 comparison summary。
- 完整 sample/QC 表和参数。
- 8-cluster composition、统一标准化 subtype evidence 与 bootstrap calls。
- Raw `score_genes` 三来源敏感性分析；每个 panel 使用独立 colorbar，并明确不可跨 panel 比较绝对值。
- 所有旧 DEG recovery tables、compact alternatives 与未恢复基因。
- DE status matrix：tested DEG / tested zero / insufficient cells。
- 完整 UpSet/Jaccard 图。
- 旧 pathway 表，并统一标注 `not recomputed`。
- DPT/root audit 和多个 root sensitivity analysis（如果后续补做）。

## 六、新 PPT 的统一措辞与页脚

### 建议固定措辞

- `Legacy QC was reconstructed exactly at the cell- and gene-count level.`
- `All downstream clusters use the current legacy_ppt_qc_leiden_res_0_5 namespace.`
- `Legacy c1–c11 and current c0–c7 are not numerically mapped.`
- `Legacy response matches are post-hoc, tissue-constrained program comparisons—not cluster identity mappings.`
- `T2–Ctrl effects are exploratory cell-level comparisons; one sample per condition within each tissue does not establish biological-replicate inference.`
- `Legacy pathway terms were recovered from the presentation/workbook but not recomputed with the original library.`

### 所有 DE/ORA 页统一页脚

`Wilcoxon cell-level exploratory comparison; FDR≤0.05 and |logFC|≥0.25 unless noted; one biological sample per tissue × condition.`

### 术语

- 第一次出现写 `T2/CML`，并在研究背景页定义；之后统一用一个名称。
- 旧群：`legacy c1–c11`。
- 新群：`current c0–c7` 或 `Legacy-QC rerun c0–c7`。
- 用 `response-program recovery/match`，不用 `cluster replication/identity`。
- 用 `exploratory expression-state ordering`，不用 `developmental lineage`。

## 七、可直接用于制图的主要输出

### QC、降维与组成

- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/qc_violin_by_sample.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/qc_scatter_counts_genes_sample.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/highly_variable_genes.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/pca_variance_ratio.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/umap_sample_condition_tissue_cluster.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/cluster_percent_by_sample.png`

### Signatures 与 markers

- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/signature_scores/cluster_subtype_evidence.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/signature_scores/signature_t2_ctrl_effect_heatmap.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/marker_validation/marker_t2_ctrl_effect_heatmap.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/marker_validation/marker_expression_condition_tissue_dotplot.png`

### DEG overlap 与 pathway

- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_tissue_vs_cluster_overlap_jaccard.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_upset_top_patterns__bone_marrow.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_upset_top_patterns__spleen.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_upset_top_patterns__thymus.png`

### Trajectory

- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/trajectory_paga_cluster_graph.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/figures/trajectory_diffmap_dpt_pseudotime.png`
- `output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/trajectory/figures/program_scores_by_cluster.png`

## 八、主要数值证据文件

- 运行状态：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/pipeline_status.txt`
- QC 验证：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/qc_validation.json`
- 完整运行审计：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/run_audit.json`
- 预处理摘要：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/summary.json`
- 完整处理对象（含 3,000-gene global DE ranking）：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/inkt_scanpy_tutorial_processed.h5ad`
- 当前 cluster 组成：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/rerun_cluster_composition.csv`
- 原 PPT 表格提取：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_de_tables_extracted.csv`
- 最大覆盖比较：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_de_best_matches.csv`
- 所有 cluster 候选：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_de_all_matches.csv`
- Legacy pathway driver retention：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_pathway_driver_retention.csv`
- Current curated ORA：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_offline_gene_set_enrichment.csv`
- Tissue-anchor overlap：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_tissue_vs_cluster_overlap_summary.csv`
- 统一 literature-list per-cell scores：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/standardized_legacy/standardized_scores_per_cell.csv.gz`
- 统一 cluster×subtype 摘要：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/standardized_legacy/standardized_score_summary.csv`
- 统一 cluster calls 与 bootstrap intervals：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/standardized_legacy/standardized_cluster_calls.csv`
- Raw `score_genes` T2−Ctrl effects（附录敏感性）：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/signature_scores/signature_t2_vs_ctrl_effects.csv`
- Raw multi-source cluster evidence（附录敏感性）：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/signature_scores/cluster_subtype_evidence.csv`
- Marker effects：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/signatures/marker_validation/marker_t2_vs_ctrl_effects.csv`
- Trajectory summary：`output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/trajectory/summary.json`

## 九、制作正式 PPT 前仍建议补做的事项

1. **定义 T2/CML。** 明确 T2 是样本代号、疾病模型、处理还是时间点，并统一全文。
2. **若必须声称复刻旧 downstream analysis，补跑旧参数。** 至少需要 100/200 PCs、neighbors=20、原聚类算法/分辨率和 t-SNE；否则继续称 reanalysis。
3. **补 full-gene global DE。** 让 global/tissue/cluster 使用同一 10,670-gene universe，才能严格定义 tissue-specific genes。
4. **决定 response matching 主口径。** 主稿建议用可审计的 maximum coverage，附录同时给 compact/Jaccard/Fisher candidates。
5. **若保留 pathway term，获取原 GMT/library 版本。** 需要物种、版本、orthology、background 和多重校正方法。
6. **若保留 trajectory，重做 root sensitivity。** 使用 c7 thymus medoid、最高 biomarker-score thymus cell 和多个候选 roots，而不是同一首细胞重算。
7. **如要做生物学推断，必须增加 replicate-aware 分析。** 当前每个 tissue×condition 只有一个样本，无法建立动物层面的可推广结论。
