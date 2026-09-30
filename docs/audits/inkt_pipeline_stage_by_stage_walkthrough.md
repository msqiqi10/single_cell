# iNKT reproduction pipeline：逐阶段详细解读

建议一边打开 [reproduction PDF](/home/zzz0054/bio3/output/iNKT_reproduction_deck/iNKT_legacy_QC_reproduction_detailed.pdf)，一边按下面顺序看。这里讲的是实际执行的分析流程，不只是 PPT 的叙事顺序。

整个流程可以压缩成：

```text
6 个 10x 样本
→ RNA/capture feature 分离
→ 精确复建 legacy QC
→ normalize + log1p
→ HVG
→ PCA
→ neighbor graph
→ UMAP + Leiden
→ marker/signature annotation
→ tissue 和 cluster×tissue DE
→ legacy DEG/program matching
→ DEG overlap + curated ORA
→ PAGA/diffusion map/DPT
→ 综合结论与证据分级
```

## 首先记住两个最重要的边界

第一，所谓“exact reproduction”只适用于 QC cohort：

- 原始数据：18,458 cells × 32,285 Gene Expression features。
- QC 后：15,532 cells × 10,670 genes。
- 六个样本的 QC 后细胞数逐个与旧 PPT 一致。

但 downstream 不是旧分析的逐位复刻：

- 旧 PPT：11 个 legacy clusters，参数记录不完整。
- 当前分析：8 个 Leiden clusters，使用当前 Scanpy 参数。
- 所以后续应该叫“在精确复建的 QC cohort 上重新分析”，不能叫整个旧分析完全复现。

第二，统计设计只有六个样本：

```text
Ctrl_BM       T2_BM
Ctrl_Spleen   T2_Spleen
Ctrl_Thymus   T2_Thymus
```

每个 tissue × condition 只有一个样本，也就是每个比较组没有独立动物重复。几千个细胞不是几千个生物学重复。因此：

- 可以描述这个数据集内部的效应方向、大小和一致性。
- 可以生成探索性的 cell-level p-value。
- 不能把这些 p-value 解释成可推广到动物总体的 T2 效应证据。

---

# Stage 0：运行、审计与复现范围

对应 PDF 的第 1、34、35 页。

实际运行顺序是：

1. 构建并执行 preprocessing notebook。
2. 验证 QC 结果是否精确匹配。
3. 检查脚本能否编译，运行 signature/trajectory 的 preflight。
4. 并行运行：
   - gene-set audit
   - signature scoring
   - marker validation
   - legacy table extraction
   - 两个 DE shards
   - extended trajectory
5. DE 完成后再做 overlap 和离线 ORA。

这次 audited run：

- Python 3.12.12
- Scanpy 1.12.1
- CPU-only
- 总运行时间约 2 分 53 秒
- 所有任务 exit code 都是 0
- QC validator 状态为 `passed`

你应该把“pipeline 成功”理解为：

> 所有预期文件都生成了，关键维度、参数和数据结构通过了自动检查。

但它不代表：

> 生物学假设已经被证明，或者旧 PPT 的每个分析参数都被恢复。

运行状态可看：[pipeline_status.txt](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/pipeline_status.txt)

---

# Stage 1：读取六个原始 10x 矩阵

对应 PDF 第 2 页的输入部分。

每个样本用 `sc.read_10x_mtx` 读取，使用 gene symbol 作为 feature name。随后：

- 给 barcode 加上样本名前缀，避免不同样本出现同名 barcode。
- 写入三个 metadata：
  - `sample`
  - `condition`
  - `tissue`
- 六个 RNA 矩阵以 outer join 合并。

每个输入包含：

- 32,285 个 Gene Expression features
- 3 个 Multiplexing Capture features：BM、Spleen、Thymus

三个 capture features 被单独记录和汇总，没有放入 RNA 表达矩阵，也没有用于重新 demultiplex 或过滤细胞。

原始样本数：

| Sample | Raw cells | QC 后 | 保留率 |
|---|---:|---:|---:|
| Ctrl_BM | 3,998 | 3,379 | 84.5% |
| Ctrl_Spleen | 3,829 | 3,422 | 89.4% |
| Ctrl_Thymus | 2,101 | 1,372 | 65.3% |
| T2_BM | 3,131 | 2,731 | 87.2% |
| T2_Spleen | 3,760 | 3,371 | 89.7% |
| T2_Thymus | 1,639 | 1,257 | 76.7% |
| Total | 18,458 | 15,532 | 84.1% |

这里值得注意的是 thymus 保留率明显较低。因此后面看到 thymus 的细胞组成、DEG 数量或 trajectory 与其他组织不同，不能完全排除 QC 选择带来的影响。

这一阶段你应该检查：

- 六个样本是否都成功读取。
- barcode 是否唯一。
- RNA 和 capture features 是否正确分离。
- raw sample counts 加总是否为 18,458。
- metadata 是否与目录名称一致。

---

# Stage 2：Legacy QC 的精确复建

这是整个 reproduction 最重要、也最容易因顺序不同而复现失败的步骤。

## 2.1 先过滤基因

第一步不是先过滤细胞，而是：

```python
sc.pp.filter_genes(adata, min_cells=100)
```

含义是：

> 一个基因必须在全部 18,458 个原始细胞中至少被 100 个细胞检测到。

结果：

```text
18,458 × 32,285
→ 18,458 × 10,670
```

共移除 21,615 个低频基因。

## 2.2 在缩小后的 10,670-gene universe 上重新计算细胞 QC

定义：

- mitochondrial：`mt-` 或 `MT-` 开头
- ribosomal：`Rps`、`Rpl`、`RPS`、`RPL`
- hemoglobin：匹配 `^(Hb[ab]|HBA|HBB)`

然后计算每个细胞的：

- `total_counts`
- `n_genes_by_counts`
- `pct_counts_mt`
- ribosomal/Hb 百分比等

关键点是：`n_genes_by_counts` 是在 10,670 个保留基因中重新计算的。若先对 32,285 genes 计算 QC，再过滤基因，会得到不同的细胞集合。

## 2.3 再过滤细胞

保留条件：

```text
n_genes_by_counts >= 200
n_genes_by_counts < 2500
pct_counts_mt < 5%
```

最终结果：

```text
15,532 cells × 10,670 genes
```

实际保留边界：

- 最低 detected genes：284
- 最高 detected genes：2,499
- 最高 mitochondrial fraction：4.99866%

## 2.4 QC 图应该怎么看

### Violin plot

你会看到每个样本的：

- detected genes
- total counts
- mitochondrial percentage

重点不是每个 violin 必须形状一样，而是：

- retained cells 不应超过 2,499 genes。
- mitochondrial percentage 不应达到或超过 5%。
- thymus 分布和保留率可能与 BM/spleen 不同。

### Counts-versus-genes scatter

一般来说：

- counts 与 detected genes 正相关。
- 极高 genes/counts 的点可能是 doublets。
- 高 mt% 的点可能是受损或濒死细胞。

不过，本流程没有进行：

- doublet detection
- ambient RNA correction
- total-count cutoff
- cell-cycle regression
- batch integration
- ribosomal/Hb 基因过滤
- mitochondrial regression

因此“通过 legacy QC”不等于完成了现代标准下所有可能的 QC。

## 2.5 这一阶段能支持什么

可以说：

> Legacy QC 在总细胞数、基因数和六个样本的细胞数层面得到了精确重建。

不能说：

> 已证明每个 cell ID、gene ID 和旧分析完全逐位一致。

旧代码和旧 cell/gene hash 不存在，所以只能证明当前本地矩阵经过这套顺序产生了完全相同的计数结果。

详细验证：[qc_validation.json](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/qc_validation.json)

---

# Stage 3：Normalization 和数据层设计

过滤完成后执行：

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()
```

近似理解为：

$$
x'_{ig} =
\log\left(
1+
\frac{x_{ig}}{\sum_g x_{ig}}
\times 10,000
\right)
$$

也就是先让每个细胞的总表达量标准化到 10,000，再做 `log(1+x)`。

最终 H5AD 中有三种重要表达表示：

- `layers["counts"]`：QC 后的原始整数 counts。
- `X`：normalized + log1p expression。
- `raw`：所有 10,670 genes 的 normalized + log1p expression。

这里有一个容易混淆的地方：

> `adata.raw` 在 Scanpy 里不一定是“raw counts”。本项目里的 `.raw` 是 log-normalized expression；真正的整数 counts 在 `layers["counts"]`。

Normalization 的作用是减弱测序深度差异，不会：

- 自动消除 batch effect。
- 使不同样本成为生物学重复。
- 消除 tissue composition。
- 证明某个基因的变化是调控效应。

---

# Stage 4：HVG selection 和 PCA

对应 PDF 第 3 页。

## 4.1 选择 3,000 个 highly variable genes

参数：

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="sample",
    flavor="seurat",
)
```

含义：

- 从 10,670 genes 中选择最能描述细胞间变化的 3,000 genes。
- 使用 `sample` 作为 batch key，让 HVG 选择考虑六个样本中的重复可变性。
- 有 355 个基因在六个样本中都被标记为 HVG intersection。

注意：

> `batch_key="sample"` 只是让 HVG 选择 sample-aware，不是 batch correction。

对象本身仍保留全部 10,670 genes。HVG 只是给 PCA 和 notebook 中的初始 DE 作为 mask。

## 4.2 PCA

参数：

```python
sc.tl.pca(
    n_comps=50,
    svd_solver="arpack",
    mask_var="highly_variable",
    random_state=0,
)
```

没有先进行 gene scaling 或 covariate regression。

结果：

- PC1：2.169%
- PC2：1.339%
- 前 2 PCs：3.508%
- 前 10 PCs：8.403%
- 前 30 PCs：12.447%
- 前 50 PCs：15.105%

## 4.3 PCA variance plot 怎么看

你应该看到 explained variance 随 PC 编号逐步降低，没有非常明显的单一大轴。

低解释率不一定是问题。scRNA-seq 数据：

- 稀疏
- 维度高
- 存在多个小的生物/技术变化来源

所以单个 PC 解释 1–2% 很常见。

但这里没有 scaling 和 integration，因此前几个 PC 可能同时携带：

- tissue identity
- sample/library differences
- mitochondrial/ribosomal variation
- activation state

不要看到 PC1/PC2 分离，就直接称为 iNKT1→iNKT2 differentiation。

---

# Stage 5：Neighbor graph、UMAP 和 Leiden clustering

对应 PDF 第 4、5、12 页。

## 5.1 Neighbor graph

参数：

```python
sc.pp.neighbors(
    n_neighbors=15,
    n_pcs=30,
    use_rep="X_pca",
)
```

实际记录：

- 前 30 PCs
- 15 neighbors
- Euclidean metric
- random state 0

含义是：对每个细胞，在 PCA 空间中找表达状态最接近的邻居，构建细胞图。

## 5.2 UMAP

```python
sc.tl.umap(random_state=0)
```

UMAP 是 neighbor graph 的二维可视化，不是新的生物学测量。

正确读法：

- 近邻细胞通常有较相似的表达状态。
- 小岛或局部结构可能表示 distinct states。
- UMAP 轴没有生物学单位。
- 两个岛之间的二维距离不是发育时间。
- UMAP 岛的面积也不直接等于细胞数量或生物学重要性。

## 5.3 Leiden clustering

运行了三个 resolution：

- resolution 0.2 → 6 clusters
- resolution 0.5 → 8 clusters
- resolution 1.0 → 16 clusters

主分析采用：

```text
leiden_res_0_5
```

参数：

- `flavor="igraph"`
- `n_iterations=2`
- `directed=False`
- random state 0

主聚类结果：

| Current cluster | Cells | Tissue composition | T2 fraction |
|---|---:|---|---:|
| c0 | 4,828 | 95.8% BM | 45.5% |
| c1 | 162 | 50.0% spleen、48.1% BM | 48.1% |
| c2 | 378 | BM/spleen/thymus mixed | 50.5% |
| c3 | 5,677 | 93.7% spleen | 47.2% |
| c4 | 887 | 52.6% BM、45.5% spleen | 48.8% |
| c5 | 1,129 | 57.9% spleen、36.2% BM | 54.4% |
| c6 | 2,379 | 97.3% thymus | 47.2% |
| c7 | 92 | 96.7% thymus | 47.8% |

## 5.4 这里真正应该看到什么

最清楚的结构是 tissue：

- c0 是 BM anchor。
- c3 是 spleen anchor。
- c6 是 thymus anchor。
- c7 是一个很小的 thymus-dominant state。
- c1、c2、c4、c5 是不同程度的 mixed states。

所有 cluster 都同时含 Ctrl 和 T2，没有 condition-exclusive cluster。说明：

> T2 更像是在已有组织/状态内部改变表达程序，而不是创造一个完全独立的新细胞群。

不过由于没有 batch integration，而且每个 tissue×condition 只有一个文库，condition effect 与 sample/library effect 无法完全拆开。

## 5.5 为什么不能把 legacy c1–c11 对应到 current c0–c7

Cluster number 是算法输出标签，没有固定生物学含义。

旧分析有 11 群，当前有 8 群，说明旧状态可能在当前分析中：

- 被合并；
- 被拆分；
- 聚类边界发生变化。

因此不能写：

```text
legacy c1 = current c6
```

正确写法是：

```text
legacy c1 thymus response program is best recovered in current c6 thymus
```

这是 response-program match，不是细胞身份的一一映射。

---

# Stage 6：Marker 和 gene-signature annotation

对应 PDF 第 6–9、31 页。

这一阶段分成三层证据。

## 6.1 手工 marker panel

四组 marker：

- iNKT/T lineage：`Trac, Trav11, Traj18, Cd3d, Cd3e, Zbtb16`
- NK/cytotoxic：`Klrb1c, Nkg7, Prf1, Gzma, Gzmb, Ifng`
- Maturation：`Cd27, Il2rb, Klrg1, Tbx21`
- Migration/activation：`Cxcr6, S1pr1, Sell, Ccr7, Il7r`

Dotplot/stacked violin/matrixplot 用来回答：

- 哪些 cluster 表达某个 marker？
- 是大量细胞低表达，还是少量细胞高表达？
- marker 是 cluster-specific，还是广泛存在？

不要用单个 marker 给 cluster 定性。尤其 activation marker 很容易受组织、应激和制样过程影响。

Cluster marker ranking 使用：

- 每个 cluster versus rest
- Wilcoxon
- 只测试 3,000 HVGs
- BH correction
- 保存每群 top 100

这是候选 marker 排序，不是已经验证的细胞类型标签。

## 6.2 审计 12 套 signatures

共有：

- Literature iNKT1/iNKT2/iNKT17
- In-house iNKT1/iNKT2/iNKT17
- Wang 2022 iNKT1/iNKT2/iNKT17
- Circulatory
- Direct TCR activation
- Residency

实际测得的基因覆盖：

| Signature | Canonical | Measured |
|---|---:|---:|
| Literature iNKT1 | 26 | 25 |
| Literature iNKT2 | 86 | 76 |
| Literature iNKT17 | 36 | 30 |
| In-house iNKT1 | 35 | 35 |
| In-house iNKT2 | 18 | 16 |
| In-house iNKT17 | 22 | 15 |
| Wang iNKT1 | 47 | 44 |
| Wang iNKT2 | 39 | 36 |
| Wang iNKT17 | 46 | 31 |
| Circulatory | 25 | 24 |
| Direct TCR activation | 14 | 11 |
| Residency | 40 | 36 |

iNKT17 的覆盖最弱。`Il17a`、`Il17f` 和 `Ccr6` 等未进入最终 gene universe，部分是严格 `min_cells=100` 的直接结果。

因此：

> 没有检测到 Il17a/Il17f 不能等价为没有 iNKT17 biology。

不同来源的同类 signatures 重叠也不高：

- iNKT1 Jaccard：0.20–0.438
- iNKT2：0.057–0.191
- iNKT17：0.216–0.386

所以不存在一个唯一、绝对的“iNKT1 score”。跨来源一致比单一列表结果更可信。

## 6.3 每个细胞计算 signature score

这一版 PDF 把两种算法分开，避免把不同统计单位混为一谈。

### 6.3.1 正文主分析：统一标准化 literature-list score

正文第 6–9 页只使用原 PPT 对应的 literature iNKT1/iNKT2/iNKT17 lists。输入是 `.raw` 中已经做过 library-size normalization 和 `log1p` 的表达值：

$$
E_{ig}
=
\log\left(
1+
\frac{C_{ig}}{\sum_h C_{ih}}
\times 10{,}000
\right)
$$

其中：

- $C_{ig}$：细胞 $i$ 中基因 $g$ 的原始 UMI count。
- $\sum_h C_{ih}$：该细胞所有保留基因的总 counts。
- $E_{ig}$：用于下游评分的 log-normalized expression。

第一步，对每个基因跨全部 15,532 个细胞做 z-score：

$$
Z_{ig}
=
\frac{E_{ig}-\mu_g}{\sigma_g}
$$

这里 $\mu_g$、$\sigma_g$ 都由全部 15,532 个细胞计算，标准差使用 `ddof=0`。因此高表达基因不会仅凭绝对表达量大就支配 gene set。

第二步，对 subtype $k$ 的所有可测成员基因等权平均：

$$
A_{ik}
=
\frac{1}{\lvert G_k^{\mathrm{obs}}\rvert}
\sum_{g\in G_k^{\mathrm{obs}}} Z_{ig}
$$

缺失基因不按 0 补入，也不进入分母。实际覆盖为：

- iNKT1：25/26；缺失 `Slam4`。
- iNKT2：76/86。
- iNKT17：30/36；缺失中包括 `Il17a`、`Ccr6`，`Il17f` 也不在最终 gene universe。

第三步，再把每一套 subtype 平均分数跨全部细胞标准化：

$$
S_{ik}
=
\frac{A_{ik}-\mu_{A_k}}{\sigma_{A_k}}
$$

最终三套 $S_{ik}$ 都具有全局均值 0、标准差 1，可以使用同一 `−3.5` 至 `+3.5 SD` 色标。

图里的每个元素因此是：

- UMAP 上每个点：一个细胞的 $S_{ik}$。
- Violin 中每个观测：一个细胞的 $S_{ik}$，小提琴外形表示分布密度。
- 3×8 heatmap 的每个格子：某个 subtype 行、某个 cluster 列内所有细胞 $S_{ik}$ 的平均值。
- 正值：该 cluster 高于全数据对该 signature 的平均水平。
- 负值：低于全数据平均水平。
- 它不是表达百分比、概率、p-value，也不是已确认的 subtype 身份。

### 6.3.2 附录敏感性分析：raw `scanpy.tl.score_genes`

旧重跑还保留了原先的 `score_genes` 分析：

- 表达来源：10,670-gene `.raw`。
- `n_bins=25`。
- control size = `max(50, signature gene count)`。
- 从全 10,670 genes 中选择表达量匹配的 control genes。

概念上：

$$
\mathrm{RawScore}
=
\operatorname{Mean}(\text{signature genes})
-
\operatorname{Mean}(\text{expression-matched controls})
$$

它适合检查“换一套基因表，空间模式是否仍存在”，但不适合直接比较不同 gene-set 的绝对分数，因为 gene 数和 control genes 均不同。新版 PDF 因此把它移到附录，并给每个来源独立 colorbar；raw argmax 不再决定正文的 subtype call。

## 6.4 从连续 score 得到保守的 cluster call

先在每个 cluster 内计算三个统一分数的均值，再把观测均值最高的 subtype 当作 candidate。随后做 10,000 次、seed 0 的 bootstrap：

1. 在每个 cluster 内按 sample 分层。
2. 每个 sample 保持原细胞数，有放回地重采样细胞。
3. 每次重采样重新计算 candidate mean。
4. 每次都计算 `candidate − max(另外两个 subtype)`；竞争项不是预先固定的 runner-up。
5. 用 2.5% 和 97.5% 分位数形成 95% interval。

Call 规则：

- candidate mean 的 CI 下界 `>0`，且 margin CI 下界 `>0`：`legacy-list iNKT*-enriched`。
- candidate 为正但没有与竞争项稳定分离：`mixed`。
- candidate 没有稳定高于 0：`unclassified`。

当前结果：

| Cluster | iNKT1 mean SD | iNKT2 mean SD | iNKT17 mean SD | 主分析 call |
|---|---:|---:|---:|---|
| c0 | +0.204 | +0.019 | −0.015 | legacy-list iNKT1-enriched |
| c1 | +0.155 | +0.019 | −0.214 | mixed |
| c2 | −0.814 | −1.978 | −0.585 | unclassified |
| c3 | −0.150 | −0.133 | −0.119 | unclassified |
| c4 | +0.725 | −0.195 | −0.256 | legacy-list iNKT1-enriched |
| c5 | −1.885 | +1.147 | +1.537 | legacy-list iNKT17-enriched；生物学叙事仍标 mixed iNKT2/iNKT17 |
| c6 | +0.812 | +0.182 | −0.194 | legacy-list iNKT1-enriched |
| c7 | −3.248 | −1.622 | −0.465 | unclassified |

c5 的 iNKT17 相对最大竞争项 margin 为 `+0.390 SD`，95% interval 约为 `+0.239` 至 `+0.542`。不过 `Il17a`、`Il17f`、`Ccr6` 不可测，且 raw 多来源结果并不稳定，所以这里的严谨表达是：

> c5 在原 PPT 对应的 literature list 上呈 iNKT17-related enrichment，同时有较高 iNKT2 signal；这支持 mixed continuous programs，不足以确认离散 iNKT17 subtype。

c7 的三个统一分数都低于 0，因此旧 raw argmax 产生的 “iNKT2 winner” 被撤回。它现在是 `unclassified`。

这些 bootstrap intervals 只说明本数据集内的 cell-resampling stability。因为每个 tissue×condition 只有一个样本，它们不是 biological-replicate confidence intervals。

## 6.5 T2 对 signatures 的影响

Global T2−Ctrl 变化很小：

- literature iNKT1：−0.0112
- literature iNKT2：+0.0018
- literature iNKT17：+0.0079
- in-house iNKT1：−0.0221
- Wang iNKT1：−0.0049

相对更清楚的是：

- Direct TCR activation：+0.121
- Residency：+0.0339

Direct TCR activation 在三个组织均增加：

- BM：+0.129
- spleen：+0.131
- thymus：+0.099

因此数据更支持：

> T2 与 activation/residency program shift 相关。

而不是：

> T2 导致整个 iNKT population 从一种 subtype 转换为另一种 subtype。

---

# Stage 7：Specific marker validation

对应 PDF 第 6–8、31 页中的 marker claim 检查。

检查的基因包括：

```text
Xcl1, Il4, Ifng, Rorc, Il17a, Il17f,
Gzma, Gzmb, Prf1, Fasl, Klf2, Zbtb16
```

其中 `Il17a` 和 `Il17f` 不在 10,670-gene universe，不能进行当前比较。

每个 marker 计算：

- mean log1p expression
- `expm1` 后的 mean normalized expression
- fraction expressing
- T2−Ctrl mean difference
- descriptive log2 fold change：

$$
\log_2
\frac{\text{T2 mean}+0.001}
{\text{Ctrl mean}+0.001}
$$

这个 log2FC 是描述性均值比，不是 DE 模型估计值。

## Xcl1 是最重要的例子

- Ctrl expressing fraction：67.15%
- T2 expressing fraction：68.11%
- global descriptive log2FC：+0.166
- BM：+0.180
- spleen：+0.241
- thymus：−0.055

所以旧稿里的“Xcl1 case-specific”不成立。正确表述是：

> Xcl1 在 Ctrl 和 T2 中均广泛表达，T2 effect 较小且具有 tissue dependence。

## 稀疏 marker 应怎么读

例如 `Rorc`：

- global log2FC：+0.350
- Ctrl expressing：1.03%
- T2 expressing：1.49%

虽然 fold change 看起来较大，但阳性细胞比例非常低。对于这种基因，必须把：

- fold change
- mean expression
- fraction expressing

放在一起看。否则极少数细胞的轻微变化会被误写为明显 population shift。

---

# Stage 8：Differential expression

对应 PDF 第 10、11、13–26 页。

这里有两套不同的 DE，必须区分。

## 8.1 Global T2 versus Ctrl

这是把三个组织合在一起比较：

```text
7,359 T2 cells vs 8,173 Ctrl cells
```

方法：

- Wilcoxon
- T2 versus Ctrl
- 只测试 3,000 HVGs
- BH correction

结果：

- 345 genes：FDR ≤ 0.05
- 132 genes：同时满足 FDR ≤ 0.05 和 |logFC| ≥ 0.25
- 其中 57 up、75 down

旧 global list 有 39 genes：

- 只有 20 个进入当前 3,000-HVG test universe。
- 20/20 均达到当前 FDR。
- 16/20 同时达到 |logFC| ≥ 0.25。
- 另外 19 个是“not tested”，不是“failed”。

Global comparison 的问题是 tissue 被混在一起。结果可能同时反映：

- 真实的 condition response
- tissue composition
- 不同组织中 effect direction 的混合
- sample/library differences

所以 tissue-stratified DE 比 pooled global DE 更重要。

另外，notebook 输出的 `top500` 是正 Wilcoxon score 端的前 500 个，并不是完整双向 DEG 表。

## 8.2 Tissue 和 cluster×tissue DE

每个组织单独比较 T2 versus Ctrl，并对每个 observed cluster×tissue 单元重复比较。

方法：

- expression：`.raw` 中的 10,670-gene log-normalized matrix
- 全部 10,670 genes
- Wilcoxon
- `tie_correct=True`
- 每个 comparison 内对 10,670 genes 做 BH
- positive logFC = T2-up
- reporting threshold：
  - FDR ≤ 0.05
  - |logFC| ≥ 0.25
- 每组至少 20 cells

一共规划 26 个 units：

- 3 个 tissue units
- 23 个 observed cluster×tissue units

结果：

- 19 completed
- 7 skipped
- 0 failed

## 8.3 Tissue-level results

| Tissue | T2/Ctrl cells | Significant DEGs | Up/down |
|---|---:|---:|---:|
| BM | 2,731 / 3,379 | 718 | 350 / 368 |
| Spleen | 3,371 / 3,422 | 864 | 356 / 508 |
| Thymus | 1,257 / 1,372 | 352 | 195 / 157 |

旧 tissue DEG 恢复：

- BM：35/39
- spleen：44/50
- thymus：31/34

在假设旧 PPT 正 logFC 代表 T2-up 的前提下，恢复的基因方向全部一致。

这是整个 downstream 中最强的 reproduction evidence 之一。

## 8.4 主要 cluster×tissue results

- c0 BM：592 DEGs，274 up / 318 down
- c3 spleen：735，275 / 460
- c6 thymus：304，140 / 164
- c4 BM：83，56 / 27
- c4 spleen：45，33 / 12
- c5 BM：39，26 / 13
- c5 spleen：73，44 / 29

## 8.5 “zero DEG”和“skipped”完全不同

通过测试但没有基因达到阈值：

- c1 BM
- c1 spleen
- c2 spleen
- c2 thymus
- c5 thymus
- c7 thymus

这表示：

> 测试执行了，但在当前效应阈值和 FDR 下没有基因通过。

因某组少于 20 cells 而跳过：

- c0 thymus
- c1 thymus
- c3 thymus
- c4 thymus
- c6 BM
- c6 spleen
- c7 spleen

这表示：

> 没有进行统计测试。

因此：

```text
0 DEG ≠ skipped
skipped ≠ no biology
0 DEG ≠ 两组完全相同
```

## 8.6 最重要的统计限制

这些 p-value 把细胞作为 observations，但 condition 和单个 sample 在每个 tissue 中完全绑定。

因此细胞数很大时，即使很小的效应也可能得到极小 p-value。阅读优先级应该是：

1. effect direction
2. effect size
3. fraction expressing
4. 是否跨组织/分析层面一致
5. 最后才是 cell-level p-value

---

# Stage 9：从旧 PPT 提取结果并进行 response-program matching

对应 PDF 第 13–29 页。

## 9.1 旧结果如何恢复

脚本读取：

- 原 PPT 的 XML tables
- gene-list workbook
- marker workbook

实际提取：

- 18 个 sheets
- 27 个 result tables
- 298 个 legacy DEG rows
- 625 个 legacy pathway rows
- 12 套 curated signatures

无法可靠识别的 scope 保留为 unresolved，没有强行猜测。

## 9.2 如何比较 legacy cluster 和 current cluster

对于一个 legacy cluster response：

1. 先固定 tissue。
2. 在同一 tissue 内遍历所有完成 DE 的 current clusters。
3. 将 legacy DEG set 与 current significant DEG set 比较。
4. 计算：

$$
Recovery =
\frac{|Legacy \cap Current|}
{|Legacy|}
$$

$$
Jaccard =
\frac{|Legacy \cap Current|}
{|Legacy \cup Current|}
$$

以及：

- shared gene count
- sign agreement
- current DEG set size

自动 best match 的优先级是：

1. shared gene count 最大
2. recovery 最大
3. Jaccard 最大
4. current DEG set 更小

因此大而广的 DEG set 更容易成为 maximum-coverage match。

如果另一个较小 cluster 的 Jaccard 或特异性更高，它被称为 compact candidate。两者不一致通常意味着：

> 一个旧 response program 在当前聚类中被拆分或合并了。

## 9.3 最强的几个匹配

- legacy c1 thymus → current c6 thymus  
  26/30 legacy genes recovered；26/26 apparent signs concordant

- legacy c5 BM → current c0 BM  
  53/61 recovered；全部共享基因方向一致

- legacy c8 spleen → current c3 spleen  
  73/84 recovered；全部共享基因方向一致

其他例子：

- legacy c2 BM → current c4：17/18
- legacy c2 spleen：
  - maximum coverage c3：19/27
  - compact candidate c4：18/27
- legacy c6 spleen：
  - compact candidate c5：22/29
- legacy c7 BM：
  - c0 和 c5 都只恢复 3/4，映射不唯一
- legacy c9 spleen → c3：19/28

正确表述是“program recovered in…”，不是“cluster identity confirmed”。

## 9.4 Legacy pathway driver retention

旧 pathway 表中的 driver genes 被拆出来，与 best-match current DEGs 比较。

结果包括：

- legacy c1 thymus：8/8
- c2 BM：4/4
- c2 spleen：7/7
- c5 BM：15/15
- c6 spleen：5/6
- c8 spleen：20/20
- c9 BM：2/2
- c9 spleen：4/4

这说明旧结果中的核心 HSP、IEG、mitochondrial 等 driver modules 被高度保留。

但这不是 pathway enrichment 的重算，只是 driver-level retention。

---

# Stage 10：DEG overlap 和 ORA

对应 PDF 第 27、32、33 页。

## 10.1 Tissue DEG 与 tissue-anchor cluster overlap

每个显著 DEG set 被分成：

- T2-up
- T2-down

然后将整个 tissue 与该 tissue 的主要 anchor cluster 比较。

结果：

| Comparison | Anchor DEGs | Shared | Up Jaccard | Down Jaccard |
|---|---:|---:|---:|---:|
| BM vs c0-BM | 592 | 489 | 0.629 | 0.566 |
| Spleen vs c3-spleen | 735 | 600 | 0.656 | 0.566 |
| Thymus vs c6-thymus | 304 | 261 | 0.667 | 0.655 |

解释：

- tissue-level response 很大一部分由 dominant cluster 承载。
- 但仍有 tissue-only 或 cluster-only genes。
- tissue-only 可能来自其他细胞状态。
- cluster-only 可能是状态特异信号，也可能来自 power/threshold 差异。

UpSet-like plots 显示不同 cluster 的 DEG 共享模式。它适合回答：

> 哪些 response genes 是组织共同核心，哪些只出现在特定 state？

但 UpSet intersection size 本身不是显著性检验。

## 10.2 当前 ORA

本轮没有提供原始 PAGER/Reactome/KEGG GMT，因此只对本地 12 套 curated signatures 做 ORA。

对于每个 completed unit、每个方向：

1. 选择 FDR ≤ 0.05 且 |logFC| ≥ 0.25 的 genes。
2. 背景为全部 10,670 个 tested genes。
3. 对每个 signature 计算 one-sided hypergeometric enrichment。
4. BH correction 在每个 unit × direction × library 内分别执行。

一共：

```text
19 units × 2 directions × 12 terms = 456 tests
```

主要结果：

- BM T2-up：
  - Direct TCR activation FDR ≈ 5.8×10⁻⁶
  - Residency ≈ 1.0×10⁻⁵
- Spleen T2-up：
  - Direct TCR activation ≈ 6.4×10⁻⁶
  - Residency ≈ 1.2×10⁻⁵
- Thymus T2-up：
  - iNKT17-list overlap
  - Direct TCR
  - Residency
- Thymus T2-down：
  - 多套 iNKT1 signatures

## 10.3 ORA 应该怎么解释

ORA 只表示：

> 当前 DEG list 与某个 gene set 的重叠，高于相同背景下随机重叠的预期。

它不表示：

- pathway 被直接测量。
- pathway 整体激活。
- pathway 是因果机制。
- 某个 subtype 的细胞数量增加。

尤其 thymus T2-up 出现 iNKT17 ORA，只能说：

> thymus T2-up DEG list 与某些 iNKT17 signatures 有重叠。

不能说：

> iNKT17 population 增加。

因为 c5 只有统一 literature-list 的连续 iNKT17-related enrichment；它同时具有较高 iNKT2 score，并缺少 Il17a/Il17f/Ccr6 等关键可测证据，尚不足以建立稳定、独立的 iNKT17 subtype 或 population expansion。

## 10.4 PaGER-scFGA 在这套 reproduction 中的位置

[PaGER-scFGA paper](/home/zzz0054/bio3/docs/references/pager-scFGA.pdf) 对应的是旧 pathway 分析的背景方法之一，但当前 reproduction 并没有完整重新运行同版本的 PaGER/Reactome/KEGG 流程。

缺少的信息包括：

- 原数据库版本
- species/orthology mapping
- background universe
- 原始 correction scope
- 部分 pathway 参数

因此当前结果能说：

> legacy pathway driver/module was retained。

不能说：

> 原 PAGER pathway analysis 已被完整复现。

旧 PPT 中诸如 “Prion disease”“Measles”“Estrogen signaling” 的名字通常来自共享的 HSP、IEG、mitochondrial genes，不能按疾病字面解释。

---

# Stage 11：PAGA、Diffusion Map 和 DPT pseudotime

对应 PDF 第 5、30 页。

这是当前结果里最需要谨慎解释的一部分。

## 11.1 PAGA

执行：

```python
sc.tl.paga(adata, groups="leiden_res_0_5")
```

PAGA 使用同一个 15-neighbor graph，计算 cluster 之间的 graph connectivity。

较强的边包括：

- c4–c6：0.393
- c0–c1：0.223
- c4–c5：0.208
- c1–c3：0.175
- c2–c7：0.167

PAGA 图可以说：

> 这些 expression states 在 neighbor graph 中连接较强。

不能说：

> 箭头表示细胞从一个 cluster 发育到另一个 cluster。

这里的 PAGA 本身是无方向 connectivity。

## 11.2 Diffusion map

计算 15 个 diffusion components，用于描述 graph 上较连续的变化结构。

它比 UMAP 更适合探索连续状态，但仍然只是表达相似性空间，不是时间测量。

## 11.3 Root selection

候选 root markers：

```text
Cd27, Itgam, Il2rb, Klrb1c, Tcrb, Zbtb16
```

实际可用：

```text
Cd27, Il2rb, Klrb1c, Zbtb16
```

由于 `Itgam` 不存在，代码 fallback 实际把 root score 简化成了 cluster mean `Cd27`。

c7 的 Cd27 mean 最高，因此选择 c7。

接着 root cell 不是：

- c7 的 medoid
- marker score 最高的细胞
- 一个预先验证的 immature thymus cell

而是代码顺序中的第一个 c7 cell：

```text
Ctrl_Spleen_ACGATGTAGCTGACCC-1
```

虽然 c7 有 89/92 cells 来自 thymus，但恰好第一个 c7 cell 是 Ctrl spleen。这使 root choice 很脆弱。

## 11.4 DPT

执行：

```python
sc.tl.dpt(n_dcs=10)
```

全部 15,532 cells 都有 finite DPT，范围 0–1。

Cluster median：

- c7：0.103
- c0：0.989
- c1：0.992
- c2：0.977
- c3：0.989
- c4：0.989
- c5：0.981
- c6：0.992

这不是漂亮的连续 early→intermediate→late 分布。它几乎是：

```text
c7 = low DPT
其余所有 clusters = 接近 1
```

所以最诚实的描述是：

> DPT primarily separates the selected c7 root state from the remaining expression states.

而不是：

> DPT reconstructs a complete developmental trajectory.

T2−Ctrl 的 mean DPT 差异也极小：

- BM：约 +0.00029
- spleen：约 +0.00041
- thymus：约 −0.00289

因此没有证据表明 T2 在这个 DPT 上引起明显的全局进展或倒退。

## 11.5 Extended trajectory analysis

扩展脚本：

- 使用 16 个 measured programs
- 共 275 个 genes
- 对每个 gene 跨所有细胞做 z-score
- program score = member-gene z-scores 的平均
- 分析 program/gene 与 DPT 的 Spearman correlation
- 用 BH 校正
- 将 DPT 分成最多 20 个 quantile bins
- 选 absolute rho 最大的 40 个 dynamic genes

问题是：

- 15,532 cells 会让非常小的相关性也得到很小 p-value。
- DPT 本身高度饱和。
- tissue structure 可同时驱动 DPT 和 gene expression。

所以应该优先看 rho 的大小和曲线形状，而不是只看 FDR。

所谓 alternative root validation 仍选择 c7，并再次使用第一个 c7 cell，因此 alternative DPT 与原 DPT 完全相同，Spearman rho=1。它只能说明两个 cluster-selection heuristic 都选择 c7，不能算独立的 root sensitivity test。

---

# Stage 12：把所有结果整合成最终生物学叙事

建议按证据强度分三层。

## 强证据

1. Legacy QC cohort 在计数层面精确复建。
2. Tissue identity 是当前数据结构的主轴。
3. BM、spleen、thymus 的 legacy DEG cores 高度恢复。
4. 主要 legacy response drivers 高度保留。
5. BM/spleen/thymus 的 dominant current clusters 承载了大部分 tissue response。

## 中等强度证据

1. 统一 literature-list 口径支持 c0、c4、c6 的相对 iNKT1 enrichment；标准化分数本身不能证明全数据的 subtype abundance。
2. c5 的 iNKT17-related score 高于 iNKT2，但两者都高；结合核心 marker 缺失，保守解释为 mixed iNKT2/iNKT17 continuous programs。
3. c7 的三个统一 scores 均低于全局均值，因此标为 unclassified；旧 raw `score_genes` 的 iNKT2 argmax 仅作为附录方法敏感性结果。
4. T2 更明显关联 Direct TCR activation、residency、IEG/HSP 和 mitochondrial/proteostasis remodeling，而不是全面 subtype conversion。
5. `Xcl1` 广泛存在于 Ctrl 和 T2，不是 case-specific marker。

## 尚未建立的结论

1. 旧 cluster 与新 cluster 的一一身份对应。
2. c5 已构成稳定、独立的 iNKT17 subtype；当前只建立了 literature-list continuous-program enrichment。
3. 没有任何 iNKT17 cells。
4. 旧 PAGER/Reactome/KEGG terms 已同库复现。
5. “Prion/Measles/Estrogen”等代表字面疾病机制。
6. PAGA/DPT 证明了发育谱系、分支方向或 ancestry。
7. Cell-level p-value 证明了可推广的动物总体 T2 effect。

---

# 阅读每一张结果图时，固定问这六个问题

1. **这张图比较的是哪些细胞？**  
   Global、tissue，还是 cluster×tissue？

2. **使用的 gene universe 是什么？**  
   3,000 HVGs，还是完整 10,670 genes？

3. **颜色或数值是什么？**  
   Raw signature score、descriptive log2FC、Wilcoxon logFC，还是 fraction expressing？

4. **“没有结果”是哪一种？**  
   Not tested、tested but nonsignificant，还是 skipped？

5. **统计单位是什么？**  
   细胞还是独立动物？这里的 DE p-value 是细胞级探索结果。

6. **Cluster 名称属于哪个 namespace？**  
   `legacy c1–c11` 还是 `current c0–c7`？

如果这六个问题能回答清楚，基本就不会误读这套结果。

## 主要文件

- [Detailed reproduction PDF](/home/zzz0054/bio3/output/iNKT_reproduction_deck/iNKT_legacy_QC_reproduction_detailed.pdf)
- [Editable PPTX](/home/zzz0054/bio3/output/iNKT_reproduction_deck/iNKT_legacy_QC_reproduction_detailed.pptx)
- [逐页结果说明](/home/zzz0054/bio3/docs/audits/inkt_legacy_ppt_vs_legacy_qc_rerun_content.md)
- [Processed H5AD](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/preprocess/inkt_scanpy_tutorial_processed.h5ad)
- [Legacy/current best matches](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/comparison_to_legacy_ppt/legacy_ppt_de_best_matches.csv)
- [Current curated ORA](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/de_pathway/overlap_pathway/current_offline_gene_set_enrichment.csv)
- [Trajectory summary](/home/zzz0054/bio3/output/iNKT_legacy_ppt_qc_runs/20260818_125231/extended/trajectory/summary.json)
