# 2026-09-19｜iNKT 细胞毒性与 GO 分析

[结果 PDF](presentation/iNKT_cytotoxicity_GO_20260919.pdf) · [PPTX（含中文讲稿）](presentation/iNKT_cytotoxicity_GO_20260919.pptx) · [图总览](presentation/contact_sheet.jpg)

## 本轮实际完成

- 最新分析输入：9月5日对象，15,532细胞 × 10,670基因；9月15日为展示更新，没有更新表达数据。保留原9个refined cluster，并补充三个组织的稳定分群敏感性。
- 新计算：PAGER核心与E1模块、背景校正评分、组织×群条件比较、五次QC匹配、组织内cluster-vs-rest差异表达、GO BP/MF/CC全库ORA、局部DPT与根细胞敏感性。
- 条件DE复用9月5日完整结果，核对基因集及各组细胞数，复制到本目录并记录SHA256。所有GO富集为本轮新计算。
- 输出 19 页论文式图表与PPTX；PPTX采用图像页，保留SVG和绘图源码可修改。代码和结果是真实文件，不是旧目录符号链接。

## 当前主要数值

**优先关注：骨髓C4杀伤核心模块T2升高（差值+0.086，BH q=0.025）；脾脏C3小幅降低（−0.023，q=0.0015）。两者五次QC匹配均同向，并在主要重叠稳定分群保持方向。**

**签名依赖性必须一并汇报：加Ctla2a后的E1代表模块在这两群均不显著（q=0.93、0.53）。杀伤相关GO条目没有通过q≤0.05，最小q=0.632。因此结论限于局部核心效应转录信号，不是整个细胞毒性通路或杀伤能力的实验证据。**

条件GO较强信号包括骨髓C0与脾脏C3的T2上调基因富集氧化磷酸化；多群另见核糖体/翻译相关条目。诸如“synapse translation”的条目可由核糖体基因驱动，不能据名称推断形成神经突触。

| 组织/原分群 | T2−Ctrl核心均值 | 细胞层面BH q | QC同向次数 | 背景校正评分同向 |
|---|---:|---:|---:|---|
| Spl/C1 | -0.131 | 0.0588 | NA（匹配后细胞不足） | 是 |
| BM/C2 | -0.100 | 0.222 | 5/5 | 是 |
| BM/C4 | +0.086 | 0.0251 | 5/5 | 是 |
| Thy/C2 | +0.067 | 0.493 | 4/5 | 是 |
| Spl/C4 | +0.061 | 0.3 | 5/5 | 是 |
| Spl/C0 | +0.059 | 0.463 | 4/5 | 是 |
| Thy/C6 | -0.049 | 0.399 | 5/5 | 是 |
| Spl/C5-2 | -0.035 | 0.335 | NA（匹配后细胞不足） | 是 |

按绝对效应排序，并非独立生物学重复验证。NA表示QC匹配后每组不足20细胞，不是0次同向。模块均值是log1p归一化表达，不是杀伤率，也不是log2FC。
总共 93 条合格模块对照，11 条通过统一BH q≤0.05；包含原分群、稳定分群和三个重叠基因集，不能当作独立发现数量。

## GO 方法与结果入口

- 官方GO版本：`releases/2026-07-26`。MGI GAF头和下载哈希见 sources/GO_provenance.json 与 download_manifest.json。
- 仅小鼠注释；排除NOT和ND，保留IEA；只沿is_a/part_of向上继承，不沿regulates继承。精确symbol与明确Ensembl→MGI别名匹配。
- 每个对照的背景为实际检出的基因∩该GO分支已注释基因；集合大小5–500；上调/下调分开。基因筛选：gene BH≤0.05，|log2FC|≥0.25。
- 每个对照×方向×GO分支在全部可检验term中BH校正，包含零命中term；另给出每种分析跨所有对照的q_analysis_global。
- 杀伤相关GO补充核查按term名称匹配cytotox、cell killing、granzyme、exocytosis、lymphocyte mediated；全部显著及不显著条目保留在GO_cytotoxicity_targeted_audit.csv。
- 图中富集不表示通路激活；不同GO条目可能由同一组基因驱动。主图取显著、有≥2命中的代表条目，Jaccard≥0.75去冗余，完整结果保留。

[完整GO结果](results/tables/GO_ORA_all.csv.gz) · [显著GO结果](results/tables/GO_ORA_significant.csv) · [背景与查询覆盖](results/tables/GO_contrast_coverage.csv) · [绘图GO条目](results/tables/displayed_GO_terms.csv)

## 拟时序方法与限制

- 每组织在9月5日的局部邻居图重新计算diffusion map，避免继承全局坐标。取最大连通分量；未纳入细胞数在root audit中逐项列出。
- 原前体面板Cd24a/Egr2/Hivep3中，Cd24a不在当前10,670基因集，实际以Egr2/Hivep3两基因评分最高的局部群选根。根的前体身份支持有限；三个根细胞做敏感性，不能把此排序当谱系。
- 两条件共享同一组织DPT的20个区间，每bin至少20细胞；仅连续有支持的区间平滑。胸腺细胞集中在不连续的DPT区间，没有可靠连续曲线；骨髓/脾脏根敏感性最低相关约0.76。跨组织不能比较速度。
- 这是转录状态排序，不是RNA velocity、分化速率或已验证谱系。

## 可复现与目录

- `code/analyze.py`：评分、GO、轨迹；`code/presentation.py`：图与报告；`code/verify.py`：数值及文件验证；`code/test_analysis.py`：统计单元测试。
- 从仓库根目录：`bash iNKT_by_date/2026-09-19/code/run.sh`。长任务在tmux运行，CPU-only，GPU全部禁用。
- `sources/` 保存GO原始注释、版本、输入哈希；`results/de/`含新marker DE和复用的条件DE；`results/figures/`含PNG/SVG；`presentation/`含PDF/PPTX和讲稿索引。
- 重跑需仓库原始输入及现有`.venv`，本日期目录未重复拷贝大体积表达矩阵。原分析文件未改写。

## 解释边界

- 每组织×条件只有1个sample标签；据9月24日会议Rob说明，每个标签由3只小鼠混样（pooled），即每组n=1个pool，仍没有可识别的独立动物重复。细胞层面P值和QC匹配不能解决生物重复缺失。
- 这是iNKT免疫杀伤相关转录分析，与图中药物肝/肾毒性项目不同。
- 采用两篇文章的方法与presentation，不声称完整重现PAGER平台、PPI、GAFA随机森林/SCORPION。

## 方法来源

- [PAGER-scFGA](https://doi.org/10.3389/fbinf.2024.1336135)：Figure 5功能模块、GO细胞组分和拟时序功能曲线。
- [Borra et al. 2026](https://doi.org/10.3390/biology15070588)：§3.5/Figure 4c亚群功能富集展示。
- [GO官方本体](https://geneontology.org/docs/download-ontology/)与[官方注释](https://geneontology.org/docs/download-go-annotations/)。

## 最终验证

5项统计单元测试与318项独立数值/文件检查通过；PDF与PPTX均19页，已完成图形目检。阶段退出码0；CPU运行，任务tmux会话已结束。详见[验证记录](results/verification.json)与[阶段状态](results/workflow_status.json)。
