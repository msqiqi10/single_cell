# 图中的富集分数：−ln p 与 −log10 q（Enrichment Scores in Figures）

> 不同图用了不同的“分数”：旧图 −ln(名义 p)，新 GO 气泡图 −log10(BH q)，GSEA 图用带符号 NES。

## 定义

为了让极小的 p 值可读，图上常用 −log 变换：值越大越显著。−ln p 是自然对数，−log10 q 是常用对数，二者差常数倍 ln(10)≈2.303；q 与 p 也不是同一个量。

## CS 类比

≈ 把跨越多个数量级的延迟画成对数刻度。

## 常见误读与注意

- 不同变换的数字不能直接横向比较；也不能把跨库、跨方法的分数并到一个色标里当作同一含义（图注：Scale does not remove method/library differences）。

## 在本项目中

- 37 页版第 23–25 页：Score = −ln(nominal enrichment P)，与 Blood Table S3 核对；星号表示通路 BH FDR≤0.05（`iNKT_by_date/2026-09-20/package/notes/presenter_guide.zh.md`）。
- 19 页版 GO 气泡图：点面积 = gene ratio（查询基因比例），颜色 = −log10(BH q)，上限 10（`iNKT_by_date/2026-09-19/code/presentation.py` 第 67–68 行）。
- Rob 在会上指出第 31 页热图的命名与显著性标记需要处理（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 6 节）。
- 9/29：Yue 建议网络节点大小 = −log10 p（越显著越大），颜色 = 功能组（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt` 约 16:28:05–16:29:02）。

## 相关概念

- [[基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）|GSEA-and-NES]] — GSEA 不需要先选“显著基因”，而是把全部基因按变化排序，看基因集是否集中在顶端或底端；NES 带方向。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
- [[交互式知识图谱 HTML（Knowledge Graph, KG; interactive HTML）|Knowledge-Graph-HTML]] — Yue 想要的最终交付：可点击的网络，节点大小 = 显著性，颜色 = 功能组，点开能看到基因。
- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
