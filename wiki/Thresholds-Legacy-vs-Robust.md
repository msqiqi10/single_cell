# 两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）

> legacy = 名义 p≤0.05 且 FC≥1.5；robust = FDR≤0.05 且 |log2FC|≥0.25。两套同时改变校正和效应门槛，不能称“更严/更松”。

## 定义

旧 PPT（第 23–28 页 KEGG 图）用 nominal P≤0.05、线性 FC≥1.5（即 T2 上调 log2FC≥log2(1.5)）；当前主分析用基因 BH FDR≤0.05、|log2FC|≥0.25，T2-up 与 Ctrl-up 分开。

## CS 类比

≈ 两套不同的告警规则：一套按 raw 阈值，一套按误报率加最小效应门槛。

## 常见误读与注意

- 两套规则得到的入选基因数不能直接比较生物学强弱。

## 在本项目中

- 输入数（T2-up）：legacy 规则骨髓 C0 166、脾脏 C5-1 286、胸腺 C6 188；robust 规则 274、42、140（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 3 节；`iNKT_by_date/2026-09-25/results/README.md` 第 2 节）。
- 旧规则下 Top 10 图只覆盖 22、24、16 个不同基因；旧六个方向的全库 KEGG ORA 均无 FDR≤0.05（同 README）。
- KEGG MAPK：paper 规则约 0.078，robust 规则约 0.627（同上）。
- 37 页版第 23–25 页为 “Paper pathways”，第 26–28 页为 “Current top pathways”（`iNKT_by_date/2026-09-20/package/slides.json`）。

## 相关概念

- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
- [[对数倍数变化（Log2 Fold Change, log2FC）|Log2-Fold-Change]] — log2FC 衡量效应大小：+1 表示 T2 约为 Ctrl 的 2 倍。
- [[KEGG 通路数据库（KEGG Pathway Database）|KEGG]] — KEGG 是人工整理的代谢、信号和疾病通路图集合；条目名常是疾病名。
- [[过表征分析（Over-Representation Analysis, ORA; hypergeometric test）|Over-Representation-Analysis]] — ORA 问：入选基因列表里落在某个基因集的数量，是否比随机预期多。
