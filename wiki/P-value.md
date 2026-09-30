# p 值（p-value）

> p 值是“如果两组其实没差别，看到这么大差异的概率”，不是效应大小，也不是结论正确的概率。

## 定义

p 值（p-value）是零假设下出现当前或更极端结果的概率。细胞数越多，即使很小的效应也能得到极小的 p 值。

## CS 类比

≈ 单元测试的 flaky 概率阈值：p 小只表示“偶然性不大”，不表示“变化很重要”。

## 常见误读与注意

- 15,532 个细胞会让极小的相关性也 p 很小；应先看方向、效应大小、阳性比例，再看 p（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 8.6）。
- 名义 p（nominal）不校正多重比较。

## 在本项目中

- 旧 PPT 规则用名义 p≤0.05（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 3 节）。
- 所有 p 值都是细胞层面的探索性结果（`iNKT_by_date/2026-09-19/README.md` “解释边界”）。

## 相关概念

- [[多重检验、BH 与 FDR / q 值（Multiple Testing, Benjamini-Hochberg, FDR, q-value）|Multiple-Testing-BH-FDR]] — 检验了一万个基因，必然有一些偶然“显著”；BH 把 p 值调整为 FDR（q 值）来控制假发现比例。
- [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] — 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。
- [[两套筛选阈值：legacy 与 robust（Thresholds: Legacy vs Robust）|Thresholds-Legacy-vs-Robust]] — legacy = 名义 p≤0.05 且 FC≥1.5；robust = FDR≤0.05 且 |log2FC|≥0.25。两套同时改变校正和效应门槛，不能称“更严/更松”。
