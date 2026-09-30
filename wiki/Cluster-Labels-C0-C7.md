# 本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）

> 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。

## 定义

当前分析的簇编号 C0–C7 来自 Leiden 分辨率 0.5；C5 又被重聚类拆为 C5-1 和 C5-2，得到 9 个 refined cluster。旧 PPT 有 11 个簇（legacy c1–c11），编号空间不同。

## CS 类比

≈ 两个不同版本的枚举值：数字相同并不代表语义相同。

## 常见误读与注意

- 不能写 legacy c1 = current c6；应写“legacy c1 的胸腺响应程序在当前 c6 中得到恢复”（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.5）。

## 在本项目中

- 簇的细胞数与组织分布（按 `iNKT_by_date/2026-09-19/results/tables/cell_metadata_scores.csv.gz` 的 `cluster_c5_split_20260830` 列，我们对行计数）：

| 簇 | 细胞数 | 骨髓 | 脾脏 | 胸腺 | 主要组织 |
|---|---:|---:|---:|---:|---|
| C0 | 4,828 | 4,623 | 185 | 20 | 骨髓 |
| C1 | 162 | 78 | 81 | 3 | 混合 |
| C2 | 378 | 145 | 119 | 114 | 混合 |
| C3 | 5,677 | 351 | 5,319 | 7 | 脾脏 |
| C4 | 887 | 467 | 404 | 16 | 骨髓/脾脏 |
| C5-1 | 812 | 209 | 584 | 19 | 脾脏 |
| C5-2 | 317 | 200 | 70 | 47 | 骨髓 |
| C6 | 2,379 | 37 | 28 | 2,314 | 胸腺 |
| C7 | 92 | 0 | 3 | 89 | 胸腺 |

- 所有簇都同时含 Ctrl 与 T2，没有 condition 独占簇（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 5.4）。
- 主要候选簇：骨髓 C0、脾脏 C3（AP-1）；骨髓 C4、脾脏 C3（细胞毒性/CCT）；C5-2（iNKT17 相关）；C5-1（T2 脾脏中占比更高）。
- 同一目录 `iNKT_by_date/2026-09-19/results/de/` 按 `cluster_tissue__C0__bone_marrow.csv.gz` 命名保存组织×簇 DE。

## 相关概念

- [[Leiden 聚类（Leiden Clustering）|Leiden-Clustering]] — Leiden 是在近邻图上做社区发现的算法；分辨率越高群越多。
- [[细胞类型与细胞状态（Cell Type vs Cell State）|Cell-Type-vs-Cell-State]] — “类型”是相对稳定的身份，“状态”是同一类细胞当前的活动模式；本项目的 C0…C7 是算法分出的表达状态簇，不是已验证的细胞类型。
- [[细胞组成与比例（Cell Composition, Proportions）|Cell-Composition-and-Proportions]] — 比例图有两种分母：条件内的簇频率，和簇内的条件占比，不能混读。
- [[目前的主要发现与证据等级（Key Findings so Far）|Key-Findings]] — 结论与证据等级汇总。
