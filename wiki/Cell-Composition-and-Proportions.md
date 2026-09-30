# 细胞组成与比例（Cell Composition, Proportions）

> 比例图有两种分母：条件内的簇频率，和簇内的条件占比，不能混读。

## 定义

细胞组成（cell composition）看每个簇占某个条件细胞的百分比。两个常混淆的量：n(cluster, condition)/N(condition)（条件内簇频率）与 n(condition, cluster)/N(cluster)（簇内条件占比，T2 share）。

## CS 类比

≈ 条件概率 P(簇|条件) 与 P(条件|簇) 的区别。

## 常见误读与注意

- 单细胞捕获比例不是组织里的绝对数量；没有动物层面的误差棒或显著性。
- 不是都以 50% 为基线：T2 总体占比约 45%（转录稿 00:43:57），骨髓和脾脏基线不同。

## 在本项目中

- 比例图 37 页版第 8–10 页（`iNKT_by_date/2026-09-20/package/slides.json`）；R12：27 个组织×簇核对（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md`）。
- C5-1 在脾脏：Ctrl 239/3,422（7.0%）→ T2 345/3,371（10.2%），簇内 T2 占比 59.1%（我们对 `iNKT_by_date/2026-09-19/results/tables/cell_metadata_scores.csv.gz` 计数得到）。
- 骨髓 C3：Ctrl 236 → T2 115（T2 占比 32.8%，同上计数）。
- 9/24 会议：ZERU 指出 C5-1 在 T2 脾脏细胞中占更大份额（转录稿 00:44:36）。

## 相关概念

- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
- [[伪重复与 n=1（Pseudoreplication, Why Cell-level p-values Are Not Animal-level Evidence）|Pseudoreplication-and-n-equals-1]] — 把同一个 pool 里的几千个细胞当作独立重复，会让 p 值虚高；本项目每组只有 n=1 个 pool。
- [[组织：骨髓、脾脏、胸腺（Bone Marrow, Spleen, Thymus）|Tissues-BM-Spleen-Thymus]] — 骨髓造血、胸腺培育 T 细胞、脾脏是血液过滤与免疫应答场所；三者是不同的“微环境”。
