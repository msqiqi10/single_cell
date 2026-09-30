# iNKT 亚型 iNKT1 / iNKT2 / iNKT17 及其标志基因（iNKT Subsets and Markers）

> iNKT1：Tbx21/Ifng；iNKT2：Gata3/Il4；iNKT17：Rorc/Il23r/Ccr6/Il17。

## 定义

iNKT 细胞按主要分泌的细胞因子和主控转录因子分为三类：iNKT1（Tbx21/T-bet，分泌 IFN-γ）、iNKT2（Gata3，分泌 IL-4）、iNKT17（Rorc/RORγt，Il23r、Ccr6，分泌 IL-17）。

## CS 类比

三种子类 ≈ 三个继承同一基类的实现，靠“主控 TF 开关 + 输出接口（细胞因子）”区分。

## 常见误读与注意

- 不同来源的同类 signature 重叠很低：iNKT1 Jaccard 0.20–0.438，iNKT2 0.057–0.191，iNKT17 0.216–0.386，不存在唯一的“iNKT1 score”（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 6.2）。
- 基因 mRNA 检出 ≠ 细胞因子分泌：Il17a 在保留细胞中 0 个检出，Il17f 仅 2 个，不能推断分泌能力。

## 在本项目中

- 37 页版第 5、6、7 页分别为 iNKT1/iNKT2/iNKT17 marker program（`iNKT_by_date/2026-09-20/package/slides.json` 中 series=old 的前 37 项）。
- C5-2 在 iNKT17 综合分、Rorc、Il23r 上突出（转录稿约 00:42:03）；Ccr6 因 `min_cells=100` 被过滤，从原始矩阵找回：原始矩阵 92 个细胞检出，保留细胞中 60 个，其中 54 个在 C5-2（`iNKT_by_date/2026-09-25/results/README.md` 第 4 节；`iNKT_by_date/2026-09-25/results/tables/full_feature_marker_status.csv`）。
- c5 的 iNKT17 分数 +1.537 SD，iNKT2 +1.147 SD，iNKT1 −1.885 SD（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md` 6.4 表；该表用旧 c5 未拆分）。
- Rob：iNKT17 特征在他们的台面实验中与 IL-17 分泌相关，是优先方向（转录稿 01:22:19、01:28:38）。

## 相关概念

- [[iNKT 细胞（invariant Natural Killer T cell）|iNKT-Cell]] — iNKT 是表达“不变” TCR、同时带 NK 样特征的 T 细胞亚群，能快速分泌大量细胞因子。
- [[细胞因子（Cytokine: IFN-γ, IL-4, IL-17）|Cytokine]] — 细胞因子是免疫细胞之间传递指令的分泌蛋白。
- [[过滤阈值与 min_cells=100（Filtering Thresholds, why Ccr6 was lost）|Filtering-Thresholds]] — 先删低频基因（min_cells=100），再删异常细胞；这一步让 Ccr6 等稀有 marker 掉出了主对象。
- [[标志基因与模块分数（Marker Genes, Module Scores）|Marker-Genes-and-Module-Scores]] — marker 用来给簇命名，模块分数把一组基因合成一个数。
- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
