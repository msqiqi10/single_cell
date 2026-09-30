# 核糖体与翻译（Ribosome and Translation）

> 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。

## 定义

翻译（translation）是按 mRNA 序列合成蛋白质的过程，由核糖体（ribosome）完成；核糖体由大、小亚基的几十种核糖体蛋白（Rpl* 大亚基，Rps* 小亚基）和 rRNA 组成。

## CS 类比

核糖体 ≈ 编译器/解释器：所有蛋白质都要经过它。它的组件基因数量多、表达高，因此在任何富集分析里都容易“统治”结果。

## 常见误读与注意

- 核糖体基因成员多，会让 “translation”“ribosome biogenesis” 以及名字很怪的条目（如 “translation at synapse”）反复出现，它们共享同一批 Rpl/Rps 基因（见 [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]]）。
- 表达升高不等于蛋白合成速率升高。

## 在本项目中

- Rob（9/24）：核糖体蛋白的大量变化可能是蛋白合成受损，也可能只是细胞处于炎症/应激（“angry”）；不是他优先的方向（转录稿 01:21:52–01:22:13）。
- Yue（9/29）：核糖体相关变化是大家已知的，“不新颖”，不是合作者感兴趣的（`/Users/zeruzhang/Downloads/Zoom会议/Zongliang Yue's Personal Meeting Room 2026-09-29 16:02(GMT-5:00).txt`，约 16:37:27–16:37:59；ASR 文本，字词有误）。
- 功能组 G02（翻译、核糖体生成、rRNA）与 G03（“translation at synapse”等，由核糖体成员带出）：`iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`（G02=9 个、G03=10 个节点）。
- 敏感性分析把 `^(Rpl|Rps|Hsp|Dnaj)` 前缀的 151 个基因同时从候选和背景剔除：`iNKT_by_date/2026-09-25/results/sensitivity_manifest.json`。
- “加速器-刹车”假设：翻译相关表达上升的同时蛋白质量控制（Hspa8、Dnaja1）下降——ZERU 在 9/24 会议提出为待检验假设，仅为表达证据（转录稿 00:47:28–00:47:56）。

## 相关概念

- [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] — 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。
- [[为什么会出现 “measles”、“synapse translation” 这样的名字（Pathway Names Are Not Diseases）|Pathway-Names-Are-Not-Diseases]] — 通路/GO 名字是历史命名，命中的其实是共享的免疫、炎症、核糖体或热休克基因。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
