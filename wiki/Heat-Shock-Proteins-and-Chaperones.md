# 热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）

> 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。

## 定义

蛋白质必须折叠成特定三维形状才有功能。分子伴侣（chaperone）帮助折叠并处理错折叠蛋白；热休克蛋白（heat shock protein, Hsp）在应激时大量诱导；Dnaj 是 Hsp70 的伴侣因子；CCT/TriC 是折叠细胞骨架蛋白（如微管蛋白）的环状复合体。

## CS 类比

分子伴侣 ≈ 内存校验/垃圾回收：保证新分配的对象格式正确，出错就重试或回收。

## 常见误读与注意

- Hsp 基因是典型的应激/解离响应基因（见 [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]]），因此很多“通路”只是同一批 Hsp 基因换名字。
- 剔除 Hsp 后信号仍在，才说明不只是热休克。

## 在本项目中

- 脾脏 C3：Cct2/3/4/6a/8、Tuba1a、Tubb4b 等表达降低，作为局部候选（`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 8；37 页版第 36 页，`iNKT_by_date/2026-09-20/package/slides.json`）。
- R10：mask heat-shock 后再找新信息，已做 Hsp/Dnaj 排除敏感性（同 H 文档要求表）。
- “Measles” 行实际 driver 为 Fos、Hspa1a、Hspa1b、Hspa2、Jun、Tlr7（`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 第 4 节）。
- 功能组 G01 含 “protein folding” 及 spermatogenesis/gamete generation 等由热休克基因带出的名称（`iNKT_by_date/2026-09-25/results/network/annotated_nodes.csv`）；README 称 G01 为“混合注释，包含热应激、折叠”。

## 相关概念

- [[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] — 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。
- [[即刻早期基因与组织解离应激（Immediate-Early Genes, Dissociation Stress）|Immediate-Early-Genes-and-Dissociation-Stress]] — Fos/Jun 等即刻早期基因几分钟内就会被诱导，组织解离/处理本身也能诱导它们——这是 AP-1 结果的主要替代解释。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
- [[功能组 G01–G07 及含义（Functional Groups G01–G07）|Functional-Groups-G01-G07]] — 语义 + GNN 融合把 97 个 GO 条目分成 7 组，G04–G06 是能量/线粒体相关，G02/G03 是核糖体/翻译，G01 最混杂。
