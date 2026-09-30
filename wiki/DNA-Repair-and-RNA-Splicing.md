# DNA 修复与 RNA 剪接（DNA Repair, RNA Splicing: Shld1, Sfpq）

> DNA 修复与 RNA 剪接是细胞维持基因组和转录本质量的两个过程；本项目仅作为次级候选。

## 定义

DNA 修复（DNA repair）纠正 DNA 损伤；Shld1 属于 shieldin 复合体，参与双链断裂修复通路选择。RNA 剪接（splicing）把前体 mRNA 的内含子剪掉，Sfpq 是参与 RNA 加工的蛋白。（背景知识。）

## CS 类比

DNA 修复 ≈ 数据库的崩溃恢复；剪接 ≈ 预处理阶段去掉注释再编译。

## 常见误读与注意

- 表达下降不等于修复系统失效；SFPQ 表达变化也不等于剪接事件改变（要检查具体剪接事件）。
- 骨髓 C0 的剪接子集 ORA 有支持，但全排序 GSEA 不显著，不能称整个 spliceosome 被激活。

## 在本项目中

- ZERU 在 9/24 会议指出 Sfpq、Dynll1 在骨髓/骨髓 C5-2 下调，Shld1 在骨髓 C5-2 下调，并强调这只是提出问题（转录稿约 00:48:19–00:48:53；00:54:33–00:55:29）。
- Sfpq、Dnaja1、Hspa8 的下调在旧 PPT 重跑中方向一致（`docs/audits/inkt_legacy_ppt_vs_legacy_qc_rerun_content.md`）。
- “次级候选：骨髓 C0 剪接”：`docs/audits/inkt_history_and_supervisor_requirements_20260915.md` 阶段 8。
- “Fanconi/DNA 修复”条目在脾脏 C5-1 的全库排序中出现，FDR=1，属弱统计支持的探索线索（转录稿 01:06:57–01:07:37）。

## 相关概念

- [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] — 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。
- [[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] — 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。
- [[基因集富集分析与 NES（Gene Set Enrichment Analysis, GSEA; Normalized Enrichment Score, NES）|GSEA-and-NES]] — GSEA 不需要先选“显著基因”，而是把全部基因按变化排序，看基因集是否集中在顶端或底端；NES 带方向。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
