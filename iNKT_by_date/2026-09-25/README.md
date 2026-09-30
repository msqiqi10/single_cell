# 2026-09-25｜9月24日会议后续：AP-1 / iNKT17 / GO功能归组

接在 [2026-09-19 细胞毒性与GO分析](../2026-09-19/README.md) 之后。数据相同（15,532细胞），按 Rob 与 Yue 的会议意见补做。

[汇报 PPTX（第20–31页，接19页细胞毒性/GO报告）](presentation/iNKT_followup_AP1_iNKT17_GO_network_20260925.pptx) · [英文简报](results/BRIEF_EN.md) · [实验讨论清单](results/EXPERIMENT_CANDIDATES.md) · [会议意见对照](results/REQUIREMENTS_CHECKLIST.csv)

## 主要结果

- **AP-1**：GO:0035976“transcription factor AP-1 complex”在骨髓C0（Fos/Fosl2/Jun/Junb）和脾脏C3（Fos/Jun/Junb/Jund）的T2上调基因中富集，global q≈0.003 / 0.002。20次QC匹配中方向全部保持。
- **MAPK**：GO“MAPK cascade”q=1，KEGG q=0.078，均不显著。不能因为AP-1显著就说MAPK激活。
- **iNKT17**：Ccr6之前被min_cells=100过滤掉，从原始矩阵找回。保留细胞中有60个表达，其中54个在C5-2。Il17a为0，Il17f只有2个细胞且都不在C5-2，因此无法判断分泌能力。
- **GO覆盖（Yue）**：BM C0共274个T2上调基因，只有66个落在显著BP条目里；脾脏C3为96/275。
- **功能网络（Yue）**：97个BP节点、577条边，分成7个探索性功能组（GOLDEN风格融合 + GraphSAGE）。这是小鼠GO适配版，不是原发表模型的完整复现。剔除Rpl/Rps/Hsp/Dnaj基因后，OXPHOS富集仍保留。
- **样本设计**：每个标签是3只小鼠pooled，每组n=1个pool，结果仍为细胞层面的探索性结论。

## 目录

- `results/`：9月25日跟进的完整结果包（表、图、网络、模型、验证、代码快照）
- `stage01_followup/`：第一阶段（完整DE候选与Rob基因核对）
- `meeting_evidence/`：会议逐条核对依据（ASR复核、截帧、候选覆盖）
- `notes/`：会议意见与下一步（中文），以及本地分析项目说明
- `code_local/`：本地运行入口代码（原项目位于 ~/Documents/Codex/2026-09-24/.../local_inkt，里面有 `.venv` 和 1.97 GiB 数据快照，没有复制到这里）
- `code/build_deck.js`：生成PPT的脚本
- `presentation/`：PPTX

## 外部依赖（未完成）

AML NK / 患者数据（Box）、用于 RNA velocity 的 FASTQ、每个pool的动物组成与处理记录、IL-17实验细节、湿实验验证。
