# 文档之间看起来不一致或需要核对的地方（Document Inconsistencies and Unsure Points）

写这个 wiki 时读到的、你可能会踩坑的地方。这里只是记录，不是判断谁对谁错。

1. **仓库 README 链接的日期目录本地不存在**：`README.md` 链接 `iNKT_by_date/2026-09-15/`、`2026-09-06/`、`2026-09-05/`，但本地只有 `iNKT_by_date/2026-09-19/`、`2026-09-20/`、`2026-09-25/`（和进行中的 `2026-09-30/`）。README 提到的 `input.zip`、`CML_NK_scRNA_TKI/` 也不在目录中。
2. **早期文档引用的 `output/` 文件大多不在本地**：`docs/inkt_history_and_supervisor_requirements_20260915.md` 与 `docs/inkt_pipeline_stage_by_stage_walkthrough.md` 指向 `output/iNKT_reproduction_deck/`、`output/iNKT_discovery_20260905/` 等（以及 `/home/zzz0054/...` Linux 路径）；本地 `output/` 只有 `iNKT_meeting_followup_20260905/` 的部分表格。所以早期数字（如 ARI=0.9868、35/39）只能引用记录它们的文档，没有回到源表逐项核对。
3. **转录稿路径**：任务书写 `/Users/zeruzhang/Downloads/zoom-9024/…`，实际在 `/Users/zeruzhang/Downloads/Zoom会议/zoom-9024/GMT20260924-133029_Recording.transcript.vtt`。
4. **会议笔记里的相对链接与实际目录名不同**：`iNKT_by_date/2026-09-25/notes/2026-09-24_Rob_Yue会议意见与下一步.md` 指向 iNKT_20260925_Rob_Yue_results/README.md（该文件夹不存在） 与 `2026-09-24_会议核对依据/`，实际为 `iNKT_by_date/2026-09-25/results/README.md` 与 `iNKT_by_date/2026-09-25/meeting_evidence/`。
5. **幻灯页码有三套编号**：37 页版（Rob 说“slide 5”“slide 31”即此编号）、19 页版（`iNKT_by_date/2026-09-19/presentation/slide_index.csv`）、58 页整合版顺序（`iNKT_by_date/2026-09-20/package/notes/deck_order.json`，`presenter_guide.zh.md` 用的是这套；另有 28 页主报告）。会议笔记还说“录像中编辑器顺序号与页面印刷号不一致”。引用页码时请写明是哪一套。本 wiki 用 37 页版页码时对应 `slides.json` 里 series=old 的前 37 项。
6. **簇数量随阶段变化**：旧 PPT 11 群；当前 8 群（分辨率 0.5）；C5 拆分后 9 个 refined cluster；单组织稳定分群另用 S0… 命名（`iNKT_by_date/2026-09-19/results/de/stable__*.csv.gz`）。
7. **两个“17”**：`docs/inkt_history_and_supervisor_requirements_20260915.md` 阶段 8 的“严格优先表有 17 条 gene×contrast 记录”与 `iNKT_by_date/2026-09-25/results/README.md` 的“17 个组织×群比较”是不同的量。
8. **“26”个 DE**：`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 8.2 规划 26 个单元（8 簇，19 完成 7 跳过）；09-25 保留 26 份既有 DE 表（含 C5 拆分后的簇；17 组织×群 + 3 组织为主分析）。数字相同，构成可能不同，未逐项核对（不确定）。
9. **PC 数**：PCA 计算 50 个成分，邻居图只用前 30 个（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 4.2、5.1），而 UMAP 比较和汇报都称 “50 PCs”；文档未解释两者关系（不确定）。9/24 转录稿写成 “500, 50, 100, 200 PCs” 是 ASR 误识。
10. **ASR 的问题**：09-24 转录稿中 “can now treat these 6 labels as 6 independent animal replicates”（00:36:45）与结论相反，疑似把 “cannot” 识别成 “can”；09-29 转录稿把 KEGG 写成 “cag”、pathway 写成 “password”、MAPK 写成 “mapness”，以及 “GNN” 有多种误识。本 wiki 对这些处标“推测”。
11. **Il17a/Il17f 的措辞**：`docs/inkt_pipeline_stage_by_stage_walkthrough.md` 说它们“不在 10,670 gene universe”；09-25 从原始矩阵重建后说“保留细胞中 0 个/2 个检出，原始输入条码中各 3 个”。二者并不冲突，但前者的说法只到“不可测”，后者才是检出计数。
12. **T2 的定义**：转录稿称 T2 是 “tumors dataset”，Rob 团队做 CML 嵌合小鼠模型，但项目文件没有写明 T2 的具体模型/时间点。
13. **Rob 的全名**：转录稿只有 “Rob”；与 Borra et al. 论文作者 Robert S. Welner（UAB）的对应是推断。
14. **09-30 目录**：已完成；`results/html` 仍为空（交互 HTML 在 `drilldown/html/`），脚本中含指向 `/Users/zeruzhang/Documents/Codex/...` 的绝对路径（`iNKT_by_date/2026-09-30/code/common.py`），本地复跑需要那套环境。
15. **README 里桥接基因的措辞**：09-30 README 早先写成“Fos、Jun（+ Junb/Jund）”，容易被读成四个基因在每一对通路都是桥；已改为“Fos + Jun，另加 Jund/Junb/Nfatc2 之一（视通路对而定）”（`iNKT_by_date/2026-09-30/README.md`）。按 `iNKT_by_date/2026-09-30/drilldown/tables/anchor_bridge_genes.csv`，Nfatc2 出现在 TCR、Th17 里且不是 DEG，所以严格说“有 DE 支持的”只有 Fos + Jun（TCR、Th17）或再加 Jund（MAPK、IL-17）/Junb（TNF）。
16. **09-30 README 的“偏离规格”一节**仍写“未做基因层面钻取与 HTML”，但同一文件后面的 Drill-down 一节已补做（追加时未更新该句）。
17. **GOLDEN 论文 KG 的 README_KG.md**：`archive/2026-09-30/` 包内 `README_KG.md` 仍是 09-29 版文字（`/Users/zeruzhang/Documents/zhou2_latest_20260929/deliverables_clean_latest/README_for_Yue_20260930.md`）。
