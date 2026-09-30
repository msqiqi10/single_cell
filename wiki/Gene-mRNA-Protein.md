# 基因、mRNA 与蛋白质（Gene / mRNA / Protein, Central Dogma）

> 基因是 DNA 上的“函数定义”，mRNA 是一次“调用的输出缓存”，蛋白质是最终干活的“进程”；scRNA-seq 只数 mRNA。

## 定义

基因（gene）是 DNA 上编码某个功能单元的一段序列；细胞把它转录（transcription）成信使 RNA（mRNA），再翻译（translation）成蛋白质（protein），这条单向流程叫中心法则（central dogma）。蛋白质才是大多数生化功能的执行者。

## CS 类比

基因 ≈ 源代码里的函数定义；mRNA ≈ 该函数被调用后写入内存的临时结果；蛋白质 ≈ 真正在运行的进程。单细胞 RNA 测序只能读到“内存里有多少份临时结果”，读不到进程本身是否在跑。

## 常见误读与注意

- “Fos 的 mRNA 变多”不等于“Fos 蛋白活性变强”。本项目所有表达量都是 mRNA 的 UMI 计数（见 [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]]），不测蛋白、磷酸化或 DNA 结合。
- 小鼠基因符号首字母大写其余小写（Fos、Prf1），人类全大写（FOS、PRF1）；跨物种数据库需要做大小写/同源映射。

## 在本项目中

- 基因符号是本项目每一张表的行名：10,670 个保留基因（`docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md`）。
- “AP-1 复合体富集不直接测量复合体形成、DNA 结合或蛋白活性”——`iNKT_by_date/2026-09-25/README.md` 与 `iNKT_by_date/2026-09-25/results/README.md` 第 3 节。
- Rob 会上指出核糖体基因变化“可能是蛋白合成受损，需要检验，也可能只是细胞状态‘angry’”（`/Users/zeruzhang/Downloads/Zoom会议/zoom-9024/GMT20260924-133029_Recording.transcript.vtt`，约 01:21:52–01:22:06）。
- KEGG 库中基因符号为大写，映射到小鼠符号时按不区分大小写匹配：`iNKT_by_date/2026-09-30/code/01_ora.py`。

## 相关概念

- [[转录因子（Transcription Factor, TF）|Transcription-Factor]] — 转录因子是决定“哪些基因被开启”的调控蛋白，相当于配置开关。
- [[核糖体与翻译（Ribosome and Translation）|Ribosome-and-Translation]] — 核糖体是把 mRNA 翻译成蛋白质的机器；Rpl/Rps 基因大量上调在这个项目里更像“细胞很愤怒”而非新发现。
- [[计数矩阵（Count Matrix, AnnData）|Count-Matrix]] — 细胞 × 基因的稀疏整数表；本项目主对象是 15,532 × 10,670。
- [[富集不等于激活（Enrichment ≠ Activation）|Enrichment-vs-Activation]] — “某通路富集”只说明列表和基因集重叠多，不说明通路被激活或起了因果作用。
