# 拟时序（Diffusion Map, DPT, PAGA）与 RNA 速度（Pseudotime vs RNA Velocity）

> 拟时序按表达相似度给细胞排序；RNA velocity 用剪接/未剪接 RNA 估计变化方向，需要额外输入。

## 定义

拟时序（pseudotime）沿细胞表达相似性构成的图给细胞排一个连续顺序：扩散图（diffusion map）+ 扩散伪时间（diffusion pseudotime, DPT）+ PAGA 图连接。RNA velocity 则用每个基因未剪接与已剪接 mRNA 的比例估计瞬时变化方向。

## CS 类比

拟时序 ≈ 在相似度图上做最短路径排序；velocity ≈ 读取“未提交的写入队列”，预测下一步状态。

## 常见误读与注意

- 拟时序不是时间、速率或已验证谱系。PAGA 无方向，DPT 只区分了根群与其余。
- RNA velocity 需要 spliced/unspliced 计数，需要 BAM 或能重新比对的 FASTQ；本项目缺输入，R04 未完成。

## 在本项目中

- 全局：root 选 c7（因代码 fallback 简化为 Cd27），root 细胞恰是第一个 c7 细胞（Ctrl 脾脏），使 root 很脆弱；DPT 几乎是 “c7 低，其余接近 1”（`docs/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 11）。
- 9/19 重新做每组织局部 diffusion map，root 取 Egr2/Hivep3 评分最高群（Cd24a 不在基因集）；骨髓/脾脏 root 敏感性最低相关约 0.76；胸腺没有可靠连续曲线（`iNKT_by_date/2026-09-19/README.md` “拟时序方法与限制”；19 页版第 13–15 页）。
- R04 velocity：缺 spliced/unspliced 或可生成它们的 BAM/FASTQ（`docs/inkt_history_and_supervisor_requirements_20260915.md` 阶段 7）；Rob 只有 FASTQ、没有 BAM（转录稿 01:34:25–01:34:42）；Yue 期望新数据可做 velocity（01:34:49）。

## 相关概念

- [[近邻图（Neighbour Graph, k-NN）|Neighbour-Graph]] — 在 PCA 空间给每个细胞找 15 个最相似的邻居，UMAP 和聚类都建立在这张图上。
- [[未决问题与下一步（Open Questions and Next Steps）|Open-Questions-and-Next-Steps]] — 下一步汇总。
- [[自然杀伤细胞（Natural Killer Cell, NK cell）|NK-Cell]] — NK 细胞是先天免疫的“杀伤者”，靠穿孔素和颗粒酶杀死异常细胞；本项目用 NK 论文做方法参照。
- [[敏感性分析（Sensitivity Analysis）|Sensitivity-Analysis]] — 换一个合理的分析选择，看结论是否还在；不是新的独立验证。
