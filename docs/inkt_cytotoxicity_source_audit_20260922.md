# 细胞毒性选取与更新数据核查（2026-09-22）

## 现有结果核对

最新实质计算为 iNKT_by_date/2026-09-19；9月20日是可编辑展示更新。输入仍为9月5日 scored_base.h5ad：15,532细胞、10,670基因。
当前三个模块：Prf1/Gzma/Gzmb；上述三基因加Ctla2a；上述三基因加Nkg7。前两个是对PAGER方法的简化/代表基因，并非经核实的完整作者补充表。

直接核对 cytotoxicity_contrasts.csv：

|组织/群|核心3基因差值 / q|加Ctla2a差值 / q|加Nkg7差值 / q|
|---|---|---|---|
|骨髓C4|+0.08620 / 0.02513|-0.00833 / 0.93258|+0.07593 / 0.01884|
|脾脏C3|-0.02297 / 0.001538|-0.00615 / 0.53344|-0.03367 / 0.0003433|

差值为T2−Ctrl的平均log1p表达差；q为细胞层面检验，不是动物重复验证。签名依赖明显，不宜仅按显著性选基因。

## 更新的外部数据：可用性与用途

1. **优先：GSE296020，小鼠iNKT1细胞毒性亚群。** GEO公开日期2026-04-30；对应论文为2026年8月《Journal of Immunology》。分选群体bulk RNA-seq，共8个文库：胸腺Fgd5阳性/阴性各2个，骨髓CD5/CD244群各2个。GEO有胸腺和骨髓原始计数CSV.gz。胸腺对照适合提取外部iNKT1状态候选签名；骨髓两群的分选定义不同，不可直接视作同一iNKT群的高低细胞毒性对照。也不能将全部差异基因等同于直接杀伤基因。
   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296020
   - https://academic.oup.com/jimmunol/article-abstract/215/8/vkag199/8768938
2. **GSE298293：胸腺小鼠iNKT单细胞参考。** 2026-01-21公开；BALB/c、B6各2个文库，提供MTX/TSV（系列压缩包156.9 MB）。适合核对亚型及成熟状态特异性；不是杀伤实验，也不是我们T2处理的重复。
   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298293
3. **GSE306154：结肠iNKT单细胞参考。** 2026-02-26公开；共9个文库，含7个人源和2个鼠源文库（EGFRi、Ctrl）。提供MTX/TSV（179.3 MB）。组织差异大、鼠源条件各1个文库，只列次选。
   - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306154
4. **Mouse MSigDB 2026.1.Mm：可下载的正式基因集。** 可选M5 GO/M2通路作独立定义；这是注释资源而非新增表达实验。不能仅凭版本名认为它比本地GO更新；本地GO本体为2026-07-26，已存MGI注释和版本哈希。
   - https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp?targetSpeciesDB=Mouse

补充基因依据：2025年NKG7研究可支持把Nkg7纳入敏感性比较，但不构成小鼠iNKT专用完整签名。
https://onlinelibrary.wiley.com/doi/10.1002/eji.202551885
现有PAGER参考的研究对象为小鼠NK：https://pmc.ncbi.nlm.nih.gov/articles/PMC11058213/

## 当前矩阵覆盖（仅查HDF5基因名，无重算）

- 存在：Prf1、Gzma、Gzmb、Nkg7、Gzmm、Ctsw、Fasl、Ccl5、Ifng、Rab27a、Unc13d、Stx11、Stxbp2。
- 不在当前过滤矩阵：Gzmk、Gnly、Gzmh、Ncr3。缺失于此矩阵不等于生物学不表达；跨物种签名需另核对同源关系和原始features。
- 上述扩展基因仅为覆盖核查候选，不是本次正式选定的新评分模块。
- 本地GO目录中，GO:0042267（NK细胞介导细胞毒性）覆盖29个基因，GO:0001906（cell killing）66个，GO:0140507（granzyme-mediated programmed cell death signaling pathway）10个。这些集合可能包含调控/靶细胞相关成员，需逐基因检查，不能全数当作同方向杀伤分数。

## 建议

保留原模块作历史对照；优先核对GSE296020胸腺对照及论文补充表，结合直接效应、脱颗粒和受体/调控的功能区分固定候选集，再评估我们数据中的覆盖与稳健性。预先固定基因选取规则，完整报告不同定义下的结果，不以能否使C4/C3显著作为选取标准。

本轮完成来源检索、样本/下载格式核对、现有结果及基因名核对；未下载外部表达矩阵，未重算评分、改写旧分析或提交代码。论文全文/补充表及外部计数内容尚待逐项核验。
