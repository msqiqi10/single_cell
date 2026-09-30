# 组织：骨髓、脾脏、胸腺（Bone Marrow, Spleen, Thymus）

> 骨髓造血、胸腺培育 T 细胞、脾脏是血液过滤与免疫应答场所；三者是不同的“微环境”。

## 定义

骨髓（bone marrow, BM）是造血和白血病细胞所在处；胸腺（thymus, Thy）是 T 细胞（含 iNKT）发育成熟的场所；脾脏（spleen, Spl）是过滤血液、集中免疫应答的外周器官。（背景知识；iNKT 在这三处均可见。）

## CS 类比

三个组织 ≈ 三个部署环境（dev / staging / prod）：同一个程序在不同环境里的配置和负载不同。

## 常见误读与注意

- 组织是数据里最强的变异轴，pooled 全局比较会把“组织差异”与“T2 效应”混在一起，因此应优先看组织内比较。

## 在本项目中

- 样本量（QC 后 T2/Ctrl）：骨髓 2,731/3,379；脾脏 3,371/3,422；胸腺 1,257/1,372；胸腺保留率最低（Ctrl 65.3%）——`docs/inkt_pipeline_stage_by_stage_walkthrough.md` Stage 1、8.3。
- 簇的组织构成：c0 95.8% 骨髓，c3 93.7% 脾脏，c6 97.3% 胸腺（同上 5.3 表）。
- 骨髓 C4 局部候选与脾脏 C3 CCT/TriC 程序（`docs/inkt_history_and_supervisor_requirements_20260915.md` 阶段 8）。
- Rob：与其它组织相比，很多信号只是细胞所处的邻域/应激（转录稿 01:22:42–01:23:10）。

## 相关概念

- [[小鼠模型与混样（Mouse Model and Pooled Samples）|Mouse-Model-and-Pooled-Samples]] — 每个样本标签由 3 只小鼠混合而成：每组只有 n=1 个 pool。
- [[本项目的簇标签 C0…C7、C5-1 / C5-2（Cluster Labels C0-C7）|Cluster-Labels-C0-C7]] — 8 个 Leiden 簇 + C5 拆成 C5-1/C5-2 = 9 个 refined cluster；C0 骨髓、C3 脾脏、C6 胸腺为主。
- [[白血病 CML / AML 与 “T2” 条件（Leukemia, CML, AML and T2）|Leukemia-CML-AML-and-T2]] — T2 = 肿瘤条件；本项目文件没有说明 T2 具体对应哪种白血病，只能说是 tumour。
- [[细胞组成与比例（Cell Composition, Proportions）|Cell-Composition-and-Proportions]] — 比例图有两种分母：条件内的簇频率，和簇内的条件占比，不能混读。
