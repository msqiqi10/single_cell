# 敏感性分析（Sensitivity Analysis）

> 换一个合理的分析选择，看结论是否还在；不是新的独立验证。

## 定义

敏感性分析（sensitivity analysis）保持主分析不变，另外系统地改变一个选择（阈值、基因集、抽样、排除项），检验结果是否稳健。

## CS 类比

≈ 回归测试 + 参数扫描：主结果是 baseline，敏感性是不同配置下的 diff。

## 常见误读与注意

- 共享同一批输入的分析不算独立验证。
- 稳健只说明不依赖该选择，不能证明生物学机制。

## 在本项目中

- 本项目做过的敏感性：QC 匹配（20 次）；Rpl/Rps/Hsp/Dnaj 排除（151 基因）；Harmony；PC 数 50/100/200；模块加/不加 Ctla2a、Nkg7；DPT root 敏感性；三组织稳定分群（见各概念页）。
- “签名依赖性必须一并汇报”：加 Ctla2a 后骨髓 C4/脾脏 C3 均不显著（q=0.93、0.53）（`iNKT_by_date/2026-09-19/README.md`）。
- R10：heat-shock mask（`docs/inkt_history_and_supervisor_requirements_20260915.md` 表格）。
- 09-30 分析同时运行“全部基因”与“排除 Rpl/Rps/Mrp/mt-/Hsp/Dnaj”两个版本（`iNKT_by_date/2026-09-30/code/common.py`）。

## 相关概念

- [[QC 匹配与标准化均值差（QC Matching, SMD）|QC-Matching-and-SMD]] — 让 Ctrl 与 T2 细胞的测序深度、检出基因数、线粒体比例匹配后再比，检查结果是否由 QC 差异造成。
- [[热休克蛋白、分子伴侣与蛋白折叠（Heat Shock Proteins, Chaperones, Protein Folding: Hsp, Dnaj, CCT）|Heat-Shock-Proteins-and-Chaperones]] — 分子伴侣帮助新合成的蛋白折叠成正确形状；Hsp/Dnaj/CCT 是最常见的伴侣家族。
- [[批次效应与 Harmony 整合（Batch Effect, Harmony, Integration）|Harmony-and-Batch-Integration]] — 批次效应是技术带来的系统偏差；Harmony 把不同批次的细胞在低维空间“对齐”。
- [[细胞毒性（Cytotoxicity: Perforin, Granzymes）|Cytotoxicity]] — 细胞毒性指杀伤靶细胞的能力，转录层面常用 Prf1/Gzma/Gzmb 三个基因的均值近似。
