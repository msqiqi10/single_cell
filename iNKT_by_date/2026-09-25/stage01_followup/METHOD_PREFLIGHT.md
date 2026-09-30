# GOLDEN 接入预检查

已读取的方法稿：Downloads 中的 GOLDEN-F_nonformatted (1).docx 和 GOLDEN-GNN-update.docx。以下只记录源代码状态，不表示已运行模型。

## GOLDEN

官方版本：94f9d367483ccb4872630ebaa3b611b11220d9a6。

- `toden_e_predict` 实際读取的 valid_nodes 用于匹配 GO/PAG 元数据和 PAG-PAG 边；README 所称的 gene-list 入口不能未经转换就接收我们的差异基因名单。
- 需要 `go_metadata/m_type_biological_process.txt`、元数据与图结构等未随此仓库发布的文件。
- CLI 使用 `txt_file_path` 调用，而函数参数名为 `pags_txt_path`，存在直接调用错误。
- 特征与邻接矩阵分别排序，且图只加入边上的节点；适配时必须按同一 GO ID 列表显式对齐，并保留孤立节点，避免样本错位或丢失。
- 不直接执行未经检查的 PAGER／总结服务调用；本次分析的输入和推断在用户本地／研究服务器处理。

## GOLDEN-GNN

官方版本：a85654b6b853ddd60e59efa7f3bb769832d687b0。

- 顶层流水线是 human-only RummaGEO 0.85 / 0.90 示例，使用随库的人类基因集，不是小鼠 GO／iNKT 输入。
- 可复用的底层 `gnn_linkpred_holdout.py` 接收边表和节点特征 NPZ。边表需含 `GS_A_ID`、`GS_B_ID`、指定权重；节点必须与特征一致。
- 新的 mouse GO 输入需要条目 ID、正式名称／描述、冻结版本基因成员、组织／群／条件证据，以及按明确规则生成的边。全条目成员与本次命中 driver 必须分开。
- 网络边重建指标只评价图表示；不能将内部 AUC、聚类稳定性当作独立生物学验证。
- 方法稿和本地 paper-ready 索引使用的模型、维数与聚类方法有版本差异。下一阶段按真实执行代码与固定配置记录，不从稿件中拼接参数。

## 尚需的实际材料

已计算 GO 全表与 coverage/catalog、冻结注释成员、相应分析代码在远程 bio3 的 2026-09-19 目录；包内 current_analysis_tables 的 22 个预期文件已登记。先取得并核验这些文件，再确定实际节点范围。没有用其他项目的 SP 标签或人类嵌入代替 iNKT 结果。
