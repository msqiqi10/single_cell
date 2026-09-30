# iNKT 本地分析项目

2026-09-25 从 `/home/zzz0054/bio3` 复制当前 iNKT 分析所需的数据、代码和注释。后续在本目录维护分析，源数据保存在 `data/bio3/`，新结果写入 `runs/`，不覆盖远端或本地输入。

## 大小与内容

- 数据快照：1,334 个文件，2,110,268,178 字节，约 1.97 GiB。
- 初始数据加基础环境约 2.7 GiB；完成模型安装与本轮分析后，当前项目约 3.7 GiB。
- 6 个样本的原始 10x 表达矩阵、主要 AnnData 对象、完整 DE 表、GO 注释及已有结果、分析代码均已复制。
- 当前主矩阵为 15,532 个细胞 × 10,670 个基因。
- 远端整个 bio3 约 25 GB；本地没有复制无关项目、旧环境和大部分历史输出。本快照不是整个远端项目的完整备份，也不包含 FASTQ。

## 已验证

全部 1,334 个文件的大小和 SHA256 与远端清单一致。主矩阵成功本地读取，6 份原始矩阵的 matrix / barcodes / features 文件齐全。见 [文件与矩阵验证](provenance/local_verification.json)。

完全离线复算骨髓 C0 的 T2-up / Ctrl-up、BP / MF / CC 共 17,912 项 GO 检验。背景数、条目成员数、输入数、命中数和命中基因与旧结果逐项一致，P 值和家族内 BH FDR 误差仅为浮点精度。见 [GO 验证](runs/BM_C0_GO_local_validation/validation.json)。这一步复用已有 DE 与冻结的 GO 成员映射，没有重跑 DE，也没有复算跨所有比较的全局 FDR。

## 使用

在本目录的终端中运行：

```bash
./run_local.sh verify
./run_local.sh go-pilot
```

运行不需要 VPN 或 SSH。独立环境位于 `.venv/`，Python 3.12；精确包版本见 `requirements.lock.txt`。如需重建环境，使用 Python 3.12 创建 venv，再安装该文件，安装依赖本身需要网络。不要复制远端 Linux 的虚拟环境到 Mac。

`src/` 保存本地入口代码；`provenance/remote_manifest.json` 保存来源和逐文件哈希。新增分析保留参数、筛选口径及验证记录，输出到新的 `runs/` 子目录。需要远端新增数据时再显式拉取新快照；当前没有配置自动双向同步。

## 接续会议要求

**更新：会议要求的本地分析已执行。** 见 [完整结果](runs/meeting_followup/README.md)、[英文简报](runs/meeting_followup/BRIEF_EN.md)、[实验讨论清单](runs/meeting_followup/EXPERIMENT_CANDIDATES.md)。已覆盖 20 个主条件比较，恢复完整 32,285 特征，运行 GOLDEN/Fusion 小鼠适配及三次 GNN，并生成可追溯功能网络和技术敏感性检查。

可用 `./run_local.sh verify-followup` 核验结果；`./run_local.sh followup` 将按冻结配置重新计算本轮各阶段。后者不会下载新数据或模型，使用已有本地环境与模型，但会更新本轮计算输出；报告文字需要在输入改变时重新核对。模型安装后更新的精确依赖仍见 requirements.lock.txt。

本地数据已用于完整候选覆盖、GO 功能归组、Rob 指定基因核对，以及从原始矩阵检查被当前特征集遗漏的基因。

GOLDEN / GNN 使用显式记录的小鼠 GO 本地适配，不声称原发表模型的完整复现。当前验证也不等于整个远端历史流程已在 Mac 上全部验证。新 AML 数据、FASTQ 和独立动物元数据不在本次快照中。
