# bio3: 小鼠 iNKT 单细胞转录组分析 (Mouse iNKT scRNA-seq analysis)

## Project summary

本仓库记录一个小鼠 iNKT 细胞单细胞 RNA-seq 项目的全部分析代码、结果表、图、幻灯片和知识库。核心问题：
**肿瘤 (T2) 条件是否改变骨髓 (BM)、脾脏 (Spleen)、胸腺 (Thymus) 中的 iNKT 细胞？** 具体到细胞群 (Leiden cluster C0-C7)、
基因（如 AP-1、iNKT17、细胞毒性模块）和功能程序 (GO / KEGG)。下游合作者 Rob、PI Yue 的会议意见驱动了 09-19 到 09-30 的几轮追加分析。
工作在服务器 `/home/zzz0054/bio3` 上完成，本仓库是 2026-09-29 下载的快照，加上 09-25 / 09-30 的后续。

## Data at a glance

- 3 个组织 (BM / Spleen / Thymus) x 2 个条件 (Ctrl / T2) = 6 个样本标签。
- 每个标签由 **3 只小鼠混合 (pooled)**，即每组 **n = 1 个 pool**，没有生物学重复；所有 p 值/q 值都是细胞层面的，存在 pseudoreplication 风险（见 wiki 的 `Pseudoreplication-and-n-equals-1`）。
- QC 后 **15,532 cells x 10,670 genes**（原始 18,458 x 32,285）。
- 8 个 Leiden 簇 C0-C7（C5 再分 C5-1 / C5-2）。

## Directory map

```text
.
├── README.md                 本文件
├── DATA_NOT_IN_GIT.md        未入库的大文件 (*.h5ad): 路径 / 大小 / SHA-256 / 来源
├── pyproject.toml, uv.lock   uv 环境 (Python 3.12, scanpy, scvi, CPU torch)
├── scripts/
│   ├── runners/              run.sh 和 run_inkt_*.sh 服务器运行脚本 (15 个)
│   ├── check_md_links.py     markdown 相对链接检查器
│   └── md_links_known_dangling.txt  快照里本来就失效的链接基线
├── notebooks/
│   ├── iNKT/                 原始 notebook
│   └── scripts/iNKT/         分析/验证/出图 Python 脚本 (run_inkt_*.py, build_*.py, test_*.py)
├── output/iNKT_meeting_followup_20260905/   09-05 的 DE/表/来源 (objects/*.h5ad 不在 git)
├── iNKT_by_date/             按日期归档的代码 + 结果 + 幻灯片 (见下方时间线)
├── docs/
│   ├── references/           参考论文 PDF
│   ├── audits/               审计/历史/流程说明 markdown 与 csv
│   └── README_download_snapshot_20260929.md   旧版 README (原样保留)
└── wiki/                     81 页 GitHub-wiki 风格知识库 (入口 wiki/Home.md)
```

## Timeline of dated analyses

| 日期 | 内容 | README |
|---|---|---|
| 2026-09-19 | iNKT 细胞毒性与 GO 分析（PAGER 核心/E1 模块、五次 QC 匹配、GO ORA、局部 DPT） | [iNKT_by_date/2026-09-19](iNKT_by_date/2026-09-19/README.md) |
| 2026-09-20 | 论文式资料包与可编辑 PPT（28 / 58 页等） | [iNKT_by_date/2026-09-20](iNKT_by_date/2026-09-20/README.md) |
| 2026-09-25 | 9/24 会议后续：AP-1 / iNKT17 / GO 功能归组 | [iNKT_by_date/2026-09-25](iNKT_by_date/2026-09-25/README.md) |
| 2026-09-30 | GO/KEGG 功能网络与锚点分组（PI Yue 的问题） | [iNKT_by_date/2026-09-30](iNKT_by_date/2026-09-30/README.md) |

更早的阶段（08-18 到 09-15）的历史与导师要求见 [docs/audits/inkt_history_and_supervisor_requirements_20260915.md](docs/audits/inkt_history_and_supervisor_requirements_20260915.md)，
流程逐步说明见 [docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md](docs/audits/inkt_pipeline_stage_by_stage_walkthrough.md)。

## Knowledge base (wiki)

从 [wiki/Home.md](wiki/Home.md) 开始：每个生物学/统计术语都有“通俗解释 + 与本项目的联系”。发布到 GitHub Wiki 的方法见 [wiki/README_PUBLISHING.md](wiki/README_PUBLISHING.md)。

## How to reproduce

- 环境：`uv sync` 创建 `.venv`（Python 3.12，Linux CPU wheels）。`.venv/` 不在 git 中。
- `scripts/runners/*.sh` **假定服务器布局**：它们先 `cd /home/zzz0054/bio3`（绝对路径），再调用 `.venv/bin/python notebooks/scripts/iNKT/...`，
  所以必须在服务器上、仓库内容位于 `/home/zzz0054/bio3`、且 `.venv` 与原始输入 (`input/`，约 7 GB 的 `input.zip`) 存在时才能运行。
  例：`bash scripts/runners/run_inkt_meeting_followup_20260905.sh scores`。
- `.h5ad` 中间对象未入库，可从服务器取回或由上面的 runner 重新生成，详见 [DATA_NOT_IN_GIT.md](DATA_NOT_IN_GIT.md)。
- 检查文档链接：`python3 scripts/check_md_links.py` 和 `python3 wiki/_tools/check_links.py --paths`。

## What is not in git

三个 >50 MB 的 `.h5ad`（scored_base / bone_marrow_stable / spleen_stable）以及其余 `.h5ad`（共 7 个，约 410 MiB）、`*.zip`、`.venv/`。
清单、SHA-256 和来源见 [DATA_NOT_IN_GIT.md](DATA_NOT_IN_GIT.md)。未使用 Git LFS。

## Portability notes

下列文件含本地或远程机器的绝对路径，属于历史记录，**未改写**；在别的机器上运行需手动替换：

- `iNKT_by_date/2026-09-30/code/common.py`：`/Users/zeruzhang/Documents/Codex/2026-09-24/...`（KEGG 输入、MiniLM 模型、`VENV_PY`）
- `iNKT_by_date/2026-09-25/stage01_followup/code/fetch_existing_GO.py`：本机 ssh control socket 路径 (`/Users/zeruzhang/.ssh/...`)，用于连接服务器
- `iNKT_by_date/2026-09-25/code_local/inventory_remote.py` 和 `iNKT_by_date/2026-09-25/results/code/inventory_remote.py`：服务器/本机绝对路径
- `iNKT_by_date/2026-09-25/results/` 与 `iNKT_by_date/2026-09-30/results/network/*/models/` 下的 manifest/log/summary JSON 记录了生成时的 `/Users/zeruzhang/...`、`/private/tmp/...` 路径（只是溯源记录）
- `iNKT_by_date/2026-09-25/results/provenance/remote_manifest.json` 记录服务器端脚本名（`run_inkt_*.sh`）的旧根目录位置
- 日期目录中的旧代码（例如 `iNKT_by_date/2026-09-19/code/analyze.py`、`iNKT_by_date/2026-09-20/code/capture_sources.py`）仍引用旧的 `docs/*.pdf` 位置；
  它们是当时运行的历史快照，未改动。论文现在位于 `docs/references/`。
- 若干 markdown 文档里的链接指向快照中不存在的文件（已在 `scripts/md_links_known_dangling.txt` 列出）。
- 文本笔记 `wiki/*.md` 中若干页面也引用了上述本地路径作为出处说明。
- `scripts/runners/run.sh` 会在当前目录生成 `run_inkt_tmux.sh`（服务器根目录），这是它原有行为。

## License

license: TBD by PI（尚未选择许可证；在 PI 决定之前请勿再分发）。
