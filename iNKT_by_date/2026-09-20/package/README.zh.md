# iNKT 可编辑汇报资料包｜2026-09-20

建议先打开 **index.html**：按整合顺序浏览每页，直接找到该页的独立源图、SVG、PNG、数值表和论文出处。

## 四个现成版本

| 版本 | 用途 | PPTX | PDF |
|---|---|---|---|
| 28页主报告 | 按论文逻辑组织的主线 | [可编辑PPT](decks/iNKT_main_talk_28_editable.pptx) | [PDF](decks/iNKT_main_talk_28_editable.pdf) |
| 58页完整版 | 主报告28页 + 详细附录30页；原37页与新增19页全部保留，加2页整合说明 | [可编辑PPT](decks/iNKT_integrated_58_editable.pptx) | [PDF](decks/iNKT_integrated_58_editable.pdf) |
| 上次37页转换版 | 9月15日综合PDF，逐页还原为可编辑元素 | [可编辑PPT](decks/iNKT_previous_37_editable.pptx) | [PDF](decks/iNKT_previous_37_editable.pdf) |
| 新增19页转换版 | 9月19日细胞毒性与GO分析，独立使用 | [可编辑PPT](decks/iNKT_cytotoxicity_GO_19_editable.pptx) | [PDF](decks/iNKT_cytotoxicity_GO_19_editable.pdf) |

## Supervisor 可以怎样编辑

- **标题、正文、图题、结论和说明**：原生PowerPoint文本框，直接点击修改。
- **数值表/方法表**：原生PowerPoint表格，逐格修改；CSV备份在data/。
- **UMAP、dotplot、GO气泡图、曲线等**：每个图形面板独立插入，可单独移动、缩放、删除或替换。程序图提供SVG矢量源和PNG回退；已有原始图片直接使用源文件。
- **没有使用整张幻灯片截图作为PPT内容。** assets/page_previews仅用于找图与质量检查，不嵌入PPT。
- 图内散点、坐标和图例仍属于图形对象，不是原生Excel数据图表；如果要改统计数值或重新绘图，请使用对应数值表与绘图代码。支持SVG的PowerPoint可尝试“转换为形状”；不承诺每个点都是原生可编辑数据点。
- 编辑PPT后，在PowerPoint使用“导出PDF”同步修改。包内PDF由同一源图布局生成，不是Office渲染导出；实际Office字体替换、换行可能略有差异。

## 主报告组织逻辑

1. 数据与QC → 细胞状态/marker → 组织分层组成（原iNKT展示和Blood Figure 3的逻辑）。
2. 差异表达 → 杀伤模块 → 群特异与条件相关GO（PAGER Figure 5与Borra Figure 4c）。
3. GO细胞组分 → 局部拟时序 → QC敏感性 → 证据边界与结论。
4. 附录保留所有详细DE、旧通路核查、GSEA、Venn、整合敏感性和既有候选发现。

## 资料来源与实际结果

- 旧部分：9月15日37页综合报告，复用其原绘图代码和数值。
- 新部分：9月19日19页细胞毒性/GO结果。未重新做生物分析，未改变统计阈值或挑选新结论。
- 细胞杀伤核心信号：骨髓C4升高、脾脏C3小幅降低；扩展E1代表模块与杀伤相关GO未显著。只支持局部转录状态候选，不能写成已证明杀伤功能改变。
- 当前数据每组织×条件只有1个sample标签，没有独立动物重复；拟时序不是velocity或谱系验证。

## 文件导航

- `assets/panels/`：可直接插入PPT的142个独立图元素及其SVG版本。
- `data/`：展示用数值表、GO完整结果、源DE表与原生表格CSV；不重复打包大型表达矩阵。
- `notes/slide_materials_index.csv`：整合页码—源页—论文—图文件—数据表索引。
- `notes/panel_source_index.csv`：每个图元素的来源与类型。
- `notes/presenter_guide.zh.md`：逐页中文讲稿提示和出处；PPT备注页也有说明。
- `references/`：PAGER、Borra与Blood三篇参考论文。
- `slides.json`：独立元素布局与文本/表格数据。
- `code/rebuild_decks.py`：从本包内素材独立重建PPTX和配套PDF。
- `code/source_builder_*.py`：原始绘图代码快照，作为来源审计；完整重绘依赖原项目数据和环境，不是独立分析入口。
- `notes/verification.json`：文本/表格/图片来源与文件结构核对。
- `notes/editability_smoke_test.json`：修改标题、表格单元格、移动图形后保存并重新打开的测试结果。

## 重建

准备含python-pptx、pypdf、matplotlib的Python环境，在解压后的package目录运行：

```bash
python code/rebuild_decks.py
```

本项目使用现有`.venv`，长任务由tmux运行，GPU全部禁用。旧数据和旧报告不覆盖。
