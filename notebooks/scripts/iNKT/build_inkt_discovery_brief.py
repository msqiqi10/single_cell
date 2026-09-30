"""Create a discovery-focused evidence brief from completed iNKT analyses.
No new DE tests: stricter post-screen prioritization is explicitly labelled.
"""
from pathlib import Path
import json, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from run_inkt_discovery import OUT, PREV, ROOT, sha256_file

FCS=['original_log2FC','stable_log2FC','QC_log2FC','stable_QC_log2FC']
QS=['original_FDR','stable_FDR','QC_FDR','stable_QC_FDR']
COLS=['Original cluster','Tissue recluster','QC matched','Recluster + QC']
FOCAL=['Tcf7','Itga4','Emb','Ly6a','S100a6','Cdca4']
BM='cluster_tissue__C4__bone_marrow'

def read(name):return pd.read_csv(OUT/'tables'/name,low_memory=False)
def save_table(d,name):d.to_csv(OUT/'tables'/name,index=False)
def text_page(title,paragraphs):
    fig=plt.figure(figsize=(12,8));fig.text(.065,.93,title,fontsize=21,weight='bold',va='top');y=.82
    for heading,body in paragraphs:
        fig.text(.065,y,heading,fontsize=14,weight='bold',va='top');y-=.043
        wrapped=textwrap.fill(body,106)
        fig.text(.065,y,wrapped,fontsize=11.5,va='top',linespacing=1.5);y-=.028*(wrapped.count('\n')+1)+.058
    fig.text(.065,.035,'iNKT discovery | 2026-09-05 | Cell-level exploratory evidence; sensitivity checks reuse the same samples.',fontsize=8,color='#555555')
    return fig

def heatmap(d,title,subtitle,labels=None):
    fig,ax=plt.subplots(figsize=(11,max(5.8,len(d)*.34+2.3)))
    arr=d[FCS].to_numpy(float);q=d[QS].to_numpy(float)
    im=ax.imshow(arr,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_xticks(range(4),COLS,fontsize=10);ax.set_yticks(range(len(d)),labels if labels is not None else d.gene,fontsize=10)
    for i in range(len(d)):
        for j in range(4):
            ax.text(j,i,f'{arr[i,j]:+.2f}'+('*' if q[i,j]<=.05 else ''),ha='center',va='center',fontsize=9,color='white' if abs(arr[i,j])>1.2 else 'black')
    fig.colorbar(im,ax=ax,pad=.03,label='log2FC: T2 versus Ctrl')
    ax.set_title(title+'\n'+subtitle,fontsize=13,pad=16)
    fig.text(.13,.045,'* Gene BH FDR <= 0.05 within each contrast. Columns are sensitivity analyses, not independent replicates.',fontsize=9)
    fig.subplots_adjust(left=.22,right=.93,top=.80,bottom=.15)
    return fig

def drivers(p):
    rows=[]
    for r in p.itertuples():
        ds=[pd.read_csv(PREV/f'de/{r.unit_id}.csv.gz').set_index('gene'),pd.read_csv(PREV/f'de/{r.stable_unit}.csv.gz').set_index('gene'),pd.read_csv(OUT/f'de_qc/{r.unit_id}.csv.gz').set_index('gene'),pd.read_csv(OUT/f'de_qc/{r.stable_unit}.csv.gz').set_index('gene')]
        for gene in r.drivers.split(';'):
            row={'unit_id':r.unit_id,'term':r.term,'gene':gene}
            for k,df in enumerate(ds):row[FCS[k]]=df.loc[gene,'logfoldchanges'];row[QS[k]]=df.loc[gene,'pvals_adj']
            rows.append(row)
    return pd.DataFrame(rows)

def main():
    all27=read('discovery_shortlist.csv')
    priority=all27[(all27[QS]<=.05).all(axis=1)&all27.same_sign_all_5].copy()
    priority['priority_rule']='Post-screen: existing local/non-dominant screen + all four gene FDR<=0.05 + same sign in five cell selections'
    save_table(priority,'discovery_priority_all4_FDR05.csv')
    p=read('pathway_candidates_all.csv')
    p=p[p.stable_and_QC_supported&(p.no_heat_FDR<=.05)&(p.FDR_across_units<=.05)&(p.n_other_drivers>=3)&(p.n_other_drivers/p.n_drivers>=.5)].copy()
    screen=read('functional_screen_all.csv.gz').set_index(['unit_id','analysis','library','sensitivity','direction','term'])
    for i,r in p.iterrows():
        key=(r.stable_unit,'QC_matched',r.library,'full',r.direction,r.term)
        p.loc[i,'stable_QC_FDR']=screen.loc[key,'fdr'] if key in screen.index else np.nan
        p.loc[i,'program_group']='CCT_TriC' if 'CCT/TriC' in r.term else ('splicing_gene_subset' if r.term=='Spliceosome' else 'shared_cytoskeletal_response_genes')
        p.loc[i,'interpretation_status']='shared driver genes; not evidence of atherosclerosis' if 'atherosclerosis' in r.term else ('ORA support; whole-set GSEA not significant' if r.term=='Spliceosome' else 'protein-folding-associated expression; not a functional activity measurement')
        p.loc[i,'legacy_PPT_context']='Slide 14: old c1_thymus; HSPA1A/HSPA1B/HSPA8 only; old and current cluster numbers not equated' if r.term=='Spliceosome' else 'No exact matching pathway title in extracted legacy PPT tables'
        gp=PREV/f'enrichment/{r.unit_id}__full__KEGG_Mouse_2019__GSEA.csv.gz'
        if r.library=='KEGG_Mouse_2019' and gp.exists():
            g=pd.read_csv(gp).set_index('Term')
            if r.term in g.index:
                p.loc[i,'original_GSEA_NES']=g.loc[r.term,'NES'];p.loc[i,'original_GSEA_FDR']=g.loc[r.term,'FDR q-val']
                g2=pd.read_csv(PREV/f'enrichment/{r.unit_id}__without_heat_shock__KEGG_Mouse_2019__GSEA.csv.gz').set_index('Term')
                p.loc[i,'no_heat_GSEA_FDR']=g2.loc[r.term,'FDR q-val']
    save_table(p,'discovery_functional_priority.csv')
    dr=drivers(p);save_table(dr,'discovery_pathway_driver_evidence.csv')
    focal=priority[priority.unit_id.eq(BM)].set_index('gene').loc[FOCAL].reset_index()
    save_table(focal,'discovery_BM_C4_focal_genes.csv')
    qa=read('QC_matching_audit.csv').set_index('unit_id');cw=read('original_stable_crosswalk.csv')
    mapped=cw[cw.tissue.eq('bone_marrow')&cw.original_cluster.eq('C4')].sort_values('overlap',ascending=False).iloc[0]
    il7=read('gene_candidates_all_sensitivities.csv');il7=il7[il7.unit_id.eq('cluster_tissue__C5-1__spleen')&il7.gene.eq('Il7r')].iloc[0]
    cct=p[p.term.str.contains('Tubulin')].iloc[0];spl=p[p.term.eq('Spliceosome')].iloc[0]
    fig1=heatmap(focal,'Priority 1: localized BM C4 expression changes',f"Ctrl {int(qa.loc[BM,'n_Ctrl'])} / T2 {int(qa.loc[BM,'n_T2'])} cells; QC matched {int(qa.loc[BM,'n_each_matched'])} per condition\n{mapped.fraction_original:.1%} of original C4 maps to BM_S1; BM_S1 also contains other original groups")
    fig1.savefig(OUT/'figures/discovery_BM_C4_focal.png',dpi=180);fig1.savefig(OUT/'figures/discovery_BM_C4_focal.pdf')
    cc=dr[dr.term.eq(cct.term)].copy();fig2=heatmap(cc,'Priority 2: lower CCT/TriC-associated expression in spleen C3',f"ORA FDR {cct.FDR:.4g}; QC {cct.QC_FDR:.4g}; recluster {cct.stable_FDR:.4g}; no Hsp/Dnaj {cct.no_heat_FDR:.4g}\nActin- and tubulin-folding term names share Cct2/3/4/6a/8 drivers")
    fig2.savefig(OUT/'figures/discovery_CCT_driver_evidence.png',dpi=180)
    sd=dr[dr.term.eq('Spliceosome')&~dr.gene.str.match(r'(?i)^(Hsp|Dnaj)')].copy()
    fig3=heatmap(sd,'Priority 3: an upregulated splicing-gene subset in BM C0',f"ORA FDR {spl.FDR:.4g}; no Hsp/Dnaj {spl.no_heat_FDR:.4g}; whole-set GSEA q = {spl.original_GSEA_FDR:.3f}\nThis does not establish a whole-pathway shift or altered splice isoforms")
    fig3.savefig(OUT/'figures/discovery_splicing_driver_evidence.png',dpi=180)
    labels=[f'{r.gene} | {r.cluster} / '+('BM' if r.tissue=='bone_marrow' else r.tissue) for r in priority.itertuples()]
    fig4=heatmap(priority,'All stricter-priority gene records',f'{len(priority)} gene x cluster records; all four FDR <= 0.05 and five matched selections keep the sign',labels)
    fig4.savefig(OUT/'figures/discovery_priority_all17.png',dpi=180)
    summary=text_page('What the current iNKT pipeline finds in this dataset',[
        ('1. BM C4: a localized expression-state candidate','Tcf7, Itga4 and Emb decrease in T2; Ly6a, S100a6 and Cdca4 increase. These changes pass all four DE sensitivities. Tcf7, Emb and Ly6a are not significant in the whole-BM comparison. This motivates a focused population-state question; it does not prove differentiation or migration.'),
        ('2. Spleen C3: CCT/TriC-associated genes decrease','Cct2, Cct3, Cct4, Cct6a, Cct8, Tuba1a and Tubb4b drive the tubulin-folding term. ORA remains significant after QC matching, tissue reclustering and removal of Hsp/Dnaj genes. Related actin-folding terms share the same core drivers.'),
        ('3. BM C0: a splicing-gene subset merits lower-priority follow-up','Spliceosome ORA survives heat-shock exclusion, driven by 11 non-Hsp/Dnaj genes. However, whole-set GSEA is not significant (q = 0.654; no heat-shock q = 0.801). The older PPT already named Spliceosome, but with three HSPA drivers in an old thymus group.'),
        ('What changed in prioritization','Spleen C5-1 Il7r keeps its positive direction, but QC-matched FDR is 0.548. It is retained as an exploratory lead rather than a first-priority candidate. None of these findings has independent animal-level validation.')])
    methods=text_page('Scope, inference and the next focused question',[
        ('What was actually run','The current iNKT workflow was extended with open DE/KEGG/Reactome screening. 35 contrasts were audited; 29 retained at least 20 cells per condition after joint QC-stratum matching. Gene tests use 10,670 measured genes and uncorrected expression. This run does not substitute the CML scVI/SCORPION workflow.'),
        ('How priority was assigned','The broader screen keeps all 806 supported gene-by-group records, including 27 local records outside the explicit Hsp/ribosomal/respiratory families. An additional post-screen rule retains 17 records with FDR <= 0.05 in all four DE analyses and consistent direction in five cell selections. The six BM C4 genes highlighted here are a readable subset, not the complete passing set.'),
        ('What these checks can and cannot establish','Original and reclustered groups overlap; matched cells and five selections come from the same samples. Only one sample label exists per tissue and condition. These are cell-level exploratory statistics, not independent treatment replicates or causal effects. A local versus tissue effect comparison is descriptive, not an interaction test.'),
        ('Next focused question within this project','First establish whether the BM C4 expression pattern is reproduced in independently identified iNKT samples and reflects a state shift or cell-mixture change. Next follow the spleen C3 folding-gene program. Keep the BM C0 splicing subset exploratory until whole-set or junction-level evidence is available. Velocity remains unavailable without spliced/unspliced inputs.')])
    with PdfPages(OUT/'discovery_brief_20260905.pdf') as pdf:
        for fig in [summary,fig1,fig2,fig3,fig4,methods]:pdf.savefig(fig);plt.close(fig)
    rows=['| 基因 | 原群 log2FC / FDR | QC 匹配 log2FC / FDR | 组织整体 FDR | Ctrl→T2 检测率 |','|---|---|---|---:|---|']
    for r in focal.itertuples():rows.append(f'| {r.gene} | {r.original_log2FC:+.3f} / {r.original_FDR:.3g} | {r.QC_log2FC:+.3f} / {r.QC_FDR:.3g} | {r.tissue_FDR:.3g} | {r.detected_Ctrl:.1%}→{r.detected_T2:.1%} |')
    md=f'''# 这批 iNKT 数据中，当前最值得追踪的是什么？

沿用最新 iNKT 分析流程，已完成全基因、完整 KEGG/Reactome 功能库筛选，并新增 QC 匹配、组织重聚类一致性和 5 次细胞抽样敏感性检查。现在优先追踪 **骨髓 C4 的局部表达变化**和**脾脏 C3 的 CCT/TriC 相关表达降低**；骨髓 C0 的剪接基因子集列为较低优先级线索。

[六页发现报告（PDF，英文图文）](discovery_brief_20260905.pdf) · [17 条严格优先记录](tables/discovery_priority_all4_FDR05.csv) · [功能候选及方法差异](tables/discovery_functional_priority.csv)

## 1. 骨髓 C4：Tcf7、Itga4、Emb 下调，Ly6a、S100a6、Cdca4 上调

原群 Ctrl={int(qa.loc[BM,'n_Ctrl'])}、T2={int(qa.loc[BM,'n_T2'])} 个细胞；QC 匹配后各 {int(qa.loc[BM,'n_each_matched'])} 个。下面六个基因在原分群、组织重聚类、原群 QC 匹配、重聚类后 QC 匹配中均同向且基因 FDR≤0.05，五次匹配抽样也同向。图和完整表保留四项效应与 FDR。

{chr(10).join(rows)}

Tcf7、Emb、Ly6a 在整个骨髓的比较中未达 FDR 0.05，说明这些信号可能在整体平均中不明显；这不是正式的 cluster×condition 交互作用检验。六个基因都未见于旧 PPT 提取的 DEG 表，但不能据此声称文献首创。

原 C4 的 {int(mapped.overlap)}/{int(mapped.overlap/mapped.fraction_original)} 个细胞（{mapped.fraction_original:.2%}）落入骨髓重聚类 S1；原 C4 只占 S1 的 {mapped.fraction_stable:.1%}，因此两种分群并非完全等同。Cd8b1 也进入严格优先表，但检测率较低，需考虑群内细胞组成，未作为主线。

当前可提出的是“骨髓这个细胞群是否出现特定表达状态变化”。Itga4 的黏附/迁移相关注释可帮助形成假设（[NCBI 小鼠 Itga4](https://www.ncbi.nlm.nih.gov/gene/16401)），现有转录变化不能直接证明迁移、分化方向或功能改变。

![骨髓 C4 六个候选的四项效应](figures/discovery_BM_C4_focal.png)

## 2. 脾脏 C3：CCT/TriC 相关表达降低

T2 下调基因中，微管蛋白折叠相关条目的 ORA FDR={cct.FDR:.4g}，跨对照校正后 FDR={cct.FDR_across_units:.4g}；QC 匹配后 {cct.QC_FDR:.4g}、重聚类后 {cct.stable_FDR:.4g}、去 Hsp/Dnaj 后 {cct.no_heat_FDR:.4g}。驱动基因为 **Cct2、Cct3、Cct4、Cct6a、Cct8、Tuba1a、Tubb4b**。

这支持优先检查蛋白折叠/细胞骨架相关的表达程序。相关 actin folding 条目共享五个 Cct 基因，应视为同一组核心证据；不能把多个通路名当作多个独立机制。通路的生物学定义见 [Reactome CCT/TriC tubulin folding](https://reactome.org/content/detail/R-HSA-389960)。本分析未测量蛋白折叠活性，也不能推出杀伤功能下降。这里使用完整 Reactome ORA；没有新增 Reactome GSEA。

![CCT/TriC 驱动基因](figures/discovery_CCT_driver_evidence.png)

## 3. 骨髓 C0：剪接基因子集上调，但整条通路证据不足

Spliceosome ORA FDR={spl.FDR:.4g}，跨对照 FDR={spl.FDR_across_units:.4g}；QC 匹配 {spl.QC_FDR:.4g}、重聚类 {spl.stable_FDR:.4g}、去 heat-shock {spl.no_heat_FDR:.4g}。去除 Hsp/Dnaj 后还有 **Lsm5、Lsm6、Lsm7、Magohb、Snrpa、Snrpb2、Snrpe、Snrpf、Snrpg、Srsf5、Zmat2** 这 11 个驱动基因。

但已有全排序 GSEA 不显著：NES={spl.original_GSEA_NES:.3f}，q={spl.original_GSEA_FDR:.3f}；去 heat-shock 后 q={spl.no_heat_GSEA_FDR:.3f}。因此只能说一组剪接相关基因值得追踪，不能说整个 spliceosome 一致激活，更不能说已发现可变剪接改变。

旧 PPT 第 14 页已有 Spliceosome，来源是旧 c1_thymus 的 HSPA1A/HSPA1B/HSPA8。这里增加的是当前骨髓 C0 中不依赖这些 heat-shock 基因的证据，不是首次看到这个通路名；旧新 cluster 编号不作身份对应。

## 需要下调的已有候选

脾脏 C5-1 的 Il7r 原始 log2FC={il7.original_log2FC:+.3f}、FDR={il7.original_FDR:.4g}，QC 匹配后 log2FC={il7.QC_log2FC:+.3f}、FDR={il7.QC_FDR:.3f}；重聚类后 QC FDR={il7.stable_QC_FDR:.3f}。方向仍为正，但支持减弱，暂不放在首要验证位置。匹配减少样本量，也可能调整掉部分真实生物差异；不能据此判定 Il7r 没变化。

## 筛选范围、优先级和限制

- 全部 35 个对照进入 QC 审计，其中 29 个匹配后仍各有至少 20 个细胞并完成全基因 DE；6 个不足的保留明确状态。比较使用 10,670 个检测基因；HVG 不作为 DE 背景。DE 使用原表达而非 Harmony 校正值。
- [原有探索层](tables/gene_candidates_supported.csv)共 806 条基因×细胞群记录；排除明确 Hsp/Dnaj、核糖体、呼吸链/线粒体家族并考虑局部效应后有 27 条。新增的**事后优先规则**要求四项基因 FDR 均≤0.05，五次细胞选择同号，得到 17 条。17 不是独立发现数量，也没有将四个相关检验合成为新的 p 值。
- 基因 FDR 在每个对照内校正。功能筛选保留整个合格库及每对照 BH，并增加按分析/库/方向/去 heat-shock 分组的跨对照 BH。功能优先表的入选依据是 ORA，不代表 GSEA 一定支持。人为的基因家族过滤只是帮助阅读，不是严密的生物过程分类。
- 功能优先表还保留 “Fluid shear stress and atherosclerosis”：它由 Actb/Actg1/Dusp1/Fos/Rac2 驱动，不能据此推断动脉粥样硬化，未将该疾病名称作为发现。
- 目前每组织×条件仅一个 sample 标签，缺少可识别的动物重复。同一数据的重聚类、匹配、抽样只能说明敏感性，不能代替独立生物重复，当前 p/FDR 不支持推广到动物总体的治疗效应。
- 当前 pipeline 使用最新 iNKT 流程；未转用 CML 的 scVI/SCORPION，也未新增 velocity。velocity 仍需 spliced/unspliced 或可生成这些计数的原始数据。

下一步应围绕骨髓 C4 的表达状态候选与脾脏 C3 的 CCT 程序安排独立样本验证；首先确认细胞身份和样本重复，区分群内状态变化与细胞组成变化。剪接基因子集作为第二梯队，避免将范围扩大到没有现有证据的机制。

## 可复核文件

- [全部候选及四项敏感性](tables/gene_candidates_all_sensitivities.csv)
- [17 条严格优先记录](tables/discovery_priority_all4_FDR05.csv)
- [功能优先证据（含 GSEA 状态）](tables/discovery_functional_priority.csv)
- [每个通路驱动基因的四项 DE 结果](tables/discovery_pathway_driver_evidence.csv)
- [完整功能检验](tables/functional_screen_all.csv.gz)
- [QC 匹配审计](tables/QC_matching_audit.csv)
- [选中细胞 ID](tables/QC_selected_cell_ids.csv.gz)
- [五次抽样结果](tables/candidate_cell_subsampling_summary.csv)
- [来源、筛选与校验记录](discovery_brief_manifest.json)
'''
    (OUT/'DISCOVERIES.zh.md').write_text(md)
    # Keep the original screen report intact: link the refined brief at the top.
    rp=OUT/'README.md';old=rp.read_text();link='> 最新发现解读与更严格的优先层：[DISCOVERIES.zh.md](DISCOVERIES.zh.md)。原有探索表和筛选规则保留如下。\n\n'
    if not old.startswith('> 最新发现解读'):rp.write_text(link+old)
    manifest={'status':'built_pending_verification','pipeline_manifest_SHA256':sha256_file(OUT/'manifest.json'),'script_SHA256':sha256_file(Path(__file__)),'strict_priority_rule':'post-screen local non-dominant candidates: all four gene FDR<=0.05 and five matched selections same sign','strict_gene_records':len(priority),'strict_unique_genes':priority.gene.nunique(),'strict_records_not_in_legacy_DEG_tables':int((~priority.present_in_legacy_DEG_tables).sum()),'functional_priority_ORA_records':len(p),'functional_priority_ORA_groups':p.program_group.nunique(),'GSEA_splicing_significant':False,'main_figures_are_selected_examples':True,'pdf_pages':6,'biological_replicates_identified':False,'files_SHA256':{str(x.relative_to(OUT)):sha256_file(x) for x in [OUT/'DISCOVERIES.zh.md',OUT/'discovery_brief_20260905.pdf',OUT/'tables/discovery_priority_all4_FDR05.csv',OUT/'tables/discovery_functional_priority.csv',OUT/'tables/discovery_pathway_driver_evidence.csv']}}
    (OUT/'discovery_brief_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
