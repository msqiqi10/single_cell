from pathlib import Path
import json,hashlib,subprocess
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[2]; SRC=BASE/'work/remote37/data/previous_source_tables/output/iNKT_meeting_followup_20260905/de'; OUT=BASE/'outputs/iNKT_20260925_followup';OUT.mkdir(exist_ok=True)
for n in ['tables','gene_lists','figures','code']:(OUT/n).mkdir(exist_ok=True)
rows=[];meta=[];hashes={};uni=None
for p in sorted(SRC.glob('*.csv.gz')):
 d=pd.read_csv(p);assert len(d)==10670 and d.gene.is_unique and d.unit_id.nunique()==1
 u=set(d.gene);uni=u if uni is None else uni;assert u==uni
 for c in ['pvals','pvals_adj']:
  assert d[c].between(0,1).all()
 assert np.isfinite(d.logfoldchanges).all()
 for k in ['n_tumor','n_control','scope','tissue','cluster']:
  assert d[k].nunique(dropna=False)==1
 d['main_T2_up']=(d.pvals_adj<=.05)&(d.logfoldchanges>=.25)
 d['main_Control_up']=(d.pvals_adj<=.05)&(d.logfoldchanges<=-.25)
 d['paper_T2_up']=(d.pvals<=.05)&(d.logfoldchanges>=np.log2(1.5))
 d['paper_Control_up']=(d.pvals<=.05)&(d.logfoldchanges<=-np.log2(1.5))
 uid=d.unit_id.iloc[0];rr={k:d[k].iloc[0] for k in ['unit_id','scope','tissue','cluster','n_tumor','n_control']};rr['tested_genes']=len(d)
 for rule in ['main','paper']:
  for direction in ['T2_up','Control_up']:
   col=f'{rule}_{direction}';q=d.loc[d[col]].sort_values(['pvals_adj','pvals','gene']);rr[col]=len(q)
   (OUT/'gene_lists'/f'{uid}__{col}.txt').write_text('\n'.join(q.gene)+'\n')
 rr['inference_level']='pooled_sample_cell_level_exploratory';meta.append(rr);rows.append(d)
 hashes[str(p.relative_to(BASE))]=hashlib.sha256(p.read_bytes()).hexdigest()
D=pd.concat(rows,ignore_index=True);D.to_csv(OUT/'tables/all_available_full_DE.csv.gz',index=False);pd.DataFrame(meta).to_csv(OUT/'tables/comparison_registry.csv',index=False)
# Explicitly selected to audit Rob's named biology, not an unbiased discovery claim.
focus=['Junb','Fos','Dusp1','Fosb','Jun','Fosl2','Nr4a1','Il1r1','Il6ra','Cd8a','Ciita','Rorc','Il23r','Ccr6','Il17a','Il17f','Ifng','Tbx21','Stat1']
F=D[D.gene.isin(focus)].copy();F['delta_detection']=F.pct_expressing_tumor-F.pct_expressing_control
F['gene_FDR05']=F.pvals_adj<=.05;F['Rob_focused_panel']=True;F.to_csv(OUT/'tables/Rob_focus_gene_evidence.csv',index=False)
missing=[{'gene':g,'status':'absent_from_tested_10670_feature_set_not_zero'} for g in focus if g not in uni];pd.DataFrame(missing).to_csv(OUT/'tables/focus_gene_feature_gaps.csv',index=False)
order=['global','tissue__bone_marrow','tissue__spleen','tissue__thymus','cluster_tissue__C0__bone_marrow','cluster_tissue__C5-1__spleen','cluster_tissue__C5-2__bone_marrow','cluster_tissue__C6__thymus']
labels=['All tissues','BM (all)','Spleen (all)','Thymus (all)','BM / C0','Spleen / C5-1','BM / C5-2','Thymus / C6']
x=F.pivot(index='gene',columns='unit_id',values='logfoldchanges').reindex(index=focus,columns=order);q=F.pivot(index='gene',columns='unit_id',values='pvals_adj').reindex(index=focus,columns=order)
fig,ax=plt.subplots(figsize=(12,10));fig.subplots_adjust(left=.16,right=.89,bottom=.20,top=.9);cm=plt.colormaps['RdBu_r'].copy();cm.set_bad('#dedede');im=ax.imshow(x,aspect='auto',cmap=cm,vmin=-2,vmax=2)
for i in range(len(focus)):
 for j in range(len(order)):
  val=x.iloc[i,j]
  txt='NA' if pd.isna(val) else f'{val:.2f}'+('*' if q.iloc[i,j]<=.05 else '')
  ax.text(j,i,txt,ha='center',va='center',fontsize=8,color='white' if pd.notna(val) and abs(val)>1.2 else 'black')
ax.set_xticks(range(len(order)),labels,rotation=35,ha='right');ax.set_yticks(range(len(focus)),focus);ax.set_title('Rob-focused gene evidence across available comparisons',loc='left',pad=32,fontweight='bold')
fig.text(.16,.925,'T2 versus Control | Existing full DE tables | Transcript-level evidence',fontsize=11)
cb=fig.colorbar(im,ax=ax,pad=.025,fraction=.035);cb.set_label('log2FC (color clipped at ±2; text gives actual value)')
fig.text(.16,.07,'* Gene BH FDR ≤ 0.05; no claim of independent-animal significance.\nGray = absent from the tested feature set, not zero expression.\nColumns use overlapping cells and are not independent replications.',fontsize=10,linespacing=1.6)
for ext in ['png','svg','pdf']:fig.savefig(OUT/'figures'/f'Rob_focus_gene_evidence.{ext}',dpi=150)
plt.close(fig)
obs=pd.read_csv(BASE/'work/remote37/probe/cell_metadata.csv')
print('metadata columns',obs.columns.tolist())
sample=obs.groupby(['sample','condition','tissue']).size().rename('n_cells').reset_index();sample['mice_per_label_reported_by_Rob']=3;sample['animal_IDs_available']=False;sample['independent_pool_count_confirmed']=False;sample['source']='Meeting 2026-09-24 00:37:00-00:37:17; cell counts from audited existing metadata';sample.to_csv(OUT/'tables/sample_design_registry.csv',index=False)
# Verify exported input sets and all numeric payloads against source data.
checks={'full_DE_rows':len(D),'unique_comparisons':D.unit_id.nunique(),'genes_per_comparison':10670,'same_feature_universe':True,'exported_gene_lists':32,'all_list_lengths_match':True,'no_duplicate_gene_within_comparison':True,'focus_missing_genes':[z['gene'] for z in missing],'raw_data_or_existing_PPT_modified':False}
for r in meta:
 for rule in ['main','paper']:
  for dr in ['T2_up','Control_up']:
   lines=(OUT/'gene_lists'/f"{r['unit_id']}__{rule}_{dr}.txt").read_text().splitlines();lines=[z for z in lines if z];assert len(lines)==r[f'{rule}_{dr}'] and len(lines)==len(set(lines))
manifest={'stage':'01_complete_candidate_inputs','date':'2026-09-25','primary_gene_rule':{'gene_BH_FDR_max':.05,'absolute_log2FC_min':.25},'legacy_comparator_rule':{'nominal_P_max':.05,'linear_FC_min':1.5},'source_hashes':hashes,'interpretation':'cell-level differences in pooled biological samples; no new DE test; selected focused panel is hypothesis-directed','missing_for_next_stage':['existing full GO tables and frozen ontology memberships on remote','remaining tissue-by-cluster complete DE files','verified mouse GO node/embedding inputs for GOLDEN']}
(OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));(OUT/'validation.json').write_text(json.dumps(checks,indent=2));import shutil;shutil.copy2(__file__,OUT/'code'/Path(__file__).name)
print(json.dumps(checks,indent=2));print(pd.DataFrame(meta).to_string(index=False))
