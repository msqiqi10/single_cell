"""Keep signature sensitivity and direct GO evidence distinct from the core result."""
import numpy as np, pandas as pd
from analyze import OUT,RES,KEY
p=RES/'tables';obs=pd.read_csv(p/'cell_metadata_scores.csv.gz',low_memory=False);contrasts=pd.read_csv(p/'cytotoxicity_contrasts.csv');go=pd.read_csv(p/'GO_ORA_all.csv.gz',low_memory=False)
kill=go[(go.analysis=='condition')&(go.definition=='original')&go.term.str.contains('cytotox|cell killing|granzyme|exocytosis|lymphocyte mediated',case=False)]
kill.to_csv(p/'GO_cytotoxicity_targeted_audit.csv',index=False)
focus=[('bone_marrow','C4'),('spleen','C3')];rows=[];cross=[]
for tissue,cluster in focus:
 d=pd.read_csv(RES/'de'/f'cluster_tissue__{cluster}__{tissue}.csv.gz')
 rows.append(d[d.gene.isin(['Prf1','Gzma','Gzmb','Ctla2a','Nkg7'])])
 m=(obs.tissue==tissue)&(obs[KEY]==cluster);vc=obs.loc[m,'stable_cluster'].value_counts();best=vc.index[0]
 r=contrasts[(contrasts.definition=='stable')&(contrasts.tissue==tissue)&(contrasts.cluster==best)&(contrasts.module=='PAGER_core')].iloc[0]
 cross.append({'tissue':tissue,'original_cluster':cluster,'best_overlap_stable_cluster':best,'overlap_cells':vc.iloc[0],'original_n_cells':m.sum(),'original_fraction_in_stable':vc.iloc[0]/m.sum(),'stable_total_n_cells':int((obs.stable_cluster==best).sum()),'stable_core_delta':r.delta_mean,'stable_core_q_all':r.q_all})
pd.concat(rows).to_csv(p/'focus_cytotoxic_gene_DE.csv',index=False);pd.DataFrame(cross).to_csv(p/'focus_stable_cluster_crosswalk.csv',index=False)
print(pd.DataFrame(cross).to_string(index=False))
print('Minimum original-cluster treatment cytotoxicity GO q:',kill.q_family.min())
