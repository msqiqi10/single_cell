"""Preserve all existing DE outputs and prior candidate evidence alongside new work."""
from pathlib import Path
import shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/bio3';OUT=ROOT/'runs/meeting_followup'
status=pd.read_csv(OUT/'tables/all_comparison_status.csv');blocks=[]
for row in status[status.status=='computed'].itertuples():
    d=pd.read_csv(DATA/'output/iNKT_meeting_followup_20260905/de'/f'{row.unit_id}.csv.gz')
    for key,value in dict(unit_id=row.unit_id,tissue=row.tissue,cluster=row.cluster,scope=row.scope,strict_min20=row.strict_min20,n_T2=row.n_tumor,n_Ctrl=row.n_control).items():d[key]=value
    blocks.append(d)
pd.concat(blocks).to_csv(OUT/'tables/all_26_existing_DE_with_eligibility.csv.gz',index=False)
dest=OUT/'existing_evidence';dest.mkdir(exist_ok=True)
for name in ['discovery_BM_C4_focal_genes.csv','gene_candidates_all_sensitivities.csv','discovery_shortlist.csv','pathway_exact_driver_groups.csv']:
    shutil.copy2(DATA/'output/iNKT_discovery_20260905/tables'/name,dest/name)
for name in ['cytotoxicity_contrasts.csv','cytotoxicity_integrated_summary.csv','cytotoxicity_QC_sensitivity.csv','GO_cytotoxicity_targeted_audit.csv','focus_stable_cluster_crosswalk.csv','module_gene_coverage.csv']:
    shutil.copy2(DATA/'iNKT_by_date/2026-09-19/results/tables'/name,dest/name)
print('Preserved',len(blocks),'full comparisons and existing evidence')
