import sys, re, gzip, functools, collections, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'code'))
from common import *   # ROOT, RES, P19, P25, KEGG, EXCL, ANCHOR_KEGG ...
DD = OUT/'drilldown'; TAB = DD/'tables'; FIGD = DD/'figures'; HTMLD = DD/'html'
KEY=[('cluster_tissue__C0__bone_marrow','BM C0'),('cluster_tissue__C3__spleen','Spleen C3'),('cluster_tissue__C4__bone_marrow','BM C4'),
     ('cluster_tissue__C5-2__bone_marrow','BM C5-2'),('cluster_tissue__C5-2__spleen','Spleen C5-2'),('cluster_tissue__C6__thymus','Thymus C6')]
KEYNAME=dict(KEY)
# drilled units: id -> short label
DRILL=[('KEGG:MAPK signaling pathway','KEGG MAPK'),('KEGG:IL-17 signaling pathway','KEGG IL-17'),
       ('KEGG:T cell receptor signaling pathway','KEGG TCR'),('KEGG:Th17 cell differentiation','KEGG Th17'),
       ('KEGG:TNF signaling pathway','KEGG TNF'),('GO:0035976','GO AP-1 complex'),('ENERGY','Energy module (G01)')]
DRILLID=[d[0] for d in DRILL]; DRILLNAME=dict(DRILL)
CAVEATS=('Exploratory, cell-level results. Each tissue x condition label is 3 mice pooled (n = 1 pool). '
         'Enrichment/annotation is not pathway activation. Co-membership edges are annotation, not physical interaction. '
         'No STRING/PPI data were used.')
def slug(s): return re.sub(r'[^A-Za-z0-9]+','_',s).strip('_')
