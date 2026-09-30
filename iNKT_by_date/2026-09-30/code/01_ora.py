"""Step 1: ORA for GO BP / GO MF / GO CC (anchor only) / KEGG_2019_Mouse, versions all (a) and exclude (b).
Same rules as 2026-09-19/code/analyze.py: hypergeom sf, term size 5-500 in universe, BH within comparison x direction x source.
Gene selection: pvals_adj<=0.05 & |logFC|>=0.25, T2_up / Ctrl_up separately. Universe = expressed genes (pct>0 either arm) & source-annotated genes."""
import sys; sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from common import *
import numpy as np, pandas as pd, collections
from scipy.stats import hypergeom
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty(len(p)); q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]); return q
de=pd.read_csv(P25/'results/tables/all_20_condition_DE.csv.gz')
cat=pd.read_csv(P19/'results/tables/GO_term_catalog.csv').set_index('go_id')
gmt={}
for l in (P19/'sources/GO_mouse_measured.gmt').read_text().splitlines():
    t,n,*g=l.split('\t'); gmt[t]=set(g)
measured=set(de.gene.unique())
libs={ns:{t:gmt[t] for t in cat.index[cat.namespace==ns] if t in gmt} for ns in ['BP','MF','CC']}
names={t:cat.loc[t,'term'] for t in gmt}
# KEGG: uppercase symbols mapped case-insensitively onto measured mouse symbols
up=collections.defaultdict(set)
for g in measured: up[g.upper()].add(g)
kegg={}; kmap=[]
for l in KEGG.read_text().splitlines():
    n,_,*g=l.split('\t'); mapped=set().union(*(up.get(x,set()) for x in g))
    kegg['KEGG:'+n]=mapped; names['KEGG:'+n]=n; kmap.append(dict(term=n,n_gmt=len(g),n_mapped_to_measured=len(mapped)))
pd.DataFrame(kmap).to_csv(RES/'tables/KEGG_gene_mapping_summary.csv',index=False)
libs['KEGG']=kegg
def ora(query,universe,sets,src):
    query=set(query)&universe; M=len(universe); n=len(query); rows=[]
    for t,gs in sets.items():
        m=gs&universe; K=len(m)
        if 5<=K<=500:
            hit=query&m; rows.append((t,names[t],src,M,K,n,len(hit),';'.join(sorted(hit))))
    d=pd.DataFrame(rows,columns=['term_id','term','source','M','K','n','k','genes'])
    if len(d):
        d['pvalue']=hypergeom.sf(d.k.to_numpy()-1,M,d.K.to_numpy(),n) if n else 1.
        d['q']=bh(d.pvalue)
    return d
out=[]; status=[]
for version in ['all','excl']:
    keep=(lambda g:True) if version=='all' else (lambda g:not EXCL.match(g))
    for uid,d in de.groupby('unit_id',sort=False):
        d=d[d.gene.map(keep)]
        expressed=set(d.loc[(d.pct_expressing_tumor>0)|(d.pct_expressing_control>0),'gene'])
        for direction in ['T2_up','Ctrl_up']:
            sel=set(d.loc[(d.pvals_adj<=.05)&((d.logfoldchanges>=.25) if direction=='T2_up' else (d.logfoldchanges<=-.25)),'gene'])
            for src,lib in libs.items():
                ann=set().union(*lib.values()); uni=expressed&ann
                r=ora(sel,uni,lib,src)
                r['version']=version; r['unit_id']=uid; r['direction']=direction; out.append(r)
                status.append(dict(version=version,unit_id=uid,direction=direction,source=src,n_selected=len(sel),n_background=len(uni),n_query_in_background=len(sel&uni),n_terms_tested=len(r),n_sig_q05=int((r.q<=.05).sum()) if len(r) else 0))
res=pd.concat(out,ignore_index=True)
res.to_csv(RES/'tables/ORA_all_terms_both_versions.csv.gz',index=False)
pd.DataFrame(status).to_csv(RES/'tables/ORA_status.csv',index=False)
# cross-check version all vs 2026-09-19 GO_ORA_all
old=pd.read_csv(P19/'results/tables/GO_ORA_all.csv.gz',low_memory=False,keep_default_na=False)
old=old[(old.analysis=='condition')&old.definition.isin(['original','tissue'])&old.namespace.isin(['BP','MF'])]
m=res[(res.version=='all')&res.source.isin(['BP','MF'])].merge(old[['go_id','unit_id','direction','namespace','q_family','k','K','M']],left_on=['term_id','unit_id','direction','source'],right_on=['go_id','unit_id','direction','namespace'],how='inner')
print('matched',len(m),'of',(res.version=='all').sum(),'; max|q diff|',(m.q-m.q_family.astype(float)).abs().max(),'; k/K/M mismatch',((m.k_x!=m.k_y)|(m.K_x!=m.K_y)|(m.M_x!=m.M_y)).sum())
pd.DataFrame([dict(n_matched=len(m),max_abs_q_diff=float((m.q-m.q_family.astype(float)).abs().max()),n_kKM_mismatch=int(((m.k_x!=m.k_y)|(m.K_x!=m.K_y)|(m.M_x!=m.M_y)).sum()))]).to_csv(RES/'tables/ORA_reuse_check_vs_0919.csv',index=False)
print(pd.DataFrame(status).groupby(['version','source']).n_sig_q05.sum())
