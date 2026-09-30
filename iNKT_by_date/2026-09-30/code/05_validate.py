import sys,json,random; sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from common import *
import numpy as np,pandas as pd
from scipy.stats import hypergeom
res=pd.read_csv(RES/'tables/ORA_all_terms_both_versions.csv.gz',keep_default_na=False,na_values=[''])
de=pd.read_csv(P25/'results/tables/all_20_condition_DE.csv.gz')
out={'exclusion_regex':EXCL_REGEX,'checks':{}}
random.seed(7)
for v in ['all','excl']:
    n=pd.read_csv(RES/'network'/v/'nodes.csv',keep_default_na=False,na_values=[''])
    sig=n[n.significant&(n.source!='GO CC anchor')].sample(3,random_state=7); rr=[]
    for t in sig.itertuples():
        u,d=t.sig_comparisons.split(';')[0].split('|'); row=res[(res.version==v)&(res.term_id==t.term_id)&(res.unit_id==u)&(res.direction==d)].iloc[0]
        x=de[de.unit_id==u]; 
        if v=='excl': x=x[~x.gene.str.match(EXCL_REGEX)]
        expr=set(x[(x.pct_expressing_tumor>0)|(x.pct_expressing_control>0)].gene)
        sel=set(x[(x.pvals_adj<=.05)&((x.logfoldchanges>=.25) if d=='T2_up' else (x.logfoldchanges<=-.25))].gene)
        # recompute term membership and full-family BH independently
        gm={}
        for l in (P19/'sources/GO_mouse_measured.gmt').read_text().splitlines():
            a,b,*g=l.split('\t'); gm[a]=set(g)
        import collections
        if t.term_id.startswith('KEGG:'):
            up=collections.defaultdict(set)
            for g in de.gene.unique(): up[g.upper()].add(g)
            for l in KEGG.read_text().splitlines():
                nm,_,*g=l.split('\t')
                if 'KEGG:'+nm==t.term_id: gm[t.term_id]=set().union(*(up.get(y,set()) for y in g))
            fam=[k for k in gm if k.startswith('KEGG:')]; src='KEGG'
        else:
            cat=pd.read_csv(P19/'results/tables/GO_term_catalog.csv').set_index('go_id'); ns=cat.loc[t.term_id,'namespace']; fam=list(cat.index[cat.namespace==ns]); src=ns
            fam=[k for k in fam if k in gm]
        if src=='KEGG':
            for l in KEGG.read_text().splitlines():
                nm,_,*g=l.split('\t'); 
                gm['KEGG:'+nm]=set().union(*(up.get(y,set()) for y in g))
        uni=expr&set().union(*(gm[k] for k in fam)); M=len(uni); q=sel&uni; n_=len(q)
        ps=[];tid=[]
        for k in fam:
            m=gm[k]&uni
            if 5<=len(m)<=500: ps.append(hypergeom.sf(len(q&m)-1,M,len(m),n_)); tid.append(k)
        ps=np.array(ps); o=np.argsort(ps); qq=np.empty(len(ps)); qq[o]=np.minimum(1,np.minimum.accumulate((ps[o]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1])
        q_re=qq[tid.index(t.term_id)]
        rr.append(dict(term=t.term_id,unit=u,direction=d,q_table=float(row.q),q_rederived=float(q_re),match=bool(np.isclose(q_re,row.q,rtol=1e-9,atol=1e-15))))
    out['checks'][f'{v}_q_rederivation']=rr
    A=[a for a in [*('KEGG:'+x for x in ANCHOR_KEGG),'GO:0035976','GO:0000165','GO:0070371','GO:0038066','GO:0007254','GO:0072538','GO:0032620']]
    out['checks'][f'{v}_all_12_anchors_present']=bool(set(A)<=set(n.term_id)) and n.is_anchor.sum()==12
    if v=='excl':
        bad=set()
        for col in ['hit_genes_union','measured_genes']:
            pass
        nn=pd.read_csv(RES/'network'/v/'nodes.csv',keep_default_na=False); 
        for col in ['hit_genes_union','measured_genes']:
            bad|={g for s in nn[col] for g in str(s).split(';') if g and EXCL.match(g)}
        E=pd.read_csv(RES/'network'/v/'edges.tsv',sep='\t'); bad|={g for s in E.shared_measured_genes for g in s.split(';') if EXCL.match(g)}
        Rr=res[res.version=='excl']; bad|={g for s in Rr.genes.astype(str) for g in s.split(';') if g and EXCL.match(g)}
        out['checks']['excl_no_excluded_genes_remaining']=len(bad)==0; out['checks']['excl_excluded_gene_hits']=sorted(bad)
        out['checks']['excl_n_genes_excluded_of_10670']=int(sum(bool(EXCL.match(g)) for g in de.gene.unique()))
    # edge rule
    E=pd.read_csv(RES/'network'/v/'edges.tsv',sep='\t'); out['checks'][f'{v}_edges_rule_ok']=bool(((E.JACCARD>=.25)&(E.n_shared>=3)).all())
    nodes_sig=n[~n.is_anchor]; out['checks'][f'{v}_nonanchor_nodes_all_sig']=bool(nodes_sig.significant.all())
out['checks']['reuse_vs_0919']=pd.read_csv(RES/'tables/ORA_reuse_check_vs_0919.csv').to_dict('records')[0]
out['network_summary']=json.load(open(RES/'network/network_summary.json'))
out['passed']=all(x['match'] for v in ['all','excl'] for x in out['checks'][f'{v}_q_rederivation']) and out['checks']['excl_no_excluded_genes_remaining'] and out['checks']['all_all_12_anchors_present'] and out['checks']['excl_all_12_anchors_present']
json.dump(out,open(RES/'validation.json','w'),indent=2,default=str); print(json.dumps({k:v for k,v in out.items() if k!='network_summary'},indent=1,default=str)[:2500])
