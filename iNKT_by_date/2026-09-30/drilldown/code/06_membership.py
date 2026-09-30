"""Step 6: full membership for drilled pathways, per-gene DE, annotation frequency, gene co-membership edges, bridge table."""
from common_dd import *
import numpy as np, pandas as pd, itertools
de=pd.read_csv(P25/'results/tables/all_20_condition_DE.csv.gz',usecols=['unit_id','gene','logfoldchanges','pvals_adj','strict_estimable_min20'])
assert de.strict_estimable_min20.all()   # all 20 (19 listed unit ids + tissue-level) eligible
de['deg']=(de.pvals_adj<=.05)&(de.logfoldchanges.abs()>=.25)
measured=set(de.gene.unique())
nodes=pd.read_csv(RES/'network/excl/nodes.csv',keep_default_na=False,na_values=[''])
nodes['mg']=nodes.measured_genes.fillna('').map(lambda s:set(s.split(';')) if s else set())
# ---- full KEGG (frozen gmt)
up=collections.defaultdict(set)
for g in measured: up[g.upper()].add(g)
kfull={}
for l in KEGG.read_text().splitlines():
    n,_,*g=l.split('\t'); s=set()
    for x in g: s|= up[x] if x in up else {x}   # unmapped keep gmt symbol (uppercase) -> not detected
    kfull['KEGG:'+n]=s
# ---- full GO AP-1 complex (GAF, same rules as 09-19 analyze.py: non-NOT, non-ND, taxon 10090, is_a/part_of propagation)
SRC=P19/'sources'
terms={};cur=None
for line in (SRC/'go-basic.obo').read_text().splitlines():
    if line=='[Term]': cur={'parents':[],'alt':[]};continue
    if line.startswith('['): cur=None;continue
    if cur is None: continue
    if line.startswith('id: '): cur['id']=line[4:];terms[line[4:]]=cur
    elif line.startswith('alt_id: '): cur['alt'].append(line[8:])
    elif line.startswith('is_a: '): cur['parents'].append(line[6:].split()[0])
    elif line.startswith('relationship: part_of '): cur['parents'].append(line.split()[2])
    elif line=='is_obsolete: true': cur['obsolete']=True
terms={k:v for k,v in terms.items() if not v.get('obsolete')}
alt={a:k for k,v in terms.items() for a in v['alt']}
@functools.lru_cache(None)
def anc(t):
    r={t}
    for p in terms[t]['parents']:
        if p in terms: r|=anc(p)
    return frozenset(r)
ap1=set(); direct=collections.defaultdict(set)
with gzip.open(SRC/'mgi.gaf.gz','rt') as f:
    for line in f:
        if line[0]=='!': continue
        x=line.rstrip('\n').split('\t')
        if len(x)<15 or 'NOT' in x[3].split('|') or x[6]=='ND' or x[12].split('|')[0]!='taxon:10090': continue
        t=alt.get(x[4],x[4])
        if t in terms and 'GO:0035976' in anc(t): ap1.add(x[2]); direct[x[2]].add(t)
gmt_ap1=None
for l in (P19/'sources/GO_mouse_measured.gmt').read_text().splitlines():
    if l.startswith('GO:0035976\t'): gmt_ap1=set(l.split('\t')[2:])
print('AP-1 GAF symbols',len(ap1),'; frozen measured gmt',len(gmt_ap1),'; gmt subset of GAF',gmt_ap1<=ap1)
extra=[g for g in ap1 if g in measured and g not in gmt_ap1]; print('GAF genes measured but not in frozen gmt (alias mapping?):',extra)
full={k:kfull[k] for k in DRILLID if k.startswith('KEGG')}
full['GO:0035976']=ap1|gmt_ap1
gsn=nodes[nodes.gnn_semantic=='G01']
_sets=[kfull[t] if t in kfull else m for t,m in zip(gsn.term_id,gsn.mg)]
_c=collections.Counter(g for m in _sets for g in m)
full['ENERGY']={g for g,c in _c.items() if c>=3}   # core = gene annotated in >=3 of the 15 G01 terms (the union is broader: disease/endocannabinoid terms)
print('energy union',len(_c),'core',len(full['ENERGY']))
energy_terms=list(gsn.term_id)
# ---- annotation frequency over significant nodes (b)
sig=nodes[nodes.significant.astype(str)=='True']
def freq(g): return sum(g in m for m in sig.mg), sum(g in m for m in nodes.mg)
gene_terms=collections.defaultdict(list)
for t,m in zip(nodes.term_id,nodes.mg):
    for g in m: gene_terms[g].append(t)
# ---- DE wide
wide=pd.DataFrame(index=sorted(measured))
for u,s in KEY:
    d=de[de.unit_id==u].set_index('gene'); wide[s+'|lfc']=d.logfoldchanges; wide[s+'|fdr']=d.pvals_adj; wide[s+'|deg']=d.deg
degany=de[de.deg].groupby('gene').unit_id.apply(lambda s:';'.join(sorted(s)))
ndegany=de[de.deg].groupby('gene').unit_id.nunique()
# ---- membership long table
drill_of=collections.defaultdict(list)
for k in DRILLID:
    for g in full[k]: drill_of[g].append(DRILLNAME[k])
rows=[]
for k in DRILLID:
    for g in sorted(full[k]):
        det=g in measured; r=dict(pathway=DRILLNAME[k],pathway_id=k,gene=g,detected_in_data='yes' if det else 'no',
            excluded_in_version_b=bool(EXCL.match(g)))
        for u,s in KEY:
            if det and g in wide.index:
                r[f'{s} log2FC']=wide.at[g,s+'|lfc']; r[f'{s} FDR']=wide.at[g,s+'|fdr']; r[f'{s} DEG']='yes' if wide.at[g,s+'|deg'] else 'no'
            else: r[f'{s} log2FC']=np.nan; r[f'{s} FDR']=np.nan; r[f'{s} DEG']='NA'
        r['n_key_comparisons_DEG']=sum(r[f'{s} DEG']=='yes' for _,s in KEY)
        r['DEG_all20_units']=degany.get(g,''); 
        n=len(drill_of[g]); r['n_drilled_pathways']=n; r['drilled_pathways']=';'.join(drill_of[g])
        a,b=freq(g); r['n_significant_nodes_containing']=a; r['n_nodes_incl_anchors_containing']=b
        r['shared_or_specific']='shared' if n>=2 else 'specific'
        rows.append(r)
mem=pd.DataFrame(rows); mem.to_csv(TAB/'pathway_gene_membership.csv',index=False)
summ=[]
for k in DRILLID:
    m=mem[mem.pathway_id==k]
    s=dict(pathway=DRILLNAME[k],n_members=len(m),n_detected=int((m.detected_in_data=='yes').sum()),n_shared=int((m.shared_or_specific=='shared').sum()))
    for _,sh in KEY: s['DEG '+sh]=int((m[f'{sh} DEG']=='yes').sum())
    s['DEG_any_key']=int((m.n_key_comparisons_DEG>0).sum()); summ.append(s)
summ=pd.DataFrame(summ); summ.to_csv(TAB/'pathway_drill_summary.csv',index=False); print(summ.to_string())
# ---- gene co-membership network
universe=sorted(set().union(*full.values()))
tm={}  # term membership used for edges
for t,m in zip(nodes.term_id,nodes.mg): tm[t]=full[t] if t in full else (kfull[t]&(set(m)|set()) if False else m)
gt=collections.defaultdict(set)
for t,m in tm.items():
    for g in m:
        if g in set(universe): gt[g].add(t)
inv=collections.defaultdict(list)
for g,ts in gt.items():
    for t in ts: inv[t].append(g)
cnt=collections.Counter()
for t,gs in inv.items():
    for a,b in itertools.combinations(sorted(gs),2): cnt[(a,b)]+=1
edges=[(a,b,w) for (a,b),w in cnt.items() if w>=2]
print('universe genes',len(universe),'edges w>=2',len(edges),'w>=3',sum(e[2]>=3 for e in edges),'w>=5',sum(e[2]>=5 for e in edges))
pd.DataFrame(edges,columns=['gene_a','gene_b','n_shared_terms']).to_csv(TAB/'gene_comembership_edges.csv',index=False)
json.dump(dict(universe=universe,gene_terms={g:sorted(v) for g,v in gt.items()},drill_of=drill_of,energy_terms=energy_terms),open(TAB/'_gene_net_aux.json','w'))
# ---- bridge analysis
def lab(k): return DRILLNAME[k]
pairs=[]
ap='GO:0035976'; ks=[k for k in DRILLID if k.startswith('KEGG')]
order=[(ap,k) for k in ks]+list(itertools.combinations(ks,2))+[('ENERGY',k) for k in [ap]+ks]
brows=[]
for a,b in order:
    sh=sorted(full[a]&full[b]); det=[g for g in sh if g in measured]
    r=dict(pathway_A=lab(a),pathway_B=lab(b),n_A=len(full[a]),n_B=len(full[b]),n_shared=len(sh),jaccard=round(len(sh)/len(full[a]|full[b]),4),
           shared_genes=';'.join(sh),n_shared_detected=len(det))
    dg=[g for g in det if any(wide.at[g,s+'|deg'] for _,s in KEY)]
    r['n_shared_DEG_any_key']=len(dg); r['shared_DEG_genes']=';'.join(dg)
    def sgn(g):
        v=[wide.at[g,s+'|lfc'] for _,s in KEY if wide.at[g,s+'|deg']]; return '+' if np.mean(v)>0 else '-'
    r['shared_DEG_genes_signed']=';'.join(g+'('+sgn(g)+')' for g in dg)
    for _,s in KEY:
        gs=[g for g in det if wide.at[g,s+'|deg']]
        r[f'shared_DEG_{s}']=';'.join(f'{g}({wide.at[g,s+"|lfc"]:+.2f})' for g in gs)
    r['n_shared_DEG_all20_units']=sum(g in degany.index for g in det)
    # DE-driven core: shared genes that are DEG in the same direction in >=2 key comparisons
    brows.append(r)
B=pd.DataFrame(brows)
def verdict(r):
    n=r.n_shared; d=r.n_shared_DEG_any_key
    if n==0: return 'No shared annotated genes: neither annotation overlap nor a DE-driven core.'
    comps=[s for _,s in KEY if r[f'shared_DEG_{s}']]
    if d==0: return f'{n} shared annotated gene(s), none a DEG in the six key comparisons: annotation overlap only, no DE-driven core.'
    genes=r.shared_DEG_genes_signed.replace(';',', ')
    return (f'{n} shared annotated gene(s), {d} of them DEG in {", ".join(comps)}: {genes} (+ = higher in T2). '
            f'The DE-supported part of the overlap is this small set ({d}/{n}); the remaining shared genes are annotation overlap only.')
B['interpretation']=B.apply(verdict,axis=1)
B.to_csv(TAB/'anchor_bridge_genes.csv',index=False)
pd.set_option('display.width',250,'display.max_colwidth',80)
print(B[['pathway_A','pathway_B','n_A','n_B','n_shared','n_shared_detected','n_shared_DEG_any_key','shared_DEG_genes']].to_string())
print(B[B.pathway_A=='GO AP-1 complex'][['pathway_B','shared_genes']].to_string())
print(mem[mem.pathway_id=='GO:0035976'].iloc[:,:14].to_string())

# ---- named candidate bridge genes
cands=['Fos','Fosb','Fosl1','Fosl2','Jun','Junb','Jund','Atf3','Nfatc2','Dusp1','Dusp2','Dusp6','Nr4a1','Nfkbia','Nfkb1','Rela','Tnfaip3','Traf2','Traf6','Chuk','Egr1','Myc','Ier2','Zfp36','Socs3','Il17ra','Il6','Cxcl2','Ccl2','Mapk8','Mapk14','Map2k7','Akt2','Grb2','Hsp90aa1','Hsp90ab1','Mmp9']
cr=[]
for g in cands:
    r=dict(gene=g,detected='yes' if g in measured else 'no')
    for k in DRILLID: r[DRILLNAME[k]]='Y' if g in full[k] else ''
    for u,s_ in KEY:
        r[s_+' log2FC(FDR)']= f'{wide.at[g,s_+"|lfc"]:+.2f} ({wide.at[g,s_+"|fdr"]:.1e})'+(' *' if wide.at[g,s_+'|deg'] else '') if g in wide.index else 'NA'
    r['DEG_all20_units']=degany.get(g,''); cr.append(r)
pd.DataFrame(cr).to_csv(TAB/'bridge_candidate_genes.csv',index=False)
print(pd.DataFrame(cr).iloc[:,:12].to_string())
