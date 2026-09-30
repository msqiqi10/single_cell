"""Meeting-requested counts, percentages and both denominators on the same figure."""
from run_inkt_meeting_followup import *

def main():
    t=pd.read_csv(OUT/'tables/frequency_reconciliation_0.csv');d=pd.read_csv(OUT/'tables/frequency_reconciliation_1.csv')
    fig,axes=plt.subplots(2,3,figsize=(17,10));fig.subplots_adjust(top=.86,bottom=.20,hspace=.30,wspace=.25)
    for j,tissue in enumerate(TISSUES):
        sub=t[t.tissue.eq(tissue)];v=d[d.tissue.eq(tissue)].set_index('cluster').reindex(REFINED_ORDER);nc=int(sub.loc[sub.condition.eq('Ctrl'),'total_tissue_condition'].iloc[0]);nt=int(sub.loc[sub.condition.eq('T2'),'total_tissue_condition'].iloc[0]);x=np.arange(len(v));co=v.control_share_within_tissue_cluster_pct;tu=v.tumor_share_within_tissue_cluster_pct
        ax=axes[0,j];ax.bar(x,co,color='#2878B5',label='Ctrl');ax.bar(x,tu,bottom=co,color='#D9534F',label='T2');baseline=100*nt/(nc+nt);ax.axhline(100-baseline,color='black',ls='--',lw=.8)
        for i,r in enumerate(v.itertuples()):
            if r.control_count+r.tumor_count:
                if r.control_share_within_tissue_cluster_pct>=10:ax.text(i,r.control_share_within_tissue_cluster_pct/2,f'{r.control_count:.0f}\n{r.control_share_within_tissue_cluster_pct:.0f}%',ha='center',va='center',fontsize=7,color='white')
                if r.tumor_share_within_tissue_cluster_pct>=10:ax.text(i,100-r.tumor_share_within_tissue_cluster_pct/2,f'{r.tumor_count:.0f}\n{r.tumor_share_within_tissue_cluster_pct:.0f}%',ha='center',va='center',fontsize=7,color='white')
        ax.set_xticks(x,REFINED_ORDER,rotation=40);ax.set_ylim(0,105);ax.set_title(f'{tissue}\nCtrl total N={nc:,}; T2 total N={nt:,}',fontsize=12);ax.set_ylabel('Condition share within cluster (%)');ax.legend(frameon=False,fontsize=8,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.04))
        ax=axes[1,j];y=v.tumor_minus_control_percentage_points;ax.bar(x,y,color=np.where(y>=0,'#D9534F','#2878B5'));ax.axhline(0,c='#888888',lw=.7);ax.set_xticks(x,[f'{c}\n{int(r.control_count)}/{int(r.tumor_count)}' for c,r in zip(REFINED_ORDER,v.itertuples())],rotation=45,ha='right',fontsize=8);ax.set_ylabel('T2 - Ctrl frequency (percentage points)');ax.set_title('Cluster labels include exact Ctrl/T2 counts',fontsize=10)
    fig.suptitle('Frequency comparison: every denominator is explicit | 2026-09-05',fontsize=17,weight='bold')
    fig.text(.05,.105,'Top: T2 share = 100 × n(T2, cluster, tissue) / [n(T2, cluster, tissue) + n(Ctrl, cluster, tissue)].',fontsize=12)
    fig.text(.05,.07,'Bottom: frequency difference = 100 × [n(T2, cluster, tissue)/N(T2, tissue) - n(Ctrl, cluster, tissue)/N(Ctrl, tissue)].',fontsize=12)
    fig.text(.05,.035,'Dashed line: expected Ctrl share if cluster frequency is unchanged; it depends on sample totals and is not generally 50%. Empty groups are unestimable.',fontsize=11)
    save(fig,'08_frequency_denominators')
if __name__=='__main__':main()
