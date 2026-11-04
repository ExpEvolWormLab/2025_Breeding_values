"""Independent preview only, with available summary inputs. Does not replace R pipeline."""
from pathlib import Path
import pandas as pd, numpy as np, matplotlib.pyplot as plt, zipfile, io
out=Path(__file__).parent/'preview';out.mkdir(exist_ok=True)
x=pd.read_csv('/mnt/data/pleiotropy_new_results/45_All_R5_clusters_GWAS_evolution_predictions.csv')
x=x[x.scenario=='primary_survivors'].copy();x['group']=x['group'].astype(str)
g=pd.read_csv('/mnt/data/11_Clusters_selection_and_GWAS_09f.csv');g['group']=g['group'].astype(str)
mer=x.merge(g[['group']+[f'beta_{tr}_NaCl' for tr in ['SF','SB','FS','FB','BS','BF','area.F']]],on='group',validate='one_to_one')
tags={'I':(5883574,7996414,[7167738,7491656,7596359]),'IV':(8650000,10600000,[8892663,10499135]),'V':(6046042,7278103,[6322080,6637695,6826421])}
tr=['SF','SB','FS','FB','BS','BF','area.F']; cols=[f'beta_{z}_NaCl' for z in tr]
sds=mer[cols].std()
traj=None
with zipfile.ZipFile('/mnt/data/chrV_deep_dive.zip') as z:
 p=[a for a in z.namelist() if a.endswith('56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv') and not a.startswith('__MACOSX/')][0]
 traj=pd.read_csv(io.BytesIO(z.read(p)))
for chr,(lo,hi,focus) in tags.items():
 d=mer[(mer.chr==chr)&mer.ps.between(lo,hi)].sort_values('ps').copy()
 fig=plt.figure(figsize=(13,7.7),layout='constrained')
 gs=fig.add_gridspec(2,2,height_ratios=(1.05,1.1))
 ax=fig.add_subplot(gs[0,:]); a2=fig.add_subplot(gs[1,0]); a3=fig.add_subplot(gs[1,1])
 ax.axhline(0,color='#888',linestyle='--',lw=1)
 for p in focus:ax.axvline(p/1e6,color='#ddd',linestyle=':',lw=1)
 good=d.dropna(subset=['primary_slope_median_phase'])
 ax.scatter(good.ps/1e6,good.primary_slope_median_phase,c=np.where(good.pred_full>0,'#B73545','#1866A3'),s=58,zorder=3)
 for _,row in good.iterrows():ax.annotate(str(row.group),(row.ps/1e6,row.primary_slope_median_phase),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
 ax.set_xlim(lo/1e6,hi/1e6);ax.set_xlabel('Position (Mb)');ax.set_ylabel('Phase-aligned evolutionary slope / generation')
 ax.set_title(f'A  |  Chromosome {chr} — current R5 clusters (n={len(d)})',loc='left',fontweight='bold')
 mm=(d[cols]/sds).T.values
 im=a2.imshow(np.clip(mm,-3,3),aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3)
 a2.set_yticks(range(7),tr);a2.set_xticks(range(len(d)),d.group,rotation=90);a2.set_xlabel('R5 group')
 a2.set_title('B  |  Trait GWAS effect (scaled within trait)',loc='left',fontweight='bold')
 fig.colorbar(im,ax=a2,shrink=.7,label='effect / genome-wide SD')
 if chr!='V':
   a3.axis('off');a3.text(.5,.57,'Exact-site Pool-seq table not available\nfor this chromosome in the uploaded archives',ha='center',va='center',transform=a3.transAxes,fontsize=12)
   a3.text(.5,.36,'Script 63 will incorporate the historical\nscript-28/30c output when supplied',ha='center',va='center',transform=a3.transAxes,fontsize=10,color='#666')
 else:
   t=traj[traj.SNP_ID.isin([f'V_{p}' for p in focus]) & traj.population.isin(['A00','GA1','GA2','GA4'])].copy()
   for id in t.SNP_ID.unique():
     sub=t[t.SNP_ID==id]; anc=sub[sub.sample_id=='A00'].ALT_frequency.iloc[0]; flip=(anc>.5)
     for pop,color in [('GA1','#2676AA'),('GA2','#D28B38'),('GA4','#298269')]:
        z=sub[(sub.population==pop)&(sub.generation>0)].sort_values('generation')
        if len(z)==0:continue
        y=np.r_[1-anc if flip else anc,1-z.ALT_frequency if flip else z.ALT_frequency]
        xx=np.r_[0,z.generation.values]
        a3.plot(xx,y,marker='o',markersize=3,alpha=.75,lw=1.5,color=color,linestyle=['-','--',':'][list(t.SNP_ID.unique()).index(id)],label=f'{pop}: {id}')
   a3.set_ylim(0,1);a3.set_xlabel('Generation');a3.set_ylabel('A00-minor allele frequency')
   a3.legend(fontsize=7,loc='upper left',ncol=2,frameon=False)
 a3.set_title('C  |  Focal exact-SNP trajectories',loc='left',fontweight='bold')
 fig.suptitle('Descriptive regional analysis — linked R5 clusters are not independent causal loci',y=1.03,fontsize=12)
 fig.savefig(out/f'Chr{chr}_regional_PREVIEW.png',dpi=180,bbox_inches='tight')
 fig.savefig(out/f'Chr{chr}_regional_PREVIEW.pdf',bbox_inches='tight')
 plt.close(fig)
print('PREVIEW OUTPUTS',*list(out.glob('*')),sep='\n')
