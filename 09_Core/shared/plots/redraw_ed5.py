"""Recreate original ED4 a/b bubble matrices and c vertical boxes with current data.
Run with /data1/liyaru/miniconda/bin/python. No model execution or metric recalculation.
"""
from pathlib import Path
import os, sys, json, hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/velocity-ed5-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
B=Path(__file__).resolve().parents[2];ROOT=B;T=B/'Reports/redrawn_ed5';T.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(B/'Code/vendor'))
import openpyxl
DATA=B/'Source_data/cellDancer_simulations'
OUT=T/'rendered';OUT.mkdir(exist_ok=True)
SRC=B/'Code/inputs/extended_v2/Extended_Data_Fig05.xlsx'
w=openpyxl.load_workbook(SRC,read_only=True,data_only=True)
tables={}
for sheet in ['Gene-level cosine','Cell-level PCA2 cosine']:
    rows=list(w[sheet].values)
    tables[sheet]=pd.DataFrame([r[1:] for r in rows[1:]],index=[r[0] for r in rows[1:]],columns=rows[0][1:],dtype=float)
assert all(x.shape==(17,42) for x in tables.values())
# Read the published unit-level CSV, rather than use values digitized from old artwork.
families=['Mono-kinetic','Transcriptional_boost','Multi-forward_branching']
unit=[]
for chunk in pd.read_csv(DATA/'per_unit_common_genes.csv.gz',chunksize=100000):
    unit.append(chunk[chunk.dataset.isin([f+'_1' for f in families])])
unit=pd.concat(unit,ignore_index=True)
summary=pd.read_csv(DATA/'dataset_summary.csv')
inv=pd.read_csv(DATA/'dataset_inventory.csv').set_index('dataset')
labels={'velocyto':'velocyto','scvelo.sto':'scVelo(sto)','scvelo.dyn':'scVelo(dyn)','dynamo':'dynamo',
        'velovae':'VeloVAE','unitvelo':'UniTVelo','deepvelo_vae':'DeepVelo(VAE)',
        'deepvelo_gcn':'DeepVelo(GCN)','velovi':'veloVI','sdevelo':'SDEvelo',
        'celldancer':'cellDancer','latentvelo':'LatentVelo','stt':'STT','bivi':'biVI'}
# Original category colors; new methods extend the same muted palette.
colors={'velocyto':'#fb8072','scvelo.sto':'#cab2d6','scvelo.dyn':'#b3a2c7','dynamo':'#fdbf6f',
        'velovae':'#b2df8a','unitvelo':'#ffff99','deepvelo_vae':'#1f78b4','deepvelo_gcn':'#80b1d3',
        'velovi':'#fccde5','sdevelo':'#8dd3c7','celldancer':'#f781bf','latentvelo':'#bc80bd',
        'stt':'#ff7f00','bivi':'#bdbdbd'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':.55,
                     'pdf.fonttype':42,'svg.fonttype':'none','xtick.major.width':.5,
                     'ytick.major.width':.5,'savefig.facecolor':'white'})
W,H=760,1080
fig=plt.figure(figsize=(W/72,H/72))
def ax_at(x,y,w,h): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,**kw):return fig.text(x/W,1-y/H,s,**kw)
cmap=LinearSegmentedColormap.from_list('original_blue',['#08306b','#2171b5','#6baed6','#c6dbef','#f7fbff'])
audit={'source_workbook_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'panels':{},
       'style_reference':str(B/'Assets/original/Extended_Figures_submission.pdf'),
       'reference_pdf_page_1based':5,'models_rerun':False,
       'display_normalization':'Per-panel global min-minus-0.01 to max for dot area; column-wise average rank for blue fill. Mean uses raw cosine across datasets and the original joint Median/Mean min-minus-0.01/max for bar width; only the top three ranks are displayed. Raw values preserved.',
       'missing':'No dot/box or summary substitution; missing results remain missing.'}
for letter,sheet,top,title in [
    ('a','Gene-level cosine',82,'Cosine Similarity of High-dim RNA Velocity in Gene Level on Dyngen Simulation Datasets (n=42)'),
    ('b','Cell-level PCA2 cosine',365,'Cosine Similarity of Low-dim RNA Velocity in Cell Level on Dyngen Simulation Datasets (n=42)')]:
    df=tables[sheet];v=df.to_numpy();m=len(df);ax=ax_at(28,top,712,190)
    ax.set(xlim=(-6.3,48.0),ylim=(m-.5,-.5));ax.axis('off')
    for i,method in enumerate(df.index):
        if i%2==0:ax.add_patch(Rectangle((-6.3,i-.5),51.3,1,facecolor='#dedede',edgecolor='none'))
        ax.text(-6.15,i,method,va='center',fontsize=6.1)
    low=float(np.nanmin(v))-.01;high=float(np.nanmax(v));sizes=(v-low)/(high-low)
    ranks=df.rank(axis=0,method='average',ascending=False)
    # Better rank is lighter blue, matching the original bubble plot.
    rankcolors=1-(ranks-1).div(df.count(axis=0)-1,axis=1)
    max_area=67
    for i in range(m):
        for j in range(42):
            if np.isfinite(v[i,j]):ax.scatter(j,i,s=max_area*sizes[i,j],c=[cmap(rankcolors.iloc[i,j])],edgecolors='#555555',linewidths=.4)
    for j,ds in enumerate(df.columns):ax.text(j-.15,-.72,ds.replace('_',' '),rotation=30,ha='left',va='bottom',fontsize=4.25)
    med=df.median(axis=1);mean=df.mean(axis=1)
    blo=min(med.min(),mean.min())-.01;bhi=max(med.max(),mean.max())
    sums=pd.DataFrame({'median_raw':med,'mean_raw':mean,'valid_datasets':df.count(axis=1)})
    for x,series,title2 in [(43,mean,'Mean')]:
        rank=series.rank(ascending=False,method='min');avgrank=series.rank(ascending=False,method='average')
        norm=(series-blo)/(bhi-blo)
        ax.text(x+.2,-.8,title2,rotation=30,ha='left',va='bottom',fontsize=5.5)
        for i,value in enumerate(norm):
            if not np.isfinite(value):continue
            ax.add_patch(Rectangle((x-.7,i-.46),1.95*value,.92,facecolor=cmap(1-(avgrank.iloc[i]-1)/(m-1)),edgecolor='#555555',linewidth=.4))
            if rank.iloc[i]<=3:ax.text(x-.7+1.95*.4,i,str(int(rank.iloc[i])),ha='center',va='center',fontsize=4.8)
        sums[title2.lower()+'_rank']=rank
    text(16,top-55,letter,fontweight='bold',fontsize=12)
    text(395,top-55,title,ha='center',fontsize=8)
    la=ax_at(590,top+198,150,22);la.axis('off');la.set(xlim=(0,15),ylim=(0,2))
    la.text(0,1.75,'Ranking',fontsize=4.8)
    for k in range(5):la.add_patch(Rectangle((k*.43,.55),.43,.5,fc=cmap(1-k/4),ec='#777',lw=.25))
    la.text(0,.12,'1',fontsize=4);la.text(1.6,.12,str(m),fontsize=4)
    la.text(5,1.75,'Score (normalized)',fontsize=4.8)
    for j,z in enumerate([0,.5,1]):
        la.scatter(5.8+j*2.5,.78,s=max_area*z,fc='white',ec='#555',lw=.4)
        la.text(5.8+j*2.5,.05,f'{z:.0%}',ha='center',fontsize=4)
    df.to_csv(OUT/f'{letter}_raw_cosine.csv',float_format='%.17g')
    pd.DataFrame(sizes,index=df.index,columns=df.columns).to_csv(OUT/f'{letter}_display_area.csv',float_format='%.17g')
    sums.to_csv(OUT/f'{letter}_summaries.csv',float_format='%.17g')
    audit['panels'][letter]={'methods':m,'datasets':42,'finite':int(np.isfinite(v).sum()),'missing':int((~np.isfinite(v)).sum()),'area_min':low,'area_max':high,'summary_min':float(blo),'summary_max':float(bhi)}

text(16,612,'c',fontweight='bold',fontsize=12)
titles=['Mono-kinetic','Transcriptional\nboost','Multi-forward\nbranching']
stats=[]
for r,level in enumerate(['cell','gene']):
    for col,family in enumerate(families):
        dataset=family+'_1';d=unit[(unit.dataset==dataset)&(unit.level==level)]
        # Original six boxes panels use descending median order independently.
        available=[m for m in labels if np.isfinite(d.loc[d.method==m,'cosine']).any()]
        order=sorted(available,key=lambda m:-d.loc[d.method==m,'cosine'].median())
        x=68+col*225;top=637+r*216
        ax=ax_at(x,top,203,137)
        for i,m in enumerate(order):
            a=d.loc[d.method==m,'cosine'].to_numpy();a=a[np.isfinite(a)]
            assert len(a)>0
            expected=summary[(summary.dataset==dataset)&(summary.method==m)&(summary.level==level)&(summary.gene_scope=='common_genes')].iloc[0]
            assert len(a)==expected.valid_units and np.isclose(np.median(a),expected['median'],atol=1e-12)
            bp=ax.boxplot([a],positions=[i],widths=.76,patch_artist=True,showfliers=True,manage_ticks=False,whis=1.5,
                          medianprops={'color':'#333','linewidth':.7},boxprops={'linewidth':.55},
                          whiskerprops={'color':'#555','linewidth':.5},capprops={'color':'#555','linewidth':.5},
                          flierprops={'marker':'d','markersize':.65,'markeredgewidth':.1,'markerfacecolor':'#555','markeredgecolor':'#555','alpha':.65})
            bp['boxes'][0].set_facecolor(colors[m]);bp['boxes'][0].set_edgecolor('#555')
            stats.append({'dataset':dataset,'level':level,'method':m,'position':i,'valid_units':len(a),'median':float(np.median(a)),'q1':float(np.quantile(a,.25)),'q3':float(np.quantile(a,.75))})
        ax.set(xlim=(-.7,len(order)-.3),ylim=(-1.05,1.05))
        ax.set_xticks(range(len(order)),[labels[m] for m in order],rotation=90,fontsize=5.5)
        ax.set_yticks([-1,-.75,-.5,-.25,0,.25,.5,.75,1]);ax.tick_params(axis='y',labelsize=5,length=2);ax.tick_params(axis='x',length=2)
        for sp in ax.spines.values():sp.set_color('#777');sp.set_linewidth(.5)
        if col==0:ax.set_ylabel(('Cell-Level' if r==0 else 'Gene-Level')+'\nCosine Similarity',fontsize=8,fontweight='bold',labelpad=9)
        if r==0:ax.set_title(titles[col],fontsize=9,fontweight='bold',pad=10)
        audit['panels'][f'c_{level}_{family}']={'method_order':order,'common_genes':int(inv.loc[dataset,'common_genes']),'available_methods':int(d.method.nunique()),'unit_rows':len(d)}
fig.savefig(OUT/'Extended_Data_Fig05.svg')
fig.savefig(OUT/'Extended_Data_Fig05.pdf')
fig.savefig(OUT/'Extended_Data_Fig05.png',dpi=210)
plt.close(fig)
pd.DataFrame(stats).to_csv(OUT/'c_box_statistics.csv',index=False,float_format='%.17g')
unit.to_csv(OUT/'c_unit_scores.csv.gz',index=False,float_format='%.17g')
audit['input_hashes']={str(p.relative_to(B)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA/'per_unit_common_genes.csv.gz',DATA/'dataset_summary.csv',DATA/'dataset_inventory.csv']}
audit['output_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file()}
(T/'plot_audit.json').write_text(json.dumps(audit,indent=2))
print('Original a/b bubbles, summaries and c vertical boxes generated.',flush=True)
