"""One-page a-d artwork for the paired spatial-context analysis; no inference."""
from pathlib import Path
import os,sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/context-ed16-mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
B=Path(__file__).resolve().parents[2];sys.path.insert(0,str(B/'Code/vendor'))
import pymupdf as fitz
D=B/'Source_data/spatial_context_consistency';OUT=B/'Assets/numbered_20261011';OUT.mkdir(exist_ok=True)
P=pd.read_csv(D/'paired_spots.csv');S=pd.read_csv(D/'summary.csv');C=pd.read_csv(D/'by_cell_type.csv')
order=S[S.metric=='spatial'].sort_values('median',ascending=False).method.tolist()
types=['Inner EPC 1','Inner EPC 2','Outer EPC 1','Outer EPC 2','P-TGC (Sec. P-TGC 1)','P-TGC (Sec. P-TGC 2)']
short=['Inner EPC1','Inner EPC2','Outer EPC1','Outer EPC2','Secondary P-TGC1','Secondary P-TGC2']
for f in ['Arial.ttf','Arial_Bold.ttf']:font_manager.fontManager.addfont('/usr/share/fonts/truetype/msttcorefonts/'+f)
plt.rcParams.update({'font.family':'Arial','font.size':6.3,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5})
W,H=595.276,841.89;fig=plt.figure(figsize=(W/72,H/72))
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def txt(x,y,s,**kw):return fig.text(x/W,1-y/H,s,**kw)
txt(24,20,'EXTENDED DATA FIGURE 16',fontsize=8,fontweight='bold')
txt(W/2,48,'Same 608 target spots: 1,292-spot full context versus 608-spot restricted context',ha='center',fontsize=7.4)
for k,(metric,x,color,title) in enumerate([('gene',99,'#3c789d','a  Shared-gene direction'),('spatial',280,'#cd783b','b  Spatial direction')]):
    a=ax(x,84,116,222);txt(25 if k==0 else 258,69,title,fontsize=8,fontweight='bold')
    for i,m in enumerate(order):
        v=P.loc[P.method==m,metric+'_cosine'].dropna()
        if len(v):
            a.boxplot([v],positions=[i],vert=False,widths=.54,showfliers=False,whis=(5,95),manage_ticks=False,patch_artist=True,boxprops={'facecolor':color,'alpha':.7,'linewidth':.4},medianprops={'color':'black','linewidth':.8},whiskerprops={'color':color,'linewidth':.5},capprops={'color':color,'linewidth':.5})
            r=S[(S.method==m)&(S.metric==metric)].iloc[0]
            a.text(1.11,i,str(int(r.shared_genes_used)) if k==0 else f'{r["median"]:.3f}',va='center',fontsize=5.7)
        else:a.text(0,i,'NA',ha='center',va='center',color='#666')
    a.set(xlim=(-1.04,1.04),ylim=(18.7,-.7),yticks=np.arange(19),yticklabels=order if k==0 else [],xticks=[-1,0,1])
    a.axvline(0,ls='--',color='#888',lw=.5);a.tick_params(axis='y',length=0,pad=4,labelsize=6.3);a.tick_params(axis='x',length=2,pad=2)
    a.set_xlabel('Cosine similarity',labelpad=2);a.spines[['top','left','right']].set_visible(False);a.grid(axis='y',color='#ededed',lw=.4);a.set_axisbelow(True)
    a.text(1.1,-1.1,'Genes' if k==0 else 'Median',fontsize=5.5)
txt(444,69,'c  Spatial median by type',fontsize=8,fontweight='bold')
a=ax(451,84,128,222);grid=C[C.metric=='spatial'].pivot(index='method',columns='cell_type',values='median').loc[order,types]
im=a.pcolormesh(np.arange(7)-.5,np.arange(20)-.5,grid,cmap='RdBu_r',vmin=-1,vmax=1,rasterized=False);a.set_ylim(18.5,-.5)
for i in range(19):
    for j in range(6):
        v=grid.iloc[i,j];a.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=5,color='white' if abs(v)>.65 else '#222')
a.set_xticks(range(6));a.set_xticklabels(short,rotation=58,ha='right',fontsize=5.4);a.set_yticks([]);a.tick_params(length=0,pad=2)
txt(25,349,'Boxes: median / IQR; whiskers: 5th–95th percentiles. All spatial rows: 608 valid pairs.',fontsize=6.4)
txt(25,359,'TopicVelo: transition-displacement proxy in a shared expression reference. VeloAE gene comparison: unavailable.',fontsize=6.1)
txt(25,378,'d  Spatial distribution of paired direction consistency',fontsize=8,fontweight='bold')
cells=pd.read_csv(B/'Source_data/spatial_completion/spatial_cells.csv');other=cells[~cells.cell_id.isin(P.cell_id)]
for i,m in enumerate(order):
    row,col=divmod(i,5);x=24+col*111;y=396+row*99
    a=ax(x,y+9,105,87);t=P[P.method==m]
    a.scatter(other.x,other.y,s=2.4,c='#d6d6d6',edgecolors='none')
    a.scatter(t.x,t.y,s=3,c=t.spatial_cosine,vmin=-1,vmax=1,cmap='RdBu_r',edgecolors='none')
    a.set_aspect('equal');a.set_xlim(cells.x.min()-60,cells.x.max()+60);a.set_ylim(cells.y.min()-60,cells.y.max()+60);a.axis('off')
    txt(x+52.5,y+3,m+f' ({t.spatial_cosine.median():.3f})',ha='center',fontsize=6.5)
txt(479,719,'Red: same direction\nBlue: opposite direction\nGrey: non-target spots\n\nMedian in parentheses',fontsize=6.5,va='top',linespacing=1.6)
cax=ax(174,805,240,7);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,-.5,0,.5,1]);cb.solids.set_rasterized(False);cb.ax.tick_params(labelsize=6,length=2,pad=2);cb.set_label('Cosine similarity (panels c and d)',fontsize=6.5,labelpad=2)
out=OUT/'Extended_Data_Fig16.pdf';fig.savefig(out);plt.close(fig)
d=fitz.open(out);assert len(d)==1 and not d[0].get_images();d[0].get_pixmap(matrix=fitz.Matrix(2,2)).save(out.with_suffix('.png'));out.with_suffix('.svg').write_text(d[0].get_svg_image(text_as_path=False))
print(out)
