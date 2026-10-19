"""Restore original ED7/8 designs as current ED13/14, with all 19 configurations.
Reads portable, barcode-aligned CSVs extracted from saved results; no model training.
The Gaussian grid follows scVelo 0.3.4 compute_velocity_on_grid(adjust_for_stream=True).
"""
from pathlib import Path
import sys, os, json
os.environ.setdefault('MPLCONFIGDIR','/tmp/velocity-reference-mpl')
B=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(B/'Code/vendor'))
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import pymupdf as fitz
D=B/'Source_data/reference_restoration'
O=B/'Assets/arial'
R=B/'Reports/reference_restoration';R.mkdir(exist_ok=True)
for font in ['Arial.ttf','Arial_Bold.ttf']:font_manager.fontManager.addfont('/usr/share/fonts/truetype/msttcorefonts/'+font)
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.titlesize':6.5,'axes.labelsize':6.5,'axes.linewidth':.5,'xtick.labelsize':5.5,'ytick.labelsize':6,'pdf.fonttype':42,'svg.fonttype':'none'})
ORDER=['velocyto','scVelo(sto)','scVelo(dyn)','VeloAE','dynamo','VeloVAE','UniTVelo','DeepVelo(VAE)','cellDancer','veloVI','LatentVelo','DeepVelo(GCN)','STT','TFvelo','TopicVelo','biVI','DeepKINET','SDEvelo','cell2fate']
COL=dict(zip(ORDER,['#D62728','#CAB2D6','#6A3D9A','#33A02C','#FFB547','#B2DF8A','#BCBD22','#1F78B4','#E377C2','#FF9896','#17BECF','#A6CEE3','#FF7F00','#7B4173','#9ACD32','#1E90FF','#FF6347','#4682B4','#CD5C5C']))
NT={'HSC':'#d62728','GMP-like':'#279e68','MEP-like':'#aa40fc','Mon':'#e377c2','Neu':'#b5bd61','Bas':'#1f77b4','Ery':'#ff7f0e','Meg':'#8c564b'}
LT={'Baso/Eos':'#279e68','Ery/Meg':'#d62728','Lym/pDC':'#1f77b4','Mono':'#aa40fc','Neutro':'#aec7e8','Progenitors':'#ff7f0e','Unknown':'#808080'}
W,H=595.276,841.89
def ax(fig,x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(fig,x,y,s,**kw):return fig.text(x/W,1-y/H,s,**kw)
def letter(fig,x,y,s):text(fig,x,y,s,fontsize=8,fontweight='bold')
def field(a,xy,v,labels,palette,title,annotation=False,point_size=None):
    for c,color in palette.items():
        ok=labels==c;a.scatter(xy[ok,0],xy[ok,1],s=point_size if point_size is not None else 4.5,c=color,alpha=.40 if annotation else .32,edgecolors='none',rasterized=False)
    lo,hi=xy.min(0),xy.max(0);span=hi-lo
    if v is not None:
        keep=np.isfinite(v).all(1);X=xy[keep];V=v[keep]
        grids=[np.linspace(lo[j]-.01*span[j],hi[j]+.0101*span[j],50) for j in range(2)]
        xx,yy=np.meshgrid(*grids);dd,idx=cKDTree(X).query(np.c_[xx.ravel(),yy.ravel()],k=max(1,int(len(X)/50)))
        scale=np.mean([g[1]-g[0] for g in grids])*.5
        weight=np.exp(-.5*(dd/scale)**2)/(scale*np.sqrt(2*np.pi));mass=weight.sum(1)
        vg=(V[idx]*weight[:,:,None]).sum(1)/np.maximum(1,mass)[:,None]
        lengths=np.linalg.norm(vg,axis=1);length=np.abs(V[idx]).mean(1).sum(1)
        bad=(lengths<min(1e-5,lengths.max()*.9))|(length<np.percentile(length,5))
        vg[bad,0]=np.nan;lw=.6*lengths.reshape(50,50)/max(lengths.max(),1e-20)
        a.streamplot(*grids,vg[:,0].reshape(50,50),vg[:,1].reshape(50,50),density=1.35,color='#252525',linewidth=lw,arrowsize=.42,arrowstyle='-|>',zorder=3)
    a.set(xlim=(lo[0]-.04*span[0],hi[0]+.04*span[0]),ylim=(lo[1]-.04*span[1],hi[1]+.04*span[1]))
    a.set_xticks([]);a.set_yticks([])
    if annotation:
        a.set(xlabel='UMAP1',ylabel='UMAP2')
        for c in palette:
            ctr=np.median(xy[labels==c],axis=0);a.text(*ctr,c,fontsize=6.5,fontweight='bold',ha='center',va='center')
    else:a.axis('off')
    a.set_title(title,fontsize=6.5,pad=4)
def boxes(a,cos):
    order=cos.median().sort_values(ascending=False).index
    for i,m in enumerate(order):
        bp=a.boxplot([cos[m].dropna().to_numpy(float)],positions=[i],widths=.72,patch_artist=True,manage_ticks=False,
          boxprops={'linewidth':.5},medianprops={'color':'#333','linewidth':.65},whiskerprops={'linewidth':.5},capprops={'linewidth':.5},flierprops={'marker':'.','markersize':.7,'markeredgewidth':.1,'alpha':.35})
        bp['boxes'][0].set_facecolor(COL[m])
    a.set_xticks(range(19),order,rotation=90,fontsize=5.5);a.tick_params(axis='x',length=0,pad=3)
    a.set(xlim=(-.7,18.7),ylim=(-1.05,1.05),ylabel='Cosine similarity');a.set_yticks([-1,-.5,0,.5,1]);a.grid(axis='y',lw=.35,color='#dddddd');a.set_axisbelow(True)
    for s in a.spines.values():s.set_color('#777777')
    return list(order)
def render(n):
    ds='scNT' if n==13 else 'scLT';palette=NT if n==13 else LT
    cells=pd.read_csv(D/(ds+'_reference_cells.csv'));proj=pd.read_csv(D/(ds+'_projected_directions.csv'));cos=pd.read_csv(D/(ds+'_per_cell_cosine.csv')).set_index('cell_id')
    xy=cells[['x','y']].to_numpy();ref=cells[['vx','vy']].to_numpy();labels=cells.cell_type.to_numpy()
    fig=plt.figure(figsize=(W/72,H/72));text(fig,24,20,f'EXTENDED DATA FIGURE {n}',fontsize=8,fontweight='bold');letter(fig,25,49,'a')
    if n==13:
        field(ax(fig,53,61,190,158),xy,None,labels,palette,'',annotation=True,point_size=32)
        letter(fig,286,49,'c');boxorder=boxes(ax(fig,327,61,235,158),cos)
        letter(fig,25,286,'b');start,ph,dy=310,98,121
    else:
        full=pd.read_csv(D/'scLT_timepoints.csv')
        for j,samples in enumerate([['d2_5'],['d5r1','d5r2']]):
            a=ax(fig,45+j*135,62,110,126)
            for sam,color in zip(samples,['#ffbf69'] if j==0 else ['#9ecae1','#3182bd']):
                q=full[full['sample']==sam];a.scatter(q.x,q.y,s=3.0,c=color,edgecolors='none',alpha=.8,label=sam)
            a.set_xticks([]);a.set_yticks([]);a.set(xlabel='CCE1',ylabel='CCE2',title='Day 2.5' if j==0 else 'Day 5');a.legend(loc='lower left',frameon=False,fontsize=5.5,markerscale=3,handletextpad=.3)
        letter(fig,304,49,'d');boxorder=boxes(ax(fig,327,62,235,126),cos)
        letter(fig,25,255,'b');letter(fig,25,385,'c');start,ph,dy=409,75,99
    # ED13b / ED14c follow their displayed reference-cosine medians (not GDC).
    stream_order=cos[ORDER].median().sort_values(ascending=False,kind='stable').index.tolist()
    for i,m in enumerate(['Reference']+stream_order):
        a=ax(fig,39+i%5*107,start+i//5*dy,94,ph)
        if m=='Reference':v=ref;title='dynamo metabolic labeling' if n==13 else 'Clone-based\nreference'
        else:
            q=proj[proj.method==m].set_index('cell_id').reindex(cells.cell_id);v=q[['vx','vy']].to_numpy();title=m+'\n'+f'{cos[m].median():.3f}'
        # Wrap the long reference title to keep every title within its column.
        if n==13 and i==0:title='dynamo metabolic\nlabeling'
        field(a,xy,v,labels,palette,title,point_size=12 if n==13 else 4.5)
    legend=ax(fig,30,798,535,20);legend.axis('off');legend.legend(handles=[Line2D([],[],marker='o',ls='',color=color,label=c,ms=3) for c,color in palette.items()],loc='center',ncol=len(palette),frameon=False,fontsize=6,handletextpad=.3,columnspacing=1.1)
    raw=R/(f'Extended_Data_Fig{n}_without_clone_strip.pdf');fig.savefig(raw);plt.close(fig)
    out=O/f'Extended_Data_Fig{n}.pdf'
    with fitz.open(raw) as doc:
        if n==14:
            with fitz.open(D/'clone_examples_arial.pdf') as strip:doc[0].show_pdf_page(fitz.Rect(42,260,552,373),strip,0)
        doc.subset_fonts();doc.save(out,garbage=4,deflate=True)
        doc[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(out.with_suffix('.png'));out.with_suffix('.svg').write_text(doc[0].get_svg_image(text_as_path=False))
    with fitz.open(out) as check:
        sizes=[s['size'] for block in check[0].get_text('dict')['blocks'] for line in block.get('lines',[]) for s in line['spans'] if s['bbox'][1]>30 and s['size']<7.5]
    record={'file':out.name,'vector_labels_replaced':0,'raster_labels_replaced':0,'new_vector_text':True,'methods':19,'stream_panels':20,'boxplot_order':boxorder,'cell_count':len(cells),'historical_clone_strip':n==14,'renderer':'scVelo Gaussian-grid equations, Matplotlib streamplot; per-cell coordinates fixed','font':'Arial','body_pt':[min(sizes),max(sizes)],'panel_pt':8,
            'stream_order':['Reference']+list(stream_order),'stream_order_rule':'Reference first; then descending reference-direction cosine median, left-to-right then top-to-bottom',
            'median_labels':{m:f'{cos[m].median():.3f}' for m in ORDER},'median_label_style':'Method name then three-decimal median, following main Fig. 4a; metric is reference-direction cosine, not GDC',
            'scatter_area_pt2':{'annotation':32,'streams':12} if n==13 else {'timepoints':3.0,'streams':4.5},'boxplot_width_pt':235}
    (R/f'Extended_Data_Fig{n}.json').write_text(json.dumps(record,indent=2))
    (B/f'Reports/extended_typography/Extended_Data_Fig{n}.json').write_text(json.dumps(record,indent=2))
    print('Rendered',out,flush=True)
if __name__=='__main__':
    for n in ([int(x) for x in sys.argv[1:]] or [13,14]):render(n)
