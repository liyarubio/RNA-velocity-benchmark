"""Restore original Extended Data encodings using saved, traceable local data.
Run with /data1/liyaru/miniconda/bin/python. No model training or velocity inference.
"""
from pathlib import Path
import os,sys,json,hashlib,re,itertools
os.environ.setdefault('MPLCONFIGDIR','/tmp/velocity-extended-original-mpl')
B=Path(__file__).resolve().parents[2];R=B
sys.path.insert(0,str(B/'Code'))
from scenario_scope import FAMILIES, CAPTION_EN
from ed6_revision import revise_text
T=B/'Reports/redrawn_extended';T.mkdir(parents=True,exist_ok=True)
D=T/'data';O=T/'rendered';D.mkdir(exist_ok=True);O.mkdir(exist_ok=True)
sys.path.insert(0,str(B/'Code/vendor'))
import numpy as np,pandas as pd,h5py,openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,PathPatch,Arc
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial import cKDTree
import pymupdf as fitz
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.linewidth':.5,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
ORDER=['velocyto','scVelo(sto)','scVelo(dyn)','VeloAE','dynamo','VeloVAE','UniTVelo','DeepVelo(VAE)','cellDancer','veloVI','LatentVelo','DeepVelo(GCN)','TFvelo','TopicVelo','STT','biVI','DeepKINET','SDEvelo','cell2fate']
COL=dict(zip(ORDER,['#D62728','#CAB2D6','#6A3D9A','#33A02C','#FFB547','#B2DF8A','#E6CF5C','#0000CD','#FFB6C1','#DB7093','#20B2AA','#87CEFA','#7B4173','#9ACD32','#FF8C00','#1E90FF','#FF6347','#4682B4','#CD5C5C']))
TOPO=['bifurcating_converging','bifurcating_cycle','bifurcating_loop','bifurcating','binary_tree','branching','consecutive_bifurcating','converging','cycle','cycle_simple','disconnected','linear','linear_simple','trifurcating']
TC=dict(zip(TOPO+['total'],['#F8766D','#E78600','#C49A00','#A3A500','#39B600','#00BA38','#00C08B','#00BFC4','#00B0F6','#9590FF','#E76BF3','#FF62BC','#F564E3','#619CFF','#000000']))
BLUE=LinearSegmentedColormap.from_list('rankblue',['#08306b','#2171b5','#6baed6','#c6dbef','#f7fbff'])
GREEN=plt.get_cmap('Greens_r')
AUDIT={'inputs':{},'figures':{},'models_rerun':False};CAPTIONS={}
def track(p):
 p=Path(p);AUDIT['inputs'][str(p.resolve())]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size};return p
def book(n):
 p=track(B/f'Code/inputs/extended_v2/Extended_Data_Fig{n:02d}.xlsx');w=openpyxl.load_workbook(p,read_only=True,data_only=True)
 out={}
 for s in w:
  it=s.iter_rows(values_only=True);h=next(it,None)
  if h:out[s.title]=pd.DataFrame(it,columns=h)
 return out
def csv(d,name):d.to_csv(D/(name+'.csv'),index=False,float_format='%.17g')
def fig(w,h):return plt.figure(figsize=(w/72,h/72))
def ax(f,x,y,w,h):W,H=f.get_size_inches()*72;return f.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(f,x,y,s,**kw):W,H=f.get_size_inches()*72;return f.text(x/W,1-y/H,s,**kw)
def panel(f,x,y,s):text(f,x,y,s,fontsize=13,fontweight='bold')
def grid(a):
 a.grid(color='#e5e5e5',lw=.4);a.set_axisbelow(True)
 for sp in a.spines.values():sp.set_color('#777')
def save(f,n):
 stem=O/f'Extended_Data_Fig{n:02d}'
 for e in ['pdf','svg']:f.savefig(str(stem)+'.'+e)
 f.savefig(str(stem)+'.png',dpi=180);plt.close(f)
def save_panel(f,name):
 p=O/(name+'.pdf');f.savefig(p);plt.close(f);return p
def category(f,k):
 n=f['obs'][k];return n['categories'].asstr()[:][n['codes'][:]] if isinstance(n,h5py.Group) else n[:].astype(str)
def ids(f):return f['obs'][f['obs'].attrs.get('_index','_index')].asstr()[:]
def labelmethod(s):return {'scVelo (stochastic)':'scVelo(sto)','scVelo (dynamic)':'scVelo(dyn)'}.get(s,s)
def bubble(a,df,cmap,summary=True,labelsize=6,area=65,label_gutter=2.9,label_rotation=48):
 v=df.to_numpy(float);m,n=v.shape;low=np.nanmin(v)-.01;high=np.nanmax(v);norm=(v-low)/(high-low)
 rank=df.rank(axis=0,ascending=False,method='average');cols=1-(rank-1).div(df.count(axis=0)-1,axis=1)
 a.set(xlim=(-label_gutter,n+(4 if summary else 0)),ylim=(m-.5,-.5));a.axis('off')
 for i,name in enumerate(df.index):
  if i%2==0:a.add_patch(Rectangle((-label_gutter,i-.5),n+7+(label_gutter-2.9),.99,fc='#dedede',ec='none'))
  a.text(-label_gutter+.05,i,name,va='center',fontsize=labelsize)
  for j in range(n):
   if np.isfinite(v[i,j]):a.scatter(j,i,s=area*norm[i,j],c=[cmap(cols.iloc[i,j])],ec='#555',lw=.35)
 for j,name in enumerate(df.columns):a.text(j,-.8,name,rotation=label_rotation,ha='left',va='bottom',fontsize=labelsize-1)
 out=[]
 if summary:
  med=df.median(axis=1);mean=df.mean(axis=1);lo=min(med.min(),mean.min())-.01;hi=max(med.max(),mean.max())
  for x,s,label in [(n+.45,med,'Median'),(n+2.4,mean,'Mean')]:
   ranks=s.rank(ascending=False,method='min');cr=s.rank(ascending=False,method='average');ss=(s-lo)/(hi-lo)
   a.text(x,-.8,label,rotation=label_rotation,fontsize=labelsize,va='bottom')
   for i,value in enumerate(ss):
    if not np.isfinite(value):continue
    a.add_patch(Rectangle((x-.3,i-.45),1.5*value,.9,fc=cmap(1-(cr.iloc[i]-1)/(m-1)),ec='#555',lw=.35))
    if ranks.iloc[i]<=6:a.text(x-.3+.75*value,i,str(int(ranks.iloc[i])),ha='center',va='center',fontsize=labelsize-1)
  out=pd.DataFrame({'method':df.index,'median':med.values,'mean':mean.values,'valid_datasets':df.count(axis=1).values})
 return out,{'display_min':float(low),'display_max':float(high),'finite_values':int(np.isfinite(v).sum()),'missing_values':int((~np.isfinite(v)).sum())}

def simulation_consistency():
 w=book(8);pair=w['Paired cell-level values'].copy();dims=w['PCA dimension values'].copy()
 pair['topology']=pair.dataset.str.replace(r'_seed\d+$','',regex=True)
 methods=[m for m in ORDER if m in set(pair.method)];f=fig(840,1190);panel(f,18,26,'a')
 summaries=[]
 for i,top in enumerate(TOPO):
  a=ax(f,48+(i%5)*157,60+(i//5)*128,130,102)
  for m in methods:
   z=pair[(pair.topology==top)&(pair.method==m)]
   if not len(z):continue
   vals=z[['high','low']].astype(float);means=vals.mean();sd=vals.std()
   a.errorbar(means.high,means.low,xerr=sd.high if np.isfinite(sd.high) else None,yerr=sd.low if np.isfinite(sd.low) else None,fmt='o',ms=2.4,lw=.5,color=COL[m])
   summaries.append({'topology':top,'method':m,'paired_datasets':len(z),'high_mean':means.high,'low_mean':means.low,'high_sd':sd.high,'low_sd':sd.low})
  a.set_title(top.replace('_',' '),fontsize=6.8,pad=4);grid(a);a.tick_params(labelsize=5)
  if i//5==2:a.set_xlabel('High-dimensional cell cosine',fontsize=5.5)
  if i%5==0:a.set_ylabel('PCA2 cell cosine',fontsize=5.5)
 la=ax(f,674,316,148,112);la.axis('off');la.legend(handles=[Line2D([],[],color=COL[m],marker='o',ls='',ms=3,label=m) for m in methods],loc='center',ncol=2,frameon=False,fontsize=5.4,columnspacing=.5,handletextpad=.3)
 panel(f,18,457,'b')
 for i,m in enumerate(methods):
  a=ax(f,48+(i%6)*128,485+(i//6)*93,111,69);z=pair[pair.method==m].sort_values('dataset');v=z[['high','low']].to_numpy(float)
  for row in v:a.plot([0,1],row,c='#bbbbbb',lw=.35,zorder=1)
  for j,c in [(0,'#08306b'),(1,'#c6dbef')]:
   data=v[:,j]
   if len(np.unique(data))>1:
    violin=a.violinplot(data,positions=[j],showextrema=False,widths=.62)
    for body in violin['bodies']:body.set(fc=c,ec='#333',lw=.6,alpha=1)
   a.scatter(np.full(len(data),j),data,s=1.3,c='#333',zorder=3)
  a.set(xlim=(-.5,1.5),ylim=(-1.05,1.05));a.set_xticks([0,1],['High','Low'],fontsize=5);a.set_yticks([-1,0,1]);a.tick_params(labelsize=5)
  a.set_title(m,fontsize=6.5,pad=3,bbox=dict(fc='white',ec='#555',lw=.5,pad=1.8));a.spines[['top','right']].set_visible(False)
  if i%6==0:a.set_ylabel('Cell cosine',fontsize=6)
 # Original R 'cosine median' block: median projected-velocity accuracy, not a second cosine.
 original_colors=['#F8766D','#E38900','#C49A00','#99A800','#53B400','#00BC56','#00C094','#00BFC4','#00B6EB','#06A4FF','#A58AFF','#DF70F8','#FB61D7','#FF66A8','#000000']
 accuracy_colors=dict(zip(TOPO+['total'],original_colors))
 dimrows=[]
 for d,g in dims.groupby('dimension'):
  g=g[np.isfinite(g.cosine.astype(float))].copy()
  g['topology']=g.Dataset.str.replace(r'_seed\d+$','',regex=True)
  for top in TOPO+['total']:
   z=g if top=='total' else g[g.topology==top]
   dimrows.append({'topology':top,'dimension':int(d),'median_velocity_cosine':float(z.cosine.median()),'valid_method_dataset_scores':len(z),'datasets':z.Dataset.nunique(),'methods':z.method.nunique()})
 curves=pd.DataFrame(dimrows);panel(f,18,802,'c');a=ax(f,60,825,490,335)
 for top in TOPO+['total']:
  z=curves[curves.topology==top].sort_values('dimension');a.plot(z.dimension,z.median_velocity_cosine,'o-',color=accuracy_colors[top],ms=2.5,lw=.85,label=top)
 a.set(xlabel='Dimension',ylabel='',title='Cosine similarity',xlim=(1.4,30.6));a.set_xticks(range(2,31));a.tick_params(labelsize=6);a.grid(False)
 for sp in a.spines.values():sp.set_color('#777')
 a.legend(title='Dataset',loc='center left',bbox_to_anchor=(1.02,.5),frameon=False,fontsize=7,title_fontsize=8)
 csv(pd.DataFrame(summaries),'Fig08_topology_scatter');csv(curves,'ExtFig06c_dimension_accuracy');csv(pair,'Fig08_paired_scores')
 save(f,8);AUDIT['figures']['8']={'current_figure':6,'style':'Original a/b panels; original topology-colored accuracy curves with points and black total','methods':len(methods),'paired_rows':len(pair),'source_metric_values_unchanged':True,'derived_curve':'median of finite method-dataset projected-velocity cosine scores within topology at each PCA dimension; total takes the median across all 42 datasets','reference':'Code/archive/figure3/Figure3df_corplot_original.R:537-579','colors':accuracy_colors}
 CAPTIONS['8']='a, Original-style scatter plots separated by the 14 dyngen topologies. Each point summarizes one method across available matched seed datasets (mean and sample SD); high-dimensional and PCA2 scores are cell-level cosine similarities. b, Paired high/low-dimensional violin plots, faceted by method; lines join the same dataset, with dataset medians shown as points. c, Original-style topology-colored dimension-sensitivity curves: cosine similarity between the vectors of high-dimensional and dimension-specific performance scores, paired by method and dataset. This is agreement between performance profiles, not mean RNA-velocity accuracy. The total curve pools all finite pairs; pair counts and derived values are supplied. Seventeen available configurations are shown; missing pairs are excluded without imputation. TopicVelo table values retain unverified native-transition provenance. No model was rerun.'
 CAPTIONS['8']=revise_text(CAPTIONS['8'])

def real_metrics():
 w=book(9);hub=w['HuBMAP raw GDC'].set_index('file_name');methods=[m for m in ORDER if m in hub];hub=hub[methods].astype(float);best=hub.max(axis=1)
 x=hub.to_numpy();filled=np.where(np.isfinite(x),x,np.nanmedian(x,axis=0));Z=linkage(filled,method='average',metric='euclidean');leaves=dendrogram(Z,no_plot=True)['leaves']
 tablefile=track(B/'Tables/Table S3 Real-world Datasets.xlsx');wb=openpyxl.load_workbook(tablefile,read_only=True,data_only=True);rows=list(wb['HuBMAP'].values);meta=pd.DataFrame(rows[1:],columns=rows[0]);lookup=meta.set_index('uuid').organ
 organs=pd.Series([lookup.get(str(s).removesuffix('.h5ad'),'Unmapped') for s in hub.index],index=hub.index)
 assert not (organs=='Unmapped').any()
 cats=list(dict.fromkeys(organs));palette={c:plt.get_cmap('tab20')(i) for i,c in enumerate(cats)}
 f=fig(920,740)
 for x,label in [(16,'a'),(483,'b')]:text(f,x,32,label,fontsize=16,fontweight='bold')
 da=ax(f,20,70,38,435);dendrogram(Z,orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#999',ax=da);da.invert_yaxis();da.axis('off')
 ha=ax(f,62,70,316,435);cm=plt.get_cmap('viridis').copy();cm.set_bad('#dedede');ha.pcolormesh(np.arange(len(methods)+1),np.arange(len(hub)+1),np.ma.masked_invalid(hub.iloc[leaves].to_numpy()),cmap=cm,vmin=-1,vmax=1,rasterized=False)
 ha.set_ylim(len(hub),0);ha.set_xticks(np.arange(len(methods))+.5,methods,rotation=90,fontsize=10);ha.set_yticks([])
 ba=ax(f,384,70,16,435);ba.pcolormesh([0,1],np.arange(len(hub)+1),best.iloc[leaves].to_numpy()[:,None],cmap=cm,vmin=-1,vmax=1,rasterized=False);ba.set_ylim(len(hub),0);ba.set_xticks([.5],['Best performance'],rotation=90,fontsize=10);ba.set_yticks([])
 oa=ax(f,406,70,7,435)
 for i,c in enumerate(organs.iloc[leaves]):oa.add_patch(Rectangle((0,i),1,1,fc=palette[c],ec='none'))
 oa.set(xlim=(0,1),ylim=(len(hub),0));oa.axis('off')
 cb=ax(f,427,70,9,130);f.colorbar(plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(-1,1),cmap=cm),cax=cb,label='GDC')
 cb.tick_params(labelsize=9.5);cb.set_ylabel('GDC',fontsize=10)
 for collection in cb.collections:collection.set_rasterized(False)
 text(f,72,53,'HuBMAP datasets (n=176)',fontsize=13)
 leg=ax(f,35,630,420,90);leg.axis('off');leg.legend(handles=[Line2D([],[],marker='s',ls='',color=palette[c],label=c,ms=6) for c in cats],loc='upper left',ncol=3,fontsize=9.5,frameon=False,columnspacing=1.2,handletextpad=.5)
 ten=w['Ten datasets raw GDC'].set_index('index').astype(float);ten=ten.reindex([m for m in ORDER if m in ten.index])
 a=ax(f,495,200,410,305);summ,audit=bubble(a,ten,GREEN,True,9.5,112,label_gutter=4.5,label_rotation=65)
 text(f,687,53,'Ten classical datasets',fontsize=13,ha='center')
 csv(pd.DataFrame({'dataset':hub.index[leaves],'organ':organs.iloc[leaves].values,'best_performance':best.iloc[leaves].values}),'Fig09_clustered_datasets');csv(summ,'Fig09_ten_summary')
 save(f,9);AUDIT['figures']['9']={'current_figure':7,'style':'clustered HuBMAP heatmap with organ strip and separate Best performance; green Dataset10 bubbles and Median/Mean bars','font_revision':'Larger labels with wider method gutter and more header/legend space; numeric data and ordering unchanged','font_sizes_before_page_scaling':{'heatmap_methods':10,'bubble_methods':9.5,'bubble_datasets_and_ranks':8.5,'organ_legend':9.5,'titles':13},'HuBMAP_methods':len(methods),'classical_methods':len(ten),'organ_annotations_added':176,'clustering_missing_rule':'column median used for row-order distance only; display and summaries retain original missing values','bubble':audit}
 CAPTIONS['9']='a, Original-style row-clustered heatmap of saved GDC values for 176 HuBMAP datasets and 18 available configurations, with organ annotations recovered by exact UUID matching to Supplementary Table 3. Best performance is the maximum finite value per dataset and is displayed separately; it is not an additional method. The shared raw GDC color scale is -1 to 1. Grey denotes missing data. Average-linkage clustering uses Euclidean distances with column-median filling only to define display order; plotted values and numerical summaries remain unfilled. b, Original-style green bubble matrix for ten classical datasets and 19 configurations, with Median/Mean bars and top-six rank labels. Dot area is globally normalized from the panel minimum minus 0.01 to its maximum; lighter color denotes better within-dataset rank. Raw values, display order, coverage and organ mapping are provided in source data. Missing entries are not scored as zero.'

def streams(a,xy,v,labels,palette,title):
 for c,color in palette.items():
  keep=labels==c;a.scatter(xy[keep,0],xy[keep,1],s=2.4,c=[color],alpha=.7,ec='none',rasterized=False)
 lo=xy.min(0);hi=xy.max(0);span=hi-lo
 if v is not None:
  valid=np.isfinite(v).all(1)&(np.linalg.norm(np.nan_to_num(v),axis=1)>0)
  if valid.sum()>5:
   gx=np.linspace(lo[0]-.04*span[0],hi[0]+.04*span[0],28);gy=np.linspace(lo[1]-.04*span[1],hi[1]+.04*span[1],28);xx,yy=np.meshgrid(gx,gy)
   dd,idx=cKDTree(xy[valid]).query(np.c_[xx.ravel(),yy.ravel()],k=min(30,valid.sum()));scale=span.mean()/28*.75;wt=np.exp(-.5*(dd/scale)**2);mass=wt.sum(1)
   field=(v[valid][idx]*wt[:,:,None]).sum(1)/np.maximum(mass[:,None],1e-10);field[mass<.2]=np.nan
   a.streamplot(gx,gy,np.ma.masked_invalid(field[:,0].reshape(xx.shape)),np.ma.masked_invalid(field[:,1].reshape(xx.shape)),density=.9,color='#303030',linewidth=.45,arrowsize=.55)
 a.set(xlim=(lo[0]-.06*span[0],hi[0]+.06*span[0]),ylim=(lo[1]-.06*span[1],hi[1]+.06*span[1]));a.set_aspect('equal');a.axis('off');a.set_title(title,fontsize=8,pad=3)
def boxes(a,cos):
 methods=cos.median().sort_values(ascending=False).index
 for i,m in enumerate(methods):
  vals=cos[m].dropna().to_numpy(float);bp=a.boxplot([vals],positions=[i],patch_artist=True,widths=.7,showfliers=True,manage_ticks=False,
     boxprops={'linewidth':.5},medianprops={'color':'#333','linewidth':.65},whiskerprops={'linewidth':.5},capprops={'linewidth':.5},flierprops={'marker':'.','markersize':.7,'markeredgewidth':.1,'alpha':.5})
  bp['boxes'][0].set_facecolor(COL.get(labelmethod(m),'#aaa'))
 a.set_xticks(range(len(methods)),[labelmethod(m) for m in methods],rotation=90,fontsize=6);a.set(ylabel='Direction cosine',ylim=(-1.05,1.05));a.tick_params(labelsize=6);a.spines[['top','right']].set_visible(False)
def refs(n,ds):
 w=book(n);cells=w['Reference cells'];cos=w['Per-cell cosine'].set_index('cell_id');proj=w['Projected directions'];methods=list(w['Method summary'].method)
 xy=cells[['x','y']].to_numpy(float);v=cells[['vx','vy']].to_numpy(float);labels=cells.cell_type.to_numpy()
 cats=list(dict.fromkeys(labels));palette={c:plt.get_cmap('tab10')(i%10) for i,c in enumerate(cats)}
 if ds=='scNT':palette=dict(zip(['HSC','GMP-like','MEP-like','Bas','Ery','Meg','Mon','Neu'],['#e41a1c','#4daf4a','#984ea3','#ffff33','#ff7f00','#f781bf','#377eb8','#999999']))
 f=fig(760,850 if ds=='scNT' else 930)
 if ds=='scNT':
  panel(f,20,24,'a');a=ax(f,40,48,265,210);streams(a,xy,None,labels,palette,'scNT-seq annotation')
  for c in cats:
   ctr=np.median(xy[labels==c],axis=0);a.text(*ctr,c,fontsize=6,fontweight='bold',ha='center')
  panel(f,335,24,'c');boxes(ax(f,377,52,342,150),cos)
  panel(f,20,304,'b');start=328;ph=180;dy=218
 else:
  p=track(B/'Source_data/derived/Fig16_timepoint_coordinates.csv')
  supplied=pd.read_csv(p);full=supplied[['x','y']].to_numpy(float);sample=supplied['sample'].to_numpy();cellid=supplied.cell_id.to_numpy()
  sample_df=pd.DataFrame({'cell_id':cellid,'x':full[:,0],'y':full[:,1],'sample':sample});csv(sample_df,'Fig16_timepoint_coordinates')
  panel(f,20,24,'a')
  for j,ss in enumerate([['d2_5'],['d5r1','d5r2']]):
   a=ax(f,44+j*173,43,155,150)
   for sam,color in zip(ss,['#ffbf69'] if j==0 else ['#9ecae1','#3182bd']):
    ok=sample==sam;a.scatter(full[ok,0],full[ok,1],s=1.1,color=color,alpha=.8,ec='none',label=sam)
   a.set(xlabel='CCE1',ylabel='CCE2' if j==0 else None,title='Day 2.5' if j==0 else 'Day 5');a.set_xticks([]);a.set_yticks([]);a.legend(frameon=False,fontsize=6,markerscale=3,loc='lower left')
  panel(f,420,24,'d');boxes(ax(f,467,45,264,145),cos)
  panel(f,20,263,'b')
  # The actual five clone identities are not in the delivered per-cell tables.
  # Preserve their original illustration and explicitly report this provenance limit.
  start=462;ph=160;dy=198;panel(f,20,442,'c')
 for i,m in enumerate(['Reference']+methods):
  a=ax(f,38+i%4*180,start+i//4*dy,162,ph)
  if m=='Reference':vel=v;name='Dynamo (metabolic labeling)' if ds=='scNT' else 'Clone-based reference'
  else:
   q=proj[proj.method==m].set_index('cell_id').reindex(cells.cell_id);assert np.allclose(q[['x','y']],xy);vel=q[['vx','vy']].to_numpy(float);name=labelmethod(m)
  streams(a,xy,vel,labels,palette,name)
 la=ax(f,40,start+2*dy-4,680,68);la.axis('off');la.legend(handles=[Line2D([],[],marker='o',ls='',color=palette[c],label=c,ms=4) for c in cats],loc='center',ncol=4,frameon=False,fontsize=7)
 if ds=='scNT':save(f,n)
 else:
  raw=save_panel(f,'Fig16_without_clone_strip');original=fitz.open(track(B/'Assets/original/Extended_Figures_submission.pdf'));d=fitz.open(raw);d[0].show_pdf_page(fitz.Rect(42,281,734,410),original,10,clip=fitz.Rect(42,189,548,282))
  out=O/'Extended_Data_Fig16.pdf';d.save(out,garbage=4,deflate=True);d[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(out.with_suffix('.png'));out.with_suffix('.svg').write_text(d[0].get_svg_image());d.close()
 AUDIT['figures'][str(n)]={'style':'original annotation / method fields / vertical boxes' if ds=='scNT' else 'original timepoint / clone examples / method fields / boxes','available_methods':7,'reference_cells':len(cells),'restored_annotation_cells':len(cells) if ds=='scNT' else 3186,'score_values_unchanged':True,'historical_clone_strip':ds=='scLT'}
 CAPTIONS[str(n)]=( ('a, Original-style UMAP annotation of all 1,947 scNT reference cells. b, Metabolic-labeling reference direction and seven saved method projections in the same UMAP coordinates. c, Original-style vertical boxplots of per-cell direction cosine relative to the metabolic-labeling reference.' if ds=='scNT' else
  'a, Original-style CCE views of Day 2.5 (1,756 cells) and Day 5 (675 and 755 cells in the two replicates), restored from the local full 3,186-cell input. b, Original five clone-transition examples retained as historical artwork: the exact clone-membership/pairing table is not present in the supplied local files, so these arrows have not been reconstructed. c, Clone-based reference and seven saved method directions for 1,756 reference cells. d, Original-style vertical boxplots of per-cell cosine relative to the clone-based reference.')+
  ' Boxes show median/IQR, 1.5-IQR whiskers and outliers; methods are sorted by median. All scores and projected directions are unchanged; only stream rendering and layout are regenerated. scVelo and SDEvelo use the preprocessed-input batch, while TopicVelo, biVI, DeepKINET and cell2fate use the aligned raw-count extension. This is not a uniform-preprocessing ranking. TopicVelo uses native transition projections; its 32 zero-direction scLT cells remain undefined, not zero scores. Additional historical score columns lack matching current projection/run evidence and are not merged into this comparison. Current coverage and missing evidence are listed in the data report.')

def chord():
 w=book(19);member=w['Top100 membership'];methods=['DeepVelo(GCN)','LatentVelo','STT','UniTVelo','dynamo','scVelo(dyn)','scVelo(sto)','veloVI'];assert set(methods)==set(member.columns[1:]);overlap=member[methods].to_numpy(int).T@member[methods].to_numpy(int)
 f=fig(345,345);a=ax(f,24,24,297,297);a.set(xlim=(-1.25,1.25),ylim=(-1.25,1.25));a.axis('off');a.set_aspect('equal')
 angles={};segments=[];links=[];arc_colors=['#ff0000','#4dbbd5','#00a65a','#9b59b6','#f28e2b','#f9e51d','#a65628','#ef77be']
 for i,m in enumerate(methods):
  genes=member.loc[member[m].astype(int)>0,'Gene'].astype(str).tolist();lo=90-i*45-40;hi=90-i*45
  # Gene slots preserve table order; every chord joins the same literal gene ID.
  for j,g in enumerate(genes):angles[m,g]=np.deg2rad(hi-(j+.5)/len(genes)*40);segments.append({'method':m,'gene':g,'position':j,'angle':float(angles[m,g])})
  a.add_patch(Arc((0,0),2.08,2.08,theta1=lo,theta2=hi,color=arc_colors[i],lw=4))
  for j,g in enumerate(genes):
   shared=int(member.loc[member.Gene.astype(str)==g,methods].sum(axis=1).iloc[0])>1
   a.add_patch(Arc((0,0),1.99,1.99,theta1=hi-(j+1)/len(genes)*40,theta2=hi-j/len(genes)*40,color='#f57c00' if shared else '#f7c48b',lw=4))
  mid=np.deg2rad((lo+hi)/2);rot=((lo+hi)/2-90)%360
  if 90<rot<270:rot=(rot+180)%360
  a.text(1.18*np.cos(mid),1.18*np.sin(mid),m,fontsize=7,ha='center',va='center',rotation=rot)
 for i,j in itertools.combinations(range(len(methods)),2):
  m,n=methods[i],methods[j];genes=member.loc[(member[m].astype(int)>0)&(member[n].astype(int)>0),'Gene'].astype(str)
  for g in genes:
   t1=angles[m,g];t2=angles[n,g];p=.96*np.array([np.cos(t1),np.sin(t1)]);q=.96*np.array([np.cos(t2),np.sin(t2)])
   path=MPath([p,p*.25,q*.25,q],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]);a.add_patch(PathPatch(path,fc='none',ec='#8064a2',lw=.34,alpha=.28))
   links.append({'source_method':m,'target_method':n,'gene':g})
 p=save_panel(f,'Fig19_gene_chord');original=fitz.open(track(B/'Assets/original/Extended_Figures_submission.pdf'));d=fitz.open();page=d.new_page(width=760,height=940)
 page.show_pdf_page(fitz.Rect(23,28,367,372),fitz.open(p),0);page.show_pdf_page(fitz.Rect(390,38,732,355),original,11,clip=fitz.Rect(306,35,582,278))
 # Keep c-e together: their labels and plot extents overlap vertically in the source.
 page.show_pdf_page(fitz.Rect(42,403,718,918),original,11,clip=fitz.Rect(25,280,579,678))
 for x,y,label in [(15,25,'a'),(378,25,'b')]:page.insert_text((x,y),label,fontname='hebo',fontsize=13)
 out=O/'Extended_Data_Fig19.pdf';d.save(out,garbage=4,deflate=True);page.get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(out.with_suffix('.png'));out.with_suffix('.svg').write_text(page.get_svg_image());d.close()
 csv(pd.DataFrame(links),'Fig19_gene_chord_links');csv(pd.DataFrame(segments),'Fig19_gene_arc_positions')
 AUDIT['figures']['19']={'style':'original gene-overlap circos; original historical enrichment circos/heatmap and young-aged fields','methods':8,'gene_chords':len(links),'memberships':len(segments),'historical_panels':['b','c','d','e']}
 CAPTIONS['19']='a, Original-style gene-overlap circos plot rebuilt from the eight available top100 membership lists. Each arc represents a method and each purple chord connects the same literal gene identifier present in two lists; the inner orange/pale-orange track marks shared/unique genes; gene slots, chord links and overlap counts are supplied. This restores the original visual encoding in place of the interim overlap heatmap. b-c, Original enriched-term overlap circos and enrichment heatmap, retained as historical illustrations because matched term tables, background genes and database/version evidence have not been located. d-e, Original young/aged ReDeeM LatentVelo fields using the reported 784 shared genes; the exact matching run and gene-ID manifest remain unresolved. These historical panels are not evidence of newly completed analyses. No new enrichment, model fitting or aging comparison was performed.'

def additional_simulations():
 source=B/'Source_data/cellDancer_simulations';summary=pd.read_csv(track(source/'dataset_summary.csv'));inv=pd.read_csv(track(source/'dataset_inventory.csv')).set_index('dataset')
 units=[]
 for chunk in pd.read_csv(track(source/'per_unit_common_genes.csv.gz'),chunksize=100000):units.append(chunk[chunk.dataset.isin(['Mono-kinetic_1','Transcriptional_boost_1','Multi-forward_branching_1'])])
 units=pd.concat(units,ignore_index=True)
 families=list(FAMILIES);names=['Mono-kinetic','Transcriptional boost','Multi-forward branching']
 summary=summary[summary.family.isin(families)].copy()
 csv(summary,'Fig06_retained_condition_summary')
 csv(inv.reset_index().query('family in @families'),'Fig06_retained_dataset_inventory')
 maps={'velocyto':'velocyto','scvelo.sto':'scVelo(sto)','scvelo.dyn':'scVelo(dyn)','dynamo':'dynamo','velovae':'VeloVAE','unitvelo':'UniTVelo','deepvelo_vae':'DeepVelo(VAE)','deepvelo_gcn':'DeepVelo(GCN)','velovi':'veloVI','sdevelo':'SDEvelo','celldancer':'cellDancer','latentvelo':'LatentVelo','stt':'STT','bivi':'biVI†'}
 f=fig(820,655)
 for r,level in enumerate(['cell','gene']):
  for c,family in enumerate(families[:3]):
   a=ax(f,57+c*255,58+r*312,224,205);z=units[(units.dataset==family+'_1')&(units.level==level)];cos=pd.DataFrame({maps[m]:pd.Series(z[z.method==m].cosine.to_numpy()) for m in maps})
   boxes(a,cos);a.set_ylabel(('Cell-level' if r==0 else 'Gene-level')+' cosine' if c==0 else '');a.set_title(names[c]+f'\n{int(inv.loc[family+"_1","common_genes"])} shared genes',fontsize=9)
   panel(f,32+c*255,32+r*312,'abcdef'[r*3+c])
 save(f,6);AUDIT['figures']['6']={'style':'original ED4c vertical boxes with separate cell/gene rows','same_values_as_Fig5c':True,'independent_new_evidence':False}
 CAPTIONS['6']='Historical V2 detail view of the same three original cellDancer scenarios and saved common-finite-gene evaluation scores shown in Extended Data Fig. 5c, using the original submission Extended Data Fig. 4c vertical-boxplot style. a-c, Cell-level cosine; d-f, gene-level cosine. Noise is 0.2; shared gene counts are printed above the panels. Boxes show median/IQR, 1.5-IQR whiskers and outliers; methods are sorted by median. This figure is a detail view, not an independent experiment. BiVI (dagger) remains exploratory on continuous U/S input. STT is unavailable for forward branching. Missing and zero-norm results are not assigned zero. Source data and exclusions are unchanged.'
 f=fig(720,670);methods=list(maps);v=summary[(summary.gene_scope=='common_genes')&summary.method.isin(methods)]['median'];lo=v.min()-.01;hi=v.max()
 for r,level in enumerate(['cell','gene']):
  for c,family in enumerate(families):
   z=summary[(summary.gene_scope=='common_genes')&(summary.level==level)&(summary.family==family)];d=z.pivot(index='method',columns='series_index',values='median').reindex(index=methods,columns=range(1,11));rank=d.rank(ascending=False)
   a=ax(f,108+c*201,67+r*290,184,225)
   for i in range(len(methods)):
    if i%2==0:a.axhspan(i-.5,i+.5,color='#e5e5e5',zorder=0)
    for j in range(10):
     val=d.iloc[i,j]
     if np.isfinite(val):a.scatter(j,i,s=94*(val-lo)/(hi-lo),c=[BLUE(1-(rank.iloc[i,j]-1)/(d.iloc[:,j].count()-1))],ec='#555',lw=.35)
     else:a.text(j,i,'×',ha='center',va='center',color='#777',fontsize=8)
   a.set(xlim=(-.5,9.5),ylim=(13.5,-.5));a.set_yticks(range(14),[maps[m] for m in methods] if c==0 else ['']*14,fontsize=6);a.set_xticks(range(10),[f'{i*.2:.1f}' for i in range(1,11)],fontsize=5.5);a.tick_params(length=0)
   a.set_title(names[c],fontsize=8,pad=15);a.set_xlabel('Noise level',fontsize=7)
   if c==0:a.set_ylabel('Cell-level' if r==0 else 'Gene-level',fontsize=8)
   for sp in a.spines.values():sp.set_visible(False)
   for j in range(10):a.text(j,-.7,str(int(inv.loc[family+f'_{j+1}','common_genes'])),fontsize=4.5,ha='center')
   panel(f,86+c*201,34+r*290,'abcdef'[r*3+c])
 text(f,108,645,'Lighter blue: better within-condition rank. Circle area: normalized median cosine. Grey ×: missing. † biVI: exploratory.',fontsize=7)
 save(f,7);AUDIT['figures']['7']={'current_figure':6,'style':'original blue bubble encoding for three simulation scenarios','retained_raw_data_unchanged':True,'display_min':float(lo),'display_max':float(hi),'datasets':30,'families':families,'panels':6}
 CAPTIONS['7']=CAPTION_EN

if __name__=='__main__':
 if '--cell-dancer-only' in sys.argv or '--simulation-only' in sys.argv or '--real-metrics-only' in sys.argv:
  if (T/'plot_audit.json').exists():AUDIT.update(json.loads((T/'plot_audit.json').read_text()))
  if (T/'captions.json').exists():CAPTIONS.update(json.loads((T/'captions.json').read_text()))
 if '--real-metrics-only' in sys.argv:
  real_metrics();print('Rendered 9 (current Extended Data Fig. 7)',flush=True)
 elif '--simulation-only' not in sys.argv:
  additional_simulations();print('Rendered historical 6–7',flush=True)
 if '--cell-dancer-only' not in sys.argv and '--real-metrics-only' not in sys.argv:
  simulation_consistency();print('Rendered 8',flush=True)
  if '--simulation-only' not in sys.argv:
   real_metrics();print('Rendered 9',flush=True)
   refs(15,'scNT');refs(16,'scLT');print('Rendered 15–16',flush=True)
   chord();print('Rendered 19',flush=True)
 (T/'captions.json').write_text(json.dumps(CAPTIONS,ensure_ascii=False,indent=2))
 (T/'plot_audit.json').write_text(json.dumps(AUDIT,ensure_ascii=False,indent=2))

if __name__=='__main__':
 for old in ([9] if '--real-metrics-only' in sys.argv else [] if '--cell-dancer-only' in sys.argv else [8] if '--simulation-only' in sys.argv else [8,9,15,16,19]):
  new=old-2
  for ext in ['pdf','svg','png']:
   import shutil
   dest=T/'current_numbering';dest.mkdir(exist_ok=True)
   shutil.copy2(O/f'Extended_Data_Fig{old:02d}.{ext}',dest/f'Extended_Data_Fig{new:02d}.{ext}')
 (T/'numbering_note.txt').write_text('rendered/ retains historical V2 identifiers for traceability. current_numbering/ uses Manuscript_1010 identifiers. The V2 Fig06 detail and Fig07 noise series are omitted from the active release. No active release asset is overwritten.\n')
