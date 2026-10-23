"""Two exports: the official plotting recipe and manuscript-comparable stream style."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
import phylovelo as pv
R=Path('/task');O=R/'results'
d=pd.read_csv(O/'phylovelo_directions.csv');palette=pd.read_csv(R/'inputs/cell_type_palette.csv')
colors=dict(zip(palette.cell_type,palette.color));xy=d[['x','y']].to_numpy();v=d[['vx','vy']].to_numpy()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
span=np.ptp(xy,axis=0);lo=xy.min(0)-span*.05;hi=xy.max(0)+span*.05
for style in ['tutorial','manuscript']:
 fig,ax=plt.subplots(figsize=(6.2,6.2));fig.subplots_adjust(left=.04,right=.97,bottom=.16,top=.91)
 for name,color in colors.items():
  m=d.cell_type==name;ax.scatter(xy[m,0],xy[m,1],s=15 if style=='tutorial' else 10,c=color,edgecolors='none',alpha=.65 if style=='tutorial' else .9)
 if style=='tutorial':pv.velocity_plot(xy,v,ax,'stream',radius=.8,lw_coef=1500,arrowsize=2)
 else:
  grid=[np.linspace(lo[i],hi[i],35) for i in range(2)];xx,yy=np.meshgrid(*grid);points=np.c_[xx.ravel(),yy.ravel()]
  dist,idx=cKDTree(xy).query(points,k=30);bw=np.mean(span)/35*.75
  weights=np.exp(-.5*(dist/bw)**2);mass=weights.sum(1)
  field=(v[idx]*weights[:,:,None]).sum(1)/np.maximum(mass[:,None],1e-10);field[mass<.2]=np.nan
  ax.streamplot(*grid,np.ma.masked_invalid(field[:,0].reshape(xx.shape)),np.ma.masked_invalid(field[:,1].reshape(xx.shape)),density=1,color='#222',linewidth=.7,arrowsize=.8)
 ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]));ax.set_aspect('equal');ax.axis('off')
 ax.set_title('PhyloVelo | KPTracer 3726_NT_T1\n754 cells; native projected directions',fontsize=11)
 fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=5,color=c,label=n) for n,c in colors.items()],loc='lower center',bbox_to_anchor=(.5,.01),ncol=3,frameon=False,fontsize=7)
 for suffix in ['pdf','svg','png']:fig.savefig(O/f'PhyloVelo_KPTracer_{style}.{suffix}',dpi=180)
 plt.close(fig)
print('Saved native tutorial and manuscript-style vector figures',flush=True)
