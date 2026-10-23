"""Compact vector panel with native backward vectors and fixed UMAP coordinates."""
from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/phylovelo-backward-panel')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parent;P=R.parent
v=pd.read_csv(R/'directions.csv');xy=v[['x','y']].to_numpy();vv=v[['vx_plot_backward','vy_plot_backward']].to_numpy()
orig=pd.read_csv(P/'results/phylovelo_directions.csv');assert np.allclose(vv,orig[['vx','vy']])
span=np.ptp(xy,axis=0);lo=xy.min(0)-span*.05;hi=xy.max(0)+span*.05
grid=[np.linspace(lo[i],hi[i],35) for i in range(2)];xx,yy=np.meshgrid(*grid);pts=np.c_[xx.ravel(),yy.ravel()]
dist,idx=cKDTree(xy).query(pts,k=30);bw=np.mean(span)/35*.75;weights=np.exp(-.5*(dist/bw)**2);mass=weights.sum(1)
f=(vv[idx]*weights[:,:,None]).sum(1)/np.maximum(mass[:,None],1e-10);f[mass<.2]=np.nan
fig=plt.figure(figsize=(122/72,102/72));ax=fig.add_axes([0,0,1,1])
for row in pd.read_csv(P/'inputs/cell_type_palette.csv').itertuples():
 m=orig.cell_type==row.cell_type;ax.scatter(xy[m,0],xy[m,1],s=4,c=row.color,edgecolors='none',alpha=.9)
ax.streamplot(*grid,np.ma.masked_invalid(f[:,0].reshape(xx.shape)),np.ma.masked_invalid(f[:,1].reshape(xx.shape)),density=1,color='#222',linewidth=.4,arrowsize=.5)
ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]));ax.set_aspect('equal');ax.axis('off')
fig.savefig(R/'PhyloVelo_backward_panel.pdf');plt.close(fig)
(R/'panel_style.json').write_text(json.dumps(dict(direction='native backward',sign=1,page_pt=[122,102],equal_aspect=True,grid=[35,35],neighbors=30,bandwidth_factor=.75,mass_threshold=.2,density=1,marker_area=4,linewidth=.4,arrowsize=.5,note='Same native vectors and smoothing as manuscript-style plot; point and line sizes adjusted for compact panel.'),indent=2)+'\n')
