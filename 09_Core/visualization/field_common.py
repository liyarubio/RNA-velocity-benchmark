"""Field plotting primitive extracted from restore_reference_figures.py."""
import numpy as np
from scipy.spatial import cKDTree
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
