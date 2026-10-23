"""Project and evaluate archived TopoVelo predictions; never fit a model.
Run with scVelo 0.2.5 in the existing velogif Docker environment. ROOT accepts
an input/output directory containing ST_{small,raw}_archived.h5ad.
"""
from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad,scvelo as scv
from scipy import sparse
from topovelo_historical_evaluation import evaluate,global_directional_correctness
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(os.environ.get('TOPO_HISTORICAL_ROOT',Path(__file__).resolve().parent));O=R/'derived';O.mkdir(exist_ok=True)
T=['Inner EPC 1','Inner EPC 2','Outer EPC 1','Outer EPC 2','P-TGC (Sec. P-TGC 1)','P-TGC (Sec. P-TGC 2)']
E=[(T[0],T[1]),(T[1],T[2]),(T[2],T[3]),(T[1],T[4]),(T[4],T[5])]
P=dict(zip(T,['#A8DADC','#FFD6A5','#B7E4C7','#FFADAD','#CDB4DB','#DDBEA9']))
def dense(v):return v.toarray() if sparse.issparse(v) else np.asarray(v)
def cos(x,y):
 x=x.astype(float);y=y.astype(float);nx=np.linalg.norm(x,axis=1);ny=np.linalg.norm(y,axis=1)
 ok=np.isfinite(x).all(1)&np.isfinite(y).all(1)&(nx>1e-12)&(ny>1e-12);z=np.full(len(x),np.nan);z[ok]=np.clip((x[ok]*y[ok]).sum(1)/(nx[ok]*ny[ok]),-1,1)
 return z,ok,nx,ny
def summary(x):
 x=np.asarray(x);x=x[np.isfinite(x)]
 return dict(valid_pairs=len(x),median=float(np.median(x)),q25=float(np.quantile(x,.25)),q75=float(np.quantile(x,.75)),mean=float(x.mean()),fraction_positive=float((x>0).mean()),fraction_negative=float((x<0).mean()),fraction_ge_0_9=float((x>=.9).mean()),median_angle_degrees=float(np.median(np.degrees(np.arccos(x)))))
objects={};records=[]
for scope in ['ST_small','ST_raw']:
 path=R/(scope+'_archived.h5ad');a=ad.read_h5ad(path);v=dense(a.layers['gat_velocity']);a.layers['velocity']=v.copy()
 a.var['velocity_genes']=np.isfinite(v).all(0)&(np.linalg.norm(v,axis=0)>0)
 np.random.seed(2024);scv.tl.velocity_graph(a,vkey='velocity',basis='spatial',n_jobs=8);scv.tl.velocity_embedding(a,basis='spatial',vkey='velocity')
 vv=a.obsm['velocity_spatial'];valid=(np.linalg.norm(v,axis=1)>0)&np.isfinite(v).all(1)
 metrics=evaluate(a,E,'clusters','velocity','X_spatial');inner={k:x for k,x in metrics['raw_inner_cluster'].items() if k in T}
 score=global_directional_correctness(inner,metrics['raw_cross_boundary'],608)
 record=dict(dataset=scope,cells=a.n_obs,genes=a.n_vars,GDC=score,denominator=608,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),velocity_key='gat_velocity',zero_velocity_ids=a.obs_names[~valid].tolist(),zero_velocity_indices=np.flatnonzero(~valid).tolist(),scvelo=scv.__version__,evaluation_neighbors='preserved uns.neighbors.indices',models_refit=False,inner={k:dict(n=len(x),mean=float(np.mean(x))) for k,x in inner.items()},edges=[dict(source=k[0],target=k[1],n=len(x),median=float(np.median(x))) for k,x in metrics['raw_cross_boundary'].items()])
 records.append(record);objects[scope]=a
 xy=a.obsm['X_spatial'];pd.DataFrame(dict(cell_id=a.obs_names,x=xy[:,0],y=xy[:,1],vx=vv[:,0],vy=vv[:,1],cell_type=a.obs.clusters.astype(str).values,velocity_valid=valid,projection_valid=np.isfinite(vv).all(1)&(np.linalg.norm(vv,axis=1)>0),method='TopoVelo')).to_csv(O/(scope+'_directions.csv'),index=False)
 if scope=='ST_small':
  fig,ax=plt.subplots(figsize=(5.94,4.2));a.uns['clusters_colors']=[P[x] for x in a.obs.clusters.cat.categories]
  scv.pl.velocity_embedding_stream(a,basis='spatial',color='clusters',palette=P,legend_loc='none',density=2,linewidth=1.35,arrow_size=1.55,arrow_color='#111111',size=371.2,alpha=.9,show=False,frameon=False,ax=ax,title='')
  for c in ax.collections:
   widths=np.asarray(c.get_linewidths(),dtype=float)
   if (~np.isfinite(widths)).any():c.set_linewidths(np.where(np.isfinite(widths),widths,0))
  for p in ax.patches:
   if not np.isfinite(p.get_linewidth()):p.set_visible(False)
  ax.set_position([0,0,1,1]);fig.savefig(O/'TopoVelo_field.pdf',bbox_inches='tight',pad_inches=0);plt.close(fig)
 print(json.dumps(record),flush=True)
a=objects['ST_raw'];b=objects['ST_small'];ids=b.obs_names;ia=a.obs_names.get_indexer(ids);assert (ia>=0).all();assert np.array_equal(a.obs.clusters.astype(str).values[ia],b.obs.clusters.astype(str).values)
translation=a.obsm['X_spatial'][ia]-b.obsm['X_spatial'];assert np.allclose(translation,translation[0])
genes=a.var_names[a.var_names.isin(b.var_names)];ja=a.var_names.get_indexer(genes);jb=b.var_names.get_indexer(genes)
x=dense(a.layers['gat_velocity'][ia][:,ja]);y=dense(b.layers['gat_velocity'][:,jb]);keep=np.isfinite(x).all(0)&np.isfinite(y).all(0)
gc,gv,gn1,gn2=cos(x[:,keep],y[:,keep]);v1=a.obsm['velocity_spatial'][ia].copy();v2=b.obsm['velocity_spatial'].copy()
# An embedding cannot define a direction for a zero gene velocity. Preserve raw
# projections above; zero such pairs here, so the ordinary cosine rule gives NA.
original_ok=(np.linalg.norm(dense(a.layers['gat_velocity'][ia]),axis=1)>0)&(np.linalg.norm(dense(b.layers['gat_velocity']),axis=1)>0)
v1[~original_ok]=0;v2[~original_ok]=0
sc,sv,sn1,sn2=cos(v1,v2)
xy=a.obsm['X_spatial'][ia];table=pd.DataFrame(dict(cell_id=ids,x=xy[:,0],y=xy[:,1],cell_type=b.obs.clusters.astype(str).values,method='TopoVelo',gene_cosine=gc,spatial_cosine=sc,gene_valid=gv,spatial_valid=sv,full_vx=v1[:,0],full_vy=v1[:,1],restricted_vx=v2[:,0],restricted_vy=v2[:,1],gene_full_norm=gn1,gene_restricted_norm=gn2,spatial_full_norm=sn1,spatial_restricted_norm=sn2))
table.to_csv(O/'paired_spots.csv',index=False);pd.DataFrame(dict(gene=genes,used=keep)).to_csv(O/'TopoVelo_shared_genes.csv',index=False)
ss=[];gs=[]
for metric in ['gene','spatial']:
 ss.append(dict(method='TopoVelo',metric=metric,paired_spots=608,shared_genes_used=int(keep.sum()),interpretation='gene_velocity_on_shared_gene_names' if metric=='gene' else 'fixed_spatial_direction',**summary(table[metric+'_cosine'])))
 for kind,part in table.groupby('cell_type'):gs.append(dict(method='TopoVelo',metric=metric,cell_type=kind,paired_spots=len(part),**summary(part[metric+'_cosine'])))
pd.DataFrame(ss).to_csv(O/'summary.csv',index=False);pd.DataFrame(gs).to_csv(O/'by_cell_type.csv',index=False)
(O/'evaluation.json').write_text(json.dumps(dict(inputs=records,coordinate_translation_full_minus_small=translation[0].tolist(),shared_genes=int(keep.sum()),paired_summary=ss,undefined_pair_ids=ids[~sv].tolist(),zero_velocity_GDC_policy='Retain archived zeros and legacy sklearn cosine zero convention; do not drop or impute them.',paired_zero_policy='Undefined directions excluded; spatial vectors set to zero for undefined archived gene velocities.'),indent=2)+'\n')
print('COMPLETE',ss,flush=True)
