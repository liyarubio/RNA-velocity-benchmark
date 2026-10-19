from pathlib import Path
import sys,importlib.util,json
import numpy as np,pandas as pd,anndata as ad
M=Path('/data1/liyaru/RNA_velocity_materials');T=M/'spatial_completion_20261010'
spec=importlib.util.spec_from_file_location('ev',M/'Manuscript_2607/Figure5_sources/runs/tfvelo_20260913/code/velogif/evaluation.py');ev=importlib.util.module_from_spec(spec);spec.loader.exec_module(ev)
dirs=pd.read_csv(T/'data/spatial_directions.csv');edges=list(pd.read_csv(T/'data/reference_edges.csv').itertuples(index=False,name=None));types=set(x for e in edges for x in e)
rows=[]
for folder,name in [('6_DeepveloSA','DeepVelo(VAE)'),('2_scvelo_stochastic','scVelo(sto)'),('7_DeepveloGB','DeepVelo(GCN)'),('5_unitvelo','UniTVelo')]:
 a=ad.read_h5ad(M/'ST_XC/result/velocity/MPSTA_raw'/folder/'anndata.h5ad');d=dirs[dirs.method==name].set_index('cell_id').loc[a.obs_names]
 a.obsm['velocity_spatial']=d[['vx','vy']].values
 metrics=ev.evaluate(a,edges,'clusters','velocity','X_spatial');n=int(a.obs.clusters.isin(types).sum())
 inner={k:v for k,v in metrics['raw_inner_cluster'].items() if k in types}
 row={'method':name,'cells':a.n_obs,'target_cells':n,'GDC_all_denom':ev.global_directional_correctness(inner,metrics['raw_cross_boundary'],a.n_obs),'GDC_target_denom':ev.global_directional_correctness(inner,metrics['raw_cross_boundary'],n),'inner':{k:float(np.mean(v)) for k,v in inner.items()},'cross':{str(k):float(np.median(v)) if len(v) else None for k,v in metrics['raw_cross_boundary'].items()}}
 rows.append(row);print(json.dumps(row),flush=True)
(T/'data/evaluation_probe.json').write_text(json.dumps(rows,indent=2))
