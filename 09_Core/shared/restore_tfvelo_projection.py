"""Rebuild only TFvelo graphs/projections from saved fits, as in the original notebook.
Run in scvelo_env. No model fitting. Source H5ADs are read-only.
"""
from pathlib import Path
import os, json, hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/velocity-reference-mpl')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/velocity-reference-numba')
os.environ.setdefault('NUMBA_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
import numpy as np
import pandas as pd
import scvelo as scv
import anndata
B=Path(__file__).resolve().parents[1]
D=B/'Source_data/reference_restoration'
records=[]
for ds,folder,basis in [('scNT','scNT_adata_graph','umap'),('scLT','scLT_multitag/adata_graph','cce')]:
    path=B.parent/folder/'TFvelo.h5ad'
    a=anndata.read_h5ad(path)
    previous=a.obsm['velocity_'+basis].copy()
    scaling=np.asarray(a.var['fit_scaling_y'],float)
    values=np.asarray(a.layers['velo_hat'])/scaling[None,:]
    # Unfitted genes are all-NaN columns; scVelo excludes these columns internally.
    fitted=np.isfinite(values).all(axis=0)
    assert (fitted | np.isnan(values).all(axis=0)).all() and fitted.any()
    a.layers['velocity']=values
    # Preserve the original notebook's scVelo graph defaults and velocity-gene selection.
    scv.tl.velocity_graph(a,vkey='velocity',n_jobs=1,basis=basis,show_progress_bar=False)
    scv.tl.velocity_embedding(a,basis=basis,vkey='velocity')
    v=a.obsm['velocity_'+basis]
    assert np.isfinite(v).all() and not np.allclose(v,previous)
    pd.DataFrame({'cell_id':a.obs_names,'vx':v[:,0],'vy':v[:,1]}).to_csv(D/(ds+'_TFvelo_restored_projection.csv'),index=False,float_format='%.17g')
    records.append({'dataset':ds,'source':str(path),'scvelo_version':scv.__version__,'anndata_version':anndata.__version__,'basis':basis,'cells':a.n_obs,'genes':a.n_vars,'velocity_definition':'velo_hat / fit_scaling_y','models_rerun':False,'velocity_graph_recomputed':True,'stale_projection_replaced':True,'fitted_genes':int(fitted.sum()),'velocity_genes_used':int((np.asarray(a.var.velocity_genes)&fitted).sum()),'nonzero_directions':int((np.linalg.norm(v,axis=1)>0).sum())})
    print(records[-1],flush=True)
(B/'Provenance/tfvelo_projection_restoration.json').write_text(json.dumps(records,indent=2))
