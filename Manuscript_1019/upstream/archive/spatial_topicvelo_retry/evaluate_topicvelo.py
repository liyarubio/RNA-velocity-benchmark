"""Evaluate native TopicVelo transitions; gene direction is T @ Ms - Ms."""
from pathlib import Path
import argparse, json, importlib.util, hashlib
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity

B=Path('/work');TASK=B/'spatial_topicvelo_retry_20261010'
parser=argparse.ArgumentParser();parser.add_argument('--attempt',default='all_topics_qc');args=parser.parse_args()
OUT=TASK/'evaluation'/args.attempt;OUT.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('ev',B/'Manuscript_2607/Figure5_sources/runs/tfvelo_20260913/code/velogif/evaluation.py')
ev=importlib.util.module_from_spec(spec);spec.loader.exec_module(ev)
edges=list(pd.read_csv(B/'Manuscript_1010/Source_data/spatial_completion/reference_edges.csv').itertuples(index=False,name=None))
types={x for e in edges for x in e}
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def evaluate_valid(a,V,v):
    labels=a.obs.clusters.to_numpy();xy=np.asarray(a.obsm['X_spatial']);neighbors=ev._neighbor_indices(a)
    high=np.isfinite(V).all(1)&(np.linalg.norm(V,axis=1)>1e-12)
    low=np.isfinite(v).all(1)&(np.linalg.norm(v,axis=1)>1e-12)
    inner={};cross={}
    for cluster in sorted(types):
        values=[]
        for i in np.flatnonzero((labels==cluster)&high):
            same=np.asarray(neighbors[i],dtype=int);same=same[(labels[same]==cluster)&high[same]]
            if len(same):values.append(float(cosine_similarity(V[i:i+1],V[same]).mean()))
        inner[cluster]=values
    for source,target in edges:
        values=[]
        for i in np.flatnonzero((labels==source)&low):
            js=np.asarray(neighbors[i],dtype=int);js=js[labels[js]==target]
            if len(js):values.append(float(cosine_similarity(xy[js]-xy[i],v[i:i+1]).mean()))
        cross[(source,target)]=values
    return inner,cross,high,low
for scope in ['ST_raw','ST_small']:
    folder=TASK/'results'/args.attempt/scope
    status=json.loads((folder/'status.json').read_text())
    if status['status']!='complete':continue
    a=sc.read_h5ad(folder/'topicvelo.h5ad');T=sparse.load_npz(folder/'native_transition.npz')
    assert np.isfinite(T.data).all() and (T.data>=0).all()
    assert np.allclose(np.asarray(T.sum(1)).ravel(),1,atol=2e-6)
    xy=np.asarray(a.obsm['X_spatial'])
    if scope=='ST_raw':
        ref=pd.read_csv(B/'Manuscript_1010/Source_data/spatial_completion/spatial_cells.csv').set_index('cell_id').loc[a.obs_names]
        assert np.array_equal(xy,ref[['x','y']].values)
    X=dense(a.layers['Ms']);V=T@X-X
    assert np.isfinite(V).all() and np.any(V)
    a.layers['velocity']=V
    a.var['velocity_genes']=np.isfinite(V).all(0)&(np.linalg.norm(V,axis=0)>0)
    # Self-only transition rows have no direction. Avoid global autoscaling,
    # which propagates their NaN into every row; retain undefined rows as NaN.
    scv.tl.velocity_embedding(a,basis='spatial',vkey='velocity',T=T.copy(),autoscale=False)
    v=np.asarray(a.obsm['velocity_spatial'])
    inner,cross,high,low=evaluate_valid(a,V,v)
    assert high.any() and low.any()
    denominator=int(a.obs.clusters.isin(types).sum());assert denominator==608
    gdc=ev.global_directional_correctness(inner,cross,denominator)
    assert np.isfinite(gdc)
    pd.DataFrame({'cell_id':a.obs_names,'x':xy[:,0],'y':xy[:,1],'vx':v[:,0],'vy':v[:,1],'cell_type':a.obs.clusters.astype(str).values}).to_csv(OUT/(scope+'_directions.csv'),index=False)
    report={'dataset':scope,'status':'complete','attempt':args.attempt,'GDC':gdc,'spots':a.n_obs,'reference_spot_denominator':denominator,
        'topics_run':status['topics_run'],'topic_gene_counts_after_qc':status['topic_gene_counts_after_qc'],
        'native_uniform_neighbor_fallback_rows':status['native_uniform_neighbor_fallback_rows'],
        'native_transition_sha256':hashlib.sha256((folder/'native_transition.npz').read_bytes()).hexdigest(),
        'gene_direction':'native T @ Ms - Ms; transition-displacement proxy, not a fitted gene velocity',
        'spatial_direction':'scVelo velocity_embedding with explicit native T; autoscale=False; undefined directions excluded from scoring',
        'uses_historical_base_scvelo_velocity':False,'scvelo':scv.__version__,
        'nonzero_direction_spots':int(low.sum()),'invalid_spatial_direction_spots':a.obs_names[~low].tolist(),
        'invalid_high_dimensional_direction_spots':a.obs_names[~high].tolist(),
        'ICCoh_spots_by_type':{k:len(x) for k,x in inner.items()},
        'ICCoh':{k:float(np.mean(x)) for k,x in inner.items()},
        'CBDir':{str(k):float(np.median(x)) if len(x) else None for k,x in cross.items()}}
    (OUT/(scope+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    a.write_h5ad(folder/'evaluated.h5ad',compression='gzip')
    print(json.dumps(report,indent=2),flush=True)
