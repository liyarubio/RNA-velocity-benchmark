"""Evaluate reruns on the four ED16 reference relations and frozen UMAP.

TopicVelo GDC explicitly uses a native-transition displacement proxy in Ms;
it never uses the historical base-scVelo velocity layer.
"""
from pathlib import Path
import argparse,json,hashlib,importlib,traceback
import numpy as np,pandas as pd
from scipy import sparse
import scanpy as sc,scvelo as scv
from evaluation_metrics import evaluate,global_directional_correctness
B=Path('/task');M=Path('/materials')
EDGES=[('AT2-like','Gastric-like'),('AT2-like','Lung progenitor-like'),('Endoderm-like','Pre-EMT'),('Pre-EMT','Early EMT-1')]
ref=sc.read_h5ad(M/'velogif_XC/real_data/anndata/KPTracer/anndata.h5ad')
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def run(method):
    np.random.seed(2024)
    if method=='control':path=M/'velogif_XC/real_data/result/velocity/KPTracer/2_scvelo_stochastic/anndata.h5ad'
    else:path=B/'results'/method/(method+'.h5ad')
    info={'method':method,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
          'scvelo':scv.__version__,'reference_edges':EDGES,'historical_protocol':'corrected eval1; 706 reference-participating cells denominator; incident-edge signs'}
    a=sc.read_h5ad(path);loc=ref.obs_names.get_indexer(a.obs_names)
    assert a.n_obs==754 and np.all(loc>=0)
    assert np.allclose(a.obsm['X_umap'],ref.obsm['X_umap'][loc])
    a.obs['Cluster-Name']=ref.obs['Cluster-Name'].iloc[loc].values
    if method=='topicvelo':
        T=sparse.load_npz(B/'results/topicvelo/native_transition.npz')
        assert T.shape==(754,754) and np.isfinite(T.data).all() and (T.data>=0).all()
        assert np.allclose(np.asarray(T.sum(1)).ravel(),1,atol=2e-6)
        X=dense(a.layers['Ms']);D=np.asarray(T@X)-X
        a.layers['velocity']=D
        a.var['velocity_genes']=np.isfinite(D).all(0)&(np.linalg.norm(D,axis=0)>0)
        # Explicit native T: no graph is inferred from the historical base velocity.
        scv.tl.velocity_embedding(a,basis='umap',vkey='velocity',T=T.copy())
        info.update(prediction_type='native_transition_displacement_proxy',
                    high_dimensional_definition='T @ Ms - Ms',embedding_definition='scVelo velocity_embedding with explicit native T',
                    transition_nnz=T.nnz,transition_rows=754,topic_cells=490,selected_topic=7,
                    uses_historical_base_scvelo_velocity=False)
    elif method=='tfvelo':
        scale=np.asarray(a.var['fit_scaling_y'],float)
        with np.errstate(divide='ignore',invalid='ignore'):V=dense(a.layers['velo_hat'])/scale[None,:]
        valid=np.isfinite(V).all(0)&(scale>0)
        assert valid.any();info['fitted_genes']=int(valid.sum());info['input_genes']=a.n_vars
        a=a[:,valid].copy();a.layers['velocity']=V[:,valid]
        for key in list(a.uns):
            if key.startswith('velocity_graph'):del a.uns[key]
        a.obsm.pop('velocity_umap',None)
        a.var['velocity_genes']=True
        vg=importlib.import_module('scvelo.tools.velocity_graph');original=vg.parallelize
        def parallelize(callback,collection,*args,**kw):
            kw['as_array']=False;return original(callback,collection,*args,**kw)
        vg.parallelize=parallelize
        scv.tl.velocity_graph(a,basis='umap',vkey='velocity',n_jobs=2)
        scv.tl.velocity_embedding(a,basis='umap',vkey='velocity')
        info.update(prediction_type='TFvelo fitted gene velocity',high_dimensional_definition='velo_hat / fit_scaling_y; finite fitted genes only')
    # Historical eval1 preserves saved neighbor indices when present and only
    # reconstructs them from connectivities when absent.
    metrics=evaluate(a,EDGES,'Cluster-Name','velocity','X_umap')
    reference_types={x for edge in EDGES for x in edge}
    denominator=int(a.obs['Cluster-Name'].isin(reference_types).sum())
    assert denominator==706
    score=global_directional_correctness(metrics['raw_inner_cluster'],metrics['raw_cross_boundary'],denominator)
    xy=np.asarray(a.obsm['X_umap']);v=np.asarray(a.obsm['velocity_umap'])
    assert np.isfinite(v).all() and np.any(v)
    info.update(status='complete',cells=a.n_obs,genes=a.n_vars,GDC=score,reference_cell_denominator=denominator,
                nonzero_direction_cells=int((np.linalg.norm(v,axis=1)>0).sum()),
                ICCoh={k:float(np.mean(x)) if len(x) else None for k,x in metrics['raw_inner_cluster'].items()},
                CBDir={s+' -> '+t:float(np.median(x)) if len(x) else None for (s,t),x in metrics['raw_cross_boundary'].items()},
                boundary_cells={s+' -> '+t:len(x) for (s,t),x in metrics['raw_cross_boundary'].items()})
    if method=='control':
        expected=.6262533588720111;info['historical_expected_GDC']=expected;info['absolute_difference']=abs(score-expected)
        assert np.isclose(score,expected,atol=1e-6),(score,expected)
    pd.DataFrame({'cell_id':a.obs_names,'x':xy[:,0],'y':xy[:,1],'vx':v[:,0],'vy':v[:,1],
                  'cell_type':a.obs['Cluster-Name'].astype(str).values}).to_csv(B/'evaluation'/(method+'_directions.csv'),index=False)
    rows=[]
    for kind,groups in [('ICCoh',metrics['raw_inner_cluster']),('CBDir',metrics['raw_cross_boundary'])]:
        for group,scores in groups.items():
            for i,value in enumerate(scores):rows.append({'metric':kind,'group':str(group),'index':i,'score':value})
    pd.DataFrame(rows).to_csv(B/'evaluation'/(method+'_unit_scores.csv'),index=False)
    if method!='control':a.write_h5ad(B/'results'/method/'evaluated.h5ad',compression='gzip')
    (B/'evaluation'/(method+'.json')).write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('method',choices=['control','topicvelo','tfvelo']);args=parser.parse_args()
    run(args.method)
