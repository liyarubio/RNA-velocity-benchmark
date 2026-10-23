"""KPTracer: official PhyloVelo tutorial on the manuscript's fixed 754 cells.

Native gene-wise slopes and native PhyloVelo embedding are retained separately.
No scVelo re-projection, tree pruning, or outcome-based parameter selection.
"""
from pathlib import Path
import sys, json, hashlib, time, traceback, copy, pickle, platform
import numpy as np, pandas as pd, scipy.sparse as sparse
import scanpy as sc, anndata as ad, phylovelo as pv
from Bio import Phylo
from evaluation_metrics import evaluate, global_directional_correctness
R=Path('/task');O=R/'results';O.mkdir(exist_ok=True)
SEED=2024;np.random.seed(SEED)
EDGES=[('AT2-like','Gastric-like'),('AT2-like','Lung progenitor-like'),('Endoderm-like','Pre-EMT'),('Pre-EMT','Early EMT-1')]
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def status(stage,**extra):
 d={'stage':stage,'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**extra}
 (R/'status.json').write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def metrics_for(a):
 m=evaluate(a,EDGES,'Cluster-Name','velocity','X_umap')
 denominator=int(a.obs['Cluster-Name'].isin({x for e in EDGES for x in e}).sum());assert denominator==706
 score=global_directional_correctness(m['raw_inner_cluster'],m['raw_cross_boundary'],denominator)
 return score,m,denominator
def main():
 started=time.monotonic();status('auditing_inputs')
 refpath=Path('/input/anndata/KPTracer/anndata.h5ad');ref=sc.read_h5ad(refpath)
 assert ref.n_obs==754 and ref.obs_names.is_unique
 treepath=R/'inputs/3726_NT_T1_tree.nwk';tree=Phylo.read(treepath,'newick')
 leaves={node.name:node for node in tree.get_terminals()}
 assert len(leaves)==len(tree.get_terminals())
 absent=[x for x in ref.obs_names if x not in leaves];assert not absent,absent
 # Topological path length, exactly as the tutorial; full tree is not pruned.
 depth=np.array([len(tree.get_path(leaves[x])) for x in ref.obs_names])
 assert np.unique(depth).size>1
 source=R/'inputs/adata_processed.nt.h5ad';full=ad.read_h5ad(source,backed='r')
 assert full.obs_names.is_unique and full.var_names.is_unique
 missing=ref.obs_names.difference(full.obs_names);assert not len(missing),missing.tolist()
 a=full[ref.obs_names].to_memory();full.file.close()
 assert a.obs_names.equals(ref.obs_names)
 xy=np.asarray(ref.obsm['X_umap'])
 original_xy=np.asarray(a.obsm['X_umap']);xy_error=float(np.abs(original_xy-xy).max())
 assert np.allclose(original_xy,xy,atol=1e-6),xy_error
 before_genes=a.n_vars;sc.pp.filter_genes(a,min_counts=10)
 X=a.to_df();assert np.isfinite(X.to_numpy()).all()
 constant=X.columns[X.nunique()==1].tolist()
 if constant:X=X.drop(columns=constant);a=a[:,X.columns].copy()
 assert (X.sum(0)>=10-1e-6).all() and a.n_vars>100
 cond=(xy[:,0]<2)&~((xy[:,1]<3)&(xy[:,0]<-5.5))&~((xy[:,1]<-5)&(xy[:,0]<1)&(xy[:,0]>-2))
 common=ref.var_names.intersection(a.var_names)
 x_delta=float(np.max(np.abs(dense(ref[:,common].X)-dense(a[:,common].X)))) if len(common) else None
 manifest={'tutorial':'https://phylovelo.readthedocs.io/en/latest/notebook/KPTracer-3726-t1.html#import-data',
  'tumor':'3726_NT_T1','cells':754,'tree_leaves':len(leaves),'matched_tree_leaves':754,'missing_tree_cells':absent,
  'tree_not_pruned':True,'depth_definition':'number of edges in original full-tree root-to-leaf path',
  'depth_min':int(depth.min()),'depth_max':int(depth.max()),'original_genes':before_genes,'filtered_genes':a.n_vars,
  'gene_filter':'sum of normalized expression >= 10, as official tutorial','constant_genes_removed':constant,
  'expression_source':'Zenodo 5847462 / KPTracer-Data/expression/adata_processed.nt.h5ad X',
  'manuscript_reference_gene_overlap':len(common),'expression_max_abs_difference_on_reference_genes':x_delta,
  'umap_max_abs_difference':xy_error,'tutorial_coordinate_filter_pass_count':int(cond.sum()),
  'cell_selection':'exact current manuscript 754 cells; no new coordinate filtering',
  'seed':SEED,'cutoff':.97,'alpha':.05,'target':'x_normed','exact':False,'projection_neighbors':754//3,
  'reference_sha256':sha(refpath),'tree_sha256':sha(treepath),'expression_sha256':sha(source)}
 (R/'records/input_audit.json').write_text(json.dumps(manifest,indent=2))
 pd.DataFrame({'cell_id':ref.obs_names,'tree_depth':depth,'cell_type':ref.obs['Cluster-Name'].astype(str).values}).to_csv(O/'cell_tree_alignment.csv',index=False)
 # Validate corrected eval1 before applying the same evaluation code to a new method.
 control=sc.read_h5ad('/input/result/velocity/KPTracer/2_scvelo_stochastic/anndata.h5ad')
 control_score,cm,denom=metrics_for(control)
 assert np.isclose(control_score,.6262533588720111,atol=1e-6),control_score
 (R/'records/control_validation.json').write_text(json.dumps({'GDC':control_score,'expected':.6262533588720111,'denominator':denom},indent=2))
 sd=pv.scData(x_normed=X,phylo_tree=tree,Xdr=pd.DataFrame(xy,index=ref.obs_names),cell_names=ref.obs_names.tolist(),cell_generation=depth)
 status('inference',cells=a.n_obs,genes=a.n_vars,tree_depth_range=[int(depth.min()),int(depth.max())])
 t=time.monotonic();pv.velocity_inference(sd,depth,cutoff=.97,alpha=.05,target='x_normed',exact=False)
 infer_seconds=time.monotonic()-t
 assert len(sd.megs)>0 and np.isfinite(sd.velocity).all() and np.count_nonzero(sd.velocity)>0
 gene=pd.DataFrame({'gene':X.columns,'native_gene_velocity':np.asarray(sd.velocity),'is_MEG':X.columns.isin(sd.megs)})
 gene.to_csv(O/'gene_velocities.csv',index=False)
 sd.latent_z.to_csv(O/'latent_expression.csv.gz')
 status('native_embedding',MEGs=len(sd.megs),inference_seconds=infer_seconds,neighbors=754//3)
 t=time.monotonic();pv.velocity_embedding(sd,target='x_normed')
 embedding_seconds=time.monotonic()-t
 v=np.asarray(sd.velocity_embeded);assert v.shape==(754,2) and np.isfinite(v).all()
 assert (np.linalg.norm(v,axis=1)>0).all()
 with (O/'phylovelo_scdata.pkl').open('wb') as f:pickle.dump(sd,f,pickle.HIGHEST_PROTOCOL)
 out=ad.AnnData(X=a.X.copy(),obs=ref.obs.copy(),var=a.var.copy())
 out.obs['phylo_depth']=depth;out.var['phylovelo_gene_velocity']=np.asarray(sd.velocity)
 out.var['phylovelo_MEG']=out.var_names.isin(sd.megs)
 out.obsm['X_umap']=xy.copy();out.obsm['velocity_umap']=v.copy()
 out.uns['phylovelo_protocol_json']=json.dumps(manifest)
 for key in ['connectivities','distances']:
  if key in ref.obsp:out.obsp[key]=ref.obsp[key].copy()
 out.uns['neighbors']=copy.deepcopy(ref.uns['neighbors'])
 out.write_h5ad(O/'phylovelo.h5ad',compression='gzip')
 pd.DataFrame({'cell_id':out.obs_names,'x':xy[:,0],'y':xy[:,1],'vx':v[:,0],'vy':v[:,1],
               'cell_type':out.obs['Cluster-Name'].astype(str).values}).to_csv(O/'phylovelo_directions.csv',index=False)
 status('evaluation',MEGs=len(sd.megs))
 # The model returns a single gene-wise slope vector. Replication here is only
 # an explicit adapter to the historical coherence metric, not cell-specific inference.
 out.layers['velocity']=np.tile(np.asarray(sd.velocity)[None,:],(out.n_obs,1))
 score,m,denom=metrics_for(out)
 del out.layers['velocity']
 scores=[]
 for metric,key in [('ICCoh','raw_inner_cluster'),('CBDir','raw_cross_boundary')]:
  for group,values in m[key].items():
   for j,value in enumerate(values):scores.append({'metric':metric,'group':str(group),'index':j,'score':value})
 pd.DataFrame(scores).to_csv(O/'evaluation_unit_scores.csv',index=False)
 result={'status':'complete','method':'PhyloVelo','cells':754,'genes':a.n_vars,'MEGs':len(sd.megs),
  'nonzero_gene_velocities':int(np.count_nonzero(sd.velocity)),'nonzero_cell_directions':754,
  'inference_seconds':infer_seconds,'embedding_seconds':embedding_seconds,'GDC':score,
  'reference_cell_denominator':denom,'reference_edges':EDGES,
  'ICCoh':{str(k):float(np.mean(x)) if len(x) else None for k,x in m['raw_inner_cluster'].items()},
  'CBDir':{s+' -> '+t:float(np.median(x)) if len(x) else None for (s,t),x in m['raw_cross_boundary'].items()},
  'boundary_cells':{s+' -> '+t:len(x) for (s,t),x in m['raw_cross_boundary'].items()},
  'evaluation_neighbor_graph':'fixed manuscript input neighbors, no recomputation',
  'native_projection':True,'velocity_semantics':'one gene-wise slope vector and native cell-specific 2D projected directions',
  'GDC_limitation':'Historical ICCoh adapter repeats the single gene-wise vector across cells, so within-cluster coherence is structurally 1; GDC is not an independent measure of cell-specific high-dimensional velocity quality.',
  'wall_seconds':time.monotonic()-started,'formal_figures_changed':False}
 (O/'evaluation.json').write_text(json.dumps(result,indent=2))
 existing=pd.read_csv(R/'inputs/current_GDC.csv')[['method','GDC']]
 combined=pd.concat([existing,pd.DataFrame([{'method':'PhyloVelo','GDC':score}])],ignore_index=True).sort_values('GDC',ascending=False,kind='stable')
 combined.insert(0,'rank',np.arange(1,len(combined)+1));combined.to_csv(O/'GDC_with_phylovelo.csv',index=False)
 status('complete',**{k:result[k] for k in ['cells','genes','MEGs','GDC','wall_seconds']})
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:
  status('failed',error=str(e),traceback=traceback.format_exc());raise
