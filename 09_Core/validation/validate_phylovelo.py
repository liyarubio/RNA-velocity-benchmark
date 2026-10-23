"""Recheck the packaged PhyloVelo direction convention and CBDir from saved H5AD."""
from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
import sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from release_paths import source_path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',type=Path,default=Path('/tmp/rna-velocity-1021'))
args=parser.parse_args();O=args.out;O.mkdir(parents=True,exist_ok=True)
P=source_path('data/kptracer_phylovelo_20261021/evaluation.json').parent
a=ad.read_h5ad(P/'results/phylovelo.h5ad');f=ad.read_h5ad(P/'presentation_backward_metrics_forward/phylovelo_forward_metrics.h5ad');e=json.loads((P/'evaluation.json').read_text())
v=np.asarray(a.obsm['velocity_umap'],float);xy=np.asarray(a.obsm['X_umap']);lab=a.obs['Cluster-Name'].astype(str).to_numpy()
assert a.obs_names.equals(f.obs_names) and np.array_equal(xy,f.obsm['X_umap'])
assert np.array_equal(v,f.obsm['velocity_umap_backward']) and np.array_equal(-v,f.obsm['velocity_umap_forward']) and np.array_equal(-v,f.obsm['velocity_umap'])
assert e['GDC_comparable'] is None and e['GDC_rank'] is None and 'velocity' not in f.layers
if 'indices' in a.uns.get('neighbors',{}):nb=a.uns['neighbors']['indices']
else:
 m=a.obsp['connectivities'].tocsr();nb=[m.indices[m.indptr[i]:m.indptr[i+1]] for i in range(a.n_obs)]
rows=[];included=set()
for edge in e['CBDir']:
 s,t=edge['source'],edge['target'];included.update([s,t]);scores=[]
 for i in np.flatnonzero(lab==s):
  js=np.asarray(nb[i],int);js=js[lab[js]==t]
  if not len(js):continue
  delta=np.asarray(xy[js]-xy[i],float);den=np.linalg.norm(delta,axis=1)*np.linalg.norm(v[i]);c=np.divide(delta@(-v[i]),den,out=np.zeros(len(js)),where=den>0)
  scores.append(float(c.mean()))
 assert len(scores)==edge['n_boundary_cells'] and abs(np.median(scores)-edge['CBDir_median'])<1e-12
 rows.append({'source':s,'target':t,'n':len(scores),'median':float(np.median(scores))})
g=np.asarray(a.var['phylovelo_gene_velocity'],float);coh=float(g@g/(np.linalg.norm(g)**2));assert abs(coh-1)<1e-12
assert all(r['median']>0 for r in rows) and np.isin(lab,list(included)).sum()==706
plots=pd.read_csv(source_path('data/figure5_1021/PhyloVelo_plot_vectors.csv'));assert np.array_equal(plots.cell_id,a.obs_names)
assert np.allclose(plots[['vx','vy']],v,rtol=0,atol=1e-16)
result={'status':'pass','cells':a.n_obs,'genes':a.n_vars,'nonzero_gene_velocities':int(np.count_nonzero(g)),'reference_cells':706,'CBDir':rows,'direction_negated_exactly_once_for_metrics':True,'plot_native_backward_preserved':True,'ICCoh_global_vector_adapter':coh,'adapted_GDC':float(e['GDC_adapted_forward']),'comparable_GDC':None,'comparable_rank':None,'model_fitting':False}
(O/'phylovelo_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
