"""Reject empty fits and verify frozen KPTracer cells and coordinates."""
from pathlib import Path
import json, hashlib
import numpy as np
import anndata as ad
from scipy import sparse

T = Path(__file__).resolve().parents[1]
P = T.parent
manifest = json.loads((T/'records/input_manifest.json').read_text())
for relative, record in manifest['inputs'].items():
    if relative.endswith('.h5ad'):
        assert hashlib.sha256((P/relative).read_bytes()).hexdigest() == record['sha256']
reference = ad.read_h5ad(P/'velogif_XC/real_data/anndata/KPTracer/anndata.h5ad')
path = T/'results/tfvelo/tfvelo.h5ad'
model = ad.read_h5ad(path)
assert model.obs_names.equals(reference.obs_names)
assert np.array_equal(model.obsm['X_umap'], reference.obsm['X_umap'])
scale = np.asarray(model.var['fit_scaling_y'], dtype=float)
raw = model.layers['velo_hat']
raw = raw.toarray() if sparse.issparse(raw) else np.asarray(raw)
with np.errstate(divide='ignore', invalid='ignore'):
    velocity = raw / scale[None, :]
valid = np.isfinite(velocity).all(axis=0) & (scale > 0)
assert valid.any()
nonzero = np.count_nonzero(np.linalg.norm(velocity[:, valid], axis=1))
assert nonzero > 0
evaluation = json.loads((T/'evaluation/tfvelo.json').read_text())
assert evaluation['status'] == 'complete'
assert evaluation['fitted_genes'] == int(valid.sum())
report = dict(status='pass', cells=model.n_obs, input_genes=model.n_vars,
              genes_with_TF_mappings=int((model.var['n_TFs'] > 0).sum()),
              finite_fitted_genes=int(valid.sum()), nonzero_velocity_cells=int(nonzero),
              input_h5ad_hashes_unchanged=True, cell_ids_unchanged=True,
              umap_coordinates_unchanged=True, GDC=evaluation['GDC'],
              model_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(T/'records/tfvelo_validation.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
