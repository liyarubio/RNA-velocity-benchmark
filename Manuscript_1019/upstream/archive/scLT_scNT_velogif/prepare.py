"""Prepare isolated inputs and YAMLs; never rewrite source data."""
import json
from pathlib import Path

import anndata as ad
import numpy as np

ROOT = Path('/task')
SOURCE = Path('/inputs')
SPECS = {
    'scLT': ('scLT/LSK_lineage.h5ad', 'cell_type2', 'cce', 3186),
    'scNT': ('scNT/scNT_HSC_scvelo_process.h5ad', 'cell_type', 'umap', 1947),
}


def main():
    records = {}
    for dataset, (relative, cluster, embedding, cells) in SPECS.items():
        obj = ad.read_h5ad(SOURCE / relative)
        assert obj.n_obs == cells and obj.n_vars == 2000
        assert obj.obs_names.is_unique and obj.var_names.is_unique
        removed = {}
        for attr in ('layers', 'obsm', 'obsp', 'uns'):
            mapping = getattr(obj, attr)
            keys = [k for k in mapping if 'velocity' in k.lower()]
            removed[attr] = keys
            for key in keys:
                del mapping[key]
        for attr in ('obs', 'var'):
            frame = getattr(obj, attr)
            keys = [k for k in frame if 'velocity' in k.lower() or k.startswith('fit_')]
            removed[attr] = keys
            frame.drop(columns=keys, inplace=True)
        # Evaluation/lineage reference information must not reach estimators.
        for key in list(obj.uns):
            if key.startswith(('cell1_', 'cell2_', 'cos_sim_')):
                removed['uns'].append(key)
                del obj.uns[key]
        if dataset == 'scNT':
            for key in list(obj.layers):
                if key not in ('spliced', 'unspliced', 'Ms', 'Mu'):
                    removed['layers'].append(key)
                    del obj.layers[key]
            for key in ('dynamics', 'pp', 'velocyto_SVR', 'grid_velocity_umap'):
                if key in obj.uns:
                    removed['uns'].append(key)
                    del obj.uns[key]
            keys = [k for k in obj.obs if k == 'time' or k == 'ntr' or 'new' in k or 'total' in k]
            removed['obs'].extend(keys)
            obj.obs.drop(columns=keys, inplace=True)
        assert cluster in obj.obs and 'X_' + embedding in obj.obsm
        assert np.isfinite(obj.obsm['X_' + embedding]).all()
        target = ROOT / 'prepared' / (dataset + '.h5ad')
        target.parent.mkdir(exist_ok=True)
        if target.exists():
            raise FileExistsError(target)
        obj.write_h5ad(target)
        config = {
            'name': dataset,
            'inputs': {'rna': '../prepared/' + dataset + '.h5ad'},
            'output_dir': '../results/' + dataset,
            'methods': ['scvelo', 'sdevelo'],
            'runtime': {'jobs': 10, 'device': 'cuda:0', 'seed': 2024,
                        'embedding': embedding, 'cluster_key': cluster,
                        'gene_number': 2000, 'velocity_layer': 'velocity'},
            'evaluation': {'edges': []}, 'visualization': {'palette': {}},
        }
        # JSON is a YAML subset; the preparation environment need not install PyYAML.
        (ROOT / 'configs' / (dataset + '.yaml')).write_text(json.dumps(config, indent=2))
        records[dataset] = {'source': str(SOURCE / relative), 'shape': list(obj.shape),
                            'removed': removed, 'cluster_key': cluster, 'embedding': embedding}
    (ROOT / 'records/preparation.json').write_text(json.dumps(records, indent=2))
    exclusions = {k: 'Missing raw spliced/unspliced counts; normalized values are not rounded.'
                  for k in ('topicvelo', 'bivi', 'deepkinet', 'cell2fate')}
    (ROOT / 'records/method_scope.json').write_text(json.dumps({
        'enabled': ['scvelo', 'sdevelo'], 'excluded_new_methods': exclusions,
        'other_methods': 'Not requested under the updated scope',
        'source_commit': '263fe2fb7b36af751a5e8d075ae3e9d9b243df24',
    }, indent=2))
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
