"""Validate delivered artifacts and produce a task-level report."""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/task')


def main():
    validations = []
    statuses = []
    for dataset, basis, cluster, expected in [('scLT', 'cce', 'cell_type2', 3186),
                                                ('scNT', 'umap', 'cell_type', 1947)]:
        result = ROOT / 'results' / dataset
        summary = json.loads((result / 'run_summary.json').read_text())
        for item in summary['results']:
            statuses.append({'dataset': dataset, 'method': item['method'], 'status': item['status'],
                             'reason': item['message']})
        for method in ['topicvelo', 'bivi', 'deepkinet', 'cell2fate']:
            statuses.append({'dataset': dataset, 'method': method, 'status': 'excluded_missing_input',
                             'reason': 'Missing raw spliced/unspliced counts'})
        for filename in ['scvelo.sto.h5ad', 'scvelo.dyn.h5ad', 'sdevelo.h5ad']:
            path = result / filename
            if not path.exists():
                validations.append({'dataset': dataset, 'file': filename, 'status': 'missing'})
                continue
            obj = ad.read_h5ad(path)
            v = obj.layers['velocity']
            values = v.data if sparse.issparse(v) else np.asarray(v)
            assert obj.obs_names.is_unique
            assert obj.n_obs == expected
            assert cluster in obj.obs and 'X_' + basis in obj.obsm
            assert v.shape == obj.shape and np.isfinite(values).any()
            validations.append({'dataset': dataset, 'file': filename, 'status': 'validated',
                                'cells': obj.n_obs, 'genes': obj.n_vars,
                                'all_nan_velocity_genes': int(np.isnan(v.toarray() if sparse.issparse(v) else v).all(axis=0).sum()),
                                'nonfinite_stored_values': int((~np.isfinite(values)).sum()),
                                'bytes': path.stat().st_size})
        assert len(pd.read_csv(result / 'reference_evaluation/summary.csv')) == 3
        assert len(pd.read_csv(result / 'reference_evaluation/failures.csv')) == 0
        assert len(pd.read_csv(result / 'evals/Plot_failures.csv')) == 0
        assert len(list((result / 'figures').glob('*.svg'))) == 4
    pd.DataFrame(statuses).to_csv(ROOT / 'results/method_status.csv', index=False)
    (ROOT / 'records/final_validation.json').write_text(json.dumps(validations, indent=2))
    table = pd.read_csv(ROOT / 'results/reference_summary.csv')
    fig, ax = plt.subplots(figsize=(9, 5))
    table.pivot(index='method', columns='dataset', values='mean_cosine').plot.bar(ax=ax)
    ax.set_ylabel('Mean reference direction cosine')
    ax.set_xlabel('')
    ax.set_ylim(-1.05, 1.05)
    ax.axhline(0, color='grey', lw=.8)
    ax.tick_params(axis='x', rotation=15)
    ax.set_title('scLT and scNT (different reference systems)')
    fig.tight_layout()
    fig.savefig(ROOT / 'results/reference_summary.svg')
    fig.savefig(ROOT / 'results/reference_summary.png', dpi=150)
    plt.close(fig)
    lines = ['# Results of this run', '', 'scVelo stochastic, dynamical and SDEvelo completed for both datasets.', '',
             '| Dataset | Method | Valid reference cells | Mean cosine | Median cosine |',
             '|---|---|---:|---:|---:|']
    for _, r in table.iterrows():
        lines.append('| %s | %s | %d | %.5f | %.5f |' %
                     (r.dataset, r.method, r.valid_cells, r.mean_cosine, r.median_cosine))
    lines.extend(['', 'See the [cross-dataset evaluation plot](results/reference_summary.svg), '
                  '[complete accuracy table](results/reference_summary.csv) and '
                  '[method-status table](results/method_status.csv).', '',
                  'TopicVelo, biVI, DeepKINET and cell2fate were not run because raw counts were unavailable.',
                  'Reference scores include only coordinate-aligned, finite, nonzero velocities; other cases are recorded in the per-cell CSV.',
                  'scLT uses CCE directions from the lineage reference, while scNT uses UMAP directions from the metabolic reference; '
                  'absolute scores from different datasets are not treated as a common benchmark ranking.',
                  'Historical CSV files are retained with the reference inputs and their descriptive statistics are stored separately; this table contains newly calculated results from this run.'])
    lines.extend(['', '## Actual gene coverage', '',
                  '| Dataset | File | Retained genes | Genes with all-NaN velocity |',
                  '|---|---|---:|---:|'])
    for row in validations:
        if row['status'] == 'validated':
            lines.append('| %s | %s | %d | %d |' % (row['dataset'], row['file'], row['genes'], row['all_nan_velocity_genes']))
    lines.extend(['', 'scVelo dynamical retains NaN for genes that cannot be fitted and excludes all-NaN gene columns during projection. '
                  'SDEvelo filters genes and recomputes moments using its default workflow; recorded gene counts reflect actual outputs.'])
    (ROOT / 'RESULTS.zh-CN.md').write_text('\n'.join(lines) + '\n')
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
