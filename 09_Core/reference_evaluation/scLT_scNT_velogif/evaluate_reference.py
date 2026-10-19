"""Project new predictions, align by ID, and evaluate reference directions."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path('/task')
OUTPUTS = [('scvelo.sto.h5ad', 'scVelo (stochastic)'),
           ('scvelo.dyn.h5ad', 'scVelo (dynamic)'), ('sdevelo.h5ad', 'SDEvelo')]
SPECS = {
    'scLT': ('/inputs/scLT/LSK_lineage_cell1_ground_truth.h5ad', 'cce', 'ground_truth_velocity_cce',
             '/inputs/scLT/Clone_info_low_dim_sim.csv'),
    'scNT': ('/inputs/scNT/scNT_HSC_scvelo_process.h5ad', 'umap', 'scNT_velocity_umap',
             '/inputs/scNT/scNT.csv'),
}


def dense(value):
    return value.toarray() if sparse.issparse(value) else np.asarray(value)


def aligned_scores(pred_ids, pred_xy, pred_v, ref_ids, ref_xy, ref_v):
    pi, ri = pd.Index(pred_ids), pd.Index(ref_ids)
    if not pi.is_unique or not ri.is_unique:
        raise ValueError('Cell IDs must be unique')
    loc = pi.get_indexer(ri)
    present = loc >= 0
    px, rx, pv, rv = map(dense, (pred_xy, ref_xy, pred_v, ref_v))
    if px.shape != pv.shape or rx.shape != rv.shape or px.shape[1] != rx.shape[1]:
        raise ValueError('Coordinate and velocity dimensions do not match')
    status = np.full(len(ri), 'missing_prediction', dtype=object)
    scores = np.full(len(ri), np.nan)
    indices = np.flatnonzero(present)
    if len(indices):
        matched = loc[present]
        coordinates_match = np.isclose(px[matched], rx[present], rtol=1e-5, atol=1e-6).all(axis=1)
        a, b = pv[matched], rv[present]
        finite = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
        na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
        nonzero = (na > 0) & (nb > 0)
        status[indices] = 'valid'
        status[indices[~nonzero]] = 'zero_vector'
        status[indices[~finite]] = 'nonfinite_velocity'
        status[indices[~coordinates_match]] = 'coordinate_mismatch'
        valid = coordinates_match & finite & nonzero
        scores[indices[valid]] = np.clip(np.sum(a[valid] * b[valid], axis=1) / (na[valid] * nb[valid]), -1, 1)
    return pd.DataFrame({'cell_id': ri, 'cosine': scores, 'status': status})


def evaluate_dataset(dataset):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scanpy as sc
    import scvelo as scv

    reference, basis, vkey, historical = SPECS[dataset]
    ref = sc.read_h5ad(reference)
    out = ROOT / 'results' / dataset / 'reference_evaluation'
    out.mkdir(parents=True, exist_ok=True)
    rows, failures, distributions = [], [], []
    for filename, method in OUTPUTS:
        path = ROOT / 'results' / dataset / filename
        if not path.is_file():
            failures.append({'method': method, 'error': 'Missing result: ' + filename})
            continue
        try:
            obj = sc.read_h5ad(path)
            if not obj.obs_names.is_unique:
                raise ValueError('Duplicate predicted cell IDs')
            velocity = dense(obj.layers['velocity'])
            if velocity.shape != obj.shape:
                raise ValueError('Velocity shape does not match expression')
            keep = ~np.isnan(velocity).all(axis=0)
            if not keep.any():
                raise ValueError('All velocity genes are NaN')
            obj = obj[:, keep].copy()
            # Rebuild on all predicted cells, using the same VeloGIF projection
            # as the regular figures, before selecting reference cells.
            for key in list(obj.uns):
                if key.startswith('velocity_graph'):
                    del obj.uns[key]
            obj.obsm.pop('velocity_' + basis, None)
            np.random.seed(2024)
            scv.pp.neighbors(obj)
            scv.tl.velocity_graph(obj, basis=basis, vkey='velocity', n_jobs=10)
            scv.tl.velocity_embedding(obj, basis=basis, vkey='velocity')
            xy, v = obj.obsm['X_' + basis], obj.obsm['velocity_' + basis]
            pd.DataFrame({'cell_id': obj.obs_names, 'x': xy[:, 0], 'y': xy[:, 1],
                          'vx': v[:, 0], 'vy': v[:, 1]}).to_csv(out / (filename + '.embedding.csv'), index=False)
            per_cell = aligned_scores(obj.obs_names, xy, v, ref.obs_names,
                                      ref.obsm['X_' + basis], ref.obsm[vkey])
            per_cell.to_csv(out / (filename + '.cosine.csv'), index=False)
            counts = per_cell.status.value_counts()
            valid = per_cell.loc[per_cell.status == 'valid', 'cosine'].to_numpy()
            row = {'dataset': dataset, 'method': method, 'reference_cells': ref.n_obs,
                   'prediction_cells': obj.n_obs, 'valid_cells': len(valid),
                   'mean_cosine': float(np.mean(valid)) if len(valid) else np.nan,
                   'median_cosine': float(np.median(valid)) if len(valid) else np.nan}
            for reason in ('missing_prediction', 'zero_vector', 'nonfinite_velocity', 'coordinate_mismatch'):
                row[reason] = int(counts.get(reason, 0))
            rows.append(row)
            distributions.append((method, valid))
            print(dataset, row, flush=True)
            if not len(valid) or row['coordinate_mismatch']:
                failures.append({'method': method, 'error': 'No valid scores or coordinate mismatch; see per-cell statuses'})
        except Exception as exc:
            failures.append({'method': method, 'error': str(exc)})
            print('EVALUATION FAILED', method, repr(exc), flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(out / 'summary.csv', index=False)
    pd.DataFrame(failures, columns=['method', 'error']).to_csv(out / 'failures.csv', index=False)
    historical_frame = pd.read_csv(historical, index_col=0)
    historical_frame.describe().T.to_csv(out / 'historical_summary_NOT_CURRENT_RUN.csv')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    usable = [(name, v) for name, v in distributions if len(v)]
    if usable:
        axes[0].boxplot([v for _, v in usable], labels=[n for n, _ in usable], showfliers=False)
        axes[0].set_ylabel('Reference direction cosine')
        axes[1].bar([n for n, _ in usable], [np.mean(v) for _, v in usable])
        axes[1].set_ylabel('Mean reference direction cosine')
    for axis in axes:
        axis.set_ylim(-1.05, 1.05)
        axis.axhline(0, color='grey', lw=0.8)
        axis.tick_params(axis='x', rotation=20)
    fig.suptitle(dataset + ': newly computed predictions')
    fig.tight_layout()
    fig.savefig(out / 'reference_comparison.svg')
    fig.savefig(out / 'reference_comparison.png', dpi=150)
    plt.close(fig)
    tables = []
    for name in SPECS:
        p = ROOT / 'results' / name / 'reference_evaluation/summary.csv'
        if p.exists():
            tables.append(pd.read_csv(p))
    if tables:
        pd.concat(tables, ignore_index=True).to_csv(ROOT / 'results/reference_summary.csv', index=False)
    return int(bool(failures))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', choices=list(SPECS), required=True)
    args = parser.parse_args()
    raise SystemExit(evaluate_dataset(args.dataset))
