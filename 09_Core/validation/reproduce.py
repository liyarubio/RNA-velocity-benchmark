"""Audit and replot the saved final-figure data; never fit models or fill missing values.

Run from any directory: python 09_Core/validation/reproduce.py --out /tmp/velocity-demo
Python >=3.9; numpy, pandas, scipy, matplotlib and openpyxl are required.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import importlib.metadata
import importlib
import numpy as np
import pandas as pd
from scipy import sparse

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from release_paths import REPO_ROOT as ROOT, source_path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'demo_output')
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    checks = []
    failures = []

    def check(name, ok, details):
        checks.append(dict(check=name, passed=bool(ok), details=details))
        if not ok:
            failures.append(name)

    def same(name, actual, expected, atol=1e-10):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        ok = a.shape == b.shape and np.allclose(a, b, atol=atol, rtol=0, equal_nan=True)
        err = float(np.nanmax(np.abs(a-b))) if a.shape == b.shape and np.isfinite(a-b).any() else None
        check(name, ok, dict(max_absolute_error=err, tolerance=atol, n=int(a.size)))

    def csv(rel):
        return pd.read_csv(source_path('data/' + rel))

    def sheet(file, name):
        return pd.read_excel(source_path('source_data/' + file), sheet_name=name, engine='openpyxl')

    manifest = ROOT / 'checksums.json'
    if manifest.exists():
        expected = json.loads(manifest.read_text())
        bad = [p for p, value in expected.items() if not (ROOT/p).is_file() or sha(ROOT/p) != value]
        check('release_file_integrity', not bad, dict(files=len(expected), mismatches=bad))

    # Figure 3: exact plotted pairing population and per-dataset envelope.
    paired = sheet('Figure03.xlsx', '3bc Paired')
    finite = paired.dropna(subset=['High', 'Low'])
    subset = finite[finite.Method != 'Dynamo_deterministic']
    check('Figure03_pairing_scope', len(finite) == 740 and len(subset) == 698,
          dict(main_pairs=len(finite), extended_pairs=len(subset), main_methods=finite.Method.nunique()))
    hub_raw = sheet('Figure03.xlsx', '3g HuBMAP')
    hub = hub_raw.set_index(hub_raw.columns[0]).drop(columns=['Dynamo_deterministic'], errors='ignore')
    # Best performance is an available-method maximum, never an executable method.
    best = sheet('Figure03.xlsx', '3g Best per dataset')
    same('Figure03_best_per_dataset', hub.max(axis=1), best.iloc[:, 1])
    check('HuBMAP_coverage', len(hub) == 176 and hub.shape[1] == 18,
          dict(datasets=len(hub), methods=hub.shape[1], missing_cells=int(hub.isna().sum().sum())))
    ten = sheet('Figure03.xlsx', '3h Dataset10').set_index('Method').drop(columns=['Mean', 'Median'], errors='ignore')
    check('Classical_dataset_coverage', ten.shape == (19, 10),
          dict(methods=ten.shape[0], datasets=ten.shape[1], missing_cells=int(ten.isna().sum().sum())))
    high = sheet('Figure03.xlsx','3a High dim')
    high = high.set_index(high.columns[0])
    check('Figure03_simulation_design', high.shape == (18,42), dict(methods=high.shape[0],datasets=high.shape[1]))
    # Current ED5 boxes are recomputed from the actual plotted unit scores.
    units = sheet('Extended_Data_Fig05.xlsx','c Unit scores')
    boxes = sheet('Extended_Data_Fig05.xlsx','c Box statistics').set_index(['dataset','level','method'])
    groups = units.groupby(['dataset','level','method']).cosine
    for name, values in [('median',groups.median()),('q1',groups.quantile(.25)),('q3',groups.quantile(.75)),('valid_units',groups.count())]:
        same('ED05_box_'+name, values.reindex(boxes.index), boxes[name])
    # ED6 topology means/SD come from paired dataset rows (not pooled cells).
    pairs6=sheet('Extended_Data_Fig06.xlsx','Paired cell-level values')
    pairs6['topology']=pairs6.dataset.str.replace(r'_seed\d+$','',regex=True)
    topology=sheet('Extended_Data_Fig06.xlsx','Topology scatter').set_index(['topology','method'])
    for val,prefix in [('high','high'),('low','low')]:
        groups6=pairs6.groupby(['topology','method'])[val]
        for stat,values in [('mean',groups6.mean()),('sd',groups6.std(ddof=1))]:
            same('ED06_'+prefix+'_'+stat,values.reindex(topology.index),topology[prefix+'_'+stat])
    ed7=sheet('Extended_Data_Fig07.xlsx','HuBMAP raw GDC')
    aliases={'scvelo':'scVelo(sto)','scvelo_dynamic':'scVelo(dyn)','DeepVelo_GB':'DeepVelo(GCN)','DeepVelo_SA':'DeepVelo(VAE)','Dynamo_stochastic':'dynamo'}
    comparable=hub.rename(columns=aliases)
    same('ED07_HuBMAP_source_identity',ed7.set_index('file_name')[comparable.columns],comparable)
    same('ED07_available_maximum',ed7['Best available'],comparable.max(axis=1))

    # Figure 6: recalculate the weighted summary, including ranks and coverage.
    summary = sheet('Figure06.xlsx', '6c Summary').set_index('Method')
    domains = ['simulation_cosine', 'scLT_cosine', 'scNT_cosine', 'HuBMAP_global',
               'Dataset10_global', 'ReDeeM_global', 'MPSTA_global', 'KPTracer_global']
    accuracy = ((summary[domains] + 1) / 2).mean(axis=1)
    speed = (1 / summary[['runtime_1000_seconds', 'runtime_10000_seconds']]).mean(axis=1)
    stability = summary[['reads_cosine_' + x for x in ['0.2', '0.4', '0.6', '0.8']]].mean(axis=1)
    score = (2*accuracy.rank(method='average') + speed.rank(method='average') + stability.rank(method='average')) / (4*len(summary))
    for col, values in [('Accuracy', accuracy), ('Scalability', speed), ('Stability', stability), ('Overall_Score', score), ('accuracy_domains_n', summary[domains].count(axis=1))]:
        same('Figure06_' + col, values, summary[col])
    summary.assign(recomputed_score=score).to_csv(out/'Figure06_recomputed.csv')
    top = score.sort_values(ascending=False).index[:3].tolist()
    check('Figure06_weighted_top3', top == ['scVelo(sto)', 'TopicVelo', 'UniTVelo'], top)
    same('Figure06_Overall_Score_rank', score.rank(ascending=False,method='min'), summary.Overall_Score_rank)

    def canonical(name):
        aliases = {'scvelo_stochastic':'scVelo(sto)', 'scvelo_dynamical':'scVelo(dyn)',
                   'Scvelo_stochastic':'scVelo(sto)', 'Scvelo_dynamical':'scVelo(dyn)',
                   'Dynamo_deterministic':'dynamo', 'Dynamo':'dynamo',
                   'Velocyto':'velocyto', 'CellDancer':'cellDancer', 'VeloVI':'veloVI',
                   'DeepVelo_GCN':'DeepVelo(GCN)', 'DeepVelo_VAE':'DeepVelo(VAE)',
                   'Cell2fate':'cell2fate', 'BiVI':'biVI'}
        return aliases.get(name, name)

    current = {}
    # ED13/14: recompute every cosine from identifier-aligned saved vectors.
    ref_error = {}
    for dataset, n in [('scNT', 1947), ('scLT', 1756)]:
        ref = csv(f'reference_restoration/{dataset}_reference_cells.csv').set_index('cell_id')
        vec = csv(f'reference_restoration/{dataset}_projected_directions.csv')
        means = csv(f'reference_restoration/{dataset}_method_summary.csv').set_index('method')
        cos = csv(f'reference_restoration/{dataset}_per_cell_cosine.csv').set_index('cell_id')
        check(dataset+'_unique_identifiers', ref.index.is_unique and not vec.duplicated(['method','cell_id']).any(), dict(reference_cells=len(ref), methods=vec.method.nunique()))
        check(dataset+'_reference_population', len(ref) == n and len(means) == 19, dict(cells=len(ref), configurations=len(means)))
        result = {}
        for method, v in vec.groupby('method', sort=False):
            v = v.set_index('cell_id').reindex(ref.index)
            a, b = v[['vx','vy']].to_numpy(), ref[['vx','vy']].to_numpy()
            den = np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)
            score_c = np.divide((a*b).sum(axis=1), den, out=np.full(len(a), np.nan), where=np.isfinite(den) & (den > 0))
            result[method] = score_c
            same(dataset+'_'+method+'_per_cell_cosine', score_c, cos.reindex(ref.index)[method], atol=2e-8)
        rec = pd.DataFrame(result, index=ref.index)
        same(dataset+'_method_means', rec.mean().reindex(means.index), means['mean'], atol=2e-8)
        same(dataset+'_method_medians', rec.median().reindex(means.index), means['median'], atol=2e-8)
        same(dataset+'_valid_denominators', rec.count().reindex(means.index), means.valid_cosines)
        rec.to_csv(out/(dataset+'_recomputed_cosines.csv'))
        current[dataset+'_cosine'] = means['mean'].rename(index=canonical)
        ref_error[dataset] = dict(cells=n, configurations=19)

    # ED16: spatial vectors provide an independent check; gene-space vectors
    # reside in external H5ADs, so only their saved summary is checked here.
    spots = csv('spatial_context_consistency/paired_spots.csv')
    aggregate = csv('spatial_context_consistency/summary.csv')
    grouped = csv('spatial_context_consistency/by_cell_type.csv')
    check('ED16_population', len(spots) == 20*608 and not spots.duplicated(['method','cell_id']).any()
          and (spots.groupby('method').cell_id.nunique() == 608).all(), dict(rows=len(spots), methods=spots.method.nunique()))
    a = spots[['full_vx','full_vy']].to_numpy()
    b = spots[['restricted_vx','restricted_vy']].to_numpy()
    den = np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)
    calc = np.divide((a*b).sum(axis=1), den, out=np.full(len(a),np.nan), where=den>0)
    same('ED16_spatial_vector_cosines', calc, spots.spatial_cosine, atol=1e-10)
    for metric in ['gene','spatial']:
        values = spots.groupby('method')[metric+'_cosine']
        target = aggregate[aggregate.metric == metric].set_index('method')
        for field, statistic in [('mean', values.mean()), ('median', values.median()), ('valid_pairs', values.count()), ('q25',values.quantile(.25)), ('q75',values.quantile(.75))]:
            same('ED16_'+metric+'_'+field, statistic.reindex(target.index), target[field])
        types = spots.groupby(['method','cell_type'])[metric+'_cosine'].median()
        target = grouped[grouped.metric==metric].set_index(['method','cell_type'])
        same('ED16_'+metric+'_by_type_median', types.reindex(target.index), target['median'])
    check('ED16_VeloAE_gene_missing', spots.loc[spots.method=='VeloAE','gene_cosine'].isna().all(), 'Unavailable gene-space result remains NA')

    # ED17: saved native transitions, final score table and scalar JSON agree.
    transition = sparse.load_npz(source_path('data/kptracer_completion/native_topicvelo_transition.npz'))
    check('ED17_native_transition', transition.shape == (754,754) and np.isfinite(transition.data).all()
          and (transition.data >= 0).all() and np.allclose(np.asarray(transition.sum(axis=1)).ravel(),1,atol=1e-6),
          dict(shape=list(transition.shape), nnz=int(transition.nnz)))
    kp = csv('kptracer_completion/GDC_display_order.csv').set_index('method')
    current['KPTracer_global'] = kp.GDC.rename(index=canonical)
    for method in ['tfvelo','topicvelo']:
        result = json.loads((source_path(f'data/kptracer_completion/{method}.json')).read_text())
        name = {'tfvelo':'TFvelo','topicvelo':'TopicVelo'}[method]
        same('ED17_'+method+'_GDC_table', [kp.loc[name,'GDC']], [result['GDC']])
        check('ED17_'+method+'_denominator', result['reference_cell_denominator'] == 706, 706)
    # Figure 5 is ST_small. Never replace it with ST_raw / ED15.
    spatial = sheet('Figure05.xlsx','5a GDC').set_index('Method').GDC.rename(index=canonical)
    topic = json.loads((source_path('data/spatial_completion/topicvelo_retry/ST_small.json')).read_text())
    tf = csv('spatial_completion/TFvelo_evaluation_summary.csv').set_index('dataset')
    spatial.loc['TopicVelo'] = topic['GDC']
    spatial.loc['TFvelo'] = tf.loc['ST_small','GDC']
    current['MPSTA_global'] = spatial
    differences = []
    for field, values in current.items():
        expected_methods = set(summary.index) | ({'TopoVelo'} if field == 'MPSTA_global' else set())
        check(field+'_current_method_names', set(values.index)==expected_methods, sorted(set(values.index)^expected_methods))
        for method, old in summary[field].items():
            new = values.get(method, np.nan)
            if not np.isclose(old, new, atol=1e-6, rtol=0, equal_nan=True):
                differences.append(dict(field=field,method=method,frozen_value=old,current_value=new,difference=new-old))
    pd.DataFrame(differences).to_csv(out/'Figure06_current_source_differences.csv',index=False)

    # ED18: server-returned gene membership, not a new enrichment analysis.
    metascape = []
    for n, expected_union in [(100,519),(200,809)]:
        membership = csv(f'metascape/top{n}/mapped_membership.csv')
        counts = csv(f'metascape/top{n}/list_counts.csv').set_index('method')
        totals = membership[counts.index].sum()
        same('ED18_top'+str(n)+'_mapped_counts', totals, counts.mapped_unique)
        check('ED18_top'+str(n)+'_union', len(membership)==expected_union and membership.Entrez_ID.is_unique,
              dict(mapped_union=len(membership), configurations=len(counts)))
        metascape.append(dict(top_n=n,mapped_union=len(membership),lists=len(counts)))

    # Representative demonstration plots are deliberately named as audit plots:
    # preserved publication artwork is in figures/, without claiming a model rerun.
    if not args.no_plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        matplotlib.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':8})
        fig, axes = plt.subplots(1,2,figsize=(10,5.7),constrained_layout=True)
        order = score.sort_values().index
        axes[0].barh(order, score.loc[order], color='#40769c')
        axes[0].set(xlabel='Weighted composite score', title='Figure 6: saved-score recalculation')
        spatial_median = spots.groupby('method').spatial_cosine.median().sort_values()
        axes[1].barh(spatial_median.index, spatial_median, color='#b56c38')
        axes[1].set(xlabel='Median cosine; 608 paired spots', title='Extended Data 16: spatial consistency', xlim=(-1,1))
        fig.savefig(out/'demo_summary.pdf'); fig.savefig(out/'demo_summary.svg');plt.close(fig)

    versions = {}
    for name in ['numpy','pandas','scipy','matplotlib','openpyxl']:
        try: versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: versions[name] = importlib.import_module(name).__version__
    report = dict(status='pass' if not failures else 'fail', checks=checks, failed_checks=failures,
                  number_of_checks=len(checks), elapsed_seconds=round(time.perf_counter()-start,3),
                  environment=dict(python=platform.python_version(), platform=platform.platform(), packages=versions),
                  weighted_composite_top3=top, composite_weights=dict(Accuracy=.5,Scalability=.25,Stability=.25), frozen_vs_current_differences=len(differences),
                  differences_are_documented_provenance_not_test_failures=True,
                  models_rerun=False, raw_inputs_recomputed=False,
                  scope='Saved numerical results, vector-derived cosines, coverage and native transition integrity. Historical artwork and external model inference are not reproduced by this demo.')
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:report[k] for k in ['status','number_of_checks','failed_checks','elapsed_seconds','frozen_vs_current_differences']},indent=2))
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()
