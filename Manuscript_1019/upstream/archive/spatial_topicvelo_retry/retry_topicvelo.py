"""Resume missing TopicVelo transition inference from preserved historical fits.

Run only in the existing Docker topicvelo environment. Historical files remain
read-only; each ratio has an independent output directory and failure record.
"""
import argparse
import hashlib
import importlib
import inspect
import json
import os
from pathlib import Path
import sys
import time
import traceback
import types

import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import scvelo as scv
from topicvelo.tm_utils import TopicGeneFiltering, remove_U, FastTopics_cluster_assign
from topicvelo.transition_matrix import Combined_Topics_Transitions

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=True, choices=['ST_raw','ST_small'])
parser.add_argument('--attempt', choices=['diagnostic','zero_signal_qc','all_topics_qc'], default='diagnostic')
args = parser.parse_args()
BASE = Path('/work')
HERE = BASE/'spatial_topicvelo_retry_20261010'
source = (BASE/'ST_XC/result/velocity/MPSTA_raw/20_TopicVelo' if args.dataset=='ST_raw' else BASE/'velogif_XC/real_data/result/velocity/MPSTA/20_TopicVelo')
out = HERE/'results'/args.attempt/args.dataset
out.mkdir(parents=True, exist_ok=True)
os.chdir(out)
np.random.seed(2024)
report = {'dataset': args.dataset, 'status': 'running', 'seed': 2024,
          'attempt': args.attempt,
          'protocol': 'Resume saved fastTopics fit; percentile 35; stochastic transitions; topics_run explicitly records selected topics',
          'python': sys.version, 'scvelo': scv.__version__, 'source_files': {}}

def record():
    (out / 'status.json').write_text(json.dumps(report, indent=2) + '\n')

record()
start = time.time()
try:
    paths = [source / 'anndata.h5ad'] + list((source / 'redeem').glob('*.csv'))
    for p in paths:
        h = hashlib.sha256()
        with p.open('rb') as handle:
            for block in iter(lambda: handle.read(8 * 1024**2), b''): h.update(block)
        report['source_files'][str(p)] = h.hexdigest()
    adata = sc.read_h5ad(source / 'anndata.h5ad')
    folder = source / 'redeem'
    genes = pd.read_csv(folder / 'filtered_SU_Genes_names.csv', index_col=0).iloc[:, 0].astype(str).to_numpy()
    cells = pd.read_csv(folder / 'filtered_SU_Cells_names.csv', index_col=0).iloc[:, 0].astype(str).to_numpy()
    assert np.array_equal(cells, adata.obs_names.to_numpy()), 'Saved topic fit cell order mismatch'
    K = int(adata.obs['clusters'].nunique())
    weights = pd.read_csv(folder / ('HH_fastTopics_fit_CellWeights_k=%d.csv' % K), index_col=0).to_numpy()
    lfc = pd.read_csv(folder / ('HH_de_postmean_k=%d.csv' % K), index_col=0).to_numpy()
    lfsr = pd.read_csv(folder / ('HH_de_lfsr_k=%d.csv' % K), index_col=0).to_numpy()
    assert lfc.shape == lfsr.shape == (len(genes), K)
    assert weights.shape == (adata.n_obs, K)
    FastTopics_cluster_assign(adata, weights, t_type='fastTopics')
    selected, up, down = TopicGeneFiltering(lfc, lfsr, lfc_up_th=.5, lfc_down_th=-.5, lfsr_up_th=.001, lfsr_down_th=.001)
    top = [sorted(g) for g in remove_U(selected, genes)]
    adata.uns['top_genes'] = top
    report.update(cells=adata.n_obs, genes=adata.n_vars, historical_selected_topic=K-1,
                  topic_gene_counts=[len(g) for g in top],
                  repair='Restore required adata.uns[top_genes], omitted by historical runner; do not swallow inference exceptions')
    print(json.dumps(report, indent=2), flush=True)
    if not top[K-1]:
        raise ValueError('Original thresholds leave no topic genes; no undocumented fallback is permitted')
    from topicvelo.transition_matrix import get_cells_indices
    tc, other = get_cells_indices(adata, [K-1], topic_weights_th_percentile=[35], topic_type='fastTopics')
    def dense(x): return x.toarray() if sparse.issparse(x) else np.asarray(x)
    topics_to_run=list(range(K)) if args.attempt=='all_topics_qc' else [K-1]
    diagnostics={};excluded={};retained=[]
    for topic in topics_to_run:
        tc,_=get_cells_indices(adata,[topic],topic_weights_th_percentile=[35],topic_type='fastTopics')
        selected_vars=adata.var_names.get_indexer(top[topic])
        assert (selected_vars>=0).all()
        sub=adata[tc[0],selected_vars]
        diagnostic=[]
        for j,g in enumerate(top[topic]):
            row={'gene':g}
            for key in ['spliced','unspliced','Ms','Mu']:
                col=dense(sub.layers[key])[:,j]
                row[key+'_std']=float(np.std(col));row[key+'_nonzero']=int(np.count_nonzero(col))
            diagnostic.append(row)
        diagnostics[str(topic)]=diagnostic
        if args.attempt!='diagnostic':
            keep=[r['gene'] for r in diagnostic if r['Ms_std']>0 and r['Mu_std']>0 and r['spliced_nonzero']>0 and r['unspliced_nonzero']>0]
            excluded[str(topic)]=[g for g in top[topic] if g not in keep]
            top[topic]=keep
            if len(keep)<2:
                report.setdefault('excluded_insufficient_dimension_topics',[]).append(topic)
                continue
        retained.append(topic)
    topics_to_run=retained
    report['topic_gene_diagnostics_by_topic']=diagnostics
    report['excluded_zero_signal_genes']=excluded
    report['topics_run']=topics_to_run
    report['topic_gene_counts_after_qc']=[len(g) for g in top]
    report['protocol'] += '; '+args.attempt
    record()
    if not topics_to_run:raise ValueError('No topic genes with variable spliced and unspliced signal')
    vg = importlib.import_module('scvelo.tools.velocity_graph')
    # scVelo 0.2.x make_dense squeezes a single gene to 1-D inside GLS,
    # breaking stacking with the still-2-D second moments. Preserve the gene
    # axis only in that function; the least-squares calculation is unchanged.
    opt = importlib.import_module('scvelo.tools.optimization')
    vel = importlib.import_module('scvelo.tools.velocity')
    vel.make_dense = dense  # Keep (spots, 1) matrices in Velocity.__init__.
    gls = opt.leastsq_generalized
    gls_globals = dict(gls.__globals__)
    gls_globals['make_dense'] = lambda x: x.toarray() if sparse.issparse(x) else np.asarray(x)
    vel.leastsq_generalized = types.FunctionType(gls.__code__, gls_globals, gls.__name__, gls.__defaults__, gls.__closure__)
    original_gls=vel.leastsq_generalized
    def diagnostic_gls(x,y,x2,y2,res_std=None,res2_std=None,*a,**kw):
        report.setdefault('GLS_diagnostics',[]).append({'shape':list(x.shape),'res_std':np.asarray(res_std).tolist(),'res2_std':np.asarray(res2_std).tolist(),'finite_inputs':all(np.isfinite(dense(v)).all() for v in [x,y,x2,y2])})
        record()
        return original_gls(x,y,x2,y2,res_std,res2_std,*a,**kw)
    vel.leastsq_generalized=diagnostic_gls
    report['compatibility_fixes'] = ['Keep singleton gene axis in scVelo Velocity and GLS', 'Keep heterogeneous parallel results as list']
    original_neighbors = scv.pp.neighbors
    def compatible_neighbors(target, *a, **kw):
        result = original_neighbors(target, *a, **kw)
        distances = target.obsp['distances'].tocsr()
        rows = np.flatnonzero(np.asarray((distances > 0).sum(axis=1)).ravel() == 0)
        if len(rows):
            conn = target.obsp['connectivities'].tocsr()
            fixed = distances.tolil(copy=True)
            epsilon = np.finfo(np.float32).eps
            for row in rows:
                neighbors = conn.indices[conn.indptr[row]:conn.indptr[row+1]]
                neighbors = neighbors[neighbors != row]
                if not len(neighbors): raise ValueError('Isolated cell without connectivity')
                col = int(neighbors[0])
                fixed[row,col] = epsilon
                fixed[col,row] = max(float(fixed[col,row]),epsilon)
            target.obsp['distances'] = fixed.tocsr()
            report.setdefault('zero_distance_rows_encoded', []).append(int(len(rows)))
            print('Encoded zero-distance duplicate neighbor rows:', len(rows), flush=True)
        return result
    scv.pp.neighbors = compatible_neighbors
    original_parallelize = vg.parallelize
    def compatible_parallelize(callback, collection, *a, **kw):
        kw['as_array'] = False
        return original_parallelize(callback, collection, *a, **kw)
    vg.parallelize = compatible_parallelize
    original_graph = scv.tl.velocity_graph
    def bounded_graph(*a, **kw):
        kw['n_jobs'] = 2
        return original_graph(*a, **kw)
    scv.tl.velocity_graph = bounded_graph
    # Equivalent sparse embedding of each topic matrix into the global cell
    # index. The upstream implementation repeatedly inserts dense length-n rows
    # into CSR. Retain its float32 accumulation and six-decimal normalization.
    code = inspect.getsource(Combined_Topics_Transitions)
    begin = code.index('        #construct the combined transition matrix')
    end = code.index('        #save transition matrix', begin)
    code = code[:begin] + '''        combined_TM = csr_matrix((n,n), dtype=np.float32)
        for x in range(len(topics)):
            cells_k_indices = np.asarray(ttc_indices[x], dtype=int)
            local = TMs[x].tocoo()
            rows = cells_k_indices[local.row]
            cols = cells_k_indices[local.col]
            values = local.data * cells_weights_for_tm[rows, x]
            contribution = csr_matrix((values,(rows,cols)),shape=(n,n))
            combined_TM = (combined_TM + contribution).astype(np.float32)
        sums = np.asarray(combined_TM.sum(axis=1)).ravel()
        empty = np.flatnonzero(sums == 0)
        adata.uns['transition_empty_rows_before_native_fallback'] = empty
        inv = np.divide(1.0,sums,out=np.zeros_like(sums),where=sums!=0)
        combined_TM = combined_TM.multiply(inv[:,None]).tocsr()
        combined_TM.data = np.round(combined_TM.data,6)
        sums = np.asarray(combined_TM.sum(axis=1)).ravel()
        inv = np.divide(1.0,sums,out=np.zeros_like(sums),where=sums!=0)
        combined_TM = combined_TM.multiply(inv[:,None]).tocsr()
        if len(empty):
            conn=adata.obsp['connectivities'].tocsr()
            fallback_rows=[];fallback_cols=[];fallback_values=[]
            for i in empty:
                neighbors=conn.indices[conn.indptr[i]:conn.indptr[i+1]]
                if not len(neighbors): raise ValueError('No fallback neighbors')
                fallback_rows.extend([i]*len(neighbors));fallback_cols.extend(neighbors)
                fallback_values.extend([1/len(neighbors)]*len(neighbors))
            combined_TM += csr_matrix((np.asarray(fallback_values,dtype=np.float32),(fallback_rows,fallback_cols)),shape=(n,n))
''' + code[end:]
    namespace = dict(Combined_Topics_Transitions.__globals__)
    exec(code, namespace)
    completed_transition = namespace['Combined_Topics_Transitions']
    (out/'effective_transition_function.py').write_text(code)
    report['transition_assembly'] = 'Sparse global embedding; historical float32 accumulation and round-to-six retained'
    T, indices = completed_transition(
        adata, topics=topics_to_run, topic_weights_th_percentile=[35]*len(topics_to_run),
        velocity_type='stochastic', recompute=True, recompute_matrix=True,
        steady_state_perc=95, subset_save_prefix=str(out) + '/')
    T = sparse.csr_matrix(T)
    assert T.shape == (adata.n_obs, adata.n_obs)
    assert np.isfinite(T.data).all() and (T.data >= 0).all()
    row_sums = np.asarray(T.sum(axis=1)).ravel()
    assert np.allclose(row_sums, 1, atol=2e-6), 'Transition rows must sum to one'
    adata.uns['top_genes'] = {str(i): np.asarray(g, dtype=str) for i,g in enumerate(top)}
    adata.uns['spatial_completion'] = {'seed': 2024, 'selected_topic': K-1,
        'historical_fit_reused': True, 'minimum_usable_genes_per_topic':2, 'base_velocity_is_not_topicvelo_output': True, 'topics_run':topics_to_run, 'attempt':args.attempt}
    sparse.save_npz(out / 'native_transition.npz', T)
    adata.write_h5ad(out / 'topicvelo.h5ad', compression='gzip')
    report.update(status='complete', elapsed_seconds=time.time()-start,
                  native_uniform_neighbor_fallback_rows=len(adata.uns['transition_empty_rows_before_native_fallback']), transition_nnz=T.nnz, row_sum_min=float(row_sums.min()),
                  row_sum_max=float(row_sums.max()), topic_cells=[len(i) for i in indices],
                  output=str(out / 'topicvelo.h5ad'))
    record()
    print('COMPLETE', json.dumps(report), flush=True)
except Exception:
    report.update(status='failed', elapsed_seconds=time.time()-start, traceback=traceback.format_exc())
    record()
    raise
