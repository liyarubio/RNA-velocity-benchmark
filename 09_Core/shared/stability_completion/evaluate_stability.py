"""Evaluate saved cell2fate velocity and completed TopicVelo native transitions.

All ratios are aligned by cell/gene identity. TopicVelo uses expected transition
displacement in ONE fixed full-depth Ms gene coordinate system, never the base
scVelo velocity layer. Native transition-row cosine is exported separately too.
"""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
RATIOS = ['0_2', '0_4', '0_6', '0_8', '1_0']
parser = argparse.ArgumentParser()
parser.add_argument('--method', choices=['cell2fate', 'topicvelo'], required=True)
args = parser.parse_args()
method = args.method
out = HERE / 'evaluation' / method
out.mkdir(parents=True, exist_ok=True)

def ids(f, axis):
    group = f[axis]
    a = group[group.attrs.get('_index', '_index')].asstr()[:]
    assert len(set(a)) == len(a), 'Duplicate identifiers'
    return a

def matrix(group):
    if isinstance(group, h5py.Dataset): return group[:]
    encoding = group.attrs.get('encoding-type', 'csr_matrix')
    ctor = sparse.csc_matrix if encoding == 'csc_matrix' else sparse.csr_matrix
    return ctor((group['data'][:], group['indices'][:], group['indptr'][:]), shape=group.attrs['shape'])

def digest(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

paths = {r: (HERE.parent / 'result' / ('young2_HSPC_Multi_'+r+'_with_cellinfo') / 'cell2fate.h5ad'
             if method == 'cell2fate' else HERE / 'results' / r / 'topicvelo.h5ad') for r in RATIOS}
metadata = {}
for r,p in paths.items():
    with h5py.File(p, 'r') as f:
        obs, var = ids(f,'obs'), ids(f,'var')
        metadata[r] = {'obs': obs, 'var': var}
cells_set = set.intersection(*(set(x['obs']) for x in metadata.values()))
genes_set = set.intersection(*(set(x['var']) for x in metadata.values()))
cells = [x for x in metadata['0_2']['obs'] if x in cells_set]
genes = [x for x in metadata['0_2']['var'] if x in genes_set]
assert len(cells) >= 50 and len(genes) >= 10
pd.Series(cells, name='cell_id').to_csv(out/'common_cells.csv', index=False)
pd.Series(genes, name='gene').to_csv(out/'common_genes.csv', index=False)
mats, transitions, sources = {}, {}, {}
if method == 'topicvelo':
    # Keeping all cells avoids silently renormalizing a truncated transition.
    assert all(len(cells) == len(x['obs']) for x in metadata.values())
    with h5py.File(paths['1_0'], 'r') as f:
        ci = pd.Index(metadata['1_0']['obs']).get_indexer(cells)
        gi = pd.Index(metadata['1_0']['var']).get_indexer(genes)
        X = matrix(f['layers/Ms'])
        X = X[ci][:,gi] if sparse.issparse(X) else X[np.ix_(ci,gi)]
        X = X.toarray() if sparse.issparse(X) else X
        X = np.asarray(X, dtype=np.float64)
        assert np.isfinite(X).all()
for r,p in paths.items():
    with h5py.File(p, 'r') as f:
        ci = pd.Index(metadata[r]['obs']).get_indexer(cells)
        gi = pd.Index(metadata[r]['var']).get_indexer(genes)
        if method == 'cell2fate':
            raw = matrix(f['layers/velocity'])
            raw = raw.toarray() if sparse.issparse(raw) else raw
            source_hash = digest(raw)
            V = np.asarray(raw[np.ix_(ci,gi)], dtype=np.float64)
        else:
            T = matrix(f['obsp/topic_modeling+scvelo_stochastic_T']).tocsr()[ci][:,ci]
            assert np.isfinite(T.data).all() and (T.data >= 0).all()
            assert np.allclose(np.asarray(T.sum(1)).ravel(),1,atol=2e-6)
            transitions[r] = T
            source_hash = digest(T.data) + ':' + digest(T.indices) + ':' + digest(T.indptr)
            V = T @ X - np.asarray(T.sum(1)) * X
        mats[r] = V
        sources[r] = {'path':str(p), 'bytes':p.stat().st_size, 'mtime_ns':p.stat().st_mtime_ns,
                      'evaluated_source_array_sha256':source_hash, 'aligned_vector_sha256':digest(V),
                      'source_cells':len(metadata[r]['obs']), 'source_genes':len(metadata[r]['var'])}

def cosine_rows(A,B):
    finite = np.isfinite(A) & np.isfinite(B)
    a = np.where(finite,A,0); b = np.where(finite,B,0)
    na = np.linalg.norm(a,axis=1); nb = np.linalg.norm(b,axis=1)
    valid = (finite.sum(1)>=10) & (na>0) & (nb>0)
    scores = np.full(len(A),np.nan)
    scores[valid] = np.clip(np.einsum('ij,ij->i',a[valid],b[valid]) / (na[valid]*nb[valid]),-1,1)
    status = np.where(valid,'valid',np.where(finite.sum(1)<10,'insufficient_finite_genes','zero_vector'))
    return scores,status

rows, per_cell, native = [], [], []
display = 'cell2fate' if method == 'cell2fate' else 'TopicVelo'
for ri in RATIOS:
    for rj in RATIOS:
        score,status = cosine_rows(mats[ri],mats[rj])
        i,j = float(ri.replace('_','.')),float(rj.replace('_','.'))
        rows.append({'Method':display,'Ratio_i':i,'Ratio_j':j,'CosineSim':float(np.nanmean(score))})
        per_cell.append(pd.DataFrame({'cell_id':cells,'Ratio_i':i,'Ratio_j':j,'cosine':score,'status':status}))
        if method == 'topicvelo':
            a,b=transitions[ri],transitions[rj]
            denom=np.sqrt(np.asarray(a.multiply(a).sum(1)).ravel()*np.asarray(b.multiply(b).sum(1)).ravel())
            z=np.asarray(a.multiply(b).sum(1)).ravel()/denom
            native.append({'Method':display,'Ratio_i':i,'Ratio_j':j,'CosineSim':float(np.mean(z))})
table=pd.DataFrame(rows)
table.to_csv(out/'pairwise_cosine.csv',index=False)
pd.concat(per_cell,ignore_index=True).to_csv(out/'per_cell_cosine.csv.gz',index=False,compression='gzip')
if native: pd.DataFrame(native).to_csv(out/'native_transition_row_cosine.csv',index=False)
summary={'method':display,'status':'complete','common_cells':len(cells),'common_genes':len(genes),
         'metric': 'gene_velocity_cosine' if method=='cell2fate' else 'integrated_transition_displacement_cosine_in_fixed_full_depth_Ms',
         'not_base_scvelo_velocity':method=='topicvelo','sources':sources,
         'mean_stability':float(table.loc[(table.Ratio_j==1)&(table.Ratio_i!=1),'CosineSim'].mean()),
         'vs_full':table.loc[table.Ratio_j==1].to_dict(orient='records'),
         'valid_cells_by_comparison':{str(i):int((f.status=='valid').sum()) for i,f in enumerate(per_cell)}}
(out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='sources'},indent=2))
