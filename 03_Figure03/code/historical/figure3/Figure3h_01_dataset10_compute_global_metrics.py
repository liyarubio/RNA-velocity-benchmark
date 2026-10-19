#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'VeloAE: five metrics and the Global composite score\n====================================================\nUses the corrected five-metric implementation from compute_veloae_metrics.py\nand the get_global composite from the historical median-metrics script, with division-by-zero protection.\n\nMetrics:\n  1. CBDir      - Cross-Boundary Direction (transition probability)\n  2. CBVelCoh   - Cross-Boundary Velocity Coherence\n  3. CBDirCorr  - Cross-Boundary Direction Correctness (embedding)\n  4. ICCoh      - In-cluster Coherence\n  5. ICConf     - In-cluster Velocity Confidence\n  6. Global     - Composite based on raw ICCoh and CBDirCorr scores\n\nLoad data as in the historical median-metrics script:\n    ../real_data_graph/{dataset}/{method}.h5ad\n\nSave wide CSV tables (rows = methods, columns = datasets):\n    veloAE_metrics_global/{Metric}_all.csv      One wide table per metric\n    veloAE_metrics_global/metrics_all_long.csv  Complete long-format table\n'

import os
import sys
import argparse
import warnings
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from sklearn.metrics.pairwise import cosine_similarity

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0
scv.settings.verbosity = 0
# --n_jobs controls CPU threads; the default uses half the cores for scvelo graph construction.
N_JOBS = max(1, (os.cpu_count() or 2) // 2)

# ==================== Dataset configuration, matching the historical median-metrics script ====================
cell_types = {
    'endocrinogenesis_day15': 'clusters',
    'DentateGyrus_velocyto': 'clusters',
    'human_cd34_bone_marrow': 'clusters',
    'Forebrain': 'clusters',
    'dentategyrus_scv': 'clusters',
    'Hindbrain_GABA_Glio': 'Celltype',
    'organogenesis_chondrocyte': 'Main_cell_type',
    'erythroid_lineage': 'celltype',
    'zebrafish_process': 'Cell_type',
    'gastrulation': 'stage',
}

edges = {
    'dentategyrus_scv': [("Neuroblast", "Granule immature"),
                         ("Granule immature", "Granule mature"),
                         ("Radial Glia-like", "Astrocytes"),
                         ("OPC", "OL")],
    'DentateGyrus_velocyto': [('Radial Glia', 'Astrocyte'),
                              ("Neuroblast", 'Immature Granule'),
                              ('Immature Granule', 'Granule'),
                              ("Neuroblast", 'CA')],
    'endocrinogenesis_day15': [('Ngn3 high EP', 'Pre-endocrine'),
                               ('Pre-endocrine', 'Delta'),
                               ('Pre-endocrine', 'Beta'),
                               ('Pre-endocrine', 'Epsilon'),
                               ('Pre-endocrine', 'Alpha')],
    'human_cd34_bone_marrow': [("HSC", 'Precursors'),
                               ("HSC", 'Ery'),
                               ("HSC", 'Mono')],
    'gastrulation': [("E6.5", 'E6.75'),
                     ('E6.75', 'E7.0'),
                     ('E7.0', 'E7.25'),
                     ('E7.25', 'E7.5'),
                     ('E7.5', 'E7.75'),
                     ('E7.75', 'E8.0'),
                     ('E8.0', 'E8.25'),
                     ('E8.25', 'E8.5')],
    'erythroid_lineage': [('Blood progenitors', 'Erythroid1'),
                          ('Erythroid1', 'Erythroid2'),
                          ('Erythroid2', 'Erythroid3')],
    'Forebrain': [('Radial Glia', 'Neuroblast'),
                  ('Neuroblast', 'Immature Neuron'),
                  ('Immature Neuron', 'Neuron')],
    'Hindbrain_GABA_Glio': [('Neural stem cells', 'Proliferating VZ progenitors'),
                            ('Proliferating VZ progenitors', 'VZ progenitors'),
                            ('VZ progenitors', 'Differentiating GABA interneurons'),
                            ('Differentiating GABA interneurons', 'GABA interneurons'),
                            ('VZ progenitors', 'Gliogenic progenitors')],
    'organogenesis_chondrocyte': [('Early mesenchyme', 'Chondroctye progenitors'),
                                  ('Chondroctye progenitors', 'Chondrocytes & osteoblasts'),
                                  ('Early mesenchyme', 'Jaw and tooth progenitors')],
    'zebrafish_process': [('Pigment Progenitor', 'Iridophore'),
                          ('Pigment Progenitor', 'Melanophore'),
                          ('Schwann Cell Precursor', 'Schwann Cell')],
}

embs = {
    'dentategyrus_scv': 'umap',
    'DentateGyrus_velocyto': 'tsne',
    'endocrinogenesis_day15': 'umap',
    'human_cd34_bone_marrow': 'tsne',
    'Forebrain': 'umap',
    'gastrulation': 'umap',
    'erythroid_lineage': 'umap',
    'Hindbrain_GABA_Glio': 'tsne',
    'organogenesis_chondrocyte': 'umap',
    'zebrafish_process': 'umap',
}

# ===== Methods and velocity keys, matching the historical median-metrics script =====
ALL_METHODS = {
    'DeepVelo_GB': 'velocity',
    'scvelo_stochastic': 'velocity',
    'Dynamo_stochastic': 'velocity_S',
    'veloVI': 'velocity',
    'scvelo_dynamical': 'velocity',
    'VeloAE': 'new_velocity',
    'UniTvelo': 'velocity',
    'Celldancer': 'cd_velocity_s',
    'DeepVelo_SA': 'DSA_velocity',
    'velocyto': 'velocity',
    'TFvelo': 'velocity',
    'STT': 'velocity',
    'veloVAE': 'vae_velocity',
    'LatentVelo': 'velo_s',
    'cell2fate': 'velocity',
    'DeepKINET': 'velocity',
    'biVI': 'velocity',
    'SDEvelo': 'sde_velocity',
    'TopicVelo': 'velocity',
}

METRICS = ['CBDir', 'CBVelCoh', 'CBDirCorr', 'ICCoh', 'ICConf', 'Global']

DATA_PATH = '../real_data_graph/'
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'veloAE_metrics_global')
os.makedirs(OUT_DIR, exist_ok=True)
_ERROR_LOG = os.path.join(OUT_DIR, 'error_log.txt')


def log_error(msg):
    'Record an error in the log file.'
    with open(_ERROR_LOG, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')
    print(f"    [Error] {msg}")


# ==================== Utility functions ====================

def _to_array(x):
    if hasattr(x, 'toarray'):
        return x.toarray()
    if hasattr(x, 'A'):
        return x.A
    return np.asarray(x, dtype=np.float32)


def keep_type(adata, nodes, target, k_cluster):
    'Select cells of the target type.'
    vals = adata.obs[k_cluster].astype(str).values
    node_vals = vals[nodes]
    mask = node_vals == str(target)
    if mask.sum() > 0:
        return nodes[mask]
    # Match case-insensitively.
    tl = str(target).lower()
    return nodes[np.array([tl in v.lower() for v in node_vals])]


def get_nb(adata, sel_idx=None):
    'Get neighbor indices, supporting indices or connectivities.'
    if 'neighbors' in adata.uns and 'indices' in adata.uns['neighbors']:
        ind = adata.uns['neighbors']['indices']
        if isinstance(ind, list):
            result = ind
        else:
            result = [ind[i] for i in range(len(ind))]
        if sel_idx is not None:
            if hasattr(sel_idx, 'values'):
                sel_idx = sel_idx.values
            return [result[i] for i in range(len(result)) if sel_idx[i]]
        return result
    # Fallback: build from connectivities
    if 'connectivities' in adata.obsp:
        conn = adata.obsp['connectivities']
        nc = adata.n_obs
        inner = []
        for i in range(nc):
            row = conn[i].tocoo()
            nz, nd = row.col, row.data
            if len(nz) >= 30:
                top = np.argsort(-np.abs(nd))[:30]
                idx = nz[top]
            elif len(nz) > 0:
                idx = nz
            else:
                idx = np.array([i])
            inner.append(idx.astype(np.int64))
        if 'neighbors' not in adata.uns:
            adata.uns['neighbors'] = {}
        adata.uns['neighbors']['indices'] = inner
        if sel_idx is not None:
            if hasattr(sel_idx, 'values'):
                sel_idx = sel_idx.values
            return [inner[i] for i in range(nc) if sel_idx[i]]
        return inner
    raise KeyError('no neighbors')


# ==================== Five metrics ====================

def cross_boundary_transition(adata, k_cluster, cluster_edges, k_velocity):
    """CBDir: Cross-Boundary Direction (transition probability)"""
    kg = f'{k_velocity}_graph'
    if kg not in adata.uns:
        return {}, np.nan
    scores, all_scores = {}, {}
    for u, v in cluster_edges:
        us, vs = str(u), str(v)
        sel = adata.obs[k_cluster].astype(str) == us
        if not sel.any():
            continue
        nbs = get_nb(adata, sel.values)
        tg = adata.uns[kg][sel.values]
        ts = []
        for tp, nodes in zip(tg, nbs):
            bn = keep_type(adata, nodes, vs, k_cluster)
            if len(bn) > 0:
                pa = tp.toarray().flatten() if hasattr(tp, 'toarray') else np.asarray(tp).flatten()
                ts.append(float(np.mean(pa[bn])))  # bn contains absolute cell indices.
        if ts:
            scores[(us, vs)] = np.mean(ts)
            all_scores[(us, vs)] = ts
    ov = np.mean([s for s in scores.values()]) if scores else np.nan
    return all_scores, ov


def cross_boundary_coh(adata, k_cluster, k_velocity, cluster_edges):
    """CBVelCoh: Cross-Boundary Velocity Coherence"""
    if k_velocity not in adata.layers:
        return {}, np.nan
    velocities = _to_array(adata.layers[k_velocity])
    scores, all_scores = {}, {}
    for u, v in cluster_edges:
        us, vs = str(u), str(v)
        sel = adata.obs[k_cluster].astype(str) == us
        if not sel.any():
            continue
        si = np.where(sel.values)[0]
        nbs = get_nb(adata, sel.values)
        vus = velocities[si]
        ts = []
        for ith, nodes in enumerate(nbs):
            bn = keep_type(adata, nodes, vs, k_cluster)
            if len(bn) > 0:
                cs = cosine_similarity(vus[[ith]], velocities[bn]).mean()
                ts.append(float(cs))
        if ts:
            scores[(us, vs)] = np.mean(ts)
            all_scores[(us, vs)] = ts
    ov = np.mean([s for s in scores.values()]) if scores else np.nan
    return all_scores, ov


def cross_boundary_correctness(adata, k_cluster, k_velocity, cluster_edges, x_emb="X_umap"):
    """CBDirCorr: Cross-Boundary Direction Correctness (embedding)"""
    if k_velocity not in adata.layers or x_emb not in adata.obsm:
        return {}, np.nan
    velocities = _to_array(adata.layers[k_velocity])
    xec = adata.obsm[x_emb]

    # Find the velocity embedding.
    vek = None
    for key in adata.obsm.keys():
        if key.startswith(k_velocity) and key != k_velocity:
            vek = key
            break
    if vek is None:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        v_emb = pca.fit_transform(velocities)
    else:
        v_emb = adata.obsm[vek]

    scores, all_scores = {}, {}
    for u, v in cluster_edges:
        us, vs = str(u), str(v)
        sel = adata.obs[k_cluster].astype(str) == us
        if not sel.any():
            continue
        si = np.where(sel.values)[0]
        nbs = get_nb(adata, sel.values)
        xp, xv = xec[si], v_emb[si]
        ts = []
        for ith, nodes in enumerate(nbs):
            bn = keep_type(adata, nodes, vs, k_cluster)
            if len(bn) == 0:
                continue
            pd_pos = xec[bn] - xp[ith]
            ds = cosine_similarity(pd_pos, xv[ith].reshape(1, -1)).flatten()
            ts.append(float(np.mean(ds)))
        if ts:
            scores[(us, vs)] = np.mean(ts)
            all_scores[(us, vs)] = ts
    ov = np.mean([s for s in scores.values()]) if scores else np.nan
    return all_scores, ov


def inner_cluster_coh(adata, k_cluster, k_velocity):
    """ICCoh: In-cluster Coherence"""
    if k_velocity not in adata.layers:
        return {}, np.nan
    velocities = _to_array(adata.layers[k_velocity])
    clusters = np.unique(adata.obs[k_cluster].astype(str))
    scores, all_scores = {}, {}
    for cat in clusters:
        sel = adata.obs[k_cluster].astype(str) == cat
        if sel.sum() < 3:
            continue
        si = np.where(sel.values)[0]
        nbs = get_nb(adata, sel.values)
        cvs = velocities[si]
        ts = []
        for ith, nodes in enumerate(nbs):
            sc2 = keep_type(adata, nodes, cat, k_cluster)
            if len(sc2) > 0:
                cs = cosine_similarity(cvs[[ith]], velocities[sc2]).mean()
                ts.append(float(cs))
        if ts:
            all_scores[cat] = ts
            scores[cat] = np.mean(ts)
    ov = np.mean([s for s in scores.values()]) if scores else np.nan
    return all_scores, ov


def in_cluster_confidence(adata, k_cluster, k_velocity):
    'ICConf: In-cluster Velocity Confidence, based on scvelo velocity_confidence.'
    conf_key = f'{k_velocity}_confidence'
    if conf_key not in adata.obs.columns:
        return {}, np.nan
    clusters = np.unique(adata.obs[k_cluster].astype(str))
    scores, all_scores = {}, {}
    for cat in clusters:
        sel = adata.obs[k_cluster].astype(str) == cat
        type_score = adata.obs[conf_key][sel].values.tolist()
        scores[cat] = np.mean(type_score)
        all_scores[cat] = type_score
    ov = np.mean([s for s in scores.values()]) if scores else np.nan
    return scores, ov


# ==================== Global composite score with division-by-zero protection ====================

def get_global(in_cluster, cross_boundary, cell_number):
    'Calculate the Global composite score.\n\n    Args:\n        in_cluster (dict): Raw ICCoh scores {cluster: [scores]}.\n        cross_boundary (dict): Raw CBDirCorr scores {edge: [scores]}.\n        cell_number (int): Total number of cells.\n\n    Returns:\n        float: Global score.\n    '
    g1 = []
    for cluster0 in in_cluster:
        data_tmp1 = in_cluster[cluster0]
        p1 = len(data_tmp1) / cell_number
        h1 = np.mean(data_tmp1)

        d1 = 0
        count_num = 0
        for edge0 in cross_boundary:
            if cluster0 in edge0:
                ss = cross_boundary[edge0]
                count_num = count_num + 1
                d1 = d1 + np.sign(np.median(ss))

        # This cluster participates in no transition edge; direction is undefined, so skip it to avoid division by zero.
        if count_num == 0:
            continue

        g1.append(p1 * h1 * d1 / count_num)

    g1 = [item for item in g1 if not np.isnan(item)]
    return float(np.sum(g1))


# ==================== Data-processing helpers ====================

def prepare_adata(adata, file_name, k_cluster):
    'Apply dataset-specific processing as in the historical median-metrics script.'
    if file_name == 'gastrulation':
        adata.obs['stage'] = pd.Categorical(
            adata.obs['stage'], ordered=True,
            categories=['E6.5', 'E6.75', 'E7.0', 'E7.25', 'E7.5', 'E7.75',
                        'E8.0', 'E8.25', 'E8.5'])
    if file_name == 'DentateGyrus_velocyto':
        replace_dict = {'RadialGlia': 'Radial Glia', 'RadialGlia2': 'Radial Glia',
                        'ImmAstro': 'Astrocyte',
                        'ImmGranule1': 'Immature Granule', 'ImmGranule2': 'Immature Granule',
                        'Nbl1': 'Neuroblast', 'Nbl2': 'Neuroblast'}
        adata.obs.replace(replace_dict, inplace=True)
    if file_name == 'human_cd34_bone_marrow':
        replace_dict = {'HSC_1': 'HSC', 'HSC_2': 'HSC',
                        'Ery_1': 'Ery', 'Ery_2': 'Ery',
                        'Mono_1': 'Mono', 'Mono_2': 'Mono'}
        adata.obs.replace(replace_dict, inplace=True)
        adata.obs[k_cluster] = pd.Categorical(
            adata.obs[k_cluster], ordered=True,
            categories=['HSC', 'Ery', 'Mono', 'Precursors', 'CLP', 'DCs', 'Mega'])
    if file_name == 'erythroid_lineage':
        replace_dict = {'Blood progenitors 1': 'Blood progenitors',
                        'Blood progenitors 2': 'Blood progenitors'}
        adata.obs.replace(replace_dict, inplace=True)
    return adata


def ensure_neighbors(adata, emb_key):
    'Ensure neighbors and UMAP/TSNE embeddings are available.'
    emb_full = f'X_{emb_key}'
    if emb_full not in adata.obsm:
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata, n_comps=min(30, adata.n_obs - 1))
        scv.pp.neighbors(adata, n_neighbors=min(30, adata.n_obs - 1))
        if emb_key == 'umap':
            sc.tl.umap(adata)
        elif emb_key == 'tsne':
            sc.tl.tsne(adata)

    if 'neighbors' not in adata.uns or 'indices' not in adata.uns.get('neighbors', {}):
        get_nb(adata)  # Build indices automatically.


# ==================== Per-method processing ====================

def process_one_method(adata, method, vkey, file_name, k_cluster, edge, emb, x_emb):
    'Calculate five metrics and Global for one method; return a results dictionary.'
    results = {'Dataset': file_name, 'Model': method}

    # Clean missing values in the velocity layer.
    arr = _to_array(adata.layers[vkey])
    all_nan_rows = np.where(np.isnan(arr).all(axis=1))[0]
    all_nan_cols = np.where(np.isnan(arr).all(axis=0))[0]
    has_row_filter = len(all_nan_rows) > 0
    has_col_filter = len(all_nan_cols) > 0

    if has_row_filter or has_col_filter:
        print(f"    Filter all-NaN values: {len(all_nan_rows)} rows, {len(all_nan_cols)} columns")
        keep_rows = np.ones(arr.shape[0], dtype=bool)
        keep_rows[all_nan_rows] = False
        keep_cols = np.ones(arr.shape[1], dtype=bool)
        keep_cols[all_nan_cols] = False

        # --- Row filtering: remap neighbor indices while preserving topology. ---
        if has_row_filter and 'neighbors' in adata.uns and 'indices' in adata.uns['neighbors']:
            old_indices = adata.uns['neighbors']['indices']
            old_to_new = np.full(adata.n_obs, -1, dtype=int)
            new_ids = np.where(keep_rows)[0]
            old_to_new[new_ids] = np.arange(len(new_ids))
            new_indices = []
            for i in new_ids:
                old_nbs = old_indices[i]
                valid = keep_rows[old_nbs]
                new_nbs = old_to_new[old_nbs[valid]]
                new_indices.append(new_nbs.astype(np.int64))
            adata.uns['neighbors']['indices'] = new_indices
            # Reorder connectivities consistently.
            for obsp_key in list(adata.obsp.keys()):
                adata.obsp[obsp_key] = adata.obsp[obsp_key][keep_rows][:, keep_rows]

        adata = adata[keep_rows, keep_cols].copy()
    adata.layers[vkey] = np.nan_to_num(adata.layers[vkey], nan=0.0)

    # Neighbors: rebuild only when row filtering invalidates the graph.
    if has_row_filter and ('neighbors' not in adata.uns or 'indices' not in adata.uns['neighbors']):
        try:
            ensure_neighbors(adata, emb)
        except Exception:
            pass

    # velocity_graph: always rebuild in the updated gene space after gene filtering.
    try:
        scv.tl.velocity_graph(adata, vkey=vkey, n_jobs=N_JOBS)
        scv.tl.velocity_confidence(adata, vkey=vkey)
        scv.tl.velocity_embedding(adata, vkey=vkey, basis=emb)
    except Exception as e:
        log_error(f"{file_name}/{method}: Graph reconstruction failed: {e}")

    # Calculate all five metrics explicitly, retaining raw scores for Global.
    metric_funcs = {
        'CBDir': lambda: cross_boundary_transition(adata, k_cluster, edge, vkey),
        'CBVelCoh': lambda: cross_boundary_coh(adata, k_cluster, vkey, edge),
        'CBDirCorr': lambda: cross_boundary_correctness(adata, k_cluster, vkey, edge, x_emb),
        'ICCoh': lambda: inner_cluster_coh(adata, k_cluster, vkey),
        'ICConf': lambda: in_cluster_confidence(adata, k_cluster, vkey),
    }

    raw_scores = {}
    for mlabel, mfunc in metric_funcs.items():
        try:
            all_scores, ov = mfunc()
            raw_scores[mlabel] = all_scores
            results[mlabel] = ov if not np.isnan(ov) else np.nan
        except Exception as e:
            raw_scores[mlabel] = {}
            results[mlabel] = np.nan
            log_error(f"{file_name}/{method}: {mlabel} failed: {e}")

    # Calculate Global using raw ICCoh and CBDirCorr scores.
    try:
        results['Global'] = get_global(
            raw_scores.get('ICCoh', {}),
            raw_scores.get('CBDirCorr', {}),
            adata.n_obs,
        )
    except Exception as e:
        results['Global'] = np.nan
        log_error(f"{file_name}/{method}: Global failed: {e}")

    return results


# ==================== Main workflow ====================

def main():
    global N_JOBS
    parser = argparse.ArgumentParser(description='Calculate five VeloAE metrics and Global')
    parser.add_argument('--n_jobs', type=int, default=N_JOBS,
                        help='CPU cores: half by default; -1 uses all cores and 1 uses one core')
    parser.add_argument('--datasets', type=str, default='',
                        help='Comma-separated dataset names; an empty value runs all datasets')
    parser.add_argument('--node', type=int, default=0,
                        help='Node identifier (1/2/3), used to distinguish output files and avoid overwriting')
    args = parser.parse_args()

    N_JOBS = args.n_jobs
    scv.settings.n_jobs = args.n_jobs

    node_suffix = f'_node{args.node}' if args.node else ''

    global _ERROR_LOG
    _ERROR_LOG = os.path.join(OUT_DIR, f'error_log{node_suffix}.txt')

    # Select datasets to process.
    if args.datasets:
        target_datasets = [d.strip() for d in args.datasets.split(',') if d.strip() in cell_types]
    else:
        target_datasets = sorted(cell_types.keys())

    # Clear the log file.
    with open(_ERROR_LOG, 'w', encoding='utf-8') as f:
        f.write('Error log for velocity evaluation (5 metrics + Global)\n')
        f.write('=' * 60 + '\n')

    print("=" * 60)
    print('VeloAE: five metrics and the Global composite score')
    print(f"Methods: {len(ALL_METHODS)}  Datasets: {len(target_datasets)}")
    print(f"Data path: {DATA_PATH}")
    print(f"Output: {OUT_DIR}")
    print(f"n_jobs: {args.n_jobs}")
    print("=" * 60)

    all_rows = []

    for file_name in sorted(target_datasets):
        k_cluster = cell_types[file_name]
        edge = edges[file_name]
        emb = embs[file_name]
        x_emb = f'X_{emb}'

        for method, vkey in ALL_METHODS.items():
            print(f"\n  {file_name} / {method}  vkey={vkey}")

            file_path = f'{DATA_PATH}{file_name}/{method}.h5ad'

            nan_row = {'Dataset': file_name, 'Model': method}
            for m in METRICS:
                nan_row[m] = np.nan

            if not os.path.exists(file_path):
                log_error(f"{file_name}/{method}: File not found: {file_path}")
                all_rows.append(nan_row)
                continue

            try:
                adata = sc.read_h5ad(file_path)
            except Exception as e:
                log_error(f"{file_name}/{method}: Read failed: {e}")
                all_rows.append(nan_row)
                continue

            if vkey not in adata.layers:
                log_error(f"{file_name}/{method}: Missing layer {vkey}")
                all_rows.append(nan_row)
                continue

            # Dataset-specific processing
            try:
                adata = prepare_adata(adata, file_name, k_cluster)
            except Exception as e:
                log_error(f"{file_name}/{method}: prepare_adata failed: {e}")

            # Neighbors and embeddings
            try:
                ensure_neighbors(adata, emb)
            except Exception as e:
                log_error(f"{file_name}/{method}: neighbors/umap failed: {e}")

            # TFvelo: remove all-NaN columns.
            if method == 'TFvelo':
                try:
                    arr = _to_array(adata.layers[vkey])
                    all_nan_cols = np.where(np.isnan(arr).all(axis=0))[0]
                    if len(all_nan_cols) > 0:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_cols))
                except Exception:
                    pass

            # scvelo_dynamical: clean NaN values in the velocity embedding.
            if method == 'scvelo_dynamical':
                try:
                    for k in list(adata.obsm.keys()):
                        if vkey in k:
                            adata.obsm[k] = np.nan_to_num(adata.obsm[k], nan=0)
                except Exception:
                    pass

            try:
                results = process_one_method(adata, method, vkey, file_name, k_cluster, edge, emb, x_emb)
                all_rows.append(results)
                parts = []
                for m in METRICS:
                    v = results.get(m, np.nan)
                    parts.append(f"{m}={v:.4f}" if not np.isnan(v) else f"{m}=NaN")
                print(f"    {'  '.join(parts)}")
            except Exception as e:
                log_error(f"{file_name}/{method}: Processing exception: {e}")
                traceback.print_exc()
                all_rows.append(nan_row)

    # ========== Save results ==========
    df_result = pd.DataFrame(all_rows)

    # Complete long-format table
    fp_all = os.path.join(OUT_DIR, 'metrics_all_long.csv')
    df_result.to_csv(fp_all, index=False, float_format='%.6f', na_rep='NaN')
    print(f"\n  Complete long-format table -> {fp_all}")

    # One wide table per metric (rows = methods, columns = datasets), as in the historical median-metrics script.
    for metric_key in METRICS:
        sub = df_result.pivot(index='Model', columns='Dataset', values=metric_key)
        fp = os.path.join(OUT_DIR, f'{metric_key}_all{node_suffix}.csv')
        sub.to_csv(fp, float_format='%.6f', na_rep='NaN')
        print(f"  {metric_key} Wide table -> {fp}")

    # Print a summary.
    print(f"\n{'=' * 50}")
    print('Summary')
    for mkey in METRICS:
        valid = df_result[mkey].dropna()
        if len(valid) > 0:
            print(f"  {mkey}: mean={valid.mean():.4f}, n={len(valid)}/{len(df_result)}")
        else:
            print(f"  {mkey}: All values are NaN")

    print(f"\nComplete! Results saved to {OUT_DIR}")


if __name__ == '__main__':
    main()
