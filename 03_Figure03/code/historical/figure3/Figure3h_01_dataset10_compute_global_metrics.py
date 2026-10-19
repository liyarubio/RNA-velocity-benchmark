#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VeloAE 五指标 + Global 综合指标计算
=====================================
基于 compute_veloae_metrics.py 的五个指标逻辑（修复后的正确版本），
并额外计算 计算指标_中位数.py 中的 get_global 复合指标（已修复除零隐患）。

计算指标:
  1. CBDir      — Cross-Boundary Direction (transition probability)
  2. CBVelCoh   — Cross-Boundary Velocity Coherence
  3. CBDirCorr  — Cross-Boundary Direction Correctness (embedding)
  4. ICCoh      — In-cluster Coherence
  5. ICConf     — In-cluster Velocity Confidence
  6. Global     — 复合指标（基于 ICCoh + CBDirCorr 的 raw 分数）

数据加载方式沿用 计算指标_中位数.py:
    ../real_data_graph/{dataset}/{method}.h5ad

保存方式沿用 计算指标_中位数.py 的宽表 CSV（行=方法, 列=数据集）:
    veloAE指标_global/{Metric}_all.csv   每个指标一个宽表
    veloAE指标_global/metrics_all_long.csv  完整长表
"""

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
# 线程数由 --n_jobs 参数控制，默认使用一半核数（多线程加速 scvelo 图构建）
N_JOBS = max(1, (os.cpu_count() or 2) // 2)

# ==================== 数据集配置（与 计算指标_中位数.py 一致） ====================
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

# ===== 所有方法及其 velocity key（与 计算指标_中位数.py 一致） =====
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
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'veloAE指标_global')
os.makedirs(OUT_DIR, exist_ok=True)
_ERROR_LOG = os.path.join(OUT_DIR, 'error_log.txt')


def log_error(msg):
    """记录错误到日志文件"""
    with open(_ERROR_LOG, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')
    print(f"    [错误] {msg}")


# ==================== 工具函数 ====================

def _to_array(x):
    if hasattr(x, 'toarray'):
        return x.toarray()
    if hasattr(x, 'A'):
        return x.A
    return np.asarray(x, dtype=np.float32)


def keep_type(adata, nodes, target, k_cluster):
    """选择目标类型细胞"""
    vals = adata.obs[k_cluster].astype(str).values
    node_vals = vals[nodes]
    mask = node_vals == str(target)
    if mask.sum() > 0:
        return nodes[mask]
    # 大小写容错
    tl = str(target).lower()
    return nodes[np.array([tl in v.lower() for v in node_vals])]


def get_nb(adata, sel_idx=None):
    """获取邻居索引（兼容 indices / connectivities）"""
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


# ==================== 五个指标 ====================

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
                ts.append(float(np.mean(pa[bn])))  # bn 是绝对细胞索引
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

    # 找 velocity embedding
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
    """ICConf: In-cluster Velocity Confidence（基于 scvelo velocity_confidence）"""
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


# ==================== Global 复合指标（已修复除零隐患） ====================

def get_global(in_cluster, cross_boundary, cell_number):
    """计算 Global 复合指标。

    Args:
        in_cluster (dict): ICCoh 的 raw 分数 {cluster: [scores]}。
        cross_boundary (dict): CBDirCorr 的 raw 分数 {edge: [scores]}。
        cell_number (int): 细胞总数。

    Returns:
        float: Global 分数。
    """
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

        # 该 cluster 不参与任何 transition edge，方向性无定义，跳过避免除零
        if count_num == 0:
            continue

        g1.append(p1 * h1 * d1 / count_num)

    g1 = [item for item in g1 if not np.isnan(item)]
    return float(np.sum(g1))


# ==================== 数据处理辅助 ====================

def prepare_adata(adata, file_name, k_cluster):
    """数据集特化处理，与 计算指标_中位数.py 保持一致"""
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
    """确保有 neighbors 和 UMAP/TSNE 嵌入"""
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
        get_nb(adata)  # 自动构建 indices


# ==================== 单方法处理 ====================

def process_one_method(adata, method, vkey, file_name, k_cluster, edge, emb, x_emb):
    """处理单个方法的五指标 + Global 计算，返回 results dict"""
    results = {'Dataset': file_name, 'Model': method}

    # 空值清洗 velocity 层
    arr = _to_array(adata.layers[vkey])
    all_nan_rows = np.where(np.isnan(arr).all(axis=1))[0]
    all_nan_cols = np.where(np.isnan(arr).all(axis=0))[0]
    has_row_filter = len(all_nan_rows) > 0
    has_col_filter = len(all_nan_cols) > 0

    if has_row_filter or has_col_filter:
        print(f"    过滤全NaN: {len(all_nan_rows)}行, {len(all_nan_cols)}列")
        keep_rows = np.ones(arr.shape[0], dtype=bool)
        keep_rows[all_nan_rows] = False
        keep_cols = np.ones(arr.shape[1], dtype=bool)
        keep_cols[all_nan_cols] = False

        # --- 行过滤：重新映射 neighbor indices，保留拓扑 ---
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
            # 同样重排 connectivities
            for obsp_key in list(adata.obsp.keys()):
                adata.obsp[obsp_key] = adata.obsp[obsp_key][keep_rows][:, keep_rows]

        adata = adata[keep_rows, keep_cols].copy()
    adata.layers[vkey] = np.nan_to_num(adata.layers[vkey], nan=0.0)

    # 邻居：只有被行过滤销毁时才重建
    if has_row_filter and ('neighbors' not in adata.uns or 'indices' not in adata.uns['neighbors']):
        try:
            ensure_neighbors(adata, emb)
        except Exception:
            pass

    # velocity_graph：始终在新基因空间重建（基因维度变了）
    try:
        scv.tl.velocity_graph(adata, vkey=vkey, n_jobs=N_JOBS)
        scv.tl.velocity_confidence(adata, vkey=vkey)
        scv.tl.velocity_embedding(adata, vkey=vkey, basis=emb)
    except Exception as e:
        log_error(f"{file_name}/{method}: 图重建失败: {e}")

    # 显式计算五个指标，保留 raw 分数供 Global 使用
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
            log_error(f"{file_name}/{method}: {mlabel} 失败: {e}")

    # 计算 Global 复合指标（基于 ICCoh 和 CBDirCorr 的 raw 分数）
    try:
        results['Global'] = get_global(
            raw_scores.get('ICCoh', {}),
            raw_scores.get('CBDirCorr', {}),
            adata.n_obs,
        )
    except Exception as e:
        results['Global'] = np.nan
        log_error(f"{file_name}/{method}: Global 失败: {e}")

    return results


# ==================== 主流程 ====================

def main():
    global N_JOBS
    parser = argparse.ArgumentParser(description='VeloAE五指标+Global计算')
    parser.add_argument('--n_jobs', type=int, default=N_JOBS,
                        help='CPU 核数（默认一半核，-1 表示全部核，1 表示单核）')
    parser.add_argument('--datasets', type=str, default='',
                        help='逗号分隔的数据集列表，为空则运行所有')
    parser.add_argument('--node', type=int, default=0,
                        help='节点编号（1/2/3），用于区分输出文件，避免覆盖')
    args = parser.parse_args()

    N_JOBS = args.n_jobs
    scv.settings.n_jobs = args.n_jobs

    node_suffix = f'_node{args.node}' if args.node else ''

    global _ERROR_LOG
    _ERROR_LOG = os.path.join(OUT_DIR, f'error_log{node_suffix}.txt')

    # 确定要处理的数据集
    if args.datasets:
        target_datasets = [d.strip() for d in args.datasets.split(',') if d.strip() in cell_types]
    else:
        target_datasets = sorted(cell_types.keys())

    # 清空日志文件
    with open(_ERROR_LOG, 'w', encoding='utf-8') as f:
        f.write('Error log for velocity evaluation (5 metrics + Global)\n')
        f.write('=' * 60 + '\n')

    print("=" * 60)
    print("VeloAE 五指标 + Global 综合指标计算")
    print(f"方法数: {len(ALL_METHODS)}  数据集数: {len(target_datasets)}")
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
                log_error(f"{file_name}/{method}: 文件不存在: {file_path}")
                all_rows.append(nan_row)
                continue

            try:
                adata = sc.read_h5ad(file_path)
            except Exception as e:
                log_error(f"{file_name}/{method}: 读取失败: {e}")
                all_rows.append(nan_row)
                continue

            if vkey not in adata.layers:
                log_error(f"{file_name}/{method}: layers 中无 {vkey}")
                all_rows.append(nan_row)
                continue

            # 数据特化处理
            try:
                adata = prepare_adata(adata, file_name, k_cluster)
            except Exception as e:
                log_error(f"{file_name}/{method}: prepare_adata 失败: {e}")

            # 邻居和嵌入
            try:
                ensure_neighbors(adata, emb)
            except Exception as e:
                log_error(f"{file_name}/{method}: neighbors/umap 失败: {e}")

            # TFvelo 特殊处理：删除全 NaN 列
            if method == 'TFvelo':
                try:
                    arr = _to_array(adata.layers[vkey])
                    all_nan_cols = np.where(np.isnan(arr).all(axis=0))[0]
                    if len(all_nan_cols) > 0:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_cols))
                except Exception:
                    pass

            # scvelo_dynamical 特殊处理：清理 velocity embedding NaN
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
                log_error(f"{file_name}/{method}: 处理异常: {e}")
                traceback.print_exc()
                all_rows.append(nan_row)

    # ========== 保存结果 ==========
    df_result = pd.DataFrame(all_rows)

    # 完整长表
    fp_all = os.path.join(OUT_DIR, 'metrics_all_long.csv')
    df_result.to_csv(fp_all, index=False, float_format='%.6f', na_rep='NaN')
    print(f"\n  完整长表 -> {fp_all}")

    # 每个指标一个宽表（行=方法, 列=数据集），沿用 计算指标_中位数.py 的保存方式
    for metric_key in METRICS:
        sub = df_result.pivot(index='Model', columns='Dataset', values=metric_key)
        fp = os.path.join(OUT_DIR, f'{metric_key}_all{node_suffix}.csv')
        sub.to_csv(fp, float_format='%.6f', na_rep='NaN')
        print(f"  {metric_key} 宽表 -> {fp}")

    # 打印摘要
    print(f"\n{'=' * 50}")
    print("摘要")
    for mkey in METRICS:
        valid = df_result[mkey].dropna()
        if len(valid) > 0:
            print(f"  {mkey}: mean={valid.mean():.4f}, n={len(valid)}/{len(df_result)}")
        else:
            print(f"  {mkey}: 全部 NaN")

    print(f"\n完成! 结果保存在 {OUT_DIR}")


if __name__ == '__main__':
    main()
