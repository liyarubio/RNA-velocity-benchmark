"""
Compare velocity matrix similarity across reads downsampling ratios.

For each method, load the velocity layer from all 5 ratio directories and
compute pairwise mean per-cell cosine similarity between ratios.

Output:
  - Result/2_reads/velocity_sim_pairwise/  (5×5 heatmap per method)
  - Plotting has been split into separate files.
  - This script computes the historical long-format metrics only.
  - Result/2_reads/velocity_cosine_sim.csv (long-format CSV)
"""

import os
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# -------------------------------------------------------
# Config
# -------------------------------------------------------
base_dir = '/media/liyaru/LYR1/Benchmark/Revised/velogif_0315/result'
out_dir  = '/media/liyaru/LYR1/Benchmark/Revised/Result/2_reads'
os.makedirs(out_dir, exist_ok=True)

pairwise_dir = os.path.join(out_dir, 'velocity_sim_pairwise')
os.makedirs(pairwise_dir, exist_ok=True)

ratios = ['0.2', '0.4', '0.6', '0.8', '1.0']
ratio_dirs = {
    '0.2': 'young2_HSPC_Multi_0_2_with_cellinfo',
    '0.4': 'young2_HSPC_Multi_0_4_with_cellinfo',
    '0.6': 'young2_HSPC_Multi_0_6_with_cellinfo',
    '0.8': 'young2_HSPC_Multi_0_8_with_cellinfo',
    '1.0': 'young2_HSPC_Multi_1_0_with_cellinfo',
}

# Method file → display name, velocity layer key
# topicvelo and cell2fate are excluded (no valid velocity)
METHODS = {
    'velocyto':     ('velocyto',        'velocity'),
    'scvelo.sto':   ('scVelo(sto)',      'velocity'),
    'scvelo.dyn':   ('scVelo(dyn)',      'velocity'),
    'veloae':       ('VeloAE',           'velocity'),
    'dynamo':       ('dynamo',           'velocity'),
    'velovae':      ('VeloVAE',          'velocity'),
    'unitvelo':     ('UniTVelo',         'velocity'),
    'deepvelo_vae': ('DeepVelo(VAE)',     'velocity'),
    'celldancer':   ('cellDancer',       'velocity'),
    'velovi':       ('veloVI',           'velocity'),
    'latentvelo':   ('LatentVelo',       'velocity'),
    'deepvelo_gcn': ('DeepVelo(GCN)',    'velocity'),
    'tfvelo':       ('TFvelo',           'velocity'),
    'stt':          ('STT',              'velocity'),
    'multivelo':    ('MultiVelo',        'velocity'),
    'phylovelo':    ('PhyloVelo',        'velocity'),
    'bivi':         ('biVI',             'velocity'),
    'sdevelo':      ('SDEvelo',          'velocity'),
    'deepkinet':    ('DeepKINET',        'velocity'),
}

# Method display order (consistent with other figures)
METHOD_ORDER = [
    'velocyto', 'scVelo(sto)', 'scVelo(dyn)', 'VeloAE', 'dynamo',
    'VeloVAE', 'MultiVelo', 'UniTVelo', 'DeepVelo(VAE)', 'cellDancer',
    'PhyloVelo', 'veloVI', 'LatentVelo', 'DeepVelo(GCN)', 'TFvelo', 'STT',
    'biVI', 'SDEvelo', 'DeepKINET',
]

# -------------------------------------------------------
# Helper: load velocity matrix with cell/gene labels
# -------------------------------------------------------
def load_velocity(fpath, vel_key='velocity'):
    """Returns (obs_names, var_names, dense_array) or None if layer missing."""
    adata = ad.read_h5ad(fpath, backed='r')
    if vel_key not in adata.layers:
        return None
    v = adata.layers[vel_key]
    if sp.issparse(v):
        v = v.toarray().astype(float)
    else:
        v = np.array(v, dtype=float)
    return list(adata.obs_names), list(adata.var_names), v


def align_matrices(mats_dict):
    """
    mats_dict: {ratio -> (obs_names, var_names, array)}
    Returns: {ratio -> aligned 2D array} over common cells × common genes.
    """
    valid = {r: m for r, m in mats_dict.items() if m is not None}
    if len(valid) < 2:
        return {r: None for r in mats_dict}

    # Intersect cells and genes
    common_cells = set(valid[list(valid.keys())[0]][0])
    common_genes = set(valid[list(valid.keys())[0]][1])
    for obs, var, _ in valid.values():
        common_cells &= set(obs)
        common_genes &= set(var)

    if len(common_cells) < 50 or len(common_genes) < 5:
        print(f"    WARNING: too few common cells={len(common_cells)} genes={len(common_genes)}, skipping")
        return {r: None for r in mats_dict}

    # Stable ordering (use first ratio's order)
    first = list(valid.keys())[0]
    obs0, var0, _ = valid[first]
    cell_order = [c for c in obs0 if c in common_cells]
    gene_order = [g for g in var0 if g in common_genes]

    result = {}
    for r, m in mats_dict.items():
        if m is None:
            result[r] = None
            continue
        obs, var, arr = m
        obs_idx = {c: i for i, c in enumerate(obs)}
        var_idx = {g: i for i, g in enumerate(var)}
        ci = [obs_idx[c] for c in cell_order]
        gi = [var_idx[g] for g in gene_order]
        result[r] = arr[np.ix_(ci, gi)]
    return result, len(cell_order), len(gene_order)

# -------------------------------------------------------
# Helper: mean per-cell cosine similarity between two matrices
# Each row = one cell's velocity vector; genes with NaN in either are masked
# -------------------------------------------------------
def mean_cosine_sim(A, B):
    """
    A, B: (n_cells, n_genes), may contain NaN
    Returns scalar: mean per-cell cosine similarity (ignoring NaN genes per cell)
    """
    sims = []
    for a, b in zip(A, B):
        mask = np.isfinite(a) & np.isfinite(b)
        if mask.sum() < 10:
            continue
        a_ = a[mask]
        b_ = b[mask]
        na = np.linalg.norm(a_)
        nb = np.linalg.norm(b_)
        if na == 0 or nb == 0:
            continue
        sims.append(np.dot(a_, b_) / (na * nb))
    if len(sims) == 0:
        return np.nan
    return float(np.nanmean(sims))

# -------------------------------------------------------
# Main computation
# -------------------------------------------------------
records = []    # long-format: method, ratio_i, ratio_j, sim
pairwise = {}   # method_name -> 5×5 array

print("Loading velocity matrices and computing similarities...")
for file_stem, (display_name, vel_key) in METHODS.items():
    print(f"  Processing {display_name} ...", flush=True)

    # Load all 5 ratios (returns tuple or None)
    raw_mats = {}
    for r in ratios:
        fpath = os.path.join(base_dir, ratio_dirs[r], f'{file_stem}.h5ad')
        if not os.path.exists(fpath):
            print(f"    {r}: file missing")
            raw_mats[r] = None
            continue
        try:
            result = load_velocity(fpath, vel_key)
            raw_mats[r] = result  # (obs_names, var_names, array) or None
            if result is not None:
                print(f"    {r}: shape={result[2].shape}, nan_frac={np.isnan(result[2]).mean():.3f}")
            else:
                print(f"    {r}: velocity layer '{vel_key}' not found")
        except Exception as e:
            print(f"    {r}: ERROR - {e}")
            raw_mats[r] = None

    # Align cells and genes across ratios
    aligned = align_matrices(raw_mats)
    if isinstance(aligned, tuple):
        aligned, n_cells, n_genes = aligned
        print(f"    Aligned: {n_cells} common cells, {n_genes} common genes")
    else:
        aligned = aligned  # all None dict

    # Compute pairwise cosine similarity
    sim_mat = np.full((5, 5), np.nan)
    for i, ri in enumerate(ratios):
        for j, rj in enumerate(ratios):
            if aligned.get(ri) is None or aligned.get(rj) is None:
                continue
            if i == j:
                sim_mat[i, j] = 1.0
                continue
            sim_mat[i, j] = mean_cosine_sim(aligned[ri], aligned[rj])

    pairwise[display_name] = sim_mat

    # Record long-format
    for i, ri in enumerate(ratios):
        for j, rj in enumerate(ratios):
            records.append({
                'Method':   display_name,
                'Ratio_i':  float(ri),
                'Ratio_j':  float(rj),
                'CosineSim': sim_mat[i, j],
            })

df = pd.DataFrame(records)
df.to_csv(os.path.join(out_dir, 'velocity_cosine_sim.csv'), index=False)
print("Saved velocity_cosine_sim.csv")
