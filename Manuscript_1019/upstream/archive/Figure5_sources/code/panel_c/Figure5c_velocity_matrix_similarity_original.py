"""
Compare velocity matrix similarity across reads downsampling ratios.

For each method, load the velocity layer from all 5 ratio directories and
compute pairwise mean per-cell cosine similarity between ratios.

Output:
  - Result/2_reads/velocity_sim_pairwise/  (5×5 heatmap per method)
  - Result/2_reads/velocity_sim_to_ref.pdf (methods × ratio, similarity to 1.0x)
  - Result/2_reads/velocity_sim_pairwise_all.pdf (grid of 5×5 per method)
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

# -------------------------------------------------------
# Figure 1: similarity to 1.0x reference (methods × ratio)
# -------------------------------------------------------
ref_idx = ratios.index('1.0')
df_ref = df[df['Ratio_j'] == 1.0].copy()
df_ref = df_ref[df_ref['Ratio_i'] != 1.0]   # exclude self-comparison

pivot = df_ref.pivot(index='Method', columns='Ratio_i', values='CosineSim')
# Reorder rows
present = [m for m in METHOD_ORDER if m in pivot.index]
pivot = pivot.loc[present]

fig, ax = plt.subplots(figsize=(6, 8))
sns.heatmap(
    pivot,
    ax=ax,
    cmap='RdYlGn',
    vmin=0, vmax=1,
    annot=True, fmt='.2f',
    linewidths=0.5,
    cbar_kws={'label': 'Mean cosine similarity to 1.0×'},
)
ax.set_title('Velocity matrix similarity to full reads (1.0×)', fontsize=13, pad=10)
ax.set_xlabel('Reads ratio', fontsize=11)
ax.set_ylabel('')
ax.set_xticklabels([f'{r}×' for r in pivot.columns], rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_to_ref.pdf'), dpi=150)
plt.close()
print("Saved velocity_sim_to_ref.pdf")

# -------------------------------------------------------
# Figure 2: per-method 5×5 pairwise heatmaps (grid)
# -------------------------------------------------------
methods_present = [m for m in METHOD_ORDER if m in pairwise]
n = len(methods_present)
ncols = 4
nrows = (n + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.2))
axes = axes.flatten()

ratio_labels = [f'{r}×' for r in ratios]

for idx, mname in enumerate(methods_present):
    ax = axes[idx]
    mat = pairwise[mname]
    sns.heatmap(
        mat,
        ax=ax,
        cmap='RdYlGn',
        vmin=-0.2, vmax=1,
        annot=True, fmt='.2f',
        linewidths=0.5,
        xticklabels=ratio_labels,
        yticklabels=ratio_labels,
        cbar=False,
    )
    ax.set_title(mname, fontsize=10, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(labelsize=8)

# Hide unused subplots
for idx in range(len(methods_present), len(axes)):
    axes[idx].axis('off')

# Add colorbar
sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(vmin=-0.2, vmax=1))
sm.set_array([])
fig.colorbar(sm, ax=axes[:len(methods_present)], shrink=0.6, label='Mean cosine similarity')

plt.suptitle('Pairwise velocity matrix similarity across reads ratios', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_pairwise_all.pdf'),
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved velocity_sim_pairwise_all.pdf")

# -------------------------------------------------------
# Figure 3: line plot — similarity to 1.0x across ratios
# -------------------------------------------------------
METHOD_COLORS = {
    'velocyto':      '#D62728FF',
    'scVelo(sto)':   '#cab2d6',
    'scVelo(dyn)':   '#6a3d9a',
    'VeloAE':        '#33a02c',
    'dynamo':        '#FFB547FF',
    'VeloVAE':       '#b2df8a',
    'MultiVelo':     '#843C39FF',
    'UniTVelo':      '#bcbd22',
    'DeepVelo(VAE)': '#1f78b4',
    'cellDancer':    '#e377c2',
    'PhyloVelo':     '#E7CB94FF',
    'veloVI':        '#ff9896',
    'LatentVelo':    '#17BECFFF',
    'DeepVelo(GCN)': '#a6cee3',
    'TFvelo':        '#7B4173FF',
    'STT':           '#ff7f00',
    'biVI':          '#8c564b',
    'SDEvelo':       '#e6ab02',
    'DeepKINET':     '#66c2a5',
}

fig, ax = plt.subplots(figsize=(8, 5))
for mname in methods_present:
    mat = pairwise[mname]
    ref_col = mat[:, ref_idx]  # similarity of each ratio to 1.0x
    xs = [float(r) for r in ratios]
    color = METHOD_COLORS.get(mname, 'grey')
    ax.plot(xs, ref_col, marker='o', markersize=4, linewidth=0.9,
            color=color, label=mname)

ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.7)
ax.set_xlabel('Reads ratio', fontsize=12)
ax.set_ylabel('Mean cosine similarity to 1.0×', fontsize=12)
ax.set_title('Velocity matrix stability across reads downsampling', fontsize=13)
ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlim(0.15, 1.05)
ax.legend(fontsize=7.5, ncol=2, loc='lower right', framealpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_to_ref_line.pdf'), dpi=150)
plt.close()
print("Saved velocity_sim_to_ref_line.pdf")

print("\nDone. All outputs in:", out_dir)
