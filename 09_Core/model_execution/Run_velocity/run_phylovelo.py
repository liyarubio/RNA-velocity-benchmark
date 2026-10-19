import phylovelo as pv
import pandas as pd
import scanpy as sc
import numpy as np
import pickle
import scvelo as scv
import sys
import os
import scipy.sparse as sp

# 1. Force flush stdout (Ensure logs are visible)
sys.stdout.reconfigure(line_buffering=True)

from config import output_path, pv_adata_path, pv_pkl_path

# ================= 1. Load Data =================
print(f"Loading PhyloVelo pickle from: {pv_pkl_path}")
if not os.path.exists(pv_pkl_path):
    raise FileNotFoundError(f"Pickle file not found: {pv_pkl_path}")

with open(pv_pkl_path, 'rb') as f:
    sd = pickle.load(f)

# Try to find cell and gene names if they exist (for alignment)
pv_cells = None
pv_genes = None
if hasattr(sd, 'x_normed') and isinstance(sd.x_normed, pd.DataFrame):
    pv_cells = sd.x_normed.index
    pv_genes = sd.x_normed.columns

# ================= 2. Run Inference =================
print("Running PhyloVelo inference (updating sd in-place)...")
# This step ensures sd.velocity is populated/updated
pv.velocity_inference(sd, sd.depths, cutoff=0.95, target='x_normed')

# Extract the raw velocity object without forcing flatten.
# sd.velocity can be scalar-like, 1D, 2D dense, or sparse depending on pipeline output.
raw_velocity = sd.velocity
if sp.issparse(raw_velocity):
    rv = raw_velocity
else:
    rv = np.asarray(raw_velocity)

rv_shape = rv.shape
rv_ndim = rv.ndim if not sp.issparse(rv) else 2
print(f"Raw velocity type: {type(raw_velocity)}")
print(f"Raw velocity shape: {rv_shape}, ndim: {rv_ndim}")

# ================= 3. Load AnnData =================
print(f"Loading AnnData from: {pv_adata_path}")
adata = sc.read_h5ad(pv_adata_path)
print(f"AnnData shape: {adata.shape} (Cells, Genes)")

# Ensure unique indices before any reindex/subset operation.
# Some input datasets have duplicated var/obs names, which breaks
# pandas Index.get_indexer used internally by AnnData slicing.
adata.var_names_make_unique()
adata.obs_names_make_unique()

# ================= 4. Adaptive Assignment Logic =================
n_cells = adata.n_obs
n_vars = adata.n_vars

print(f"Checking dimensions: Velocity shape={rv_shape}, Cells={n_cells}, Genes={n_vars}")

velocity_matrix = None

if rv_ndim == 1:
    v = np.asarray(rv).reshape(-1)
    v_len = v.shape[0]
    print(f"1D velocity length: {v_len}")

    # Gene-wise vector
    if v_len == n_vars:
        print(">>> Detected GENE-WISE velocity matching adata genes.")
        velocity_matrix = np.tile(v.reshape(1, -1), (n_cells, 1))
    elif pv_genes is not None and v_len == len(pv_genes):
        print(">>> Detected GENE-WISE velocity in x_normed gene space; subsetting to common genes.")
        pv_genes = pd.Index(pv_genes)
        if not pv_genes.is_unique:
            print("   -> pv_genes is not unique; keeping first occurrence for each gene.")
            pv_genes = pv_genes[~pv_genes.duplicated(keep="first")]
            v = v[~pd.Index(sd.x_normed.columns).duplicated(keep="first")]

        common_genes = adata.var_names.intersection(pv_genes)
        if len(common_genes) == 0:
            raise ValueError("No overlapping genes between adata.var_names and x_normed.columns.")
        if len(common_genes) < adata.n_vars:
            print(f"   -> Subsetting adata genes to {len(common_genes)} common genes.")
        adata = adata[:, common_genes].copy()
        vel_series = pd.Series(v, index=pv_genes)
        aligned_gene_vel = vel_series.loc[adata.var_names].values
        velocity_matrix = np.tile(aligned_gene_vel.reshape(1, -1), (n_cells, 1))

    # Cell-wise vector
    elif v_len == n_cells:
        print(">>> Detected CELL-WISE velocity matching adata cells.")
        velocity_matrix = np.tile(v.reshape(-1, 1), (1, n_vars))
    elif pv_cells is not None and v_len == len(pv_cells):
        print(">>> Detected CELL-WISE velocity in x_normed cell space; subsetting to common cells.")
        pv_cells = pd.Index(pv_cells)
        if not pv_cells.is_unique:
            print("   -> pv_cells is not unique; keeping first occurrence for each cell.")
            pv_cells = pv_cells[~pv_cells.duplicated(keep="first")]
            v = v[~pd.Index(sd.x_normed.index).duplicated(keep="first")]

        common_cells = adata.obs_names.intersection(pv_cells)
        if len(common_cells) == 0:
            raise ValueError("No overlapping cells between adata.obs_names and x_normed.index.")
        if len(common_cells) < adata.n_obs:
            print(f"   -> Subsetting adata cells to {len(common_cells)} common cells.")
        adata = adata[common_cells].copy()
        vel_series = pd.Series(v, index=pv_cells)
        aligned_cell_vel = vel_series.loc[adata.obs_names].values.reshape(-1, 1)
        velocity_matrix = np.tile(aligned_cell_vel, (1, n_vars))
    else:
        raise ValueError(
            f"Unsupported 1D velocity length ({v_len}). "
            f"Expected one of: n_obs ({n_cells}), n_vars ({n_vars}), "
            f"len(x_normed.index) ({len(pv_cells) if pv_cells is not None else 'N/A'}), "
            f"len(x_normed.columns) ({len(pv_genes) if pv_genes is not None else 'N/A'})."
        )

elif rv_ndim == 2:
    if sp.issparse(rv):
        dense_v = rv.toarray()
    else:
        dense_v = np.asarray(rv)

    r, c = dense_v.shape
    print(f"2D velocity matrix shape: ({r}, {c})")

    # Direct shape match
    if (r, c) == (n_cells, n_vars):
        print(">>> Velocity matrix already matches adata shape.")
        velocity_matrix = dense_v
    # Transposed direct shape match
    elif (r, c) == (n_vars, n_cells):
        print(">>> Velocity matrix appears transposed; transposing to (cells, genes).")
        velocity_matrix = dense_v.T
    # x_normed-aligned matrix: map by names into adata space
    elif pv_cells is not None and pv_genes is not None and (r, c) == (len(pv_cells), len(pv_genes)):
        print(">>> Velocity matrix in x_normed space; subsetting to common cells/genes.")
        pv_cells = pd.Index(pv_cells)
        pv_genes = pd.Index(pv_genes)
        if not pv_cells.is_unique:
            keep_rows = ~pv_cells.duplicated(keep="first")
            pv_cells = pv_cells[keep_rows]
            dense_v = dense_v[keep_rows, :]
            print("   -> Removed duplicated pv_cells by first occurrence.")
        if not pv_genes.is_unique:
            keep_cols = ~pv_genes.duplicated(keep="first")
            pv_genes = pv_genes[keep_cols]
            dense_v = dense_v[:, keep_cols]
            print("   -> Removed duplicated pv_genes by first occurrence.")

        common_cells = adata.obs_names.intersection(pv_cells)
        common_genes = adata.var_names.intersection(pv_genes)
        if len(common_cells) == 0 or len(common_genes) == 0:
            raise ValueError("No overlapping cells or genes between adata and x_normed space.")
        if len(common_cells) < adata.n_obs or len(common_genes) < adata.n_vars:
            print(f"   -> Subsetting adata to common shape ({len(common_cells)}, {len(common_genes)}).")
        adata = adata[common_cells, common_genes].copy()
        vel_df = pd.DataFrame(dense_v, index=pv_cells, columns=pv_genes)
        vel_df = vel_df.loc[adata.obs_names, adata.var_names]
        velocity_matrix = vel_df.values
    elif pv_cells is not None and pv_genes is not None and (r, c) == (len(pv_genes), len(pv_cells)):
        print(">>> Velocity matrix transposed x_normed space; transposing then subsetting to common cells/genes.")
        pv_cells = pd.Index(pv_cells)
        pv_genes = pd.Index(pv_genes)
        if not pv_cells.is_unique:
            keep_cols = ~pv_cells.duplicated(keep="first")
            pv_cells = pv_cells[keep_cols]
            dense_v = dense_v[:, keep_cols]
            print("   -> Removed duplicated pv_cells by first occurrence.")
        if not pv_genes.is_unique:
            keep_rows = ~pv_genes.duplicated(keep="first")
            pv_genes = pv_genes[keep_rows]
            dense_v = dense_v[keep_rows, :]
            print("   -> Removed duplicated pv_genes by first occurrence.")

        common_cells = adata.obs_names.intersection(pv_cells)
        common_genes = adata.var_names.intersection(pv_genes)
        if len(common_cells) == 0 or len(common_genes) == 0:
            raise ValueError("No overlapping cells or genes between adata and x_normed space.")
        if len(common_cells) < adata.n_obs or len(common_genes) < adata.n_vars:
            print(f"   -> Subsetting adata to common shape ({len(common_cells)}, {len(common_genes)}).")
        adata = adata[common_cells, common_genes].copy()
        vel_df = pd.DataFrame(dense_v.T, index=pv_cells, columns=pv_genes)
        vel_df = vel_df.loc[adata.obs_names, adata.var_names]
        velocity_matrix = vel_df.values
    else:
        raise ValueError(
            f"Unsupported 2D velocity shape ({r}, {c}). "
            f"Expected one of {(n_cells, n_vars)}, {(n_vars, n_cells)}, "
            f"(len(x_normed.index), len(x_normed.columns))="
            f"({len(pv_cells) if pv_cells is not None else 'N/A'}, {len(pv_genes) if pv_genes is not None else 'N/A'})."
        )
else:
    raise ValueError(f"Unsupported velocity ndim: {rv_ndim}")

# ================= 5. Save =================
# Calculate moments (required for streamplot visualization)
print("Calculating moments...")
scv.pp.moments(adata)

print(f"Assigning velocity layer with shape {velocity_matrix.shape}...")
adata.layers['velocity'] = velocity_matrix

save_path = output_path + 'phylovelo.h5ad'
print(f"Saving result to: {save_path}")
adata.write_h5ad(save_path)
print("PhyloVelo finished successfully.")