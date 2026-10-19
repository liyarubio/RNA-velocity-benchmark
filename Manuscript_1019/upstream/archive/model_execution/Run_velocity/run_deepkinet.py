import sys
import os
# Import input_path to get reference embedding
from config import input_path, output_path, seed, count_path

sys.path.append('tools/deepkinet')
# import utils as utils
import scanpy as sc
import scvelo as scv
import numpy as np
import torch

np.random.seed(seed)
torch.manual_seed(seed)
import workflow

h5ad_path = count_path

print(f"read h5ad: {h5ad_path}\n")

# Load raw count data
adata = sc.read_h5ad(h5ad_path)
scv.pp.filter_and_normalize(adata, min_shared_counts=0, n_top_genes=2000)
raw_adata = sc.read_h5ad(h5ad_path)

print("------------------- DeepKINET ------------------------\n")

# Process
# sc.tl.leiden(adata) # need label information
adata.var['highly_variable'] = True # use all genes
adata.layers['spliced'] = raw_adata[:, adata.var_names].layers['spliced']
adata.layers['unspliced'] = raw_adata[:, adata.var_names].layers['unspliced']

# Run DeepKINET estimation
adata, deepkinet_exp = workflow.estimate_kinetics(adata)
adata.layers['velocity'] = adata.layers['DeepKINET_velocity'].copy()

# ==================== FIX: Transfer Embedding ====================
# Read the processed reference data to get embeddings (e.g., UMAP)
print(f"Reading reference embedding from: {input_path}")
adata_ref = sc.read_h5ad(input_path)

# Transfer X_umap coordinates
# Ensure cell numbers match before assignment
if adata.n_obs == adata_ref.n_obs:
    # Check for standard UMAP key
    if 'X_umap' in adata_ref.obsm.keys():
        adata.obsm['X_umap'] = adata_ref.obsm['X_umap']
        print("Transferred X_umap from reference data.")
    # Check for RNA-specific UMAP (common in Seurat objects)
    elif 'X_umap.rna' in adata_ref.obsm.keys():
        adata.obsm['X_umap'] = adata_ref.obsm['X_umap.rna']
        print("Transferred X_umap.rna as X_umap.")
    # Check for WNN UMAP
    elif 'X_wnn.umap' in adata_ref.obsm.keys():
        adata.obsm['X_umap'] = adata_ref.obsm['X_wnn.umap']
        print("Transferred X_wnn.umap as X_umap.")
    
    # Transfer CellType annotations if missing in count data
    if 'CellType' in adata_ref.obs.columns:
        adata.obs['CellType'] = adata_ref.obs['CellType'].values
        print("Transferred CellType annotations.")
else:
    print("Warning: Cell count mismatch between DeepKINET result and reference input. Skipping embedding transfer.")
# =================================================================

print(f"Saving DeepKINET results to {output_path}deepkinet.h5ad")
adata.write_h5ad(output_path + "deepkinet.h5ad")