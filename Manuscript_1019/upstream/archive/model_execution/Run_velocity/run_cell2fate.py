import torch
import numpy as np
import anndata as ad
import cell2fate as c2f
import scanpy as sc
import os
import scvelo as scv
from config import input_path, output_path, seed, data_cluster, n_job, count_path

# Set random seeds
np.random.seed(seed)
torch.manual_seed(seed) # Add torch seed
torch.cuda.empty_cache() 

# ========== Step 1. Load Data ==========
h5ad_path = count_path
print(f"read h5ad: {h5ad_path}\n")
adata = sc.read_h5ad(h5ad_path)

# [Suggestion 3]: Remove preprocessing here to avoid conflict with get_training_data
# scv.pp.filter_genes_dispersion(adata,n_top_genes=2000) 
print("Original adata:", adata)

cluster_column = data_cluster
# Check use_gpu parameter compatibility
use_gpu = torch.cuda.is_available()
accelerator = "gpu" if use_gpu else "cpu"

# ========== Step 2. Get Training Data (includes normalization and filtering) ==========
adata_train = c2f.utils.get_training_data(
    adata,
    cells_per_cluster=10**5,
    cluster_column=cluster_column,
    remove_clusters=[],
    min_shared_counts=0,
    n_var_genes=2000 # Highly variable gene selection happens here
)

# ========== Step 3. Determine Number of Modules ==========
suggested_modules = c2f.utils.get_max_modules(adata_train)
print('-----------------------------------------------')
print(f"Suggested modules: {suggested_modules}")
n_modules = 6 # You chose to manually override to 6

# ========== Step 4. Setup AnnData ==========
c2f.Cell2fate_DynamicalModel.setup_anndata(
    adata_train,
    spliced_label='spliced',
    unspliced_label='unspliced',
    cluster_label=cluster_column # [Suggestion 4]: Register cluster labels
)

# ========== Step 5. Initialize Model and Train ==========
mod = c2f.Cell2fate_DynamicalModel(adata_train, n_modules=n_modules)
torch.autograd.set_detect_anomaly(True)

try:
    mod.train(
        batch_size=64,
        lr=5e-4,
        # precision=32, # Can be kept; usually default for float32 data
        # [Suggestion 5]: For compatibility, older scvi-tools versions might prefer use_gpu
        # If error "unknown argument 'accelerator'" occurs, change to use_gpu=use_gpu
        accelerator=accelerator, 
        max_epochs=300 # Recommended to explicitly set epochs, though source has defaults
    )
except Exception as e:
    print("❌ Training interrupted:", e)
    for name, p in mod.module.named_parameters():
        if torch.isnan(p).any():
            print(f" NaN detected in parameter: {name}")
            break
    raise e

# ========== Step 6. Export Results ==========
# [Suggestion 2]: Removed duplicate export calls
adata_train = mod.export_posterior(adata_train)
adata_train = mod.compute_module_summary_statistics(adata_train)

# ========== Step 7. Post-processing and Plotting ==========
adata_raw = sc.read_h5ad(input_path)

# [Critical Fix 1]: Use index alignment to assign UMAP to prevent dimension errors
common_cells = adata_train.obs_names.intersection(adata_raw.obs_names)
if len(common_cells) < len(adata_train):
    print("Warning: Some training cells not found in raw data for UMAP mapping.")
# Safe assignment
adata_train.obsm['X_umap'] = adata_raw[adata_train.obs_names].obsm['X_umap'].copy()

mod.plot_module_summary_statistics(adata_train)

# delete=False is not needed for plotting unless you want to use this Velocity layer for something else
mod.compute_and_plot_total_velocity(
    adata_train,
    delete=False,
    plotting_kwargs={
        "color": cluster_column,
        "legend_fontsize": 10, 
        "legend_loc": "right_margin"
    }
)

# Rename layer
if "Velocity" in adata_train.layers:
    adata_train.layers["velocity"] = np.array(adata_train.layers["Velocity"])
    del adata_train.layers["Velocity"]

print(f"Saving to {output_path}cell2fate.h5ad")
adata_train.write_h5ad(output_path + "cell2fate.h5ad")