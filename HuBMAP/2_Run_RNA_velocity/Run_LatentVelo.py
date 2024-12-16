# %%
import latentvelo as ltv
import numpy as np
import scanpy as sc
import scvelo as scv
import pandas as pd
import os

SEED = 2024
np.random.seed(SEED)

# %%
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)

out_path = '/home/liyr/HuBMAP/RNA_velocity_result/LatentVelo/'
out_path2 = '/home/liyr/HuBMAP/RNA_velocity_result/LatentVelo_latent/'

cluster = 'CytoTRACE2_Potency'

spliced_key = 'spliced'
unspliced_key = 'unspliced'

# %%
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        print(adata)

        adata.var['velocity_genes'] = True

        adata = ltv.utils.standard_clean_recipe(adata)

        spliced_library_sizes = adata.layers[spliced_key].sum(1)
        unspliced_library_sizes = adata.layers[unspliced_key].sum(1)

        if len(spliced_library_sizes.shape) == 1:
            spliced_library_sizes = spliced_library_sizes[:,None]
        if len(unspliced_library_sizes.shape) == 1:
            unspliced_library_sizes = unspliced_library_sizes[:,None]

        adata.obs['spliced_size_factor'] = spliced_library_sizes #spliced_all_size_factors
        adata.obs['unspliced_size_factor'] = unspliced_library_sizes #unspliced_all_size_factors

        model = ltv.models.VAE(observed = 2000) # observed: number of genes
        epochs, val_ae, val_traj = ltv.train(model,adata,name="veloVAE")

        latent_adata, adata = ltv.output_results(model, adata, gene_velocity=True,embedding='umap')
        adata.write_h5ad(out_path + file)
        latent_adata.write_h5ad(out_path2 + file)
        
    except Exception as e:
        print('find error {}'.format(e))



