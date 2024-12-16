# %%
import latentvelo as ltv
import numpy as np
import scanpy as sc
import scvelo as scv
import pandas as pd
import os
import logging

SEED = 2024
np.random.seed(SEED)

# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/LatentVelo/'
out_path2 = './RNA_velocity_result/LatentVelo_latent/'

if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(out_path2):
    os.makedirs(out_path2)

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

spliced_key = 'spliced'
unspliced_key = 'unspliced'

logging.basicConfig(filename='logs/LatentVelo.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


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

        model = ltv.models.VAE(observed = 1000) # observed: number of genes
        epochs, val_ae, val_traj = ltv.train(model,adata,name="veloVAE")

        latent_adata, adata = ltv.output_results(model, adata, gene_velocity=True,embedding='umap')
        adata.write_h5ad(out_path + file)
        latent_adata.write_h5ad(out_path2 + file)
        
    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)



