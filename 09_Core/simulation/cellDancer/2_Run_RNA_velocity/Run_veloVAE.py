# %%
import scanpy as sc
import numpy as np
import sys
import torch
import os
import scvelo as scv
from scipy.sparse import csr_matrix

SEED = 2024
np.random.seed(SEED)
torch.manual_seed(SEED)
import logging

sys.path.append('/home/liyr/Benchmark/VeloVAE-master')
import velovae as vv

# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/veloVAE/'

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)


if not os.path.exists(out_path):
    os.makedirs(out_path)

logging.basicConfig(filename='logs/veloVAE.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


# %%
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        print(adata)
        
        adata.layers['spliced'] = csr_matrix(adata.layers['spliced'])
        adata.layers['unspliced'] = csr_matrix(adata.layers['unspliced'])


        # process
        vv.preprocess(adata, n_gene=2000)
        vae = vv.VAE(adata, 
                    tmax=20, 
                    dim_z=len(set(adata.obs['leiden'])), 
                    device='cuda:0')
        
        config = {
            # You can change any hyperparameters here!
            # 'learning_rate': 1e-3,
            # 'learning_rate_ode': 2e-3,
            # 'learning_rate_post': 1e-3
        }
        vae.train(adata,
                config=config,
                plot=False,
                #gene_plot=gene_plot,
                #figure_path=figure_path,
                embed='umap')
        
        vae.save_model("Sup/velovae", 'encoder_vae', 'decoder_vae')
        vae.save_anndata(adata, 'vae',out_path, file_name=file)
        
    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)


