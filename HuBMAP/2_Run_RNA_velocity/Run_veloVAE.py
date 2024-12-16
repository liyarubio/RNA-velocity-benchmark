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

sys.path.append('/home/liyr/Benchmark/VeloVAE-master')
import velovae as vv

# %%
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)

out_path = '/home/liyr/HuBMAP/RNA_velocity_result/VeloVAE/'

cluster = 'CytoTRACE2_Potency'


# %%
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        print(adata)

        # process
        vv.preprocess(adata, n_gene=2000)
        vae = vv.VAE(adata, 
                    tmax=20, 
                    dim_z=len(set(adata.obs[cluster])), 
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
        
        vae.save_model("/home/liyr/HuBMAP/sup_result/velovae", 'encoder_vae', 'decoder_vae')
        vae.save_anndata(adata, 'vae',out_path, file_name=file)
        
    except Exception as e:
        print('find error {}'.format(e))


