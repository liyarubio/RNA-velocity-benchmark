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


for file  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try :
        folder_path = "/home/liyr/dygen/dyngen_new/velocity/" + file+"/13_veloVAE/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        
        h5ad_path = "/home/liyr/dygen/dyngen_new/anndata/"+file+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        print(f"-----------------Simulation: {file}-------------------")
        print("------------------- veloVAE ------------------------\n")
        # read data
        adata = sc.read_h5ad(h5ad_path)
        
        # process
        vv.preprocess(adata, n_gene=2000)
        vae = vv.VAE(adata, 
                    tmax=20, 
                    dim_z=len(set(adata.obs["clusters"])), 
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
                embed='pca')
        
        vae.save_model(folder_path, 'encoder_vae', 'decoder_vae')
        vae.save_anndata(adata, 'vae',folder_path, file_name="anndata.h5ad")
        
        adata = scv.read(folder_path+"anndata.h5ad")
        adata.layers["velocity"] = adata.layers["vae_velocity"]
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata,basis = "pca")
        scv.tl.velocity_embedding(adata, basis = "umap")
        adata.write_h5ad(folder_path+"/anndata.h5ad")
    except Exception as e:
        print('find error {}'.format(e))


