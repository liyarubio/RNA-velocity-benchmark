import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import anndata as ad
import time
SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)


def run_velovi(adata,out_path):
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    latent_time = vae.get_latent_time()
    t = latent_time
    scaling = 20 / t.max(0)
    velocity_u = vae.get_velocity(velo_mode='unspliced')
    adata.layers['velocity_u'] = velocity_u / scaling
    velocity = vae.get_velocity(velo_mode='spliced')
    adata.layers['velocity'] = velocity / scaling
    
    train_loss =vae.history["elbo_train"]
    val_loss =vae.history["elbo_validation"]
    train_loss = pd.DataFrame(train_loss)
    val_loss = pd.DataFrame(val_loss)

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata,basis = "pca")
    scv.tl.velocity_embedding(adata, basis = "umap")

    adata.write_h5ad(out_path+"/anndata.h5ad")
    
    train_loss.to_csv(out_path+'train_loss.csv')
    val_loss.to_csv(out_path+'val_loss.csv')
    return
    


for simuID  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    print(f"-----------------Simulation: {simuID}-------------------")
    print("------------------- veloVI ------------------------\n")
        
        
    folder_path = "/home/liyr/dygen/dyngen_new/velocity/" + simuID+"/8_veloVI/"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    h5ad_path = "/home/liyr/dygen/dyngen_new/anndata/" +simuID+"/anndata.h5ad"
    print(f"read h5ad: {h5ad_path}\n")
    from velovi import preprocess_data, VELOVI
    adata = sc.read_h5ad(h5ad_path)
    cached_memory = torch.cuda.memory_cached()
    torch.cuda.empty_cache()
    run_velovi(adata,folder_path)

