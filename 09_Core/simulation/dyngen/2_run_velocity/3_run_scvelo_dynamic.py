import os
import sys
import glob
import pandas as pd
import numpy as np
import math
import random
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt



import scvelo as scv

SEED = 2024
np.random.seed(SEED)
random.seed(SEED)

os.chdir("/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/tmp")

for simuID  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try:

        folder_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/" + simuID+"/3_scvelo_dynamical/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        
        h5ad_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/anndata/"+simuID+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        print(f"-----------------Simulation: {simuID}-------------------")
        print("------------------- scvelo ------------------------\n")
        
        adata = sc.read_h5ad(h5ad_path)
        
        scv.tl.recover_dynamics(adata, n_jobs=10,var_names= "all")
        scv.tl.velocity(adata, mode='dynamical')
        
        
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata,basis = "pca")
        scv.tl.velocity_embedding(adata, basis = "umap")
        adata.write_h5ad(folder_path+"/anndata.h5ad")
        
    except Exception as e:
        print('find error {}'.format(e))