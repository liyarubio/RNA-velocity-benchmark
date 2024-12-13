import sys
sys.path.append("/home/liyr/RNAvelocity/7_STT/STT-release/example_notebooks")

import sctt as st

import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import os

np.random.seed(2024)



# 'cycle_simple_seed3' # error
# 'disconnected_seed2' # error

#for f  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
for f  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    print(f)
    
    h5ad_path = "/home/liyr/dygen/dyngen_new/anndata/"+f+"/anndata.h5ad"
    print(f"read h5ad: {h5ad_path}\n")
    
    adata = sc.read_h5ad(h5ad_path)

    folder_path = "/home/liyr/dygen/dyngen_new/velocity/" + f+"/4_STT/"
    if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    # scv.pp.moments(adata)
    adata.obs['attractor'] = adata.obs['clusters'].values

    try:
        print(f"-----------------Simulation: {f}-------------------")
        print("------------------- STT ------------------------\n")
        n_states = len(adata.obs['clusters'].values.unique())
        if n_states <= 2:
            n_states = 2

        adata_aggr = st.dynamical_iteration(adata,
                                        return_aggr_obj=True,
                                        n_states = n_states)
        scv.tl.velocity_graph(adata_aggr)
        scv.tl.velocity_embedding(adata_aggr,basis = "pca")
#       scv.tl.velocity_embedding(adata_aggr, basis = "umap")
        
        adata_aggr.write_h5ad(folder_path +"anndata.h5ad")
    
    except Exception as e:
        print(f"An error occurred: {e}")
        #n_states = 2
        #adata_aggr = st.dynamical_iteration(adata,
        #                                return_aggr_obj=True,
        #                                n_states = n_states)
        #scv.tl.velocity_graph(adata_aggr)
        #scv.tl.velocity_embedding(adata_aggr,basis = "pca")
        #scv.tl.velocity_embedding(adata_aggr, basis = "umap")
        #adata_aggr.write_h5ad(folder_path +"anndata.h5ad")
        continue
    
