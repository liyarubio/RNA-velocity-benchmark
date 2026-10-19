import velocyto as vcy
import scanpy as sc
import scvelo as scv
import numpy as np
import anndata as ad
import os
import sys
import pandas as pd


SEED = 2024
np.random.seed(SEED)

os.chdir("/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/tmp")

for simuID  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try:
        folder_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/" + simuID+"/12_velocyto/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        h5ad_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/anndata/"+simuID+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        adata = ad.read_h5ad(h5ad_path)
        
        
        print(f"-----------------Simulation: {simuID}-------------------")
        
        
        print("------------------- Velocyto ------------------------\n")
        
        
        adata.layers["velocity"] = adata.layers["velocity"].todense()
        
        vlm = scv.utils.convert_to_loom(adata)
        
        #vlm._normalize_S(relative_size=vlm.S.sum(0),
        #             target_size=vlm.S.sum(0).mean())
        #vlm._normalize_U(relative_size=vlm.U.sum(0),
        #             target_size=vlm.U.sum(0).mean())
        
        #vlm.perform_PCA()
        #vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=2000, b_maxl=1500, n_jobs=16)
        #vlm.knn_imputation(n_pca_dims=30, balanced=False, n_jobs=16)
        
        
        vlm.fit_gammas(fit_offset=False)
        vlm.predict_U()
        vlm.calculate_velocity()
        vlm.calculate_shift()
        vlm.extrapolate_cell_at_t()
        
        velocity_s = np.transpose(vlm.velocity)
        
        
        #result_path = folder_path+"/Velocyto.npy"
        #np.save(result_path,velocity_s)
        #print(f"Save velocity matrix in : {result_path}\n\n\n\n")
        
        adata.layers["velocity"] = velocity_s
        
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata,basis = "pca")
        scv.tl.velocity_embedding(adata, basis = "umap")
        adata.write_h5ad(folder_path+"/anndata.h5ad")

    except Exception as e:
        print('find error {}'.format(e))