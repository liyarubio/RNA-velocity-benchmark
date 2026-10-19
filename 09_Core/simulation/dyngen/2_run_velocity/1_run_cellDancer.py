import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os
import celldancer as cd
import celldancer.utilities as cdutil
from celldancer.utilities import export_velocity_to_dynamo

SEED = 2024
np.random.seed(SEED)



os.chdir("/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/tmp")

for simuID  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try:
        folder_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/" + simuID+"/1_celldancer/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        h5ad_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/anndata/"+simuID+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        adata = sc.read_h5ad(h5ad_path)
        
        print(f"-----------------Simulation: {simuID}-------------------")
        print("------------------- cellDancer ------------------------\n")
        
        cdutil.adata_to_df_with_embed(adata,
                                      us_para=['Mu','Ms'],
                                      cell_type_para='clusters',
                                      embed_para='X_pca',
                                      save_path= folder_path + 'cell_type_u_s.csv')
        
        
        path = folder_path + 'cell_type_u_s.csv'
        print(f"read gene X cell unsplcie and splice dataframe: {path}\n")
        
        cell_type_u_s=pd.read_csv(path)
        
        
        print("Run cellDaner\n")
        loss_df, cellDancer_df=cd.velocity(cell_type_u_s,n_jobs=10)
        
        # save
        celldancer_result_path = folder_path + "cellDancer.csv"
        cellDancer_df.to_csv(celldancer_result_path)
        print("Save cellDacner dataframe in : {celldancer_result_path}\n")
        
        adata_cd = export_velocity_to_dynamo(cellDancer_df,adata)
        print(adata_cd)
        adata_cd.layers["velocity"] = adata_cd.layers["velocity_S"].toarray()
        
        
        scv.tl.velocity_graph(adata_cd)
        scv.tl.velocity_embedding(adata_cd,basis = "pca")
        scv.tl.velocity_embedding(adata_cd, basis = "umap")
        adata_cd.write_h5ad(folder_path+"/anndata.h5ad")
        
    except Exception as e:
        print('find error {}'.format(e))
