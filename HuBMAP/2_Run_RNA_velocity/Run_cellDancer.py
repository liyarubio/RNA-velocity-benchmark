# %%
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

os.chdir('/home/liyr/HuBMAP/sup_result')

# %%
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)

out_path1 = '/home/liyr/HuBMAP/RNA_velocity_result/cellDancer_input/'
out_path2 = '/home/liyr/HuBMAP/RNA_velocity_result/cellDancer_output/'
out_path3 = '/home/liyr/HuBMAP/RNA_velocity_result/cellDancer/'

cluster = 'CytoTRACE2_Potency'

n_jobs = 28

# %%
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        print(adata)

        # adata to dataframe
        cdutil.adata_to_df_with_embed(adata,
                              us_para=['Mu','Ms'],
                              cell_type_para=cluster,
                              embed_para='X_umap',
                              save_path=out_path1 + file.strip('.h5ad') + '.csv'
                             )
        
        # run cellDancer
        df = pd.read_csv(out_path1 + file.strip('.h5ad') + '.csv')
        loss_df, cellDancer_df=cd.velocity(df,n_jobs=n_jobs,speed_up = False)
        # save cellDancer dataframe 
        cellDancer_df.to_csv(out_path2 + file.strip('.h5ad') + '.csv')

        # dataframe to adata
        adata_cd = export_velocity_to_dynamo(cellDancer_df,adata)
        print(adata_cd)
        adata_cd.layers["velocity_S"] = adata_cd.layers["velocity_S"].toarray()
        adata_cd.write_h5ad(out_path3 + file)

    except Exception as e:
        print('find error {}'.format(e))
            


