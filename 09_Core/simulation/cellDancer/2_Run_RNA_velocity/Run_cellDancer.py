# %%
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os
import celldancer as cd
import celldancer.utilities as cdutil
from celldancer.utilities import export_velocity_to_dynamo
import logging

SEED = 2024
np.random.seed(SEED)

# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path1 = './Sup/cellDancer_input/'
out_path2 = './Sup/cellDancer_output/'
out_path3 = './RNA_velocity_result/cellDancer/'

out_file = os.listdir(out_path3)

file_list = [item for item in file_list if item not in out_file]
print(file_list)

n_jobs = 26

logging.basicConfig(filename='logs/cellDancer.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# %%
for file in file_list:
    try :
        # read data
        print(file)
        adata = sc.read_h5ad(file_path+ file)
        print(adata)

        # adata to dataframe
        cdutil.adata_to_df_with_embed(adata,
                              us_para=['Mu','Ms'],
                              cell_type_para='leiden',
                              embed_para='X_umap',
                              save_path=out_path1 + file.strip('.h5ad') + '.csv'
                             )
        
        # run cellDancer
        df = pd.read_csv(out_path1 + file.strip('.h5ad') + '.csv')
        loss_df, cellDancer_df=cd.velocity(df,n_jobs=n_jobs,speed_up = False)
        # save cellDancer dataframe 
        cellDancer_df.to_csv(out_path2 + file.strip('.h5ad') + '.csv')
        # cellDancer_df = pd.read_csv(out_path2 + file.strip('.h5ad') + '.csv')
        
        # dataframe to adata
        cellDancer_df['cellID'] = cellDancer_df['cellID'].astype(str) ### LYR edit
        adata_cd = export_velocity_to_dynamo(cellDancer_df,adata)

        print(adata_cd)
        adata_cd.layers["velocity_S"] = adata_cd.layers["velocity_S"].toarray()
        adata_cd.write_h5ad(out_path3 + file)

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)
            