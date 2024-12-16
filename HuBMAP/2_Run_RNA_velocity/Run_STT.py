# %%
import sys
sys.path.append("/home/liyr/RNAvelocity/7_STT/STT-release/example_notebooks")

import sctt as st

import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import os

np.random.seed(2024)

# %%
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)

out_path = '/home/liyr/HuBMAP/RNA_velocity_result/STT/'
out_file = os.listdir(out_path)
cluster = 'CytoTRACE2_Potency'

file_list = [item for item in file_list if item not in out_file]

# %%
for file in file_list:
    adata = sc.read_h5ad(file_path + file)
    print(adata)

    adata.obs['attractor'] = adata.obs[cluster].values
    # n_states = len(adata.obs['attractor'].values.unique())
    n_states = 2

    try :
        adata_aggr = st.dynamical_iteration(adata,
                                        return_aggr_obj=True,
                                        n_states = n_states)
        
        adata_aggr.write_h5ad(out_path + file)

    except Exception as e:
        print('find error {}'.format(e))

    



