# %%
import sys
sys.path.append("/home/liyr/RNAvelocity/7_STT/STT-release/example_notebooks")

import sctt as st

import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import os
import logging

np.random.seed(2024)

# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/STT/'

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

n_jobs = 26

logging.basicConfig(filename='logs/STT.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


# %%
for file in file_list:
    adata = sc.read_h5ad(file_path + file)
    print(adata)

    adata.obs['attractor'] = adata.obs['leiden'].values
    # n_states = len(adata.obs['attractor'].values.unique())
    n_states = 2

    try :
        adata_aggr = st.dynamical_iteration(adata,
                                        return_aggr_obj=True,
                                        n_states = n_states)
        
        adata_aggr.obs_names = adata.obs_names
        overlap_var = adata_aggr.var.index[adata_aggr.var.index.isin(adata.var.index)]

        adata_sub = adata[:,overlap_var]
        adata_aggr = adata_aggr[:,overlap_var]

        adata_aggr.layers['ground_truth_velocity'] = adata_sub.layers['ground_truth_velocity']

        adata_aggr.write_h5ad(out_path + file)

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)

    


