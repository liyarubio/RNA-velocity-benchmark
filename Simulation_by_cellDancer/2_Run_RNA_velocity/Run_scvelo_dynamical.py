# %%
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os
import logging

# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/scvelo_dynamical/'

n_jobs = 26

# %%
logging.basicConfig(filename='logs/scvelo_dynamical.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

for file in file_list :
    adata = sc.read_h5ad(file_path+file)
    try :
        scv.tl.recover_dynamics(adata, n_jobs=n_jobs,var_names= "all")
        scv.tl.velocity(adata, mode='dynamical')
        adata.write_h5ad(out_path + file)

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)
