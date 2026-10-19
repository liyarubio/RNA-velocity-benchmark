import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import dynamo as dyn
import scanpy as sc
import scvelo as scv
import time
import logging

SEED = 2024
np.random.seed(SEED)
# %%
file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/dynamo/'

n_jobs = 26

# %%
logging.basicConfig(filename='logs/dynamo.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

for file in file_list :
    adata = sc.read_h5ad(file_path+file)
    try :
        dyn.pp.recipe_monocle(adata)
        dyn.tl.dynamics(adata)
        adata.write_h5ad(out_path + file)
 

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)

