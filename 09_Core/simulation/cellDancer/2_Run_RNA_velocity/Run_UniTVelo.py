import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import unitvelo as utv
import tensorflow as tf
import logging

os.environ['TF_USE_LEGACY_KERAS'] = 'True'

SEED = 2024
np.random.seed(SEED)
tf.random.set_seed(SEED)

file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/UniTVelo/'

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

n_jobs = 26

logging.basicConfig(filename='logs/UniTVelo.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

velo_config = utv.config.Configuration()

for file in file_list:
    print(file)

    try :
        adata = sc.read_h5ad(file_path+ file)
        adata.var['highly_variable'] = True # use all genes
        adata = utv.run_model(adata,label='leiden',config_file=velo_config)
        print(adata)

        adata.write_h5ad(out_path + file)

    except Exception as e:
            logging.exception('Error processing file %s: %s', file, e)