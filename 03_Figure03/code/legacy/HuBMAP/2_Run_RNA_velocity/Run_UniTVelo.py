import pandas as pd
import numpy as np
import os
import scanpy as sc
import scvelo as scv
import unitvelo as utv
import tensorflow as tf
os.environ['TF_USE_LEGACY_KERAS'] = 'True'

os.chdir('/home/liyr/HuBMAP/sup_result')

SEED = 2024
np.random.seed(SEED)
tf.random.set_seed(SEED)

file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)
out_path = '/home/liyr/HuBMAP/RNA_velocity_result/UniTVelo/'
cluster = 'CytoTRACE2_Potency'

velo_config = utv.config.Configuration()

for file in file_list:
    print(file)

    try :
        adata = sc.read_h5ad(file_path+ file)
        adata.var['highly_variable'] = True # use all genes
        adata = utv.run_model(adata,label=cluster,config_file=velo_config)
        print(adata)

        adata.write_h5ad(out_path + file)

    except Exception as e:
            print('find error {}'.format(e))