import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import dynamo as dyn
import scvelo as scv
import time
import logging

SEED = 2024
np.random.seed(SEED)
def run_dynamo(adata,out_path,file_name):
    dyn.pp.recipe_monocle(adata)

    # dyn.tl.dynamics(adata,model='deterministic')
    dyn.tl.dynamics(adata) # auto : stochastic

    folder_path = out_path
    result_path = folder_path+file_name
    adata.write_h5ad(result_path)

out_path1 = '/home/liyr/HuBMAP/RNA_velocity_result/Dynamo_stochastic/'
if not os.path.exists(out_path1):
    os.makedirs(out_path1)
print('has create dir')
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir('/home/liyr/HuBMAP/RNA_velocity_result/scvelo/')

file2 = os.listdir('/home/liyr/HuBMAP/RNA_velocity_result/Dynamo_stochastic')
file_list = [item for item in file_list if item not in file2]

key = 'CytoTRACE2_Potency'

logging.basicConfig(filename='dynamo_stochastic.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

for file in file_list:
    emb='umap'
    type= 'CytoTRACE2_Potency'
    adata = ad.read_h5ad(file_path+file)

    try :
        run_dynamo(adata, out_path1, file)
    except Exception as e:
            logging.exception('Error processing file %s: %s', file, e)

