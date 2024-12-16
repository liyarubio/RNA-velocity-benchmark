import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import dynamo as dyn
import scvelo as scv
import time

SEED = 2024
np.random.seed(SEED)
def run_dynamo(adata,out_path,file_name):

    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata,model='deterministic')
    folder_path = out_path
    result_path = folder_path+file_name
    adata.write_h5ad(result_path)
# def run_dynamo_dtm(adata,out_path,file_name,emb):
#     time1 = time.time()
#
#     dyn.pp.recipe_monocle(adata)
#     dyn.tl.dynamics(adata, model='deterministic', cores=10)
#     time2 = time.time()
#
#     folder_path = out_path
#     if not os.path.exists(folder_path):
#         os.makedirs(folder_path)
#     result_path = folder_path+"{}".format(file_name)
#
#     adata.write_h5ad(result_path)
#     return time2-time1

# out_path = '../big_data/result/Dynamo/'
out_path1 = '/home/liyr/HuBMAP/RNA_velocity_result/Dynamo/'
if not os.path.exists(out_path1):
    os.makedirs(out_path1)
print('has create dir')
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir('/home/liyr/HuBMAP/RNA_velocity_result/scvelo/')
key = 'CytoTRACE2_Potency'
kk = 0
total_file = 151
for file in file_list:
    kk=kk+1
    emb='umap'
    type= 'CytoTRACE2_Potency'
    print('run {} / {}'.format(kk,total_file))
    adata = ad.read_h5ad(file_path+file)
    run_dynamo(adata, out_path1, file)


