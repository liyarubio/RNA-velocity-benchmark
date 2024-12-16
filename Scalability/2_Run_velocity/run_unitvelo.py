import os
import sys
import glob
import pandas as pd
import numpy as np
import math
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import scvelo as scv
import unitvelo as utv
import time
import os
os.environ['TF_USE_LEGACY_KERAS'] = 'True'

SEED = 2024
np.random.seed(SEED)


def run_unitvelo(adata,out_path,file_name,label):
    file_name = str(file_name)
    velo_config = utv.config.Configuration()
    velo_config.MAX_ITER = 10000
    time1 = time.time()
    adata = utv.run_model(adata,label,config_file=velo_config)
    time2 = time.time()
    adata.write_h5ad(out_path+file_name+'.h5ad')
    return time2-time1




size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []
key =  'celltype'

out_path = '../time_analyse1/unitvelo/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
for i in size:
    print('run:',i)
    h5ad_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(i)
    adata = sc.read_h5ad(h5ad_path)
    adata.uns['losses'] = []
    paste_time = run_unitvelo(adata,out_path,i,key)
    time_list.append(paste_time)
    size_list.append(i)
    df = pd.DataFrame()
    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/unitvelo/{}_time.csv'.format(i))
