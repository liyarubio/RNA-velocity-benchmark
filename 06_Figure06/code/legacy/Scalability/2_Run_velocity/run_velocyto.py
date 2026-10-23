import velocyto as vcy
# import scanpy as sc
import scvelo as scv
import numpy as np
import anndata as ad
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import time
SEED = 2024
np.random.seed(SEED)

def run_velocyto(adata,out_path,file_name):
    print(file_name)

    vlm = scv.utils.convert_to_loom(adata)
    time1 = time.time()
    vlm.fit_gammas(fit_offset=False)
    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift()
    vlm.extrapolate_cell_at_t()
    time2 = time.time()
    
    paste_time = time2 - time1
    folder_path = out_path+file_name
    adata.write(folder_path+'.h5ad')
    return paste_time

size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []
i= 1000
for i in size:
    print(f'run {i}')
    file_path = '../time_analyse/scvelo_deterministic/'+'{}.h5ad'.format(str(i))
    adata = ad.read_h5ad(file_path)
    out_path = '../time_analyse1/velocyto/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    paste_time = run_velocyto(adata,out_path,str(i))
    time_list.append(paste_time)
    size_list.append(i)
   
df = pd.DataFrame()
df['cellnumber'] = size_list
df['time'] = time_list
df.to_csv('../time_analyse1/velocyto/time.csv', index=False)
