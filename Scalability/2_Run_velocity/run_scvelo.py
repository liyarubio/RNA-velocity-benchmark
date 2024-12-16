
import os
import tracemalloc
tracemalloc.start()
import time
# from memory_profiler import profile
import sys
import glob
import pandas as pd
import numpy as np
import math
import random
import scanpy as sc
import anndata as ad
import scvelo as scv
# import warnings
# warnings.filterwarnings("ignore")
SEED = 2024
np.random.seed(SEED)
random.seed(SEED)
size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
def run_scvelo_dynamic(adata):
    scv.tl.recover_dynamics(adata, n_jobs=26,var_names= "all")
    scv.tl.velocity(adata,  n_jobs=26,mode='dynamical')
    return adata
def run_scvelo_sto(adata):
    scv.tl.velocity(adata, mode='stochastic')
    return adata
sto_time_list = []
sto_size_list = []
dym_time_list = []
dym_size_list = []
out_path1 = '../time_analyse1/scvelo_stochastic/'
out_path2 = '../time_analyse1/scvelo_deterministic/'
if not os.path.exists(out_path1):
    os.makedirs(out_path1)
if not os.path.exists(out_path2):
    os.makedirs(out_path2)
for i in size:
    print(i)
    file = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(i)
    adata = sc.read_h5ad(file)
    time1 = time.time()
    # scv.tl.velocity(adata, mode='stochastic')
    adata1 = run_scvelo_sto(adata)
    time2 = time.time()
    adata2 = run_scvelo_dynamic(adata)
    time4 = time.time()
    spaste_time = time2 - time1
    spaste_time1 = time4 - time2
    sto_time_list.append(spaste_time)
    sto_size_list.append(i)
    dym_time_list.append(spaste_time1)
    dym_size_list.append(i)
    adata1.write('../time_analyse1/scvelo_stochastic/{}.h5ad'.format(i))
    adata2.write('../time_analyse1/scvelo_deterministic/{}.h5ad'.format(i))

df1 = pd.DataFrame()
df1['cellnumber'] = sto_size_list
df1['time'] = sto_time_list
df1.to_csv('../time_analyse1/scvelo_stochastic/time.csv')
df2 = pd.DataFrame()
df2['cellnumber'] = dym_size_list
df2['time'] = dym_time_list
df2.to_csv('../time_analyse1/scvelo_deterministic/time.csv')