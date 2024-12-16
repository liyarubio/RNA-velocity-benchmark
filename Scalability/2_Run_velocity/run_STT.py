# %%
import sys
sys.path.append("/home/liyr/RNAvelocity/7_STT/STT-release/example_notebooks")

import sctt as st

import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import os
import time
np.random.seed(2024)


size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []
out_path = '../time_analyse1/STT/'
os.makedirs(out_path, exist_ok=True)
# %%
df = pd.DataFrame()
for file in size:
    print(file)
    file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(str(file))
    adata = sc.read_h5ad(file_path)
    scv.pp.moments(adata)
    print(adata)
    cluster = 'celltype'
    adata.obs['attractor'] = adata.obs[cluster].values
    n_states = 4
    # print(n_states )

    # try :
    time1 = time.time()
    adata_aggr = st.dynamical_iteration(adata,
                                        return_aggr_obj=True,
                                        n_states = n_states)
    time2 = time.time()
    paste_time = time2 - time1
    time_list.append(paste_time)
    size_list.append(file)
    adata_aggr.write_h5ad(out_path + str(file)+'.h5ad')
    #
    # except Exception as e:
    #     print('find error {}'.format(e))

    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/STT/{}_time.csv'.format(file), index=False)


