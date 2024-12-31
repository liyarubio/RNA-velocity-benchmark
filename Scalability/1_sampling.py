import scanpy as sc
import pandas as pd
import numpy as np

# load data
adata = sc.read_h5ad('../data_2000gene/gastrulation.h5ad')
size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
# i = 1000
time_list = pd.DataFrame()
for i in size:
    random_indices = np.random.choice(adata.obs_names, i, replace=False)
    adata1 = adata[random_indices]
    result_path = './cell_splite/gastrulation_'+str(i)+'cell.h5ad'
    adata1.write_h5ad(result_path)
