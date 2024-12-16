import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import anndata as ad
import time
SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)
#
size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []

def run_velovi(adata,out_path,file_name):
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    time1 = time.time()
    vae = VELOVI(adata)
    vae.train()

    latent_time = vae.get_latent_time()
    t = latent_time
    scaling = 20 / t.max(0)

    velocity_u = vae.get_velocity(velo_mode='unspliced')
    adata.layers['velocity_u'] = velocity_u / scaling

    velocity = vae.get_velocity(velo_mode='spliced')
    adata.layers['velocity'] = velocity / scaling

    time2 = time.time()

    plt.legend()
    paste_time = time2 - time1
    folder_path = out_path+file_name
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    adata.write_h5ad(folder_path+'/'+str(file_name)+'.h5ad')
    return paste_time
out_path = '../result_new/veloVI/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
print('has create dir')
for i in size:

    print(f'run {i}')
    from velovi import preprocess_data, VELOVI

    file_path = '../time_analyse/scvelo_deterministic/'+'{}.h5ad'.format(str(i))
    adata = sc.read_h5ad(file_path)
    out_path = '../time_analyse1/veloVI/'
    cached_memory = torch.cuda.memory_cached()
    print(cached_memory)
    torch.cuda.empty_cache()
    i = str(i)
    paste_time = run_velovi(adata,out_path,i)
    time_list.append(paste_time)
    size_list.append(i)
df = pd.DataFrame()
df['cellnumber'] = size_list
df['time'] = time_list
df.to_csv('../time_analyse1/veloVI/time.csv')