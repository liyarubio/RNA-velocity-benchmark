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
def run_velovi(adata,out_path,file_name):
    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    vae = VELOVI(adata)
    vae.train()
    latent_time = vae.get_latent_time()
    t = latent_time
    scaling = 20 / t.max(0)
    velocity_u = vae.get_velocity(velo_mode='unspliced')
    adata.layers['velocity_u'] = velocity_u / scaling
    velocity = vae.get_velocity(velo_mode='spliced')
    adata.layers['velocity'] = velocity / scaling
    train_loss =vae.history["elbo_train"]
    val_loss =vae.history["elbo_validation"]
    train_loss = pd.DataFrame(train_loss)
    val_loss = pd.DataFrame(val_loss)
    plt.legend()
    adata.write_h5ad(out_path+file_name)
    file_name = file_name.split('.')[0]
    train_loss.to_csv(out_path+'{}_train_loss.csv'.format(file_name))
    val_loss.to_csv(out_path+'{}_val_loss.csv'.format(file_name))
    return
out_path = '/home/liyr/HuBMAP/RNA_velocity_result/veloVI/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
print('has create dir')
file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir('/home/liyr/HuBMAP/RNA_velocity_result/scvelo/')
size_list = []
key = 'CytoTRACE2_Potency'
kk = 0
total_file = len(file_list)

for i in range(total_file):
    kk=kk+1
    emb='umap'
    file = file_list[i]
    from velovi import preprocess_data, VELOVI
    print('run {} / {}'.format(kk,total_file))
    adata = sc.read_h5ad(file_path+file)
    cached_memory = torch.cuda.memory_cached()
    torch.cuda.empty_cache()
    run_velovi(adata,out_path,file)

