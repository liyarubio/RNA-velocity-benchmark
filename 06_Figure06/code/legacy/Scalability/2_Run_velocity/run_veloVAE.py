import anndata
import numpy as np
import scvelo as scv
import scanpy as sc
import sys
import pandas as pd
import torch
import os.path
sys.path.insert(1, '../VeloVAE-master/')
import velovae as vv
import time
# %load_ext autoreload
# %autoreload 2

n_gene = 2000
torch.manual_seed(2024)
np.random.seed(2024)

rate_prior = {
    'alpha': (0.0, 1.0),
    'beta': (0.0, 0.5),
    'gamma': (0.0, 0.5)
}


data = {
    "clusters" : ['dentategyrus_scv','DentateGyrus_velocyto','endocrinogenesis_day15','human_cd34_bone_marrow','Forebrain'],
    "celltype": ['gastrulation','erythroid_lineage'],
    "Celltype" :['Hindbrain_GABA_Glio'],
    'Main_cell_type' : ['organogenesis_chondrocyte'],
    "Cell_type" :['zebrafish_process']
}

size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]

time_list = []
size_list = []
for file in size:
    # file = 5000
    file = str(file)
    type= 'celltype'
    print(file,'begain')
    file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(str(file))
    adata = anndata.read(file_path)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


    figure_path = '../time_analyse1/veloVAE/'
    data_path = '../time_analyse/veloVAE/'
    name = '{}.h5ad'.format(file)
    out_path = '../time_analyse1/veloVAE/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    gene_plot = ["Gng12", "Smoc1"]
    torch.manual_seed(2024)
    np.random.seed(2024)
    time1 =time.time()
    vae = vv.VAE(adata,
                 tmax=20,
                 dim_z=5,
                 device='cuda:0')
    config = {

        # You can change any hyperparameters here!
    }
    loss = vae.train(adata,
          config=config,
          plot=False,
          gene_plot=gene_plot,
          figure_path=figure_path,
          cluster_key=type,
          embed='umap')
    time2 = time.time()
    time_paste =time2-time1
    time_list.append(time_paste)
    size_list.append(file)
    file_name = file + '.h5ad'
    vae.save_anndata(adata, 'vae', data_path, file_name=name)
    df = pd.DataFrame()
    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/veloVAE/{}_time.csv'.format(str(file)))