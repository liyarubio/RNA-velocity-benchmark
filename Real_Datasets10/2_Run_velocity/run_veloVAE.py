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
import random


n_gene = 2000
torch.manual_seed(2024)
np.random.seed(2024)

rate_prior = {
    'alpha': (0.0, 1.0),
    'beta': (0.0, 0.5),
    'gamma': (0.0, 0.5)
}


cell_types = {
    'endocrinogenesis_day15':'clusters',
    'DentateGyrus_velocyto':'clusters',
    'human_cd34_bone_marrow':'clusters',
    'Forebrain':'clusters',
    'dentategyrus_scv':'clusters',
    'Hindbrain_GABA_Glio':'Celltype',
    'organogenesis_chondrocyte':'Main_cell_type',
    'erythroid_lineage':'celltype',
    'zebrafish_process':'Cell_type',
    'gastrulation':'stage',
    'MultiVelo_10X_multiome_mouse_brain':'celltype'

}
embs = {
    'dentategyrus_scv':'umap',
        'DentateGyrus_velocyto':'tsne',
        'endocrinogenesis_day15':'umap',
        'human_cd34_bone_marrow':'tsne',
        'Forebrain':'umap',
        'gastrulation':'umap',
        'erythroid_lineage':'umap',
        'Hindbrain_GABA_Glio':'tsne',
    'organogenesis_chondrocyte':'umap',
    'zebrafish_process':'umap',
    'MultiVelo_10X_multiome_mouse_brain':'umap'
}
SEED = 2024
np.random.seed(SEED)
random.seed(SEED)
data_path = '../real_data/'
file_name_list = os.listdir(data_path)
os.makedirs('../all_result/',exist_ok=True)
out_path = '../all_result/'
for file in cell_types:

    emb=embs[file]
    emb = 'X_' + emb
    type= cell_types[file]
    print('run {}  /   10'.format(file))
    file_path = '../real_data/'+ file + '.h5ad'
    print(file_path)
    adata = anndata.read(file_path)
    print(adata)

    figure_path = out_path+file+'/veloVAE'
    data_path = out_path+file +'/'
    name = file
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
                     # gene_plot=gene_plot,
                     figure_path=figure_path,
                     cluster_key=type,
                     embed=emb)

    file_name = file
    vae.save_anndata(adata, 'vae', data_path, file_name='veloVAE.h5ad')


