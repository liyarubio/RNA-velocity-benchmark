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

    velo_config = utv.config.Configuration()
    velo_config.MAX_ITER = 20000

    time1 = time.time()
    adata = utv.run_model(adata,label,config_file=velo_config)
    time2 = time.time()
    print(out_path+file_name+'/UniTvelo.h5ad')
    adata.write_h5ad(out_path+file_name+'/UniTvelo.h5ad')
    return time2-time1

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

file_path = '../real_data/'
out_path = '../all_result/'
for file in cell_types:
    key = cell_types[file]
    print('run: {}'.format(file))
    h5ad_path = '../real_data/{}.h5ad'.format(file)
    adata = sc.read_h5ad(h5ad_path)
    sc.pp.highly_variable_genes(adata,n_top_genes=2000)
    print(adata.var)
    adata.uns['losses'] = []
    paste_time = run_unitvelo(adata,out_path,file,key)
    
# file = 'organogenesis_chondrocyte'
# key = cell_types[file]
# adata = out_path+'organogenesis_chondrocyte/Dynamo_stochastic.h5ad'
# adata = ad.read_h5ad(adata)
# sc.pp.highly_variable_genes(adata,n_top_genes=2000)
# print(adata.var)
# adata.uns['losses'] = []
# paste_time = run_unitvelo(adata,out_path,file,key)