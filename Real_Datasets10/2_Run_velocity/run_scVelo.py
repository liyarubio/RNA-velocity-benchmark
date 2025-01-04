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
import os
SEED = 2024
np.random.seed(SEED)
random.seed(SEED)
data_path = '../real_data/'
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
def run_scvelo_dynamic(adata):
    # print(file_name)
    scv.tl.recover_dynamics(adata, n_jobs=26,var_names= "all")
    scv.tl.velocity(adata, mode='dynamical')
    return adata
def run_scvelo_sto(adata):
    scv.tl.velocity(adata, mode='stochastic')
    return adata

sto_time_list = []
sto_size_list = []
dym_time_list = []
dym_size_list = []


m_data = []
for file in cell_types :
    cluster = cell_types[file]
    adata = sc.read_h5ad(data_path+file+'.h5ad')
    series = adata.obs[cluster].value_counts()
    count_positive = (series > 0).sum()
    out_path1 = '../all_result/{}/'.format(file)
    if not os.path.exists(out_path1):
        os.makedirs(out_path1)

    if count_positive > 1:
        print('file_name is {}'.format(file))
        try :
            sc.pp.neighbors(adata)
            scv.tl.umap(adata)
            scv.pp.filter_and_normalize(adata, n_top_genes=2000)
            time1 = time.time()
            adata1 = run_scvelo_sto(adata)
            adata2 = run_scvelo_dynamic(adata)
            adata1.write('../all_result/{}/scvelo_stochastic.h5ad'.format(file))
            adata2.write('../all_result/{}/scvelo_dynamical.h5ad'.format(file))
        except Exception as e:
            print('find error {}'.format(e))
