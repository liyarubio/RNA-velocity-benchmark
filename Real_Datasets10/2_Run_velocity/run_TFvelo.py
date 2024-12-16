# '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('AGG')
import os, sys
sys.path.insert(1, './home/liyr/hpz/处理cy2文件/TFvelo-main/')
# import TFvelo as TFv
import random

n_jobs = 26

def check_data_type(adata):
    for key in list(adata.var):
        if adata.var[key][0] in ['True', 'False']:
            adata.var[key] = adata.var[key].map({'True': True, 'False': False})
    return
def data_type_tostr(adata, key):
    if key in adata.var.keys():
        if adata.var[key][0] in [True, False]:
            adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    return
def run_tfvelo(adata,key):
    adata.obs['clusters'] = adata.obs[key].copy()
    adata.layers["total"] = adata.layers["spliced"].todense() + adata.layers["unspliced"].todense()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    n_cells, n_genes = adata.X.shape

    TFv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) #include the following steps
    adata.X = adata.layers["total"].copy()
    gene_names = []
    for tmp in adata.var_names:
        gene_names.append(tmp.upper())
    TFv.pp.moments(adata, n_pcs=30, n_neighbors=15)
    TFv.pp.get_TFs(adata, databases='ENCODE ChEA')
    print(adata)
    adata.uns['genes_pp'] = np.array(adata.var_names)
    n_jobs = 28
    flag = TFv.tl.recover_dynamics(adata, n_jobs=n_jobs, max_iter=20, var_names='all',
                                   WX_method = "lsq_linear", WX_thres=20, max_n_TF=99, n_top_genes=2000,
                                   fit_scaling=True, use_raw=0, init_weight_method= "correlation",
                                   n_time_points=1000)
    return



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
    'gastrulation':'stage'
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
for file  in cell_types:

    print('run: {}'.format(file))
    adata = sc.read_h5ad(data_path+file+'.h5ad')
    type = cell_types[file]
    run_tfvelo(adata,type)

    adata = sc.read_h5ad('../all_result/{}/TFvelo.h5ad'.format(file))
    n_cells = adata.shape[0]
    expanded_scaling_y = np.expand_dims(np.array(adata.var['fit_scaling_y']),0).repeat(n_cells,axis=0)
    print(expanded_scaling_y)
    adata.layers['velocity'] = adata.layers['velo_hat'] / expanded_scaling_y
    print(adata.layers['velocity'])
    folder_path = '../all_result/{}/TFvelo.h5ad'.format(file)

    adata.write(folder_path)
