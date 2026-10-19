import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import dynamo as dyn
import scvelo as scv


SEED = 2024
np.random.seed(SEED)
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

def run_dynamo(adata,out_path,file_name,emb):

    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='stochastic')

    folder_path = out_path+file_name
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    result_path = folder_path+"/Dynamo_stochastic.npy"
    # np.save(result_path,velocity_s)
    cluster = 'velocity_S'
    scv.tl.velocity_graph(adata, basis='umap', vkey=cluster)
    adata.write_h5ad(folder_path+'/Dynamo_stochastic.h5ad')

def run_dynamo_dtm(adata,out_path,file_name,emb):
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='deterministic', cores=10)

    folder_path = out_path+file_name
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    velocity_s = adata.layers['velocity_S'].toarray()
    cluster = 'velocity_S'
    scv.tl.velocity_graph(adata, basis='umap', vkey=cluster)
    print(adata)
    adata.write_h5ad(folder_path + '/Dynamo_deterministic.h5ad')

data_path = '../real_data/'
file_name_list = os.listdir(data_path)
os.makedirs('../all_result/',exist_ok=True)
out_path = '../all_result/'

for file_name in embs:
    print('run:{}'.format(file_name))
    # try:
    file_path = data_path + file_name + '.h5ad'
    adata = ad.read_h5ad(file_path)
    print(adata.var_names)
    emb = embs[file_name]
    run_dynamo(adata, out_path, file_name,emb)
    run_dynamo_dtm(adata,out_path, file_name,emb)


