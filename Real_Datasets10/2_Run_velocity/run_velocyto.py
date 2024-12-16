import time

import velocyto as vcy
import scanpy as sc
import scvelo as scv
import numpy as np
import anndata as ad
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import time
SEED = 2024
np.random.seed(SEED)


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
def run_velocyto(adata,out_path,file_name,emb):
    print(file_name)
    vlm = scv.utils.convert_to_loom(adata)
    vlm.fit_gammas(fit_offset=False)
    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift()
    vlm.extrapolate_cell_at_t()
    velocity_s = np.transpose(vlm.velocity)

    adata.layers['velocity'] = velocity_s
    folder_path = out_path+file_name + '/'
    os.makedirs(folder_path,exist_ok=True)
    result_path = folder_path+'Velocyto.h5ad'
    adata.write(result_path)
    print(adata)
    print(f"Save velocity matrix in : {result_path}\n\n\n\n")
    return

data_path = '../real_data/'
file_name_list = os.listdir(data_path)
os.makedirs('../all_result/',exist_ok=True)
out_path = '../all_result/'
for file_name in cell_types:
    print('run:{}'.format(file_name))
    file_path = data_path + file_name + '.h5ad'
    adata = ad.read_h5ad(file_path)
    emb = embs[file_name]
    out_path = out_path
    run_velocyto(adata, out_path, file_name,emb)



