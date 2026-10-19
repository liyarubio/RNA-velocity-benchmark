# %%
import sys
sys.path.append("/home/liyr/RNAvelocity/7_STT/STT-release/example_notebooks")

import sctt as st

import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import os

np.random.seed(2024)
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
# %%
file_path = '../real_data/'
file_list = os.listdir(file_path)

out_path = '../all_result/'
out_file = os.listdir(out_path)



# %%
for file in cell_types:
    print(file)
    cluster = cell_types[file]
    adata = sc.read_h5ad(file_path + file+'.h5ad')
    print(adata)

    adata.obs['attractor'] = adata.obs[cluster].values
    n_states = 4
    # print(n_states )

    try :
        adata_aggr = st.dynamical_iteration(adata,
                                            return_aggr_obj=True,
                                            n_states = n_states)

        adata_aggr.write_h5ad(out_path + file+'/STT.h5ad')

    except Exception as e:
        print('find error {}'.format(e))




