import latentvelo as ltv
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import anndata as ad
import time
import random

# warnings.filterwarnings("ignore")

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
    type = cell_types[file]
    output_path = out_path +file + '/LatentVelo.h5ad'
    print('run  {}'.format(file))
    adata = ad.read('../real_data/'+file+'.h5ad')
    adata = ltv.utils.standard_clean_recipe(adata,normalize_library = False)
    adata.var['velocity_genes'] = True
    try :
        spliced_key = 'spliced'
        unspliced_key = 'unspliced'

        spliced_library_sizes = adata.layers[spliced_key].sum(1)
        unspliced_library_sizes = adata.layers[unspliced_key].sum(1)

        if len(spliced_library_sizes.shape) == 1:
            spliced_library_sizes = spliced_library_sizes[:,None]
        if len(unspliced_library_sizes.shape) == 1:
            unspliced_library_sizes = unspliced_library_sizes[:,None]

        adata.obs['spliced_size_factor'] = spliced_library_sizes #spliced_all_size_factors
        adata.obs['unspliced_size_factor'] = unspliced_library_sizes #unspliced_all_size_factors
        n = adata.X.shape[1]
        print(n)
        model = ltv.models.VAE(observed = n) # observed: number of genes
        time1 = time.time()
        epochs, val_ae, val_traj = ltv.train(model, adata,name=file)
        latent_adata, adata = ltv.output_results(model, adata, gene_velocity=True,
                                                 embedding=emb)
        adata.write(output_path)
        print(adata)
    except Exception as e:
        print(e)
    # size_list.append(file)


