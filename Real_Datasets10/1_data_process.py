import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import scanpy as sc
import scvelo as scv

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

file_path = '/home/liyr/hpz/real_data/'
for file in cell_types:
    cluster = cell_types[file]
    adata = sc.read_h5ad(file_path+file +'.h5ad')
    scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
    scv.pp.normalize_per_cell(adata)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=50)
    scv.tl.umap(adata)
    adata.write(file_path+file+'.h5ad')
    print(adata)
    # scv.pl.umap(adata, color=cluster)