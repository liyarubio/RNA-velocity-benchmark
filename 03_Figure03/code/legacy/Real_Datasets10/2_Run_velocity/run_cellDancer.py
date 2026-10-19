import scanpy as sc
import scvelo as scv

import os
import sys
import glob
import pandas as pd
import math
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import seaborn as sns
import celldancer as cd
import celldancer.simulation as cdsim
import celldancer.utilities as cdutil
import celldancer.cdplt as cdplt
from celldancer.cdplt import colormap
from celldancer.utilities import export_velocity_to_dynamo
import time
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity

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
for file in cell_types:
    os.makedirs('../all_result/{}/CellDancer/'.format(file),exist_ok=True)
    emb = 'X_'+embs[file]
    type = cell_types[file]
    file_path = '../real_data/'+file + '.h5ad'
    adata = ad.read(file_path)
    print(adata.uns['neighbors'])
    
    print('run {} '.format(file))
    x ='../all_result/{}/CellDancer/'.format(file)+file+'.csv'
    cdutil.adata_to_df_with_embed(adata,
                           us_para=['unspliced','spliced'],
                           cell_type_para=type,
                           embed_para='X_umap',
                           save_path=x
                           )
    df = pd.read_csv(x)
    loss_df, cellDancer_df=cd.velocity(df,n_jobs=26,speed_up = False)
    cellDancer_df['cellID'] = cellDancer_df['cellID'].astype(str) ### LYR edit
    adata_cd = export_velocity_to_dynamo(cellDancer_df,adata)
    print(adata_cd)
    adata_cd.layers["velocity_S"] = adata_cd.layers["velocity_S"].toarray()
    adata_cd.write_h5ad('../all_result/{}/CellDancer.h5ad'.format(file))

