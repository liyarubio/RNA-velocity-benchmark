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

SEED = 2024
np.random.seed(SEED)
random.seed(SEED)
size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]

time_list = []
size_list = []
output_path = '../time_analyse1/latentvelo/'

if not os.path.exists(output_path):
    os.makedirs(output_path)
for file in size:
    print('run:',file)
    file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(file)
    output_path = '../time_analyse1/latentvelo/'
    print('RESULT IN :',output_path)
    adata = ad.read(file_path)
    emb_list= []
    for emb in adata.obsm:
        emb_list.append(emb)
    if 'X_tsne' in emb_list:
        emb = 'tsne'
    else :
        emb = 'umap'
    adata = ltv.utils.standard_clean_recipe(adata,normalize_library = False)
    adata.var['velocity_genes'] = True
    file = str(file)
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
    n = len(adata.var)
    model = ltv.models.VAE(observed = n) # observed: number of genes
    try :
        time1 = time.time()
        epochs, val_ae, val_traj = ltv.train(model, adata,name=file)
        # try:
        latent_adata, adata = ltv.output_results(model, adata, gene_velocity=True,
                                                 embedding=emb)
        time2 = time.time()
        adata.write(output_path + file +'.h5ad')
        spase_time = time2-time1
        time_list.append(spase_time)
        size_list.append(file)
    except Exception as e:
        print(f"在迭代{file}时发生异常：{e}")
        continue

    df = pd.DataFrame()
    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/latentvelo/{}time1.csv'.format(str(file)), index=False)
