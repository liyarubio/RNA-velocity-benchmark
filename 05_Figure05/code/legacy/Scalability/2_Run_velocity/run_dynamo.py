import pandas as pd
import anndata as ad
import numpy  as np
import os
import sys
import dynamo as dyn
import scvelo as scv
import time

SEED = 2024
np.random.seed(SEED)
embs = {'dentategyrus_scv':'umap',
        'DentateGyrus_velocyto':'tsne',
        'endocrinogenesis_day15':'umap',
        'human_cd34_bone_marrow':'tsne',
        'Forebrain':'umap',
        'gastrulation':'umap',
        'erythroid_lineage':'umap',
        'Hindbrain_GABA_Glio':'tsne',
        'organogenesis_chondrocyte':'pca',
        'zebrafish_process':'umap'}
def run_dynamo(adata,out_path,file_name,emb):

    time1 = time.time()
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='stochastic', cores=26)
    time2 = time.time()
    folder_path = out_path
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    result_path = folder_path+"{}.h5ad".format(file_name)
    adata.write_h5ad(result_path)
    return time2-time1


def run_dynamo_dtm(adata,out_path,file_name,emb):
    time1 = time.time()

    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='deterministic', cores=26)
    time2 = time.time()

    folder_path = out_path
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    result_path = folder_path+"{}.h5ad".format(file_name)
    adata.write_h5ad(result_path)
    return time2-time1

i = 1000
size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
sto_time_list = []
sto_size_list = []
dym_time_list = []
dym_size_list = []
os.makedirs('../time_analyse1/dynamo_deterministic/',exist_ok=True)
os.makedirs('../time_analyse1/dynamo_stochastic/', exist_ok=True)
for i in size:
    print('run:',i)
    h5ad_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(i)
    adata = ad.read_h5ad(h5ad_path)
    emb = 'umap'
    out_path1 = '../time_analyse1/dynamo_stochastic/'
    out_path2 = '../time_analyse1/dynamo_deterministic/'
    sto_time = run_dynamo(adata, out_path1, i,emb)
    dym_time = run_dynamo_dtm(adata, out_path2, i,emb)
    sto_time_list.append(sto_time)
    sto_size_list.append(i)
    dym_time_list.append(dym_time)
    dym_size_list.append(i)
df1 = pd.DataFrame()
df1['cellnumber'] = sto_size_list
df1['time'] = sto_time_list
df1.to_csv('../time_analyse1/dynamo_stochastic/time.csv')
df2 = pd.DataFrame()
df2['cellnumber'] = dym_size_list
df2['time'] = dym_time_list
df2.to_csv('../time_analyse1/dynamo_deterministic/time.csv')