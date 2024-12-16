import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import anndata as ad
import time
SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
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
def run_velovi(adata,out_path,file_name):

    adata = preprocess_data(adata)
    VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
    time1 = time.time()
    vae = VELOVI(adata)
    vae.train()

    latent_time = vae.get_latent_time()
    t = latent_time
    scaling = 20 / t.max(0)

    velocity_u = vae.get_velocity(velo_mode='unspliced')
    adata.layers['velocity_u'] = velocity_u / scaling

    velocity = vae.get_velocity(velo_mode='spliced')
    adata.layers['velocity'] = velocity / scaling

    train_loss =vae.history["elbo_train"]
    val_loss =vae.history["elbo_validation"]
    train_loss = pd.DataFrame(train_loss)
    val_loss = pd.DataFrame(val_loss)
    plt.legend()
    folder_path = out_path+file_name
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    print(folder_path+'/veloVI.h5ad')
    adata.write_h5ad(folder_path+'/veloVI.h5ad')
    train_loss.to_csv(folder_path+'/veloVI_train_loss.csv')
    val_loss.to_csv(folder_path+'/veloVI_val_loss.csv')
    return
out_path = '../all_result/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
file_path = '../real_data/'
size_list = []
for file in cell_types:
        print(f'run {file}')
        emb = embs[file]
        from velovi import preprocess_data, VELOVI
        adata = sc.read_h5ad(file_path+'{}.h5ad'.format(file))
        cached_memory = torch.cuda.memory_cached()
        torch.cuda.empty_cache()
        run_velovi(adata,out_path,file)
        plt_path = out_path +file + '/velovi_stream'+ '.png'
        plt_name = file + '_VeloVI'
        cluster = 'velocity'
        type =  cell_types[file]
        scv.tl.velocity_graph(adata)
        scv.pl.velocity_embedding(adata, basis=emb, title='p1_name',vkey = cluster,color=[type],dpi=150, smooth=0.5,show=False)
        scv.tl.velocity_confidence(adata, vkey=cluster)
        scv.pl.velocity_embedding_grid(adata, basis=emb,title='p2_name' ,vkey = cluster,color=[type],dpi=150,smooth=0.5,show=False)
        scv.pl.velocity_embedding_stream(adata, basis=emb, title='p3_name',vkey =  cluster, color=[type],dpi=150,smooth=0.5,show=False)

