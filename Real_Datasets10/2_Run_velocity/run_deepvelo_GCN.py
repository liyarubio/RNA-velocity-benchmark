import numpy as np
import scvelo as scv
import torch
from umap import UMAP
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu
import scanpy as sc
import pandas as pd
import os
import sys
import time
import anndata as ad
from deepvelo.utils.scatter import scatter
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo.utils.plot import statplot, compare_plot
from deepvelo import train, Constants
from deepvelo.utils import (
    velocity,
    velocity_confidence,
    continuity_confidence,
    update_dict,
    cross_boundary_correctness,
)
import os

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
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

SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

#
def run_deepvelogb(adata,file_name):
    configs = {
        "name": "DeepVelo",  # name of the experiment
        'n_gpu': 0,
        "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
        "arch": {'args': {'pred_unspliced': True}},
        "trainer": {
            "verbosity": 0},  # increase verbosity to show training progress
        'loss_dir':'/home/liyr/hpz/all_result/MultiVelo_10X_multiome_mouse_brain/{}_loss.csv'.format(file_name)
    }
    configs = update_dict(Constants.default_configs, configs)
    print(configs)
    time1 = time.time()
    velocity(adata, mask_zero=False)
    trainer = train(adata, configs)
    out = '/home/liyr/hpz/all_result/{}/DeepVelo_GB.h5ad'.format(file_name)
    adata.write_h5ad(out)

file_path='../real_data/'


for file in cell_types:

    print('run : {}'.format(file))
    adata = sc.read_h5ad(file_path+file+'.h5ad')
    try:
       run_deepvelogb(adata,file)
    except Exception as e:
        print(e)