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
        'loss_dir':'/home/liyr/HuBMAP/RNA_velocity_result/DeepVelo_GB/{}_loss.csv'.format(file_name)
    }
    configs = update_dict(Constants.default_configs, configs)
    print(configs)
    time1 = time.time()
    velocity(adata, mask_zero=False)
    trainer = train(adata, configs)
    out = '/home/liyr/HuBMAP/RNA_velocity_result/DeepVelo_GB/' +file_name
    adata.write_h5ad(out)

file_path= '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)
os.makedirs( '/home/liyr/HuBMAP/RNA_velocity_result/DeepVelo_GB/', exist_ok=True )
emb= "umap"
cluster = 'velocity'
k = 0
file_list2 = os.listdir('/home/liyr/HuBMAP/RNA_velocity_result/DeepVelo_GB/')
list_set = [i for i in file_list if i not in file_list2]
print(len(list_set))

for file in list_set:
    k = k + 1
    # file_name = file_list[file]
    print('run : {}/{}'.format(k,len(list_set)))
    adata = sc.read_h5ad(file_path+file)
    # out_path = './time_analyse/deepvelo_GB/'
    try:
        time_paste = run_deepvelogb(adata,file)
    except Exception as e:
        print(e)