

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

def run_deepvelogb(adata,file_name,outpath):
    configs = {
        "name": "DeepVelo",  # name of the experiment
        'n_gpu': 1,
        "loss": {"args": {"coeff_s": autoset_coeff_s(adata),
                          'inner_batch_size': 100}},
        "arch": {'args': {'pred_unspliced': True}},
        "trainer": {
            "verbosity": 0},  # increase verbosity to show training progress
        'loss_dir':'../gene_samples_result/DeepVelo_GB/{}_loss.csv'.format(file_name)
    }
    configs = update_dict(Constants.default_configs, configs)
    print(configs)
    time1 = time.time()
    velocity(adata, mask_zero=False)
    trainer = train(adata, configs)
    time2 = time.time()
    out = outpath + str(file_name) + '.h5ad'
    adata.write_h5ad(out)
    return time2 - time1

size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []
kk = 0
out_path = '../time_analyse1/Deepvelo_GB/'
os.makedirs(out_path, exist_ok=True)
for i in size:
    kk=kk+1
    print('run: {}/{}'.format(kk,i))
    file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(str(i))
    adata = sc.read_h5ad(file_path)
    try:
        time_paste = run_deepvelogb(adata,i,out_path)
        time_list.append(time_paste)
        size_list.append(i)
    except Exception as e:
        print(e)

    df = pd.DataFrame()
    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/Deepvelo_GB/{}_time.csv'.format(str(i)))

