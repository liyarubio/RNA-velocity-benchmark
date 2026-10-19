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
import logging

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
# conda activate Deepvelo_GB

SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/DeepVeloGB/'

# check finished file 
out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

n_jobs = 26

logging.basicConfig(filename='logs/DeepVeloGB.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def run_deepvelogb(adata,file):
    configs = {
        "name": "DeepVelo",  # name of the experiment
        'n_gpu': 0,
        "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
        "arch": {'args': {'pred_unspliced': True}},
        "trainer": {"verbosity": 0}  # increase verbosity to show training progress
    }

    configs = update_dict(Constants.default_configs, configs)
    print(configs)
    time1 = time.time()
    velocity(adata, mask_zero=False)
    trainer = train(adata, configs)
    out = out_path +file
    adata.write_h5ad(out)


for file in file_list:
    adata = sc.read_h5ad(file_path + file)
    print(adata)

    try :
        run_deepvelogb(adata,file)

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)
