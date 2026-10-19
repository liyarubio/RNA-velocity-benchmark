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

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
# fix random seeds for reproducibility
SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

os.chdir("/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/tmp")


def run_deepvelogb(adata,file_name):
    configs = {
        "name": "DeepVelo",  # name of the experiment
        'n_gpu': 0,
        "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
        "arch": {'args': {'pred_unspliced': True}},
        "trainer": {"verbosity": 0},  # increase verbosity to show training progress
        'loss_dir':'/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/{}/7_DeepveloGB/loss.csv'.format(file_name)
    }
    configs = update_dict(Constants.default_configs, configs)
    print(configs)

    velocity(adata, mask_zero=False)
    trainer = train(adata, configs)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata,basis = "pca")
    scv.tl.velocity_embedding(adata, basis = "umap")
    out = '/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/{}/7_DeepveloGB/'.format(file_name) + 'anndata.h5ad'
    adata.write_h5ad(out)


for simuID  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try:
        
        folder_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/velocity/" + simuID+"/7_DeepveloGB/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        
        h5ad_path = "/lustre/home/zhangxiaochang/project/liyr/data/dyngen/result_0921/anndata/"+simuID+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        adata = sc.read_h5ad(h5ad_path)
        
        print(f"-----------------Simulation: {simuID}-------------------")
        print("------------------- DeepVelo (GB) ------------------------\n")
        
        time_paste = run_deepvelogb(adata,simuID)
        
    except Exception as e:
        print('find error {}'.format(e))

