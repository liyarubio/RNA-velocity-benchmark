import numpy as np
import scvelo as scv
import torch
import scanpy as sc
import os
from deepvelo.utils.scatter import scatter
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo.utils.plot import statplot, compare_plot
from deepvelo import train, Constants
from deepvelo.utils import (
    velocity,
    velocity_confidence,
    continuity_confidence,
    update_dict,
    cross_boundary_correctness,)

os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

SEED = 2024
torch.manual_seed(SEED)
np.random.seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

print(torch.cuda.is_available())

adata = sc.read_h5ad("adata/redeem_young.h5ad")
print(adata)

configs = {
    "name": "DeepVelo_GB", # name of the experiment
    'n_gpu': 1, # whether use gpu
    "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
    "arch":{'args': {'pred_unspliced': True}},
    "trainer": {"verbosity": 0}, # increase verbosity to show training progress
}
configs = update_dict(Constants.default_configs, configs)

# initial velocity
# velocity(adata, mask_zero=False) # use scvelo stochastic
trainer = train(adata, configs)

print(adata)
adata.write_h5ad("adata/DeepVelo_GB.h5ad")