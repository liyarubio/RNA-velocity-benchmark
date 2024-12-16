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
import logging
from velovi import preprocess_data, VELOVI

SEED = 2024
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/veloVI/'
if not os.path.exists(out_path):
    os.makedirs(out_path)

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

logging.basicConfig(filename='logs/veloVI.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

for file in file_list:
    try :
        adata = sc.read_h5ad(file_path+file)

        adata = preprocess_data(adata)
        VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
        vae = VELOVI(adata)
        vae.train()

        latent_time = vae.get_latent_time()
        t = latent_time
        scaling = 20 / t.max(0)

        velocity_u = vae.get_velocity(velo_mode='unspliced')
        adata.layers['velocity_u'] = velocity_u / scaling
        
        velocity = vae.get_velocity(velo_mode='spliced')
        adata.layers['velocity'] = velocity / scaling

        adata.write_h5ad(out_path+file)
        
    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)

