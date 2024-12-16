import velocyto as vcy
import scanpy as sc
import scvelo as scv
import numpy as np
import anndata as ad
import os
import sys
import pandas as pd
import logging

SEED = 2024
np.random.seed(SEED)

file_path = './1_simulation_datasets/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/velocyto/'
if not os.path.exists(out_path):
    os.makedirs(out_path)

n_jobs = 26

logging.basicConfig(filename='logs/velocyto.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

for file in file_list:
    print(file)

    try :
        adata = sc.read_h5ad(file_path+ file)

        vlm = scv.utils.convert_to_loom(adata)

        vlm.fit_gammas(fit_offset=False)
        vlm.predict_U()
        vlm.calculate_velocity()
        vlm.calculate_shift()
        vlm.extrapolate_cell_at_t()
        velocity_s = np.transpose(vlm.velocity)

        adata.layers['velocity'] = velocity_s
        adata.write_h5ad(out_path + file)


    except Exception as e:
            logging.exception('Error processing file %s: %s', file, e)