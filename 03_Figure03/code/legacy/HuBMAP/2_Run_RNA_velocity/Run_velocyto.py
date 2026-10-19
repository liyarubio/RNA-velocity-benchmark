import velocyto as vcy
import scanpy as sc
import scvelo as scv
import numpy as np
import anndata as ad
import os
import sys
import pandas as pd

SEED = 2024
np.random.seed(SEED)

file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)
out_path = '/home/liyr/HuBMAP/RNA_velocity_result/velocyto/'

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
            print('find error {}'.format(e))