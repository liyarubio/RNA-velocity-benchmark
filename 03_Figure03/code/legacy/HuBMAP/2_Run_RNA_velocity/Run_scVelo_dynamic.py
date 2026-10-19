import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os

file_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
file_list = os.listdir(file_path)

out_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo_dynamic/'

cluster = 'CytoTRACE2_Potency'
n_jobs = 28

for file in file_list:
    print(file)
    
    try :
        adata = sc.read_h5ad(file_path+ file)
        scv.tl.recover_dynamics(adata, n_jobs=n_jobs,var_names= "all")
        scv.tl.velocity(adata, mode='dynamical')
        adata.write_h5ad(out_path + file)
        
    except Exception as e:
            print('find error {}'.format(e))
