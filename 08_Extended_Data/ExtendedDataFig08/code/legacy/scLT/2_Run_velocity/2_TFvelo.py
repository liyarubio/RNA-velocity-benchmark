# '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('AGG')
import os, sys
sys.path.insert(1, './TFvelo-master/')
import TFvelo as TFv


n_jobs = 26
#
def check_data_type(adata):
    for key in list(adata.var):
        if adata.var[key][0] in ['True', 'False']:
            adata.var[key] = adata.var[key].map({'True': True, 'False': False})
    return
def data_type_tostr(adata, key):
    if key in adata.var.keys():
        if adata.var[key][0] in [True, False]:
            adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    return
def run_tfvelo(adata,key):
    adata.obs['clusters'] = adata.obs[key].copy()
    adata.layers["total"] = adata.layers["spliced"].todense() + adata.layers["unspliced"].todense()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    n_cells, n_genes = adata.X.shape
    sc.pp.filter_genes(adata, min_cells=int(n_cells/50))
    sc.pp.filter_cells(adata, min_genes=int(n_genes/50))
    TFv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000) #include the following steps
    adata.X = adata.layers["total"].copy()
    gene_names = []
    for tmp in adata.var_names:
        gene_names.append(tmp.upper())
    adata.var_names = gene_names
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    TFv.pp.moments(adata, n_pcs=30, n_neighbors=15)
    TFv.pp.get_TFs(adata, databases='ENCODE ChEA')
    print(adata)
    adata.uns['genes_pp'] = np.array(adata.var_names)
    n_jobs = 28
    flag = TFv.tl.recover_dynamics(adata, n_jobs=n_jobs, max_iter=20, var_names='all',
                                   WX_method = "lsq_linear", WX_thres=20, max_n_TF=99, n_top_genes=2000,
                                   fit_scaling=True, use_raw=0, init_weight_method= "correlation",
                                   n_time_points=1000)
    return

file_list = {'multitag':'celltype','scNT':'cell_type'}

for i in file_list:
    file_path1= '/home/liyr/{}/adata_graph/scvelo_stochastic.h5ad'.format(i)
    file_result = '/home/liyr/{}/adata_graph/TFvelo.h5ad'.format(i)
    adata = sc.read_h5ad(file_result)
    print(adata)
   # key = file_list[i]
  #  run_tfvelo(adata,key)
   # adata.write(file_result)
