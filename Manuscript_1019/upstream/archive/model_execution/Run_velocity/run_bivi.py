import os, sys
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scvi
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from BIVI import biVI
import torch
from config import input_path, output_path,seed, data_cluster ,embed, n_job
np.random.seed(seed)

# from config import input_path , output_path
def run_biVI(adata):
    if 'feature_id' in adata.var.columns:
        adata.var_names = adata.var['feature_id'].to_list()
    adata.var_names_make_unique()
    # Subset to highly variable genes
    adata = adata[:, adata.var.highly_variable]
    adata_raw = adata
    # separate spliced and unspliced layers
    adata_spliced   = anndata.AnnData(adata.layers['spliced'])
    adata_unspliced = anndata.AnnData(adata.layers['unspliced'])
    adata_spliced.var = adata.var.copy()
    adata_unspliced.var = adata.var.copy()
    # annotate spliced vs unspliced genes
    adata_spliced.var['Spliced']   = True
    adata_unspliced.var['Spliced'] = False
    # add a -u to unspliced gene names
    adata_unspliced.var_names = adata_unspliced.var_names + '-u'
    # NOTICE -- ORDER IS IMPORTANT -- concatenate so that matrices are ordered unspliced, spliced
    adata = anndata.concat([adata_unspliced,adata_spliced],axis=1)
    ## Create a counts layer to use in training
    adata.layers["counts"] = adata.X.copy() # preserve counts
    # Update obs,var
    adata.obs = adata_raw.obs.copy()
    # set up anndata with biVI, specifying training layer
    biVI.biVI.setup_anndata(adata,layer="counts")
    # model arguments
    model_args = {    'n_latent'     : 10,
                      'n_layers'     : 3,
                      'dispersion'   : 'gene',
                      'THETA_IS' : 'NAS_SHAPE',
                      'n_hidden'     : 128,
                      'dropout_rate' :  0.1,
                      'log_variational'    :  True,
                      'latent_distribution':  'normal',
                      }
    # define biVI model
    model = biVI.biVI(adata,mode='Bursty',**model_args)

    max_epochs = 100
    plan_kwargs = {'lr' : 0.001,
                   'n_epochs_kl_warmup' : max_epochs/2,
                   }

    model.train(max_epochs = max_epochs,
                train_size = 0.9,
                check_val_every_n_epoch  = 1,
                plan_kwargs = plan_kwargs)

    latent_representation = model.get_latent_representation(adata)

    adata.obsm['latent_representation'] = latent_representation

    parameters = model.get_likelihood_parameters(adata)
    velocity = parameters['rel_splicing_rate'] * adata_unspliced.X.toarray() - parameters['rel_degradation_rate'] * adata_spliced.X.toarray()
    adata_raw.layers['velocity'] = velocity
    adata_raw.obsm['latent_representation'] = latent_representation

    return adata_raw


print(f"read h5ad: {input_path}\n")

adata = sc.read_h5ad(input_path)

print("------------------- biVI ------------------------\n")

# Process
#sc.tl.leiden(adata) # need label information
adata.var['highly_variable'] = True 
adata = run_biVI(adata)

adata.write_h5ad(output_path +"bivi.h5ad")
