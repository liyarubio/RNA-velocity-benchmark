import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
import unitvelo as utv
import tensorflow as tf
import os
os.environ['TF_USE_LEGACY_KERAS'] = 'True'


SEED = 2024
np.random.seed(SEED)
tf.random.set_seed(SEED)

print(tf.config.list_physical_devices('GPU'))

adata = sc.read_h5ad("adata/redeem_young.h5ad")
print(adata)

adata.var['highly_variable'] = True # use all genes
velo_config = utv.config.Configuration()
adata = utv.run_model(adata,label="CellType",config_file=velo_config)

adata.write_h5ad("adata/uniTvelo.h5ad")