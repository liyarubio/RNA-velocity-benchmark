import scanpy as sc
import pandas as pd


size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]

adata = sc.read_h5ad('../data_2000gene/gastrulation.h5ad')

# stratified sampling
strata = adata.obs['cell_type']
strata_groups = strata.groupby(strata.columns.tolist())

for i in size:
    sampled_adata = []
    for group_name, group_indices in strata_groups.indices.items():
        sample_indices = group_indices[:i]
        sampled_adata.append(adata[sample_indices, :])
    sampled_adata = sc.concat(sampled_adata)
    sampled_adata.write_h5ad('../runtime_test/cell_splite/gastrulation_{}cell.h5ad'.format(str(i)))