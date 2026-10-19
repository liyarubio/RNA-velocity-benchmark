import sdevelo as sv
import numpy as np
import scanpy as sc
import os
import torch
from config import input_path , output_path ,seed, data_cluster ,embed , n_job, device as cfg_device


h5ad_path = input_path

print(f"read h5ad: {h5ad_path}\n")

adata = sc.read_h5ad(h5ad_path)

print("------------------- SDEvelo ------------------------\n")

# Process
#sc.tl.leiden(adata) # need label information
adata.var['highly_variable'] = True # use all genes

np.random.seed(seed)
torch.manual_seed(seed)


emb = embed
args = sv.Config()
args.vis_key = 'x_' + emb
args.seed = seed

# Use device from top-level config.
# sdevelo expects args.cuda_device as an integer CUDA index and uses args.no_cuda to disable CUDA.
cfg_device = str(cfg_device).lower().strip()
if cfg_device.startswith("cuda"):
    # Accept formats like "cuda:3" or "cuda"
    parts = cfg_device.split(":")
    cuda_idx = int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else 0
    if torch.cuda.is_available():
        if torch.cuda.device_count() > 0 and cuda_idx >= torch.cuda.device_count():
            print(
                f"[Warning] cfg_device={cfg_device} out of range (device_count={torch.cuda.device_count()}). "
                f"Falling back to cuda:0.",
                flush=True,
            )
            cuda_idx = 0
        args.cuda_device = cuda_idx
        args.no_cuda = False
    else:
        print("[Warning] CUDA not available, falling back to CPU for SDEvelo.", flush=True)
        args.no_cuda = True
        args.cuda_device = 0
else:
    # Treat anything else (e.g. "cpu") as CPU mode
    args.no_cuda = True
    args.cuda_device = 0

args.scv_n_jobs = n_job
cell_num = adata.X.shape[0]
if hasattr(adata.X, "toarray"):
    adata.X = adata.X.toarray()
args.n_cell = cell_num
args.n_count = 0
model = sv.SDENN(args, adata)
adata = model.train(300)
adata.layers['velocity'] = adata.layers['sde_velocity']

adata.write_h5ad(output_path + "sdevelo.h5ad")


