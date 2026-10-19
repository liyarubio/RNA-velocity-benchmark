import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('AGG') 
import os, sys
import random
import traceback
import time
import faulthandler
import psutil  # Requires psutil: pip install psutil
from scipy import sparse

# 1. Enable low-level error handling (for segfaults/C crashes)
faulthandler.enable()

# Force stdout flushing
sys.stdout.reconfigure(line_buffering=True)

from config import input_path, output_path, seed, data_cluster, embed, n_job, tfvelo_path

original_cwd = os.getcwd()
# Resolve TFvelo path robustly (works even if script is started outside project root).
tfvelo_abs_path = tfvelo_path if os.path.isabs(tfvelo_path) else os.path.abspath(
    os.path.join(original_cwd, tfvelo_path)
)
sys.path.insert(1, tfvelo_abs_path)
import TFvelo as TFv

np.random.seed(seed)
random.seed(seed)

def print_memory_usage(step_name):
    """Print memory usage of the current process"""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"  [Memory] {step_name}: {mem_info.rss / 1024 / 1024:.2f} MB", flush=True)

def check_database_files(base_path, databases):
    if isinstance(databases, str):
        databases = databases.split()
    missing = []
    if 'all' in databases:
        p = os.path.join(base_path, "data/TF_names_v_1.01.txt")
        if not os.path.exists(p): missing.append(p)
    if 'ENCODE' in databases:
        p = os.path.join(base_path, "ENCODE/processed/")
        if not os.path.exists(p): missing.append(p)
    if 'ChEA' in databases:
        p = os.path.join(base_path, "ChEA/ChEA_2016.txt")
        if not os.path.exists(p): missing.append(p)
    
    if missing:
        print("\n[Error] Missing database files! TFvelo will fail.", flush=True)
        return False
    return True

def run_tfvelo_robust(adata, key):
    # --- Preprocessing ---
    print("  [Step 1] Preprocessing layers...", flush=True)
    print_memory_usage("Start")

    if key not in adata.obs:
        raise KeyError(f"Cluster key '{key}' not found in adata.obs")
    adata.obs['clusters'] = adata.obs[key].copy()

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise KeyError("adata.layers must contain both 'spliced' and 'unspliced'")

    spliced = adata.layers["spliced"]
    unspliced = adata.layers["unspliced"]
    if sparse.issparse(spliced) and sparse.issparse(unspliced):
        # Keep sparse to avoid large dense copies and OOM.
        adata.layers["total"] = (spliced + unspliced).tocsr()
    else:
        adata.layers["total"] = np.asarray(spliced) + np.asarray(unspliced)
        
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    TFv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    adata.X = adata.layers["total"].copy()
    TFv.pp.moments(adata, n_pcs=30, n_neighbors=15)
    
    print_memory_usage("After Moments")
    
    # --- Critical Step: Load Database ---
    print("  [Step 2] Calculating TF-Gene Correlations...", flush=True)
    t0 = time.time()
    databases_to_use = 'ENCODE ChEA' 
    TFv.pp.get_TFs(adata, databases=databases_to_use)
    print(f"           -> Finished in {(time.time()-t0)/60:.1f} min.", flush=True)
    
    adata.uns['genes_pp'] = np.array(adata.var_names)
    print_memory_usage("Before Recover Dynamics")
    
    # --- Recover Dynamics ---
    # Dynamically adjust concurrency: if memory usage is high, force reduce concurrency
    current_mem_gb = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
    safe_n_jobs = n_job
    
    # Simple protection strategy: if current main process exceeds 10GB and n_job > 4, force downgrade
    if current_mem_gb > 10 and n_job > 4:
        print(f"  [Warning] High memory usage ({current_mem_gb:.1f} GB). Reducing n_jobs to 4 to prevent crash.", flush=True)
        safe_n_jobs = 4
    
    print(f"  [Step 3] Recovering Dynamics (n_jobs={safe_n_jobs})...", flush=True)
    
    try:
        TFv.tl.recover_dynamics(adata, 
                                n_jobs=safe_n_jobs, 
                                backend="loky",  
                                max_iter=20, 
                                var_names='all',
                                WX_method="lsq_linear", WX_thres=20, max_n_TF=99, 
                                n_top_genes=2000, fit_scaling=True, use_raw=0, 
                                init_weight_method="correlation", n_time_points=1000)
    except Exception as e:
        print(f"  [Error] Multiprocessing failed: {e}", flush=True)
        print("  -> Retrying with SINGLE CORE (n_jobs=1)...", flush=True)
        TFv.tl.recover_dynamics(adata, n_jobs=1, max_iter=20, var_names='all',
                               WX_method="lsq_linear", WX_thres=20, max_n_TF=99, 
                               n_top_genes=2000, fit_scaling=True, use_raw=0, 
                               init_weight_method="correlation", n_time_points=1000)
    return

# ================= Main Program =================

print(f"Reading input: {input_path}", flush=True)
adata = sc.read_h5ad(input_path)

if not os.path.exists(output_path):
    os.makedirs(output_path, exist_ok=True)

print(f"Switching working directory to: {tfvelo_abs_path}", flush=True)
os.chdir(tfvelo_abs_path) 

databases_to_use = "ENCODE ChEA"
if not check_database_files(os.getcwd(), databases_to_use):
    os.chdir(original_cwd)
    sys.exit(1)

try:
    run_tfvelo_robust(adata, data_cluster)
    
    n_cells = adata.shape[0]
    if 'fit_scaling_y' in adata.var.keys():
        expanded_scaling_y = np.expand_dims(np.array(adata.var['fit_scaling_y']), 0).repeat(n_cells, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            adata.layers['velocity'] = np.nan_to_num(adata.layers['velo_hat'] / expanded_scaling_y)
    elif "velo_hat" in adata.layers:
        # Fallback if scaling is unavailable in this TFvelo run.
        adata.layers["velocity"] = np.array(adata.layers["velo_hat"])
    else:
        raise KeyError("Neither 'fit_scaling_y' nor 'velo_hat' found after TFvelo run.")
    
    print(f"Restoring working directory to: {original_cwd}", flush=True)
    os.chdir(original_cwd)

    save_file = output_path + 'tfvelo.h5ad'
    print(f"Saving result to: {save_file}", flush=True)
    adata.write_h5ad(save_file)
    print("TFvelo finished successfully!", flush=True)

except Exception as e:
    print("\n!!! TFvelo Failed !!!", flush=True)
    traceback.print_exc()
    os.chdir(original_cwd)
    sys.exit(1)