import os, sys
# Get the directory where the current script (run_stt.py) is located
current_script_dir = os.path.dirname(os.path.abspath(__file__))
# Get the project root directory (parent directory of Run_velocity)
project_root_dir = os.path.dirname(current_script_dir)
# Add the project root directory to Python search path
sys.path.append(project_root_dir)

from tools import sctt as st
from tools import networks
import scanpy as sc
import numpy as np
import traceback # Used for printing detailed error stack traces

from config import input_path, output_path, seed, data_cluster, embed

np.random.seed(seed)

input_file = input_path

print(f"Reading file: {input_file}")
adata = sc.read_h5ad(input_file)
# Ensure unique var names to avoid AnnData warnings and downstream errors
adata.var_names_make_unique()
adata_orig = adata.copy()

adata.obs['attractor'] = adata.obs[data_cluster].values

# Initialize n_states (use nunique for pandas compatibility)
initial_n_states = int(adata.obs[data_cluster].nunique())
n_states = initial_n_states
min_states = 2  # STT needs at least 2 states; ensure loop runs even with 1 cluster
n_states = max(n_states, min_states)

adata_aggr = None
success = False

print(f"Initial n_states based on clusters: {n_states}")

# --- Optimization start: Add automatic retry loop ---
while n_states >= min_states:
    try:
        print(f"--- Attempting STT with n_states={n_states} ---")
        # Attempt to run the core function
        adata_aggr = st.dynamical_iteration(adata, return_aggr_obj=True, n_states=n_states)
        success = True
        print(f"Successfully finished STT with n_states={n_states}")
        break # If successful, break the loop

    except Exception as e:
        # Catch errors (including ValueError: Discretizing leads to a cluster...)
        print(f"Warning: Failed with n_states={n_states}. Error: {e}")
        
        # Only attempt to reduce n_states if it is greater than the minimum
        if n_states > min_states:
            n_states -= 1
            print(f"Reducing n_states to {n_states} and retrying...")
        else:
            print("Reached minimum n_states. Cannot reduce further.")
            # Print full error stack trace for debugging
            traceback.print_exc()
            raise e # Re-raise exception, terminate program

if not success:
    print("STT failed to converge.")
    sys.exit(1)
# --- Optimization end ---

# Subsequent processing logic remains unchanged
if adata_aggr is not None:
    # Ensure index alignment (adding a check is recommended)
    if len(adata_aggr.obs_names) == len(adata.obs_names):
        adata_aggr.obs_names = adata.obs_names
    else:
        print("Warning: Resulting adata_aggr has different cell number than input.")

    gene1 = set(adata_aggr.var_names)
    gene2 = set(adata_orig.var_names)
    gene = list(gene1 & gene2)

    adata_aggr = adata_aggr[:, gene]
    adata_orig = adata_orig[:, gene]

    # Ensure index consistency again
    adata_aggr.obs_names = adata_orig.obs_names
    adata_aggr.obs[data_cluster] = adata_orig.obs[data_cluster]
    
    emb_key = 'X_' + embed
    if emb_key in adata_orig.obsm.keys():
        adata_aggr.obsm[emb_key] = adata_orig.obsm[emb_key]
    
    save_path = output_path + 'stt.h5ad'
    print(f"Saving results to {save_path}")
    adata_aggr.write(save_path)