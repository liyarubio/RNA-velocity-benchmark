import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.io import mmwrite 
from scipy.sparse import csr_matrix, hstack

import tomotopy as tp
from tomotopy.coherence import Coherence
from tomotopy.utils import Corpus

import scanpy as sc
import scvelo as scv
from scvelo.utils import get_transition_matrix

from topicvelo.tm_utils import cells_to_documents, FastTopics_cluster_assign, TopicGeneFiltering, remove_U, TopicGenesQC, aggregate_clusters
from topicvelo.embed_pl import Plot_Topics, pl_TopTopicGenes, Plot_Genes, Plot_Velocity, mfpt_plot, comparision_stacked_bar_plot, comparision_violin_plot_v2, relative_flux_plot 
from topicvelo.dist_pl import Experimental_JD_Plot, Burst_Simulation_JD_Plot, OS_Analytical_JD_Plot, ExpJD_Cluster_Focus_HeatMap
from topicvelo.transcription_simulation import GeometricBurstTranscription, JointDistributionAnalysis, JointDistributionAnalysis_exp
from topicvelo.inference_tools import Burst_Inference_Gene, KLdivergence, topic_threshold_heuristic_plot, Burst_Inference
from topicvelo.transition_matrix import Combined_Topics_Transitions, get_cells_indices, velocity_graph
from topicvelo.vel_eval_utils import *

from config import input_path , output_path ,seed ,n_job

r_script_path = "Run_velocity/run_topicvelo.R" 

np.random.seed(seed)

adata = sc.read_h5ad(input_path)
scv.pp.filter_and_normalize(adata, n_top_genes=2000)
scv.pp.moments(adata)
adata.var_names_make_unique()
print(adata)

# Ensure output directory exists (many downstream steps assume it does).
os.makedirs(output_path, exist_ok=True)

adata.var['gene_name'] = adata.var.index

genes_S = adata.var['gene_name'].to_list()
genes_U = [g+'_U' for g in genes_S]
gene_names = np.hstack((genes_S, genes_U))

# Construct Spliced/Unspliced matrix
S_U_all = csr_matrix(hstack([adata.layers['spliced'], adata.layers['unspliced']]), dtype=np.int32)
t = pd.DataFrame(S_U_all.toarray())
t.columns = gene_names
t = t.loc[:, ~(t == 0).all()]  
gene_names = t.columns

dense_matrix = t.to_numpy()
S_U = csr_matrix(dense_matrix)

folder_path = output_path
mmwrite(folder_path+'filtered_SU_Counts.mtx', S_U)
pd.DataFrame(gene_names).to_csv(folder_path+"filtered_SU_Genes_names.csv")
pd.DataFrame(adata.obs_names.to_list()).to_csv(folder_path+"filtered_SU_Cells_names.csv")


K = len(set(adata.obs['CellType']))

# Execute the R script for topic modeling
# os.system("Rscript run_topicvelo.R "+folder_path+" "+str(K))
exit_code = os.system(f"Rscript {r_script_path} {output_path} {str(K)}")

if exit_code != 0:
    raise RuntimeError(
        f"Rscript topicvelo failed (exit_code={exit_code}). "
        f"Check {r_script_path} logs and that R + dependencies are installed."
    )


# Load the cell weights from fastTopics, skipping the header and the first column (junk/index)
expected_cellweights = folder_path + f'HH_fastTopics_fit_CellWeights_k={K}.csv'
expected_de_postmean = folder_path + f'HH_de_postmean_k={K}.csv'
expected_de_lfsr = folder_path + f'HH_de_lfsr_k={K}.csv'
expected_de_z = folder_path + f'HH_de_z_k={K}.csv'

missing = [p for p in [expected_cellweights, expected_de_postmean, expected_de_lfsr, expected_de_z] if not os.path.exists(p)]
if missing:
    raise FileNotFoundError(
        "Expected topicvelo output files are missing after successful Rscript run:\n"
        + "\n".join(missing)
    )

cellWeights = np.genfromtxt(folder_path+'HH_fastTopics_fit_CellWeights_k='+str(K)+'.csv', delimiter=",", skip_header=1)[:,1:]

# Add topic weights
FastTopics_cluster_assign(adata, cellWeights, t_type = 'fastTopics')


# Load DE statistics (postmean, lfsr, z-score) from fastTopics, skipping header and first column
de_postmean = np.genfromtxt(folder_path+'HH_de_postmean_k='+str(K)+'.csv', delimiter=",", skip_header=1)[:,1:]
de_lfsr = np.genfromtxt(folder_path+'HH_de_lfsr_k='+str(K)+'.csv', delimiter=",", skip_header=1)[:,1:]
de_z = np.genfromtxt(folder_path+'HH_de_z_k='+str(K)+'.csv', delimiter=",", skip_header=1)[:,1:]


# Select the topic-specific genes
ttgs, ttg_up, ttg_down = TopicGeneFiltering(de_postmean, de_lfsr, 
                                            lfc_up_th = 0.5, lfc_down_th=-0.5, 
                                            lfsr_up_th = 0.001, lfsr_down_th = 0.001)


t = pd.DataFrame(S_U.toarray())
t.columns = gene_names
t = t.loc[:, ~(t == 0).all()]  
gene_names = t.columns

top_genes = remove_U(ttgs, gene_names)
top_genes_up = remove_U(ttg_up, gene_names)
top_genes_down = remove_U(ttg_down, gene_names)

# adata.uns['top_genes'] = top_genes
# scv.pp.neighbors(adata)


topics = [ K-1 ]
topics_to_use_percentile_thresholds = [35]

try:
    combined_TM, ttc_indices = Combined_Topics_Transitions(adata,
                                                            # topics=None,  topic_weights_th_percentile=None,
                                                            topics = topics,
                                                            topic_weights_th_percentile = topics_to_use_percentile_thresholds,
                                                            velocity_type = 'stochastic',
                                                            recompute = True, recompute_matrix = True, steady_state_perc = 95)
except Exception as e:
    print(e)

adata.write(output_path + "topicvelo.h5ad")


### Move files
import os
import shutil
import glob

def organize_topicvelo_results():
    # 1. Define paths
    base_result_dir = output_path
    target_dir_name = 'topicvelo'
    target_path = os.path.join(base_result_dir, target_dir_name)

    # 2. Create ./result/topicvelo directory if it does not exist
    if not os.path.exists(target_path):
        try:
            os.makedirs(target_path)
            print(f"[Info] Created directory: {target_path}")
        except OSError as e:
            print(f"[Error] Failed to create directory {target_path}: {e}")
            return

    # 3. Define file extensions to move
    extensions = ['*.rds', '*.csv', '*.mtx']

    print(f"[Info] Moving .rds, .csv, .mtx files to {target_path}...")

    # 4. Iterate and move files
    for ext in extensions:
        # Search for matching files in the base_result_dir
        # Note: This assumes files are generated in the root of base_result_dir.
        # If files are generated in the current script execution directory, change base_result_dir below to '.'
        files_to_move = glob.glob(os.path.join(base_result_dir, ext))

        for file_path in files_to_move:
            file_name = os.path.basename(file_path)
            dst_path = os.path.join(target_path, file_name)

            # Prevent moving the file to itself (safety check, although glob usually doesn't recurse)
            if os.path.abspath(file_path) == os.path.abspath(dst_path):
                continue

            try:
                shutil.move(file_path, dst_path)
                print(f"   -> Moved: {file_name}")
            except Exception as e:
                print(f"   [Error] Could not move {file_name}: {e}")

# Execute the organization function
organize_topicvelo_results()