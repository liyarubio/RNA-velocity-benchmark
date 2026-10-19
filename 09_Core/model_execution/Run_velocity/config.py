# Path for conventional scRNA-seq

# raw count 
count_path = 'Input_data/process_10_count.h5ad'
# processed data
input_path = 'Input_data/process_10.h5ad'
output_path = 'result/young2_HSPC_Multi_1_0_with_cellinfo/'
figures_path = output_path + 'figures/'
evals_path = output_path + 'evals/'

# # tfvelo database path
tfvelo_path = 'tools/TFvelo-main'

# multivelo path
mv_rna_path = 'Input_data/process_10_RNA.h5ad'
mv_atac_path = 'Input_data/process_10_ATAC.h5ad'
# phylovelo path
pv_adata_path = 'Input_data/process_10_lineage.h5ad'
pv_pkl_path = 'Input_data/process_10_lineage.pkl'
# Detailed parameters
n_job = 50 # Number of parallel jobs.
device = 'cuda:0' # GPU
seed = 2024 # random seed
embed = 'umap' # Key for embedding
data_cluster = 'CellType' # Key for annotations of observations/cells, a column included in adata.obs
gene_number = 2000 # Gene number
velocity_layer = 'velocity'# Key for velocity to be visualized

# Select algorithms for visualization
Methods =[
    # conventional scRNA-seq
    'velocyto',
    'scvelo',
    'veloae',
    'dynamo',
    'velovae',
    'unitvelo',
    'deepvelo_vae',
    'celldancer',
    'velovi',
    'latentvelo',
    'deepVelo_gcn',
    'stt',
    'tfvelo', # with GRN
    'multivelo', # with scATAC

    'phylovelo',  # with lineage info. # Check

    'topicvelo',
    'bivi', # with count matrix

    'deepkinet', # Check
    'sdevelo',
    'cell2fate', # Check
    # 'topovelo'

]


# methods and title for visualiazation
Methods_name ={
    # conventional scRNA-seq
    'velocyto':'velocyto',
    'scvelo.sto':"scVelo (stochastic)",
    'scvelo.dyn':"scVelo (dynamic)",
    'veloae':"veloAE",
    'dynamo':"Dynamo",
    'velovae':"veloVAE",
    'unitvelo':"UniTVelo",
    'deepvelo_vae':"DeepVelo (VAE-based)",
    'celldancer':"cellDancer",
    'velovi':"veloVI",
    'latentvelo':"LatentVelo",
    'deepvelo_gcn':'DeepVelo (GCN-based)',

    'stt':"STT",
    'tfvelo':"TFvelo", 
    'multivelo':'MultiVelo',
    'phylovelo':"PhyloVelo",

    'topicvelo':'TopicVelo',
    'bivi':'biVI', 
    'deepkinet':'DeepKINET',
    'sdevelo':'SDEvelo',
    'cell2fate':'cell2fate',
    'topovelo':'TopoVelo'
}



# Define transition direction for evaluation
edges = [("HSC", "MPP"),
         ('MPP','CMP'),
         ('MPP','LMPP'),
         ('MPP','GMP'),
         ('CMP','MEP'),
         ('MEP','MKP'),
         ('MEP','EryP'),
         ('GMP','MDP'),
         ('LMPP','CLP'),
         ('CLP','ProB'),
         ('MDP','Mono')]

# palette for visualization
celltype_palette = {
    'HSC':'#d62728',
    'MPP':'#ad494a',
    'CMP':'#1f77b4',
    'GMP':'#aec7e8',
    'MEP':'#ff7f0e',
    'MKP':'#ff9896',
    'EryP':'#ffbb78',
    'MDP':'#8c6d31',
    'Mono':'#c49c94',
    'LMPP':'#aa40fc',
    'CLP':'#c5b0d5',
    'ProB':'#98df8a'}
