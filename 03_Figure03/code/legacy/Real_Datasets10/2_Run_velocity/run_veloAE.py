import pickle, os
import numpy as np
import scvelo as scv
import scanpy as sc
import scipy
import torch
import time
import tqdm
from veloproj import *
import pandas as pd
scv.settings.verbosity = 1
parser = get_parser()
args = parser.parse_args(args=['--lr', '1e-6',
                               '--n-epochs', '20000',
                               '--g-rep-dim', '100',
                               '--k-dim', '100',
                               '--model-name', 'humanbonemarrow_model.cpt',
                               '--exp-name', 'CohAE_humanbonemarrow',
                               '--device', 'cuda:0',
                               '--nb_g_src', 'X',
                               '--ld_nb_g_src', 'X',
                               '--gumbsoft_tau', '1',
                               '--n_raw_gene', '2000',
                               '--n_conn_nb', '30',
                               '--n_nb_newadata', '30',
                               '--aux_weight', '1',
                               '--fit_offset_train', 'false',
                               '--fit_offset_pred', 'true',
                               '--use_offset_pred', 'true',
                               '--gnn_layer', 'GAT',
                               '--vis-key', 'umap',
                               '--vis_type_col', 'clusters',
                               '--scv_n_jobs', '26',
                               ])
torch.manual_seed(args.seed)
torch.cuda.manual_seed(args.seed)
np.random.seed(args.seed)
torch.backends.cudnn.deterministic = True

device = torch.device(args.device if args.device.startswith('cuda') and torch.cuda.is_available() else "cpu")
print('device:{}'.format(device))
exp_metrics = {}

def main_AE(args, adata):
    spliced = adata.layers['Ms']
    unspliced = adata.layers['Mu']
    tensor_s = torch.FloatTensor(spliced).to(device)
    tensor_u = torch.FloatTensor(unspliced).to(device)
    tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)
    tensor_v = torch.FloatTensor(adata.layers['stc_velocity']).to(device)

    model = init_model(adata, args, device)

    inputs = [tensor_s, tensor_u]
    xyids = [0, 1]
    if args.use_x:
        inputs.append(tensor_x)

    model ,loss = fit_model(args, adata, model, inputs, tensor_v, xyids, device)
    return tensor_s, tensor_u, tensor_x,loss



def exp(adata, exp_metrics,path):
    model = init_model(adata, args, device)
    model.load_state_dict(torch.load(args.model_name))
    model = model.to(device)
    model.eval()
    with torch.no_grad():
        spliced = adata.layers['Ms']
        unspliced = adata.layers['Mu']
        tensor_s = torch.FloatTensor(spliced).to(device)
        tensor_u = torch.FloatTensor(unspliced).to(device)
        tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)

        x = model.encoder(tensor_x)
        s = model.encoder(tensor_s)
        u = model.encoder(tensor_u)

        v = estimate_ld_velocity(s, u, device=device, perc=[5, 95],
                                 norm=args.use_norm, fit_offset=args.fit_offset_pred,
                                 use_offset=args.use_offset_pred).cpu().numpy()
        x = x.cpu().numpy()
        s = s.cpu().numpy()
        u = u.cpu().numpy()

    adata = new_adata(adata, x, s, u, v, g_basis=args.ld_nb_g_src, n_nb_newadata=args.n_nb_newadata,
                      X_emb_key=args.vis_key)
    # scv.tl.velocity_graph(adata, vkey='new_velocity')
    # scv.pl.velocity_embedding_stream(adata,
    #                                  vkey="new_velocity", basis=args.vis_key, color=[args.vis_type_col],
    #                                  title="Project Original Velocity into Low-Dim Space", smooth=0.5,
    #                                  dpi=150,
    #                                  save=path)
    # scv.tl.velocity_confidence(adata, vkey='new_velocity')
    # exp_metrics['Cohort AutoEncoder'] = evaluate(adata, cluster_edges, args.vis_type_col, "new_velocity",
    #                                              x_emb=args.vis_key)
    return(adata)
def run_veloAE(adata,out_path,file_name):
    tensor_s, tensor_u, tensor_x , loss = main_AE(args, adata)
    adata = exp(adata, exp_metrics,out_path)
    folder_path = out_path+file_name
    ad_path = folder_path + '/VeloAE.h5ad'
    adata.write_h5ad(ad_path)
    return loss

cell_types = {
    'endocrinogenesis_day15':'clusters',
    'DentateGyrus_velocyto':'clusters',
    'human_cd34_bone_marrow':'clusters',
    'Forebrain':'clusters',
    'dentategyrus_scv':'clusters',
    'Hindbrain_GABA_Glio':'Celltype',
    'organogenesis_chondrocyte':'Main_cell_type',
    'erythroid_lineage':'celltype',
    'zebrafish_process':'Cell_type',
    'gastrulation':'stage',
    'MultiVelo_10X_multiome_mouse_brain':'celltype'

}
embs = {
    'dentategyrus_scv':'umap',
        'DentateGyrus_velocyto':'tsne',
        'endocrinogenesis_day15':'umap',
        'human_cd34_bone_marrow':'tsne',
        'Forebrain':'umap',
        'gastrulation':'umap',
        'erythroid_lineage':'umap',
        'Hindbrain_GABA_Glio':'tsne',
    'organogenesis_chondrocyte':'umap',
    'zebrafish_process':'umap',
    'MultiVelo_10X_multiome_mouse_brain':'umap'
}

SEED = 2024
np.random.seed(SEED)
torch.manual_seed(SEED)
kk = 0
file_path = '../real_data/'
out_path = '../all_result/'
for file in cell_types:

    print('--------run veloAE for : {}--------'.format(file))
    args.vis_type_col = cell_types[file]
    args.model_name = file + '_VeloAE.cpt'
    args.exp_name = 'CohAE_'+file
    args.vis_key = embs[file]
    EXP_NAME = args.exp_name
    adata = sc.read_h5ad(file_path+file+'.h5ad')
    scv.pp.remove_duplicate_cells(adata)
    scv.utils.show_proportions(adata)
    scv.pl.proportions(adata, dpi=350)
    scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=args.n_raw_gene)
    scv.pp.moments(adata)
    args.vis_key = "X_" + embs[file]
    scv.tl.velocity(adata, vkey='stc_velocity', mode="stochastic")
    scv.tl.velocity_graph(adata, vkey='stc_velocity', n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey='stc_velocity')

    if adata:
        print('run : {}'.format(file))
    # try :
    print(adata)
    loss= run_veloAE(adata,out_path,file)
