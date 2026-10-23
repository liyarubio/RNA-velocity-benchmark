import pickle, os
import numpy as np
import scvelo as scv
import scanpy as sc
import scipy
import torch
import time
from veloproj import *
import pandas as pd
scv.settings.verbosity = 1
parser = get_parser()
args = parser.parse_args(args=['--lr', '1e-6',
                               '--n-epochs', '10000',
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
                               '--vis-key', 'X_tsne',
                               '--vis_type_col', 'clusters',
                               '--scv_n_jobs', '26',
                               ])
torch.manual_seed(2024)
torch.cuda.manual_seed(2024)
np.random.seed(2024)
torch.backends.cudnn.deterministic = True

device = torch.device(args.device if args.device.startswith('cuda') and torch.cuda.is_available() else "cpu")
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
    time1 = time.time()
    tensor_s, tensor_u, tensor_x , loss = main_AE(args, adata)
    adata = exp(adata, exp_metrics,out_path)
    time2 = time.time()
    folder_path = out_path+str(file_name) + '.h5ad'
    # i_path = out_path+file_name+'/'+str(file_name)+'Project.svg'
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)
    ad_path = folder_path
    adata.write(ad_path)
    # velocity_new = adata.layers['new_velocity']
    # result_path = folder_path+'/veloAE.npy'
    # np.save(result_path,velocity_new)
    return time2-time1 , loss

data = {
    "clusters" : ['dentategyrus_scv','DentateGyrus_velocyto','endocrinogenesis_day15','human_cd34_bone_marrow','Forebrain'],
    "celltype": ['gastrulation','erythroid_lineage'],
    "Celltype" :['Hindbrain_GABA_Glio'],
    'Main_cell_type' : ['organogenesis_chondrocyte'],
    "Cell_type" :['zebrafish_process']
}
size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
# size = [10000]
time_list = []
size_list = []
out_path = '../time_analyse1/veloAE/'
os.makedirs(out_path, exist_ok=True)
for i in size :
    # i = 1000
    key = 'celltype'
    file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(i)
    args.vis_type_col = key
    args.model_name = str(i) + '_VeloAE'
    args.exp_name = 'CohAE_'+str(i)
    EXP_NAME = args.exp_name
    adata = sc.read_h5ad(file_path)

    # for j in adata.obsm:
    args.vis_key = 'X_umap'
    scv.tl.velocity(adata, vkey='stc_velocity', mode="stochastic")
    scv.tl.velocity_graph(adata, vkey='stc_velocity', n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey='stc_velocity')
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    if adata:
        print('read {}'.format(i))
    # try :
    time_paste , loss = run_veloAE(adata,out_path,i)
    time_list.append(time_paste)
    size_list.append(i)
df = pd.DataFrame()
df['cellnumber'] = size_list
df['time'] = time_list
df.to_csv('../time_analyse1/veloAE/{}time.csv'.format(i))
