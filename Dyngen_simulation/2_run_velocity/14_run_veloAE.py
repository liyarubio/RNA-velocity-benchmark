import pickle, os
import numpy as np
import scvelo as scv
import pandas as pd
import anndata as ad
import torch
from veloproj import *
import scanpy as sc

SEED = 2024
np.random.seed(SEED)

os.chdir('/home/liyr/dygen/dyngen_new/velocity/')


cluster = 'clusters'


torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
np.random.seed(SEED)
torch.backends.cudnn.deterministic = True

parser = get_parser()
args = parser.parse_args(args=['--device', 'cuda:0',
                               '--g-rep-dim', '10', # dimentionality of gene representation (default: 100)
                               '--k-dim', '10']) # dimentionality of the hidden representation Z (default: 100)'


device = torch.device(args.device if args.device.startswith('cuda') and torch.cuda.is_available() else "cpu")
print(device)

print(args)


# %%
def main_AE(args, adata):
    spliced = adata.layers['Ms']
    unspliced = adata.layers['Mu']
    tensor_s = torch.FloatTensor(spliced).to(device)
    tensor_u = torch.FloatTensor(unspliced).to(device)
    tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)
    tensor_v = torch.FloatTensor(adata.layers['velocity']).to(device)

    model = init_model(adata, args, device)

    inputs = [tensor_s, tensor_u]
    xyids = [0, 1]
    if args.use_x:
        inputs.append(tensor_x)

    model = fit_model(args, adata, model, inputs, tensor_v, xyids, device)
    return tensor_s, tensor_u, tensor_x  

# %%
for file  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    try :
        # read data
        folder_path = "/home/liyr/dygen/dyngen_new/velocity/" + file+"/14_veloAE/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        
        h5ad_path = "/home/liyr/dygen/dyngen_new/anndata/"+file+"/anndata.h5ad"
        print(f"read h5ad: {h5ad_path}\n")
        
        print(f"-----------------Simulation: {file}-------------------")
        print("------------------- veloAE ------------------------\n")
        # read data
        adata = sc.read_h5ad(h5ad_path)

        print(adata)
        scv.tl.velocity(adata, mode='stochastic')
        
        # Run VeloAE
        tensor_s, tensor_u, tensor_x = main_AE(args, adata)

        model = init_model(adata, args, device)
        model.load_state_dict(torch.load(args.model_name))
        model = model.to(device)
        model.eval()
        with torch.no_grad():
            x = model.encoder(tensor_x)
            s = model.encoder(tensor_s)
            u = model.encoder(tensor_u)
            
            v = estimate_ld_velocity(s, u, device=device, perc=[5, 95], 
                                        norm=args.use_norm, fit_offset=args.fit_offset_pred, 
                                        use_offset=args.use_offset_pred).cpu().numpy()
            x = x.cpu().numpy()
            s = s.cpu().numpy()
            u = u.cpu().numpy()

        
        adata_new = new_adata(adata, 
                  x, s, u, v, 
                  n_nb_newadata=args.n_nb_newadata, # default 30 
                  X_emb_key='X_pca')
        
        adata_new.layers["velocity"] = adata_new.layers["new_velocity"]
        #scv.tl.velocity_graph(adata_new)
        #scv.tl.velocity_embedding(adata_new,basis = "pca")
        #scv.tl.velocity_embedding(adata_new, basis = "umap")
        
        adata_new.write_h5ad(folder_path + "anndata.h5ad")

    except Exception as e:
        print('find error {}'.format(e))
            


