# %%
import pickle, os
import numpy as np
import scvelo as scv
import pandas as pd
import anndata as ad
import torch
from veloproj import *
import scanpy as sc
import logging

SEED = 2024
np.random.seed(SEED)

# %%
file_path = './RNA_velocity_result/scvelo_stochastic/'
file_list = os.listdir(file_path)
file_list = [file for file in file_list if file.endswith('.h5ad')]

out_path = './RNA_velocity_result/veloAE/'

out_file = os.listdir(out_path)
file_list = [item for item in file_list if item not in out_file]
print(file_list)

n_jobs = 26

logging.basicConfig(filename='logs/veloAE.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# %%
parser = get_parser()
args = parser.parse_args(args=['--device', 'cuda:0'])

torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
np.random.seed(SEED)
torch.backends.cudnn.deterministic = True

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
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        print(adata)

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
                  X_emb_key='X_umap')
        
        adata_new.write_h5ad(out_path + file)

    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)
            


