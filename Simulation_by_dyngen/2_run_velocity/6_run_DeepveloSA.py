import scvelo as scv
import scanpy as sc
import numpy as np
import sklearn
import pandas as pd
import seaborn as sns
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import seaborn as sns
import os

import vae_from_deepvelo_SA as dv

SEED = 2024
np.random.seed(SEED)
tf.random.set_seed(SEED)




for file  in ['bifurcating_converging_seed1', 'bifurcating_converging_seed2', 'bifurcating_converging_seed3' ,'bifurcating_cycle_seed1' ,'bifurcating_cycle_seed2' ,'bifurcating_cycle_seed3', 'bifurcating_loop_seed1', 'bifurcating_loop_seed2' ,'bifurcating_loop_seed3', 'bifurcating_seed1', 'bifurcating_seed2' ,'bifurcating_seed3' ,'binary_tree_seed1', 'binary_tree_seed2' ,'binary_tree_seed3', 'branching_seed1' ,'branching_seed2', 'branching_seed3' ,'consecutive_bifurcating_seed1', 'consecutive_bifurcating_seed2' ,'consecutive_bifurcating_seed3' ,'converging_seed1', 'converging_seed2', 'converging_seed3', 'cycle_seed1' ,'cycle_seed2', 'cycle_seed3', 'cycle_simple_seed1', 'cycle_simple_seed2' ,'cycle_simple_seed3', 'disconnected_seed1' ,'disconnected_seed2' ,'disconnected_seed3', 'linear_seed1' ,'linear_seed2' ,'linear_seed3' ,'linear_simple_seed1' ,'linear_simple_seed2' ,'linear_simple_seed3', 'trifurcating_seed1', 'trifurcating_seed2' ,'trifurcating_seed3']:
    
    try :
        # read data
        print(f"-----------------Simulation: {file}-------------------")
        print("------------------- DeepVelo_SA ------------------------\n")
        
        h5ad_path = "/home/liyr/dygen/dyngen_new/anndata/"+file+"/anndata.h5ad"
        folder_path = "/home/liyr/dygen/dyngen_new/velocity/" + file+"/6_DeepveloSA/"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            
        adata = sc.read_h5ad(h5ad_path)
        scv.tl.velocity(adata, mode='stochastic')
        adata.var["velocity_genes"] = True
        # print(adata)

        adata_raw = adata.copy()

        # Add noise to data
        X = np.tile(adata_raw.X.A[:, adata.var["velocity_genes"]], (5, 1))
        Y = np.tile(adata.layers["velocity"][:, adata.var["velocity_genes"]], (5, 1))
        noise_sigma = (adata_raw.X.A.std()/70)**2
        X[adata_raw.shape[0]:, :] += np.random.normal(0, noise_sigma, X[adata_raw.shape[0]:, :].shape)

        XYpath = folder_path+"XY.npz"
        np.savez(XYpath, X, Y)

        X = np.load(XYpath)["arr_0"]
        Y = np.load(XYpath)["arr_1"]
        #Xpath = folder_path + "/DeepVelo_SA_prepropcessed_X.np"
        #Ypath = folder_path + "/DeepVelo_SA_prepropcessed_Y.np"
        #np.save(Xpath, X)
        #np.save(Ypath, Y)

        #X = np.load(Xpath)
        #Y = np.load(Ypath)

        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, 
                                                            Y, 
                                                            test_size=0.1, 
                                                            random_state=2024 # set 2024
                                                            )
        encoder = dv.create_encoder(X.shape[1])
        decoder = dv.create_decoder(X.shape[1])

        autoencoder = dv.VAE(encoder, decoder)
        opt = keras.optimizers.Adam(learning_rate = 0.00005) # default: learning_rate 0.001; in deepvelo tutorial 0.00005
        autoencoder.compile(optimizer=opt)

        es = keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)  # as tutorial figure2 set
        autoencoder.fit(X_train, y_train,
                epochs=100, # as tutorial figure2 set
                batch_size=2, # as tutorial figure2 set
                shuffle=True, # as tutorial figure2 set
                validation_data=(X_test, y_test),
                callbacks=[es])
        
        X = adata_raw.X.A[:, adata.var["velocity_genes"]]
        velocity_deepvelo = autoencoder.predict(X)
        print(velocity_deepvelo.shape)
        adata.layers['velocity'] = velocity_deepvelo
        
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata,basis = "pca")
        scv.tl.velocity_embedding(adata, basis = "umap")
        adata.write_h5ad(folder_path+"/anndata.h5ad")

        
    except Exception as e:
        print('find error {}'.format(e))
