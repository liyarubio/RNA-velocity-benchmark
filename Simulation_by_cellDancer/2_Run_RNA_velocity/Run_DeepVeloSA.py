# %%
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
import logging
import vae_from_deepvelo_SA as dv

SEED = 2024
np.random.seed(SEED)
tf.random.set_seed(SEED)


# %%
file_path = './RNA_velocity_result/scvelo_stochastic/'
file_list = os.listdir(file_path)

out_path = './RNA_velocity_result/DeepVeloSA/'
out_file = os.listdir(out_path)

file_list = [item for item in file_list if item not in out_file]
print(file_list)

logging.basicConfig(filename='logs/deepveloSA.log', filemode='a', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# %%
for file in file_list:
    try :
        # read data
        adata = sc.read_h5ad(file_path+ file)
        
        adata.var["velocity_genes"] = True
        print(adata)

        adata_raw = adata.copy()

        # Add noise to data
        X = np.tile(adata_raw.X.A[:, adata.var["velocity_genes"]], (5, 1))
        Y = np.tile(adata.layers["velocity"][:, adata.var["velocity_genes"]], (5, 1))
        noise_sigma = (adata_raw.X.A.std()/70)**2
        X[adata_raw.shape[0]:, :] += \
            np.random.normal(0, noise_sigma, X[adata_raw.shape[0]:, :].shape)

        XYpath = "Sup/" + file + ".npz"
        np.savez(XYpath, X, Y)

        X = np.load(XYpath)["arr_0"]
        Y = np.load(XYpath)["arr_1"]

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
        adata.layers['velocity_dv'] = velocity_deepvelo

        adata.write_h5ad(out_path+ file)

        
    except Exception as e:
        logging.exception('Error processing file %s: %s', file, e)


