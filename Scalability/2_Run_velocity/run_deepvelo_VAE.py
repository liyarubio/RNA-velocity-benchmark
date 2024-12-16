import os

import scvelo as scv
import scanpy as sc
import numpy as np
import sklearn
import scipy
import pandas as pd
import seaborn as sns
import pickle
scv.set_figure_params()
import anndata as ad
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import time
import scipy
#from vae import create_encoder, create_decoder, VAE
tf.config.list_physical_devices('GPU')
# adata = ad.read('../data_2000gene/zebrafish_process.h5ad')
n_gene = 2000
# torch.manual_seed(2024)
np.random.seed(2024)
size = [1000,2000]
#1000, 3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000
class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

def create_encoder(input_size = 1000, latent_dim = 16, verbose=1):
    encoder_inputs = keras.Input(shape=(input_size,))
    x = layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(encoder_inputs)
    z_mean = layers.Dense(latent_dim, name="z_mean", activity_regularizer=keras.regularizers.l1(1e-6))(x)
    z_log_var = layers.Dense(latent_dim, name="z_log_var", activity_regularizer=keras.regularizers.l1(1e-6))(x)
    z = Sampling()([z_mean, z_log_var])
    encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z], name="encoder")
    if verbose == 1:
        encoder.summary()
    return encoder

def create_decoder(output_size = 1000, latent_dim = 16, verbose=1):
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(16, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(latent_inputs)
    x = layers.Dense(64, activation="relu", activity_regularizer=keras.regularizers.l1(1e-6))(x)
    decoder_outputs = layers.Dense(output_size)(x)
    decoder = keras.Model(latent_inputs, decoder_outputs, name="decoder")
    if verbose == 1:
        decoder.summary()
    return decoder
"""
## Define the VAE as a `Model` with a custom `train_step`
"""

class VAE(keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.kl_loss_tracker = keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        x, y = data
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(x)
            reconstruction = self.decoder(z)
            reconstruction_loss = tf.reduce_mean(
                keras.losses.MSE(y, reconstruction)
            )
            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
            total_loss = reconstruction_loss + kl_loss
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
        }

    def call(self, x):
        z_mean, z_log_var, z = self.encoder(x)
        reconstruction = self.decoder(z)
        return reconstruction

    def test_step(self, data):
        x, y = data
        z_mean, z_log_var, z = self.encoder(x)
        reconstruction = self.decoder(z)
        reconstruction_loss = tf.reduce_mean(
            keras.losses.MSE(y, reconstruction)
        )
        kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
        kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
        total_loss = reconstruction_loss + kl_loss
        return {
            "loss": total_loss,
            "reconstruction_loss": reconstruction_loss,
            "kl_loss": kl_loss,
        }

# for type  in data:
time_list = []
size_list = []
os.makedirs('../time_analyse1/deepveloSA/',exist_ok=True)
for file in size:
    try:
        print(file)
        file_path = '../time_analyse/data/'+'gastrulation_{}cell.h5ad'.format(str(file))
        adata = ad.read(file_path)
        adata_raw = adata.copy()
        scv.tl.recover_dynamics(adata, n_jobs=26)
        # adata.var["velocity_genes"] = True
        scv.tl.latent_time(adata)
        scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)
        X = np.tile(adata_raw.X.A[:, adata.var["velocity_genes"]], (5, 1))
        Y = np.tile(adata.layers["velocity"][:, adata.var["velocity_genes"]], (5, 1))
        noise_sigma = (adata_raw.X.A.std()/70)**2
        X[adata_raw.shape[0]:, :] += \
            np.random.normal(0, noise_sigma, X[adata_raw.shape[0]:, :].shape)
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X,
                                                                                    Y,
                                                                                    test_size=0.1,
                                                                                    random_state=42)
        encoder = create_encoder(X.shape[1])
        decoder = create_decoder(X.shape[1])
        tf.config.list_physical_devices('GPU')
        autoencoder = VAE(encoder, decoder)
        opt = keras.optimizers.Adam(learning_rate=0.00001)
        autoencoder.compile(optimizer=opt)
        time1 = time.time()
        autoencoder.fit(X_train, y_train,
                        epochs=100,
                        batch_size=10,
                        shuffle=True,
                        validation_data=(X_test, y_test))


        velocity_deepvelo = autoencoder.predict(adata_raw.X.A[:, adata.var["velocity_genes"]])
        time2 = time.time()
        paste_time = time2 - time1
        time_list.append(paste_time)
        size_list.append(file)
        adata = adata[:, adata.var["velocity_genes"]]
        adata.layers['DSA_velocity'] = velocity_deepvelo
        adata.write_h5ad('../time_analyse1/deepveloSA/{}.h5ad'.format(str(file)))
    except Exception as e:
        print('find error {}'.format(e))
    df = pd.DataFrame()
    df['cellnumber'] = size_list
    df['time'] = time_list
    df.to_csv('../time_analyse1/deepveloSA/{}_time.csv'.format(file), index=False)
