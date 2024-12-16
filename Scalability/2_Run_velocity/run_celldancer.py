import scanpy as sc
import scvelo as scv

import os
import sys
import glob
import pandas as pd
import math
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import seaborn as sns
import celldancer as cd
import celldancer.simulation as cdsim
import celldancer.utilities as cdutil
import celldancer.cdplt as cdplt
from celldancer.cdplt import colormap
from celldancer.utilities import export_velocity_to_dynamo
import time
from scipy.spatial.distance import cosine
from sklearn.metrics.pairwise import cosine_similarity


SEED = 2024
np.random.seed(SEED)
size = [1000,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
time_list = []
size_list = []



result_path='../time_analyse1/cellDancer/'
if not os.path.exists(result_path):
    os.makedirs(result_path)
df = pd.DataFrame()

for i in size:
    print('-------run:{}--------'.format(i))
    file_path = '../celldacner_input/'+'{}.csv'.format(str(i))
    df = pd.read_csv(file_path)
    print(df)
    time1 = time.time()
    loss_df, cellDancer_df=cd.velocity(df,n_jobs=26, speed_up = False)
    time2 = time.time()
    paste_time = time2-time1
    cellDancer_df.to_csv(result_path+str(i)+'.csv')
    time_list.append(paste_time)
    size_list.append(i)
df['cellnumber'] = size_list
df['time'] = time_list
df.to_csv('../time_analyse1/celldancer/time.csv'.format(str(i)), index=False)

