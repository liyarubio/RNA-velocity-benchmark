#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:35:08 2024

@author: liyaru
"""

import pandas as pd
import numpy as np
from itertools import product


clone = pd.read_csv("/media/liyaru/LYR/Diff_change/7_lineage_tracing_multitag/CellTag-multi-2023-main/clone_tables/hsc.rna&atac.r1&2_master_v2.csv")
clone = clone[clone['assay'] == 'rna']
cells_d2 = clone.loc[clone['day'] == 'd2',['cell.bc','clone.id']]
cells_d5 = clone.loc[clone['day'] == 'd5',['cell.bc','clone.id']]


# all combinations
combinations = list(product(cells_d2['cell.bc'], cells_d5['cell.bc']))
df_combinations = pd.DataFrame(combinations, columns=['Cell1', 'Cell2'])

cell_pair = pd.merge(df_combinations,cells_d2,how="left",left_on="Cell1",right_on='cell.bc')
cell_pair = pd.merge(cell_pair,cells_d5,how="left",left_on="Cell2",right_on='cell.bc',suffixes=('', '_cell2'))

cell_pair_true = cell_pair.loc[cell_pair["clone.id"] == cell_pair["clone.id_cell2"],['Cell1', 'Cell2', 'clone.id']]

cell_pair_true.to_csv("/home/liyaru/DATA/8_transition/1_DATA/4_NBT_2023/potential_cell_transfer_pair.csv",index=False)
