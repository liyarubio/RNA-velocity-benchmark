import os
import scanpy as sc
import numpy as np
import scvelo as scv
# import unitvelo as utv
import pandas as pd
import math
import pickle
methods = {
    # 'cellDancer':'velocity_S',
    'DeepVelo_GB':'velocity',
    # 'DeepVelo_SA':'DSA_velocity',
    # 'Dynamo_deterministic':'velocity_S',
    # 'Dynamo_stochastic':'velocity_S',
    # 'LatentVelo':'velo_s',
    # 'scvelo':'velocity',
    # 'scvelo_dynamic':'velocity',
    # 'UniTVelo':'velocity',
    # 'VeloAE':'new_velocity',
    # 'velocyto':'velocity',
    # 'VeloVAE':'vae_velocity',
    # 'veloVI':'velocity_s',
    # 'TFvelo':'velocity',
    # 'STT':'velocity'
}

import numpy as np
from sklearn.metrics.pairwise import cosine_similarity


def summary_scores(all_scores):
    """Summarize group scores.

    Args:
        all_scores (dict{str,list}): {group name: score list of individual cells}.

    Returns:
        dict{str,float}: Group-wise aggregation scores.
        float: score aggregated on all samples

    """
    sep_scores = {k:np.mean(s) for k, s in all_scores.items() if s }
    overal_agg = np.mean([s for k, s in sep_scores.items() if s])
    return sep_scores, overal_agg


def keep_type(adata, nodes, target, k_cluster):
    """Select cells of targeted type

    Args:
        adata (Anndata): Anndata object.
        nodes (list): Indexes for cells
        target (str): Cluster name.
        k_cluster (str): Cluster key in adata.obs dataframe

    Returns:
        list: Selected cells.

    """
    return nodes[adata.obs[k_cluster][nodes].values == target]


def cross_boundary_scvelo_probs(adata, k_cluster, cluster_edges, k_trans_g, return_raw=False):
    """Compute Cross-Boundary Confidence Score (A->B).

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        cluster_edges (list of tuples("A", "B")): pairs of clusters has transition direction A->B
        k_trans_g (str): key to the transition graph computed using velocity.
        return_raw (bool): return aggregated or raw scores.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """

    scores = {}
    all_scores = {}
    for u, v in cluster_edges:
        # print(adata.obs)
        # print(k_cluster)
        sel = adata.obs[k_cluster] == u
        # sel = sel.values.tolist()
        # sel = sum(sel, [])
        # print(type(sel))
        nbs = adata.uns['neighbors']['indices'][sel]
        boundary_nodes = map(lambda nodes:keep_type(adata, nodes, v, k_cluster), nbs)
        type_score = [trans_probs.toarray()[:, nodes].mean()
                      for trans_probs, nodes in zip(adata.uns[k_trans_g][sel], boundary_nodes)
                      if len(nodes) > 0]
        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score
    if return_raw:
        return all_scores
    return scores, np.mean([sc for sc in scores.values()])


def cross_boundary_coh(adata, k_cluster, k_velocity, cluster_edges, return_raw=False):
    """Cross-Boundary Velocity Coherence Score (A->B).

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_velocity (str): key to the velocity matrix in adata.obsm.
        cluster_edges (list of tuples("A", "B")): pairs of clusters has transition direction A->B
        return_raw (bool): return aggregated or raw scores.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """
    scores = {}
    all_scores = {}
    for u, v in cluster_edges:
        sel = adata.obs[k_cluster] == u
        nbs = adata.uns['neighbors']['indices'][sel]
        boundary_nodes = map(lambda nodes:keep_type(adata, nodes, v, k_cluster), nbs)

        velocities = adata.layers[k_velocity]
        v_us = velocities[sel]
        type_score = [cosine_similarity(v_us[[ith]], velocities[nodes]).mean()
                      for ith, nodes in enumerate(boundary_nodes)
                      if len(nodes) > 0]
        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values()])

def in_cluster_scvelo_coh(adata, k_cluster, k_confidence, return_raw=False):
    """In-Cluster Confidence Score.

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_confidence (str): key to the column of cell velocity confidence in adata.obs.
        return_raw (bool): return aggregated or raw scores.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """
    clusters = np.unique(adata.obs[k_cluster])
    scores = {}
    all_scores = {}

    for cat in clusters:
        sel = adata.obs[k_cluster] == cat
        type_score = adata.obs[k_confidence][sel].values.tolist()
        scores[cat] = np.mean(type_score)
        all_scores[cat] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([s for _, s in scores.items()])

def cross_boundary_correctness(adata,
                               k_cluster,
                               k_velocity,
                               cluster_edges,
                               return_raw=False,
                               x_emb="X_umap"):
    """Cross-Boundary Direction Correctness Score (A->B)

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_velocity (str): key to the velocity matrix in adata.obsm.
        cluster_edges (list of tuples("A", "B")): pairs of clusters has transition direction A->B
        return_raw (bool): return aggregated or raw scores.
        x_emb (str): key to x embedding for visualization.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """
    scores = {}
    all_scores = {}

    x_emb = adata.obsm[x_emb]
    if x_emb == "X_umap":
        v_emb = adata.obsm['{}_umap'.format(k_velocity)]
    else:
        v_emb = adata.obsm[[key for key in adata.obsm if key.startswith(k_velocity)][0]]

    for u, v in cluster_edges:
        sel = adata.obs[k_cluster] == u
        nbs = adata.uns['neighbors']['indices'][sel] # [n * 30]

        boundary_nodes = map(lambda nodes:keep_type(adata, nodes, v, k_cluster), nbs)
        x_points = x_emb[sel]
        x_velocities = v_emb[sel]

        type_score = []
        for x_pos, x_vel, nodes in zip(x_points, x_velocities, boundary_nodes):
            if len(nodes) == 0: continue

            position_dif = x_emb[nodes] - x_pos
            dir_scores = cosine_similarity(position_dif, x_vel.reshape(1,-1)).flatten()
            type_score.append(np.mean(dir_scores))

        scores[(u, v)] = np.mean(type_score)
        all_scores[(u, v)] = type_score

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values()])

def inner_cluster_coh(adata, k_cluster, k_velocity, return_raw=False):
    """In-cluster Coherence Score.

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_velocity (str): key to the velocity matrix in adata.obsm.
        return_raw (bool): return aggregated or raw scores.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """
    clusters = np.unique(adata.obs[k_cluster])
    scores = {}
    all_scores = {}
    for cat in clusters:
        sel = adata.obs[k_cluster] == cat
        nbs = adata.uns['neighbors']['indices'][sel]
        same_cat_nodes = map(lambda nodes:keep_type(adata, nodes, cat, k_cluster), nbs)
        velocities = adata.layers[k_velocity]
        # nan_mask = np.isnan(velocities)
        # nan_rows, nan_cols = np.where(nan_mask)

        # velocities = velocities[:, ~np.isnan(velocities).any(axis=0)]
        cat_vels = velocities[sel]
        # print(cat_vels)
        # print(velocities[ nodes  for ith, nodes in enumerate(same_cat_nodes) if len(nodes) > 0] )
        # for ith, nodes in enumerate(same_cat_nodes):
        #     cosine_similarity(cat_vels[[ith]], velocities[nodes])
        #     print(cosine_similarity(cat_vels[[ith]], velocities[nodes]))
        cat_score = [cosine_similarity(cat_vels[[ith]], velocities[nodes]).mean()
                     for ith, nodes in enumerate(same_cat_nodes)
                     if len(nodes) > 0]
        all_scores[cat] = cat_score
        scores[cat] = np.mean(cat_score)

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values()])
def evaluate(adata, cluster_edges, k_cluster, k_velocity, x_emb="X_umap", verbose=True):
    """Evaluate velocity estimation results using 5 metrics.

    Args:
        adata (Anndata): Anndata object.
        cluster_edges (list of tuples("A", "B")): pairs of clusters has transition direction A->B
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_velocity (str): key to the velocity matrix in adata.obsm.
        x_emb (str): key to x embedding for visualization.

    Returns:
        dict: aggregated metric scores.

    """
    df_all = {}
    list_mean = []

    # trans_probs = cross_boundary_scvelo_probs(adata, k_cluster, cluster_edges, "{}_graph".format(k_velocity), True)
    # crs_bdr_coh = cross_boundary_coh(adata, k_cluster, k_velocity, cluster_edges, True)
    crs_bdr_crc = cross_boundary_correctness(adata, k_cluster, k_velocity, cluster_edges, True, x_emb)
    ic_coh = inner_cluster_coh(adata, k_cluster, k_velocity, True)
    # ic_scvelo_coh = in_cluster_scvelo_coh(adata, k_cluster, "{}_confidence".format(k_velocity), True)
    # list_mean.append(summary_scores(trans_probs)[1])
    # list_mean.append(summary_scores(crs_bdr_coh)[1])
    list_mean.append(summary_scores(crs_bdr_crc)[1])
    list_mean.append(summary_scores(ic_coh)[1])
    # list_mean.append(summary_scores(ic_scvelo_coh)[1])
    # df_all['trans_probs'] = trans_probs
    # df_all['crs_bdr_coh'] = crs_bdr_coh
    df_all['crs_bdr_crc'] = crs_bdr_crc
    df_all['ic_coh'] = ic_coh
    # df_all['ic_scvelo_coh'] = ic_scvelo_coh
    return list_mean, df_all


def get_global(in_cluster, cross_boundary,cell_number):
    g1 = []
    for cluster0 in in_cluster:
        data_tmp1 = in_cluster[cluster0]
        p1 = len(data_tmp1)/cell_number
        h1 = np.mean(data_tmp1)

        d1 = 0
        count_num = 0
        for edge0 in cross_boundary:
            if cluster0 in edge0:
                ss = cross_boundary[edge0]
                count_num = count_num + 1
                d1 = d1 + np.sign(np.median(ss)) ########## sign
                #d1 <- d1 + np.median(ss) ########## median
        g1.append( p1 * h1 * d1 / count_num)
    g1 = [item for item in g1 if not math.isnan(item)]
    g1 = sum(g1)
    return g1

edges_direct = ['Totipotent','Pluripotent','Multipotent','Oligopotent','Unipotent','Differentiated']
type = 'CytoTRACE2_Potency'
kk = 0
method_eval = {}
emb= "umap"
scvelo_path = '/home/liyr/HuBMAP/RNA_velocity_result/scvelo/'

eval_df = pd.DataFrame()
df_crs_bdr_crc = pd.DataFrame( )
df_ic_coh = pd.DataFrame( )
df_global = pd.DataFrame( )
global_all = pd.read_csv('/home/liyr/HuBMAP/new_hubmap_global_all.csv',sep='\t',index_col=0)
for method in methods:
    # method = 'STT'
    file_path  = '/home/liyr/HuBMAP/RNA_velocity_result/{}/'.format('scvelo')
    file_list = os.listdir(file_path)
    # file_list = ['b23c6f176280ad4c9fe8f74509cb6cf6.h5ad','fb94e862fe33f171de6e07d2d207d0de.h5ad','d0f164331f53f745cfa3e67bf3feb22d.h5ad','27db7c3986c56c8ce03e71779e7ce698.h5ad','a48ab0bf5d8084da24859c4e64336e9c.h5ad','aaedb5272190e61bf557d3b0a1bc591f.h5ad','27db7c3986c56c8ce03e71779e7ce698.h5ad','a48ab0bf5d8084da24859c4e64336e9c.h5ad','aaedb5272190e61bf557d3b0a1bc591f.h5ad']
    file_nolist = pd.read_csv('/home/liyr/HuBMAP/{}_eval.csv'.format(method),sep='\t')
    file_nolist = file_nolist['file_name'].values.tolist()
    xx = [x for x in  file_list if x not in file_nolist]
    vkey = methods[method]
    global_list = []
    ic_coh_list = []
    crs_bdr_crc = []

    file_name = []

    i = 0
    for i in range(len(xx)):
        if xx:
            try :
                print('run {} : {}/{}'.format(method,i,len(xx)))
                file = xx[i]
                adata =sc.read_h5ad(file_path+file)
                adata1 = sc.read_h5ad(scvelo_path+file)
                print(adata)

                if method == 'STT':
                    gene1 = set(adata.var_names)
                    gene2 = set(adata1.var_names)
                    gene = list(gene1 & gene2)
                    adata= adata[:,gene ]
                    adata1 = adata1[:,gene]
                    velocity_adata1 = adata.layers['velocity']
                    adata1.layers['velocity'] = np.zeros(adata1.shape)  # 先为 adata2 创建一个全0的 velocity 层
                    adata1.layers['velocity'] = velocity_adata1
                    adata = adata1.copy()
                if method == 'scvelo_dynamic':
                    arr = adata.layers['velocity']
                    all_nan_columns = np.where(np.isnan(arr).all(axis=0))[0].tolist()
                    if all_nan_columns:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_columns))
                if 'Dynamo' in method:
                    print('yes')
                    adata.layers[vkey] = adata.layers[vkey].toarray()#dynamo
                    arr = adata.layers[vkey]
                    all_nan_columns = np.where(np.isnan(arr).all(axis=0))[0].tolist()
                    if all_nan_columns:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_columns))
                if 'TFvelo' in method:
                    print('---------------yes this is TFvelo------------------')
                    arr = adata.layers[vkey]
                    all_nan_columns = np.where(np.isnan(arr).all(axis=0))[0].tolist()
                    print(all_nan_columns)
                    if all_nan_columns:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_columns))
                cell_number = adata.obs.shape[0]
                # adata.obs['CytoTRACE2_Potency'] = adata1.obs['CytoTRACE2_Potency'].values.tolist()
                print(adata)
                development_type = adata.obs['CytoTRACE2_Potency'].values.tolist()
                development_type = list(set(development_type))
                series = adata.obs['CytoTRACE2_Potency'].value_counts()
                edges = []
                for i in range(len(development_type)):
                    period = development_type[i]
                    index = edges_direct.index(period)
                    new_list = edges_direct[index + 1:]
                    for x in development_type :
                        if x in new_list:
                            a = (period, x)
                            edges.append(a)
                print(edges)
                del adata.uns['neighbors']

                # scv.tl.umap(adata)
                scv.pp.neighbors(adata)
                scv.tl.velocity_graph(adata, basis=emb,vkey=vkey)
                scv.pl.velocity_embedding(adata, basis='umap',vkey = vkey,color='CytoTRACE2_Potency',dpi=150, smooth=0.5,show=False)
                # scv.tl.velocity_confidence(adata, vkey=cluster)
                # scv.pl.velocity_embedding_grid(adata, basis='umap',vkey = 'velocity',color='CytoTRACE2_Potency',dpi=150,smooth=0.5,show=False)
                scv.pl.velocity_embedding_stream(adata, basis='umap',vkey = vkey,color='CytoTRACE2_Potency',dpi=150,smooth=0.5,show=False)
                cluster_edges, k_cluster, k_velocity, x_emb=edges,'CytoTRACE2_Potency',vkey,'X_'+emb
                # print(adata.uns['neighbors'])
                list1, df   = evaluate(adata, cluster_edges, k_cluster, k_velocity, x_emb=x_emb, verbose=True)
                in_cluster = df["ic_coh"]
                cross_boundary =df["crs_bdr_crc"]
                eval_global = get_global(in_cluster,cross_boundary,cell_number)
                global_list.append(eval_global)
                ic_coh_list.append(list1[1])
                crs_bdr_crc.append(list1[0])
                file_name.append(file)
                print(global_list,ic_coh_list,crs_bdr_crc)
            except Exception as e:
                global_list.append(np.nan)
                ic_coh_list.append(np.nan)
                crs_bdr_crc.append(np.nan)
                file_name.append(file)
                print(f"{e}")
                continue
            for index, value in zip(file_name, global_list):
                global_all.loc[index, method] = value
    eval_df = pd.DataFrame()
    eval_df['global'],eval_df['ic_coh'],eval_df['crs_bdr_crc'] ,eval_df['file_name']= global_list , ic_coh_list , crs_bdr_crc,file_name
    # method_eval[method] = eval_df
    eval_df.to_csv('/home/liyr/HuBMAP/{}_eval.csv'.format(method),sep='\t')
    df_global[method] = global_list
    df_global.index = file_name
    df_ic_coh[method] = ic_coh_list
    df_ic_coh.index = file_name
    df_crs_bdr_crc[method] = crs_bdr_crc
    df_crs_bdr_crc.index= file_name
    print(df_global)
global_all.to_csv('/home/liyr/HuBMAP/new_hubmap_global_all.csv',sep='\t')

df_global = pd.DataFrame(df_global).to_csv('/home/liyr/HuBMAP/hubmap_global_all.csv',sep='\t')
df_ic_coh = pd.DataFrame(df_ic_coh).to_csv('/home/liyr/HuBMAP/ic_coh.csv',sep='\t')
df_crs_bdr_crc = pd.DataFrame(df_crs_bdr_crc).to_csv('/home/liyr/HuBMAP/crs_bdr_crc.csv',sep='\t')



