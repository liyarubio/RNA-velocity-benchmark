import os
import scanpy as sc
import numpy as np
import scvelo as scv
import unitvelo as utv
import pandas as pd
import math
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

cell_types = {
    'gastrulation':'stage',

}
methods = {
    'DeepVelo_GB':'velocity',
    'scvelo_stochastic':'velocity',
    'Dynamo_deterministic':'velocity_S',
    'Dynamo_stochastic':'velocity_S',
    'veloVI':'velocity_s',
    'scvelo_dynamical':'velocity',
    'veloAE':'new_velocity',
    'UniTvelo':'velocity',
    'VeloVAE':'vae_velocity_s',
    'latentvelo':'velo_s',
    'Celldancer':'cd_velocity_s',
    'DeepVelo_SA':'DSA_velocity',
    'velocyto':'velocity',
    'TFvelo':'velocity',
    'STT':'velocity_stt'
}

# New
edges = {'dentategyrus_scv':[("Neuroblast", "Granule immature"),
                             ("Granule immature","Granule mature"),
                             ("Radial Glia-like", "Astrocytes"),
                             ("OPC", "OL")],
         'DentateGyrus_velocyto':[('Radial Glia','Astrocyte'),
                                  ("Neuroblast",'Immature Granule'),
                                  ('Immature Granule','Granule'),
                                  ("Neuroblast",'CA')],
         'endocrinogenesis_day15':[('Ngn3 high EP', 'Pre-endocrine'),
                                   ('Pre-endocrine', 'Delta'),
                                   ('Pre-endocrine', 'Beta'),
                                   ('Pre-endocrine', 'Epsilon'),
                                   ('Pre-endocrine', 'Alpha')],
         'human_cd34_bone_marrow':[("HSC",'Precursors'),
                                   ("HSC",'Ery'),
                                   ("HSC",'Mono')],
         'gastrulation':[("E6.5",'E6.75'),
                         ('E6.75','E7.0'),
                         ('E7.0','E7.25'),
                         ('E7.25','E7.5'),
                         ('E7.5','E7.75'),
                         ('E7.75','E8.0'),
                         ('E8.0','E8.25'),
                         ('E8.25','E8.5')],
         'erythroid_lineage':[('Blood progenitors','Erythroid1'),
                              ('Erythroid1','Erythroid2'),
                              ('Erythroid2','Erythroid3')],
         'Forebrain':[('Radial Glia','Neuroblast'),
                      ('Neuroblast','Immature Neuron'),
                      ('Immature Neuron','Neuron')],
         'Hindbrain_GABA_Glio':[('Neural stem cells','Proliferating VZ progenitors'),
                                ('Proliferating VZ progenitors','VZ progenitors'),
                                ('VZ progenitors','Differentiating GABA interneurons'),
                                ('Differentiating GABA interneurons','GABA interneurons'),
                                ('VZ progenitors','Gliogenic progenitors')],
         'organogenesis_chondrocyte':[('Early mesenchyme','Chondroctye progenitors'),
                                      ('Chondroctye progenitors','Chondrocytes & osteoblasts'),
                                      ('Early mesenchyme','Jaw and tooth progenitors')],
         'zebrafish_process':[('Pigment Progenitor','Iridophore'),
                              ('Pigment Progenitor','Melanophore'),
                              ('Schwann Cell Precursor','Schwann Cell')]

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
    'zebrafish_process':'umap'
}

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


    if x_emb == "X_umap":
        v_emb = adata.obsm['{}_umap'.format(k_velocity)]
        print('{}_umap'.format(k_velocity))
    else:
        v_emb = adata.obsm[[key for key in adata.obsm if key.startswith(k_velocity)][0]]
        print('else')
        # print(adata.uns['neighbors']['indices'])
    x_emb = adata.obsm[x_emb]
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
            # print('pos:',position_dif)
            # print('vel:',x_vel)
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
        cat_vels = velocities[sel]
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

    crs_bdr_crc = cross_boundary_correctness(adata, k_cluster, k_velocity, cluster_edges, True, x_emb)
    ic_coh = inner_cluster_coh(adata, k_cluster, k_velocity, True)
    list_mean.append(summary_scores(crs_bdr_crc)[1])
    list_mean.append(summary_scores(ic_coh)[1])
    df_all['crs_bdr_crc'] = crs_bdr_crc
    df_all['ic_coh'] = ic_coh
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


data_path = '../time_analyse1/'
false_list = []
os.makedirs('../time_analyse1/eval/', exist_ok=True)

x = [i for i in methods]
global_df = pd.DataFrame(x, columns=['method'])
print(global_df)

size = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000]
for method in methods:
    df_all = {}
    false_list = []
    global_list = []
    ic_coh_list = []
    crs_bdr_crc = []
    method_list = []
    size_list = []
    df_mean = pd.DataFrame()
    edge =edges['gastrulation']
    type = cell_types['gastrulation']

    for file_name in size:

        try:
            file = str(file_name)

            print(file,method)
            file_path = data_path + method  + '/' + file+ '.h5ad'
            adata = sc.read_h5ad(file_path)
            print(adata)
            vkey = methods[method]
            cluster = methods[method]
            emb = embs['gastrulation']
            type = cell_types['gastrulation']
            # if method == 'STT' :
            if 'scvelo' in method  :
                arr = adata.layers[vkey]
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
            if 'velocyto' in method:
                    print('yes')
                    adata.layers[vkey] = adata.layers[vkey]
                    arr = adata.layers[vkey]
                    all_nan_columns = np.where(np.isnan(arr).all(axis=0))[0].tolist()
                    if all_nan_columns:
                        adata._inplace_subset_var(np.delete(adata.var_names, all_nan_columns))
            if file_name == 'gastrulation':
                adata.obs['stage'] = pd.Categorical(adata.obs['stage'], ordered=True,
                                                    categories=['E6.5','E6.75',
                                                                'E7.0','E7.25','E7.5','E7.75',
                                                                'E8.0','E8.25','E8.5'])
                type =  'stage'
            if file_name == 'DentateGyrus_velocyto':
                replace_dict = {'RadialGlia': 'Radial Glia',
                                'RadialGlia2': 'Radial Glia',
                                'ImmAstro': 'Astrocyte',
                                'ImmGranule1':'Immature Granule',
                                'ImmGranule2':'Immature Granule',
                                'Nbl1':'Neuroblast',
                                'Nbl2':'Neuroblast'
                                }
                adata.obs.replace(replace_dict, inplace=True)
                type = cell_types[file_name]
            if file_name == 'human_cd34_bone_marrow':
                replace_dict = {'HSC_1': 'HSC',
                                'HSC_2': 'HSC',
                                'Ery_1': 'Ery',
                                'Ery_2': 'Ery',
                                'Mono_1': 'Mono',
                                'Mono_2': 'Mono'}
                adata.obs.replace(replace_dict, inplace=True)
                adata.obs['clusters'] = pd.Categorical(adata.obs['clusters'], ordered=True,
                                                       categories=['HSC', 'Ery', 'Mono','Precursors','CLP','DCs','Mega'])
                adata.uns['clusters_colors'] = ['#f781bf', '#a6d854', '#e41a1c', '#984ea3', '#e5c494','#e41a1c', '#a6cee3', '#ff7f00']
                type = cell_types[file_name]

            if file_name == 'erythroid_lineage':
                replace_dict = {'Blood progenitors 1': 'Blood progenitors',
                                'Blood progenitors 2': 'Blood progenitors',
                                }
                adata.obs.replace(replace_dict, inplace=True)
                adata.uns['celltype_colors'] = ['#c9a997', '#C72228', '#f79083', '#EF4E22']
            cell_number = adata.obs.shape[0]
            scv.pp.neighbors(adata)

            scv.tl.umap(adata)
            scv.tl.velocity_graph(adata, basis=emb,vkey=vkey)
            scv.pl.velocity_embedding(adata, basis=emb,vkey = vkey,color=type,dpi=150, smooth=0.5,show=False)
            scv.pl.velocity_embedding_stream(adata, basis=emb,vkey = vkey,color=type,dpi=150,smooth=0.5,show=False)
            cluster_edges, k_cluster, k_velocity, x_emb=edge,type,cluster,'X_'+emb
            print(k_cluster, k_velocity, x_emb)
            list, df   = evaluate(adata, cluster_edges, k_cluster, k_velocity, x_emb=x_emb, verbose=True)
            in_cluster = df["ic_coh"]
            cross_boundary =df["crs_bdr_crc"]
            eval_global = get_global(in_cluster,cross_boundary,cell_number)
            global_list.append(eval_global)
            ic_coh_list.append(list[1])
            crs_bdr_crc.append(list[0])
            size_list.append(file)
            print(global_list)
        except Exception as e:
            # false_list.append(method)
            print(f"An exception occurred while iterating {file_name}:{e}")
            continue
    eval_df = pd.DataFrame()
    eval_df['global'],eval_df['ic_coh'],eval_df['crs_bdr_crc'] ,eval_df['size']= global_list , ic_coh_list , crs_bdr_crc,size_list
    eval_df.to_csv('../time_analyse1/eval/eval_{}.csv'.format(method),index=False)




