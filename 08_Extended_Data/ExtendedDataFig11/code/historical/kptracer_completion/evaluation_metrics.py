"""Robust RNA velocity evaluation metrics used by VeloGIF."""

import math

import numpy as np
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity


def _finite_mean(values):
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    return float(array.mean()) if array.size else math.nan


def summary_scores(all_scores):
    separated = {key: _finite_mean(scores) for key, scores in all_scores.items() if len(scores)}
    return separated, _finite_mean(separated.values())


def _neighbor_indices(adata):
    neighbors = adata.uns.get("neighbors", {})
    if "indices" in neighbors:
        return np.asarray(neighbors["indices"])
    connectivities = adata.obsp.get("connectivities")
    if connectivities is None:
        raise KeyError("No neighbor indices or connectivities found")
    matrix = connectivities.tocsr() if sparse.issparse(connectivities) else sparse.csr_matrix(connectivities)
    return [matrix.indices[matrix.indptr[i]:matrix.indptr[i + 1]] for i in range(matrix.shape[0])]


def _dense_rows(matrix, rows):
    selected = matrix[rows]
    return selected.toarray() if sparse.issparse(selected) else np.asarray(selected)


def cross_boundary_correctness(adata, cluster_key, velocity_key, edges, embedding_key="X_umap"):
    labels = np.asarray(adata.obs[cluster_key])
    coordinates = np.asarray(adata.obsm[embedding_key])
    embedding_name = embedding_key[2:] if embedding_key.startswith("X_") else embedding_key
    velocity_embedding_key = "%s_%s" % (velocity_key, embedding_name)
    if velocity_embedding_key not in adata.obsm:
        candidates = [key for key in adata.obsm if key.startswith(velocity_key)]
        if not candidates:
            raise KeyError("No velocity embedding found for %s" % velocity_key)
        velocity_embedding_key = candidates[0]
    velocities = np.asarray(adata.obsm[velocity_embedding_key])
    neighbors = _neighbor_indices(adata)
    scores = {}
    for source, target in edges:
        per_cell = []
        for index in np.flatnonzero(labels == source):
            target_neighbors = np.asarray(neighbors[index], dtype=int)
            target_neighbors = target_neighbors[labels[target_neighbors] == target]
            if not target_neighbors.size:
                continue
            direction = coordinates[target_neighbors] - coordinates[index]
            per_cell.append(float(cosine_similarity(direction, velocities[index].reshape(1, -1)).mean()))
        scores[(source, target)] = per_cell
    return scores


def inner_cluster_coherence(adata, cluster_key, velocity_key):
    labels = np.asarray(adata.obs[cluster_key])
    velocities = adata.layers[velocity_key]
    neighbors = _neighbor_indices(adata)
    scores = {}
    for cluster in np.unique(labels):
        per_cell = []
        for index in np.flatnonzero(labels == cluster):
            same = np.asarray(neighbors[index], dtype=int)
            same = same[labels[same] == cluster]
            if same.size:
                per_cell.append(float(cosine_similarity(_dense_rows(velocities, [index]), _dense_rows(velocities, same)).mean()))
        scores[cluster] = per_cell
    return scores


def evaluate(adata, edges, cluster_key, velocity_key, embedding_key="X_umap"):
    cross = cross_boundary_correctness(adata, cluster_key, velocity_key, edges, embedding_key)
    inner = inner_cluster_coherence(adata, cluster_key, velocity_key)
    return {
        "CBDir": summary_scores(cross)[1],
        "ICCoh": summary_scores(inner)[1],
        "raw_cross_boundary": cross,
        "raw_inner_cluster": inner,
    }


def global_directional_correctness(inner_scores, cross_scores, cell_count):
    if not cell_count:
        return math.nan
    terms = []
    for cluster, values in inner_scores.items():
        if not values:
            continue
        relevant = [scores for edge, scores in cross_scores.items() if cluster in edge and len(scores)]
        if not relevant:
            continue
        direction = sum(np.sign(np.median(scores)) for scores in relevant) / len(relevant)
        terms.append((len(values) / cell_count) * _finite_mean(values) * direction)
    return float(sum(value for value in terms if np.isfinite(value))) if terms else math.nan
