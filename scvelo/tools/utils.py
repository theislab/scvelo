from scipy.sparse import csr_matrix, issparse
import numpy as np
import warnings


def mean(x, axis=0):
    return x.mean(axis).A1 if issparse(x) else x.mean(axis)


def make_dense(X):
    return X.A if issparse(X) and X.ndim==2 else X.A1 if issparse(X) else X


def prod_sum_obs(A, B):
    """dot product and sum over axis 0 (obs) equivalent to np.sum(A * B, 0)
    """
    return A.multiply(B).sum(0).A1 if issparse(A) else np.einsum('ij, ij -> j', A, B)


def prod_sum_var(A, B):
    """dot product and sum over axis 1 (var) equivalent to np.sum(A * B, 1)
    """
    return A.multiply(B).sum(1).A1 if issparse(A) else np.einsum('ij, ij -> i', A, B)


def norm(A):
    """computes the L2-norm along axis 1 (e.g. genes or embedding dimensions) equivalent to np.linalg.norm(A, axis=1)
    """
    return np.sqrt(A.multiply(A).sum(1).A1) if issparse(A) else np.sqrt(np.einsum('ij, ij -> i', A, A))


def vector_norm(x):
    """computes the L2-norm along axis 1 (e.g. genes or embedding dimensions) equivalent to np.linalg.norm(A, axis=1)
    """
    return np.sqrt(np.einsum('i, i -> ', x, x))


def R_squared(residual, total):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = np.ones(residual.shape[1]) - prod_sum_obs(residual, residual) / prod_sum_obs(total, total)
    r2[np.isnan(r2)] = 0
    return r2


def cosine_correlation(dX, Vi):
    dX -= dX.mean(-1)[:, None]
    Vi_norm = vector_norm(Vi)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.zeros(dX.shape[0]) if Vi_norm == 0 else np.einsum('ij, j', dX, Vi) / (norm(dX) * Vi_norm)[None, :]
    return result


def normalize(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X.multiply(csr_matrix(1. / np.abs(X).sum(1))) if issparse(X) else X / X.sum(1)
    return X


def scale(X, min=0, max=1):
    X = X - X.min() + min
    X = X / X.max() * max
    return X


def get_indices(dist):
    n_neighbors = (dist > 0).sum(1).min()
    rows_idx = np.where((dist > 0).sum(1) > n_neighbors)[0]

    for row_idx in rows_idx:
        col_idx = dist[row_idx].indices[n_neighbors:]
        dist[row_idx, col_idx] = 0

    dist.eliminate_zeros()

    indices = dist.indices.reshape((-1, n_neighbors))
    return indices, dist


def get_iterative_indices(indices, index, n_recurse_neighbors, max_neighs=None):
    def iterate_indices(indices, index, n_recurse_neighbors):
        return indices[iterate_indices(indices, index, n_recurse_neighbors - 1)] \
            if n_recurse_neighbors > 1 else indices[index]
    indices = np.unique(iterate_indices(indices, index, n_recurse_neighbors))
    if max_neighs is not None and len(indices) > max_neighs:
        indices = np.random.choice(indices, max_neighs, replace=False)
    return indices


def groups_to_bool(adata, groups, groupby=None):
    groups = [groups] if isinstance(groups, str) else groups
    if isinstance(groups, (list, tuple, np.ndarray, np.record)):
        groupby = groupby if groupby in adata.obs.keys() else 'clusters' if 'clusters' in adata.obs.keys() \
            else 'louvain' if 'louvain' in adata.obs.keys() else None
        if groupby is not None:
            groups = np.array([key in groups for key in adata.obs[groupby]])
        else: raise ValueError('groupby attribute not valid.')
    return groups


def most_common_in_list(lst):
    lst = list(lst)
    return max(set(lst), key=lst.count)


def randomized_velocity(adata, vkey='velocity', add_key='velocity_random'):
    V_rnd = adata.layers[vkey].copy()
    for i in range(V_rnd.shape[1]):
        np.random.shuffle(V_rnd[:, i])
        V_rnd[:, i] = V_rnd[:, i] * np.random.choice(np.array([+1, -1]), size=V_rnd.shape[0])
    adata.layers[add_key] = V_rnd

    from .velocity_graph import velocity_graph
    from .velocity_embedding import velocity_embedding

    velocity_graph(adata, vkey=add_key)
    velocity_embedding(adata, vkey=add_key, autoscale=False)
