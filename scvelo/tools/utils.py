from scipy.sparse import csr_matrix
import numpy as np
import warnings


def prod_sum_obs(A, B):
    """dot product and sum over axis 0 (obs) equivalent to np.sum(A * B, 0)
    """
    return np.einsum('ij, ij -> j', A, B)


def prod_sum_var(A, B):
    """dot product and sum over axis 1 (var) equivalent to np.sum(A * B, 1)
    """
    return np.einsum('ij, ij -> i', A, B)


def norm(A):
    """computes the L2-norm along axis 1 (e.g. genes or embedding dimensions) equivalent to np.linalg.norm(A, axis=1)
    """
    return np.sqrt(np.einsum('ij, ij -> i', A, A))


def cosine(dX, v, v_norm):
    dX -= dX.mean(-1)[:, None]
    return np.einsum('ij, j', dX, v) / (norm(dX) * v_norm)[None, :]


def normalize_sparse(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X.multiply(csr_matrix(1. / np.abs(X).sum(1)))
    return X


def normalize_dense(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X / X.sum(1)
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


def get_iterative_indices(indices, index, n_recurse_neighbors):
    return indices[get_iterative_indices(indices, index, n_recurse_neighbors-1)] \
        if n_recurse_neighbors > 1 else indices[index]
