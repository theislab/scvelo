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


def normalize_sparse(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X.multiply(csr_matrix(1. / X.sum(1)))
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
