import numpy as np


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
