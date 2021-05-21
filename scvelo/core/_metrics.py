from typing import Union

import numpy as np
from numpy import ndarray
from scipy.sparse import issparse, spmatrix


# TODO: Add case `axis == None`
def l2_norm(x: Union[ndarray, spmatrix], axis: int = 1) -> Union[float, ndarray]:
    """Calculate l2 norm along a given axis.

    Arguments
    ---------
    x
        Array to calculate l2 norm of.
    axis
        Axis along which to calculate l2 norm.

    Returns
    -------
    Union[float, ndarray]
        L2 norm along a given axis.
    """

    if issparse(x):
        return np.sqrt(x.multiply(x).sum(axis=axis).A1)
    elif x.ndim == 1:
        return np.sqrt(np.einsum("i, i -> ", x, x))
    elif axis == 0:
        return np.sqrt(np.einsum("ij, ij -> j", x, x))
    elif axis == 1:
        return np.sqrt(np.einsum("ij, ij -> i", x, x))
