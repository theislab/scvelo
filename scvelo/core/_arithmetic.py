import warnings
from typing import Optional, Union

import numpy as np
from numpy import ndarray
from scipy.sparse import issparse, spmatrix


def clipped_log(x: ndarray, lb: float = 0, ub: float = 1, eps: float = 1e-6) -> ndarray:
    """Logarithmize between [lb + epsilon, ub - epsilon].

    Arguments
    ---------
    x
        Array to invert.
    lb
        Lower bound of interval to which array entries are clipped.
    ub
        Upper bound of interval to which array entries are clipped.
    eps
        Offset of boundaries of clipping interval.

    Returns
    -------
    ndarray
        Logarithm of clipped array.
    """

    return np.log(np.clip(x, lb + eps, ub - eps))


def invert(x: ndarray) -> ndarray:
    """Invert array and set infinity to NaN.

    Arguments
    ---------
    x
        Array to invert.

    Returns
    -------
    ndarray
        Inverted array.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x_inv = 1 / x * (x != 0)
    return x_inv


def prod_sum(
    a1: Union[ndarray, spmatrix], a2: Union[ndarray, spmatrix], axis: Optional[int]
) -> ndarray:
    """Take sum of product of two arrays along given axis.

    Arguments
    ---------
    a1
        First array.
    a2
        Second array.
    axis
        Axis along which to sum elements. If `None`, all elements will be summed.
        Defaults to `None`.

    Returns
    -------
    ndarray
        Sum of product of arrays along given axis.
    """

    if issparse(a1):
        return a1.multiply(a2).sum(axis=axis).A1
    elif axis == 0:
        return np.einsum("ij, ij -> j", a1, a2) if a1.ndim > 1 else (a1 * a2).sum()
    elif axis == 1:
        return np.einsum("ij, ij -> i", a1, a2) if a1.ndim > 1 else (a1 * a2).sum()


def sum(a: Union[ndarray, spmatrix], axis: Optional[int] = None) -> ndarray:
    """Sum array elements over a given axis.

    Arguments
    ---------
    a
        Elements to sum.
    axis
        Axis along which to sum elements. If `None`, all elements will be summed.
        Defaults to `None`.

    Returns
    -------
    ndarray
        Sum of array along given axis.
    """

    if a.ndim == 1:
        axis = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return a.sum(axis=axis).A1 if issparse(a) else a.sum(axis=axis)
