from typing import Optional
import warnings

import numpy as np
from scipy.sparse import issparse


def clipped_log(x, lb=0, ub=1, eps=1e-6):
    """Logarithmize between [lb + epsilon, ub - epsilon]."""

    return np.log(np.clip(x, lb + eps, ub - eps))


def invert(x):
    """Invert array and set infinity to NaN.

    Args:
        x: Array to invert.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x_inv = 1 / x * (x != 0)
    return x_inv


# TODO: Add docstrings
def prod_sum(a1, a2, axis=0):
    """Take sum of product along given axis."""

    if issparse(a1):
        return a1.multiply(a2).sum(axis=axis).A1
    elif axis == 0:
        return np.einsum("ij, ij -> j", a1, a2) if a1.ndim > 1 else (a1 * a2).sum()
    elif axis == 1:
        return np.einsum("ij, ij -> i", a1, a2) if a1.ndim > 1 else (a1 * a2).sum()


# TODO: Finish type hints.
# TODO: Finish docstrings
def sum(a, axis: Optional[int] = None):
    """Sum array elements over a given axis.

    Args:
        a (): Elements to sum.
        axis (`None` or `int`): Axis along which to sum elements. If `None`, all
            elements will be summed. Defaults to `None`.

    Returns:
        sum_along_axis (): Sum of array along given axis.
    """

    if a.ndim == 1:
        axis = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return a.sum(axis=axis).A1 if issparse(a) else a.sum(axis=axis)
