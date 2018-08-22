from scipy.sparse import csr_matrix
import numpy as np
from ..logging import *


def moments(adata, renormalize=False, mode='connectivities', copy=False):
    """Computes first order moments for velocity estimation.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    renormalize: `bool` (default: `True`)
        Renormalize the moments by total counts per cell to its median.
    mode: `'connectivities'` or `'distances'`  (default: `'connectivities'`)
        Distance metric to use for moment computation.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    Ms: `.layers`
        dense matrix with first order moments of spliced counts.
    Mu: `.layers`
        dense matrix with first order moments of unspliced counts.
    """
    if 'neighbors' not in adata.uns:
        raise ValueError('You need to run `pp.neighbors` first to compute a neighborhood graph.')

    logg.info('computing moments', r=True)

    if mode in adata.uns['neighbors']: connectivities = adata.uns['neighbors'][mode]
    else: raise ValueError('mode can only be  \'connectivities\' or \'distances\'')

    connectivities.setdiag(1)
    connectivities = connectivities.multiply(1. / connectivities.sum(1))

    connectivities += connectivities.dot(connectivities*.5)

    adata.layers['Ms'] = csr_matrix.dot(connectivities, adata.layers['spliced']).toarray()
    adata.layers['Mu'] = csr_matrix.dot(connectivities, adata.layers['unspliced']).toarray()

    if renormalize:
        adata.layers['Ms'] = adata.layers['Ms'] / adata.layers['Ms'].sum(1)[:, None] * np.median(adata.layers['Ms'].sum(1))
        adata.layers['Mu'] = adata.layers['Mu'] / adata.layers['Mu'].sum(1)[:, None] * np.median(adata.layers['Mu'].sum(1))

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.layers`\n'
        '    \'Ms\', moments of spliced abundances\n'
        '    \'Mu\', moments of unspliced abundances')

    return adata if copy else None


def second_order_moments(adata):
    """Computes second order moments for stochastic velocity estimation.

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    Returns
    -------
    Mss: Second order moments for spliced abundances
    Mus: Second order moments for spliced with unspliced abundances
    """
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')

    connectivities = adata.uns['neighbors']['connectivities'] > 0
    connectivities = connectivities.multiply(1. / connectivities.sum(1))

    Mss = csr_matrix.dot(connectivities, adata.layers['spliced'].multiply(adata.layers['spliced'])).toarray()
    Mus = csr_matrix.dot(connectivities, adata.layers['spliced'].multiply(adata.layers['unspliced'])).toarray()

    return Mss, Mus


def second_order_moments_u(adata):
    """Computes second order moments for stochastic velocity estimation.

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    Returns
    -------
    Muu: Second order moments for unspliced abundances
    """
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')

    connectivities = adata.uns['neighbors']['connectivities'] > 0
    connectivities = connectivities.multiply(1. / connectivities.sum(1))

    Muu = csr_matrix.dot(connectivities, adata.layers['unspliced'].multiply(adata.layers['unspliced'])).toarray()

    return Muu