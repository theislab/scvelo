from .. import settings
from .. import logging as logg
from .utils import normalize_layers
from .neighbors import neighbors, get_connectivities

from scipy.sparse import csr_matrix
import numpy as np


def moments(data, n_neighbors=30, n_pcs=30, mode='connectivities', use_rep=None, recurse_neighbors=False,
            renormalize=False, copy=False):
    """Computes moments for velocity estimation.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    n_neighbors: `int` (default: 30)
        Number of neighbors to use.
    n_pcs: `int` (default: 30)
        Number of principal components to use.
    mode: `'connectivities'` or `'distances'`  (default: `'connectivities'`)
        Distance metric to use for moment computation.
    renormalize: `bool` (default: `False`)
        Renormalize the moments by total counts per cell to its median.
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
    adata = data.copy() if copy else data

    if 'spliced' not in adata.layers.keys() or 'unspliced' not in adata.layers.keys():
        raise ValueError('Could not find spliced / unspliced counts.')
    if 'neighbors' not in adata.uns.keys() or n_neighbors > adata.uns['neighbors']['params']['n_neighbors']:
        neighbors(adata, n_neighbors=n_neighbors, use_rep=('X_pca' if use_rep is None else use_rep), n_pcs=n_pcs)
    if mode not in adata.uns['neighbors']:
        raise ValueError('mode can only be  \'connectivities\' or \'distances\'')

    logg.info('computing moments based on ' + str(mode), r=True)
    normalize_layers(adata)

    connectivities = get_connectivities(adata, mode, n_neighbors=n_neighbors, recurse_neighbors=recurse_neighbors)

    adata.layers['Ms'] = csr_matrix.dot(connectivities, csr_matrix(adata.layers['spliced'])).astype(np.float32).A
    adata.layers['Mu'] = csr_matrix.dot(connectivities, csr_matrix(adata.layers['unspliced'])).astype(np.float32).A
    if renormalize: normalize_layers(adata, layers={'Ms', 'Mu'}, enforce=True)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'Ms\' and \'Mu\', moments of spliced/unspliced abundances (adata.layers)')
    return adata if copy else None


def second_order_moments(adata, adjusted=False):
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
        raise ValueError('You need to run `pp.neighbors` first to compute a neighborhood graph.')

    connectivities = get_connectivities(adata)
    s, u = csr_matrix(adata.layers['spliced']), csr_matrix(adata.layers['unspliced'])
    Mss = csr_matrix.dot(connectivities, s.multiply(s)).astype(np.float32).A
    Mus = csr_matrix.dot(connectivities, s.multiply(u)).astype(np.float32).A
    if adjusted:
        Mss = 2 * Mss - adata.layers['Ms'].reshape(Mss.shape)
        Mus = 2 * Mus - adata.layers['Mu'].reshape(Mus.shape)
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
        raise ValueError('You need to run `pp.neighbors` first to compute a neighborhood graph.')

    connectivities = get_connectivities(adata)
    u = csr_matrix(adata.layers['unspliced'])
    Muu = csr_matrix.dot(connectivities, u.multiply(u)).astype(np.float32).A

    return Muu
