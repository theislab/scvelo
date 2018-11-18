from .. import settings
from .. import logging as logg
from .utils import normalize_layers, X_is_logarithmized

from scanpy.api import Neighbors
from scanpy.api.pp import pca
from scipy.sparse import csr_matrix
import numpy as np


def neighbors(adata, n_neighbors=30, n_pcs=30, use_rep=None, knn=True, random_state=0, method='umap',
              metric='euclidean', metric_kwds={}, copy=False):
    """
    Compute a neighborhood graph of observations [McInnes18]_.
    The neighbor search efficiency of this heavily relies on UMAP [McInnes18]_,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='diffmap'`,
    connectivities are computed according to [Coifman05]_, in the adaption of
    [Haghverdi16]_.
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    n_pcs : `int` or `None` (default: None)
        Use this many PCs. If n_pcs==0 use .X if use_rep is None.

    use_rep : `None`, `'X'` or any key for `.obsm` (default: None)
        Use the indicated representation. If `None`, the representation is chosen automatically:
        for .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    random_state
        A numpy random seed.
    method : {{'umap', 'gauss', `None`}}  (default: `'umap'`)
        Use 'umap' [McInnes18]_ or 'gauss' (Gauss kernel following [Coifman05]_
        with adaptive width [Haghverdi16]_) for computing connectivities.
    metric
        A known metric’s name or a callable that returns a distance.
    metric_kwds
        Options for the metric.
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:
    connectivities : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    distances : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """
    logg.info('computing neighbors', r=True)
    adata = adata.copy() if copy else adata
    if adata.isview: adata._init_as_actual(adata.copy())

    neighbors = Neighbors(adata)
    neighbors.compute_neighbors(n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs, use_rep=use_rep,
                                method=method, metric=metric, metric_kwds=metric_kwds, random_state=random_state)
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': method}
    adata.uns['neighbors']['distances'] = neighbors.distances
    adata.uns['neighbors']['connectivities'] = neighbors.connectivities
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns[\'neighbors\']`\n'
        '    \'distances\', weighted adjacency matrix\n'
        '    \'connectivities\', weighted adjacency matrix')
    return adata if copy else None


def get_connectivities(adata, mode='connectivities', recurse_neighbors=False):
    connectivities = adata.uns['neighbors'][mode] > 0
    connectivities.setdiag(1)
    if recurse_neighbors:
        connectivities += connectivities.dot(connectivities * .5)
        connectivities.data = np.clip(connectivities.data, 0, 1)
    connectivities = connectivities.multiply(1. / connectivities.sum(1))
    return connectivities.tocsr().astype(np.float32)


def moments(data, n_neighbors=30, n_pcs=30, mode='connectivities', use_rep=None, recurse_neighbors=False,
            renormalize=False, plot=False, copy=False):
    """Computes first order moments for velocity estimation.

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

    if not X_is_logarithmized(adata):
        logg.info('Consider to logarithmize adata.X with `scv.pp.log1p` for better results.')

    if 'neighbors' not in adata.uns.keys() or n_neighbors > adata.uns['neighbors']['params']['n_neighbors']:
        if 'X_pca' not in adata.obsm.keys() or n_pcs > adata.obsm['X_pca'].shape[1]:
            pca(adata, n_comps=n_pcs, svd_solver='arpack')
            if plot:
                from scanpy.plotting.tools import pca_variance_ratio
                pca_variance_ratio(adata)
        neighbors(adata, n_neighbors=n_neighbors, use_rep=('X_pca' if use_rep is None else use_rep), n_pcs=n_pcs)

    if mode not in adata.uns['neighbors']:
        raise ValueError('mode can only be  \'connectivities\' or \'distances\'')

    logg.info('computing moments based on ' + str(mode), r=True)
    normalize_layers(adata)

    connectivities = get_connectivities(adata, mode, recurse_neighbors=recurse_neighbors)

    adata.layers['Ms'] = csr_matrix.dot(connectivities, csr_matrix(adata.layers['spliced'])).astype(np.float32).A
    adata.layers['Mu'] = csr_matrix.dot(connectivities, csr_matrix(adata.layers['unspliced'])).astype(np.float32).A
    if renormalize: normalize_layers(adata, layers={'Ms', 'Mu'})

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
