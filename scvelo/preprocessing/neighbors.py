from .. import settings
from .. import logging as logg

from scanpy.api import Neighbors
from scanpy.api.pp import pca
from scipy.sparse import issparse
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

    if (use_rep is None or use_rep is 'X_pca') \
            and ('X_pca' not in adata.obsm.keys() or n_pcs > adata.obsm['X_pca'].shape[1]):
        pca(adata, n_comps=n_pcs, svd_solver='arpack')

    neighbors = Neighbors(adata)
    neighbors.compute_neighbors(n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs, use_rep=use_rep, method=method,
                                metric=metric, metric_kwds=metric_kwds, random_state=random_state, write_knn_indices=True)
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': method}
    adata.uns['neighbors']['distances'] = neighbors.distances
    adata.uns['neighbors']['connectivities'] = neighbors.connectivities
    adata.uns['neighbors']['indices'] = neighbors.knn_indices
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns[\'neighbors\']`\n'
        '    \'distances\', weighted adjacency matrix\n'
        '    \'connectivities\', weighted adjacency matrix')
    return adata if copy else None


def select_distances(dist, n_neighbors=None):
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if issparse(D) else (D > 0).sum(1)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()
    return D


def select_connectivities(connectivities, n_neighbors=None):
    C = connectivities.copy()
    n_counts = (C > 0).sum(1).A1 if issparse(C) else (C > 0).sum(1)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = C.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[::-1][n_neighbors:]
        dat[rm_idx] = 0
    C.eliminate_zeros()
    return C


def get_connectivities(adata, mode='connectivities', n_neighbors=None, recurse_neighbors=False):
    C = adata.uns['neighbors'][mode]
    if n_neighbors is not None and n_neighbors < adata.uns['neighbors']['params']['n_neighbors']:
        C = select_connectivities(C, n_neighbors) if mode == 'connectivities' else select_distances(C, n_neighbors)
    connectivities = C > 0
    connectivities.setdiag(1)
    if recurse_neighbors:
        connectivities += connectivities.dot(connectivities * .5)
        connectivities.data = np.clip(connectivities.data, 0, 1)
    connectivities = connectivities.multiply(1. / connectivities.sum(1))
    return connectivities.tocsr().astype(np.float32)