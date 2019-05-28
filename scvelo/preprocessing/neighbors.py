from .. import settings
from .. import logging as logg

from scanpy.neighbors import compute_connectivities_umap
from scanpy.api import Neighbors
from scanpy.api.pp import pca
from scipy.sparse import issparse
import numpy as np


def neighbors(adata, n_neighbors=30, n_pcs=30, use_rep=None, knn=True, random_state=0, method='umap',
              metric='euclidean', metric_kwds={}, num_threads=-1, copy=False):
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
    method : {{'umap', 'gauss', 'hnsw', 'sklearn', `None`}}  (default: `'umap'`)
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

    if method is 'sklearn':
        from sklearn.neighbors import NearestNeighbors
        X = adata.obsm['X_pca'] if use_rep is None else adata.obsm[use_rep]
        neighbors = NearestNeighbors(n_neighbors=n_neighbors, metric=metric, metric_params=metric_kwds, n_jobs=num_threads)
        neighbors.fit(X)
        knn_distances, neighbors.knn_indices = neighbors.kneighbors()
        neighbors.distances, neighbors.connectivities = \
            compute_connectivities_umap(neighbors.knn_indices, knn_distances, X.shape[0], n_neighbors=30)

    elif method is 'hnsw':
        X = adata.obsm['X_pca'] if use_rep is None else adata.obsm[use_rep]
        neighbors = FastNeighbors(n_neighbors=n_neighbors, num_threads=num_threads)
        neighbors.fit(X, metric=metric, random_state=random_state, **metric_kwds)

    else:
        neighbors = Neighbors(adata)
        neighbors.compute_neighbors(n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs, use_rep=use_rep, method=method,
                                    metric=metric, metric_kwds=metric_kwds, random_state=random_state, write_knn_indices=True)

    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': method}

    adata.uns['neighbors']['distances'] = neighbors.distances
    adata.uns['neighbors']['connectivities'] = neighbors.connectivities
    if hasattr(neighbors, 'knn_indices'):
        adata.uns['neighbors']['indices'] = neighbors.knn_indices

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns[\'neighbors\']`\n'
        '    \'distances\', weighted adjacency matrix\n'
        '    \'connectivities\', weighted adjacency matrix')
    return adata if copy else None


class FastNeighbors:
    def __init__(self, n_neighbors=30, num_threads=-1):
        self.n_neighbors = n_neighbors
        self.num_threads = num_threads
        self.knn_indices, self.knn_distances = None, None
        self.distances, self.connectivities = None, None

    def fit(self, X, metric='l2', M=16, ef=100, ef_construction=100, random_state=0):
        try:
            import hnswlib
        except ImportError:
            print("In order to use fast approx neighbor search, you need to install hnswlib via \n \n"
                  "pip install -U pybind11 \n"
                  "pip install -U git+https://github.com/nmslib/hnswlib#subdirectory=python_bindings")

        ef_c, ef = max(ef_construction, self.n_neighbors), max(self.n_neighbors, ef)
        metric = 'l2' if metric is 'euclidean' else metric
        ns, dim = X.shape

        knn = hnswlib.Index(space=metric, dim=dim)
        knn.init_index(max_elements=X.shape[0], ef_construction=ef_c, M=M, random_seed=random_state)
        knn.add_items(X)
        knn.set_ef(ef)

        knn_indices, knn_distances = knn.knn_query(X, k=self.n_neighbors, num_threads=self.num_threads)

        n_neighbors = self.n_neighbors
        if knn_distances[0, 0] == 0:
            knn_distances = knn_distances[:, 1:]
            knn_indices = knn_indices[:, 1:].astype(int)
            n_neighbors -= 1

        if metric is 'l2':
            knn_distances = np.sqrt(knn_distances)

        self.distances, self.connectivities = compute_connectivities_umap(knn_indices, knn_distances, ns, n_neighbors)
        self.knn_indices = knn_indices


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


def neighbors_to_be_recomputed(adata, n_neighbors=None):
    # check if neighbors graph is disrupted
    if 'neighbors' not in adata.uns.keys() \
            or 'distances' not in adata.uns['neighbors'] or 'params' not in adata.uns['neighbors']:
        return True
    else:
        n_neighs = (adata.uns['neighbors']['distances'] > 0).sum(1)
        result = n_neighs.max() - n_neighs.min() >= 2
        # check if neighbors graph has sufficient number of neighbors
        if n_neighbors is not None:
            result = result or n_neighbors > adata.uns['neighbors']['params']['n_neighbors']
    return result


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
