from .. import settings
from .. import logging as logg

from scipy.sparse import issparse, coo_matrix
import numpy as np
import warnings
from scanpy.preprocessing import pca
from scanpy import Neighbors


def neighbors(adata, n_neighbors=30, n_pcs=None, use_rep=None, knn=True, random_state=0, method='umap',
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
        Number of principal components to use.
        If not specified, the full space is used of a pre-computed PCA,
        or 30 components are used when PCA is computed internally.
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
    method : {{'umap', 'hnsw', 'sklearn'}}  (default: `'umap'`)
        The methods only differ in runtime. Connectivities are computed with adaptive width [Haghverdi16]_).
        The 'hnsw' method is most efficient and requires to `pip install hnswlib`.
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
    adata = adata.copy() if copy else adata

    if use_rep is None:
        use_rep = 'X' if adata.n_vars < 50 or n_pcs is 0 else 'X_pca'
        n_pcs = None if use_rep is 'X' else n_pcs
    elif use_rep not in adata.obsm.keys() and 'X_' + use_rep in adata.obsm.keys():
        use_rep = 'X_' + use_rep

    if use_rep is 'X_pca':
        if 'X_pca' not in adata.obsm.keys() or n_pcs is not None and n_pcs > adata.obsm['X_pca'].shape[1]:
            pca(adata, n_comps=30 if n_pcs is None else n_pcs, svd_solver='arpack')
        elif n_pcs is None and adata.obsm['X_pca'].shape[1] < 10:
            logg.warn('Neighbors are computed on ', adata.obsm['X_pca'].shape[1], ' principal components only.')

        if len(set(np.sum(adata.obsm['X_pca'], 1))) < adata.n_obs:
            logg.warn('You seem to have duplicate cells in your data. '
                      'Consider removing these via pp.remove_duplicate_cells.')

    logg.info('computing neighbors', r=True)

    if method is 'sklearn':
        from sklearn.neighbors import NearestNeighbors
        X = adata.X if use_rep is 'X' else adata.obsm[use_rep]
        neighbors = NearestNeighbors(n_neighbors=n_neighbors, metric=metric, metric_params=metric_kwds, n_jobs=num_threads)
        neighbors.fit(X if n_pcs is None else X[:, :n_pcs])
        knn_distances, neighbors.knn_indices = neighbors.kneighbors()
        neighbors.distances, neighbors.connectivities = \
            compute_connectivities_umap(neighbors.knn_indices, knn_distances, X.shape[0], n_neighbors=30)

    elif method is 'hnsw':
        X = adata.X if use_rep is 'X' else adata.obsm[use_rep]
        neighbors = FastNeighbors(n_neighbors=n_neighbors, num_threads=num_threads)
        neighbors.fit(X if n_pcs is None else X[:, :n_pcs], metric=metric, random_state=random_state, **metric_kwds)

    else:
        logg.switch_verbosity('off', module='scanpy')
        with warnings.catch_warnings():  # ignore numba warning (reported in umap/issues/252)
            warnings.simplefilter("ignore")
            neighbors = Neighbors(adata)
            neighbors.compute_neighbors(n_neighbors=n_neighbors, knn=knn, n_pcs=n_pcs, use_rep=use_rep, method=method,
                                        metric=metric, metric_kwds=metric_kwds, random_state=random_state,
                                        write_knn_indices=True)
        logg.switch_verbosity('on', module='scanpy')

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
            print("In order to use fast approx neighbor search, you need to `pip install hnswlib`\n")

        ef_c, ef = max(ef_construction, self.n_neighbors), max(self.n_neighbors, ef)
        metric = 'l2' if metric is 'euclidean' else metric

        X = X.A if issparse(X) else X
        ns, dim = X.shape

        knn = hnswlib.Index(space=metric, dim=dim)
        knn.init_index(max_elements=ns, ef_construction=ef_c, M=M, random_seed=random_state)
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
    # check whether neighbors graph is disrupted or whether graph has insufficient number of neighbors
    invalid_neighs = 'neighbors' not in adata.uns.keys() \
                     or 'distances' not in adata.uns['neighbors'] \
                     or 'params' not in adata.uns['neighbors'] \
                     or (n_neighbors is not None and n_neighbors > adata.uns['neighbors']['params']['n_neighbors'] )
    if invalid_neighs:
        return True
    else:
        n_neighs = (adata.uns['neighbors']['distances'] > 0).sum(1)
        return n_neighs.max() * .1 > n_neighs.min()


def get_connectivities(adata, mode='connectivities', n_neighbors=None, recurse_neighbors=False):
    if 'neighbors' in adata.uns.keys():
        C = adata.uns['neighbors'][mode]
        if n_neighbors is not None and n_neighbors < adata.uns['neighbors']['params']['n_neighbors']:
            C = select_connectivities(C, n_neighbors) if mode == 'connectivities' else select_distances(C, n_neighbors)
        connectivities = C > 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            connectivities.setdiag(1)
            if recurse_neighbors:
                connectivities += connectivities.dot(connectivities * .5)
                connectivities.data = np.clip(connectivities.data, 0, 1)
            connectivities = connectivities.multiply(1. / connectivities.sum(1))
        return connectivities.tocsr().astype(np.float32)
    else:
        return None


def get_csr_from_indices(knn_indices, knn_dists, n_obs, n_neighbors):
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()


def compute_connectivities_umap(knn_indices, knn_dists, n_obs, n_neighbors, set_op_mix_ratio=1.0,local_connectivity=1.0):
    """\
    This is from umap.fuzzy_simplicial_set [McInnes18]_.
    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    from umap.umap_ import fuzzy_simplicial_set

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(X, n_neighbors, None, None,
                                          knn_indices=knn_indices, knn_dists=knn_dists,
                                          set_op_mix_ratio=set_op_mix_ratio,
                                          local_connectivity=local_connectivity)

    if isinstance(connectivities, tuple):  # returns (result, sigmas, rhos) in umap-learn 0.4
        connectivities = connectivities[0]

    distances = get_csr_from_indices(knn_indices, knn_dists, n_obs, n_neighbors)

    return distances, connectivities.tocsr()


def remove_duplicate_cells(adata):
    if 'X_pca' not in adata.obsm.keys(): pca(adata)
    l = list(np.sum(adata.obsm['X_pca'], 1) + adata.obs['n_counts'])
    n_unique_obs = len(set(l))
    if n_unique_obs < adata.n_obs:
        idx = [l.index(x) for x in set(l)]
        logg.info('Removed ', adata.n_obs - n_unique_obs, ' duplicate cells.')
        adata._inplace_subset_obs(idx)
    neighbors(adata)
