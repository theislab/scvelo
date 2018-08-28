from ..logging import logg, settings
from .utils import norm
from .velocity import velocity
from scanpy.api import Neighbors
from scipy.sparse import coo_matrix
import numpy as np
import warnings


def get_indices(dist):
    n_neighbors = (dist > 0).sum(1).min()
    rows_idx = np.where((dist > 0).sum(1) > n_neighbors)[0]

    for row_idx in rows_idx:
        col_idx = dist[row_idx].indices[n_neighbors:]
        dist[row_idx, col_idx] = 0

    dist.eliminate_zeros()

    indices = dist.indices.reshape((-1, n_neighbors))
    return indices, dist


def get_iterative_indices(indices, index, n_recurse_neighbors):
    return indices[get_iterative_indices(indices, index, n_recurse_neighbors-1)] \
        if n_recurse_neighbors > 1 else indices[index]


class Cosines:
    def __init__(self, X, V, indices, n_recurse_neighbors, sqrt_transform=False):
        self.X = X.copy()
        self.indices = indices
        self.n_recurse_neighbors = n_recurse_neighbors
        self.sqrt_transform = sqrt_transform

        self.V = V.copy()
        if sqrt_transform: self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V = self.V - self.V.mean(1)[:, None]
        self.V_norm = norm(self.V)

    def compute(self, ixs):
        vals, rows, cols, async_id = [], [], [], []

        if self.sqrt_transform:
            for i in ixs:
                knn_ixs = np.unique(get_iterative_indices(self.indices, i, self.n_recurse_neighbors))
                dX = self.X[knn_ixs] - self.X[i, None]
                dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                dX -= dX.mean(1)[:, None]

                vals.extend(np.einsum('ij, j', dX, self.V[i]) / (norm(dX) * self.V_norm[i])[None, :])
                rows.extend(np.ones(len(knn_ixs)) * i)
                cols.extend(knn_ixs)

        else:
            for i in ixs:
                knn_ixs = np.unique(get_iterative_indices(self.indices, i, self.n_recurse_neighbors))
                dX = self.X[knn_ixs] - self.X[i, None]
                dX -= dX.mean(1)[:, None]

                vals.extend(np.einsum('ij, j', dX, self.V[i]) / (norm(dX) * self.V_norm[i])[None, :])
                rows.extend(np.ones(len(knn_ixs)) * i)
                cols.extend(knn_ixs)
        async_id = int(ixs[0] / len(ixs))  # keep track of order

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 0
        vals = np.clip(vals, 1e-10, 1)

        return [async_id, vals, rows, cols]


def velocity_graph(adata, vkey='velocity', n_recurse_neighbors=2, n_neighbors=None, sqrt_transform=False, copy=False):
    """Computes a velocity graph based on cosine similarities.

    The cosine similarities are computed between velocities and potential cell state transitions.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    n_neighbors: `int` or `None` (default: None)
        Use fixed number of neighbors or do recursive neighbor search (if `None`).
    n_recurse_neighbors: `int` (default: 2)
        Number of recursions to be done for neighbors search.
    sqrt_transform: `bool` (default: `False`)
        Whether to variance-transform the cell states and velocities before computing cosine similarities.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_graph: `.uns`
        sparse matrix with transition probabilities
    """
    if vkey not in adata.layers.keys(): velocity(adata)

    logg.info('computing velocity graph', r=True)

    if n_neighbors is not None:
        keys = [key for key in ['X_pca', 'X_tsne', 'X_umap'] if key in adata.obsm.keys()]
        neighs = Neighbors(adata)
        neighs.compute_neighbors(n_neighbors=n_neighbors, use_rep=keys[-1], n_pcs=10)
        indices = get_indices(dist=neighs.distances)[0]
        n_recurse_neighbors = 1
    else:
        indices = get_indices(dist=adata.uns['neighbors']['distances'])[0]

    cos = Cosines(adata.layers['Ms'][:, adata.var['velocity_genes']],
                  adata.layers[vkey][:, adata.var['velocity_genes']], indices, n_recurse_neighbors, sqrt_transform)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if True:
            _, vals, rows, cols = cos.compute(range(adata.n_obs))  # list(map(cos.compute, range(adata.n_obs)))
        else:
            from multiprocessing import Pool
            pool = Pool(n_jobs)
            ranges = np.array_split(range(adata.n_obs), n_jobs)
            result = pool.map_async(cos.compute, ranges).get()  # async will probably destroy order
            result.sort(key=lambda tup: tup[0])  # restore order
            vals, rows, cols = [], [], []
            for r in result:
                vals.extend(r[1])
                rows.extend(r[2])
                cols.extend(r[3])

    graph = coo_matrix((vals, (rows, cols)), shape=(adata.n_obs, adata.n_obs))
    graph.eliminate_zeros()

    adata.uns[vkey+'_graph'] = graph.tocsr()

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns`\n'
        '    \'' + vkey + '_graph\', sparse matrix with cosine correlations')

    return adata if copy else None
