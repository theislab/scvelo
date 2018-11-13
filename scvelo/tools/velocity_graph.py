from .. import settings
from .. import logging as logg
from .utils import cosine_correlation, get_indices, get_iterative_indices
from .velocity import velocity
from .rank_velocity_genes import rank_velocity_genes

from scipy.sparse import coo_matrix, issparse
import numpy as np


def vals_to_csr(vals, rows, cols, shape):
    graph = coo_matrix((vals, (rows, cols)), shape=shape)
    graph_neg = graph.copy()

    graph.data = np.clip(graph.data, 0, 1)
    graph_neg.data = np.clip(graph_neg.data, -1, 0)

    graph.eliminate_zeros()
    graph_neg.eliminate_zeros()

    return graph.tocsr(), graph_neg.tocsr()


class VelocityGraph:
    def __init__(self, adata, vkey='velocity', xkey='Ms', basis=None, n_neighbors=None, n_recurse_neighbors=None,
                 random_neighbors_at_max=None, n_top_genes=None, sqrt_transform=False):

        subset = np.ones(adata.n_vars, dtype=bool)
        if n_top_genes is not None and 'velocity_score' in adata.var.keys():
            subset = adata.var['velocity_score'][::-1][:n_top_genes]
        elif 'velocity_genes' in adata.var.keys(): subset = adata.var['velocity_genes']

        X = adata[:, subset].layers[xkey].A if issparse(adata.layers[xkey]) else adata[:, subset].layers[xkey]
        V = adata[:, subset].layers[vkey].A if issparse(adata.layers[vkey]) else adata[:, subset].layers[vkey]

        self.X = np.array(X.copy(), dtype=np.float32)
        self.V = np.array(V.copy(), dtype=np.float32)

        self.sqrt_transform = sqrt_transform
        if sqrt_transform: self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= self.V.mean(1)[:, None]

        self.n_recurse_neighbors = 1 if n_neighbors is not None \
            else 2 if n_recurse_neighbors is None else n_recurse_neighbors

        if n_neighbors is None:
            if 'neighbors' not in adata.uns.keys():
                from ..preprocessing.moments import neighbors
                neighbors(adata)
            self.indices = get_indices(dist=adata.uns['neighbors']['distances'])[0]
        else:
            from .. import Neighbors
            neighs = Neighbors(adata)
            if basis is None: basis = [key for key in ['X_pca', 'X_tsne', 'X_umap'] if key in adata.obsm.keys()][-1]
            neighs.compute_neighbors(n_neighbors=n_neighbors, use_rep=basis, n_pcs=10)
            self.indices = get_indices(dist=neighs.distances)[0]

        self.max_neighs = random_neighbors_at_max

        self.graph = adata.uns[vkey + '_graph'] if vkey + '_graph' in adata.uns.keys() else []
        self.graph_neg = adata.uns[vkey + '_graph_neg'] if vkey + '_graph_neg' in adata.uns.keys() else []

    def compute_cosines(self):
        vals, rows, cols, n_obs = [], [], [], self.X.shape[0]
        progress = logg.ProgressReporter(n_obs)
        for i in range(n_obs):
            neighs_idx = get_iterative_indices(self.indices, i, self.n_recurse_neighbors, self.max_neighs)
            if self.V[i].max() != 0 or self.V[i].min() != 0:
                dX = self.X[neighs_idx] - self.X[i, None]  # 60% of runtime
                if self.sqrt_transform: dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                val = cosine_correlation(dX, self.V[i])  # 40% of runtime
            else:
                val = np.zeros(len(neighs_idx))
            vals.extend(val)
            rows.extend(np.ones(len(neighs_idx)) * i)
            cols.extend(neighs_idx)
            progress.update()
        progress.finish()

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 1e-10  # actually zero; just to store these entries in sparse matrix.

        self.graph, self.graph_neg = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs))


def velocity_graph(data, vkey='velocity', xkey='Ms', basis=None, n_neighbors=None, n_recurse_neighbors=None,
                   random_neighbors_at_max=None, n_top_genes=None, sqrt_transform=False, copy=False):
    """Computes a velocity graph based on cosine similarities.

    The cosine similarities are computed between velocities and potential cell state transitions.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    n_neighbors: `int` or `None` (default: None)
        Use fixed number of neighbors or do recursive neighbor search (if `None`).
    n_recurse_neighbors: `int` (default: 2)
        Number of recursions to be done for neighbors search.
    random_neighbors_at_max: `int` or `None` (default: `None`)
        If number of iterative neighbors for an individual cell is higher than this threshold,
        a random selection of such are chosen as reference neighbors.
    sqrt_transform: `bool` (default: `False`)
        Whether to variance-transform the cell states changes and velocities before computing cosine similarities.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_graph: `.uns`
        sparse matrix with transition probabilities
    """
    adata = data.copy() if copy else data
    if vkey not in adata.layers.keys(): velocity(adata, vkey=vkey)
    if n_top_genes is not None and 'velocity_score' not in adata.var.keys(): rank_velocity_genes(adata, n_genes=100)

    vgraph = VelocityGraph(adata, vkey=vkey, xkey=xkey, basis=basis, n_neighbors=n_neighbors,
                           n_recurse_neighbors=n_recurse_neighbors, random_neighbors_at_max=random_neighbors_at_max,
                           n_top_genes=n_top_genes, sqrt_transform=sqrt_transform)

    logg.info('computing velocity graph', r=True)
    vgraph.compute_cosines()

    adata.uns[vkey+'_graph'] = vgraph.graph
    adata.uns[vkey+'_graph_neg'] = vgraph.graph_neg

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + vkey + '_graph\', sparse matrix with cosine correlations (adata.uns)')

    return adata if copy else None
