from .. import settings
from .. import logging as logg
from ..preprocessing.neighbors import pca, neighbors
from .utils import cosine_correlation, get_indices, get_iterative_indices
from .velocity import velocity

from scipy.sparse import coo_matrix, csr_matrix, issparse
import numpy as np


def vals_to_csr(vals, rows, cols, shape, split_negative=False):
    graph = coo_matrix((vals, (rows, cols)), shape=shape)

    if split_negative:
        graph_neg = graph.copy()

        graph.data = np.clip(graph.data, 0, 1)
        graph_neg.data = np.clip(graph_neg.data, -1, 0)

        graph.eliminate_zeros()
        graph_neg.eliminate_zeros()

        return graph.tocsr(), graph_neg.tocsr()

    else:
        return graph.tocsr()


class VelocityGraph:
    def __init__(self, adata, vkey='velocity', xkey='Ms', tkey=None, basis=None, n_neighbors=None, n_recurse_neighbors=None,
                 random_neighbors_at_max=None, sqrt_transform=False, approx=False, report=False):
        subset = np.ones(adata.n_vars, dtype=bool)
        if 'velocity_genes' in adata.var.keys(): subset = adata.var['velocity_genes']

        X = adata[:, subset].layers[xkey].A if issparse(adata.layers[xkey]) else adata[:, subset].layers[xkey]
        V = adata[:, subset].layers[vkey].A if issparse(adata.layers[vkey]) else adata[:, subset].layers[vkey]
        if approx and X.shape[1] > 100:
            X_pca, PCs, _, _ = pca(X,  n_comps=50, svd_solver='arpack', return_info=True)
            self.X = np.array(X_pca, dtype=np.float32)
            self.V = (V - V.mean(0)).dot(PCs.T)
            self.V[V.sum(1) == 0] = 0
        else:
            self.X = np.array(X.copy(), dtype=np.float32)
            self.V = np.array(V.copy(), dtype=np.float32)

        self.sqrt_transform = sqrt_transform
        if sqrt_transform: self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= self.V.mean(1)[:, None]

        self.n_recurse_neighbors = 1 if n_neighbors is not None \
            else 2 if n_recurse_neighbors is None else n_recurse_neighbors

        if 'neighbors' not in adata.uns.keys(): neighbors(adata)
        if n_neighbors is None or n_neighbors < adata.uns['neighbors']['params']['n_neighbors']:
            self.indices = get_indices(dist=adata.uns['neighbors']['distances'], n_neighbors=n_neighbors)[0]
        else:
            from .. import Neighbors
            neighs = Neighbors(adata)
            if basis is None: basis = [key for key in ['X_pca', 'X_tsne', 'X_umap'] if key in adata.obsm.keys()][-1]
            neighs.compute_neighbors(n_neighbors=n_neighbors, use_rep=basis, n_pcs=10)
            self.indices = get_indices(dist=neighs.distances)[0]

        self.max_neighs = random_neighbors_at_max

        self.graph = adata.uns[vkey + '_graph'] if vkey + '_graph' in adata.uns.keys() else []
        self.graph_neg = adata.uns[vkey + '_graph_neg'] if vkey + '_graph_neg' in adata.uns.keys() else []

        if tkey in adata.obs.keys():
            self.t0 = adata.obs[tkey].copy()
            init = min(self.t0) if isinstance(min(self.t0), int) else 0
            self.t0.cat.categories = np.arange(init, len(self.t0.cat.categories))
            self.t1 = self.t0 + 1
        else: self.t0 = None

        self.report = report

    def compute_cosines(self):
        vals, rows, cols, n_obs = [], [], [], self.X.shape[0]
        progress = logg.ProgressReporter(n_obs)
        for i in range(n_obs):
            neighs_idx = get_iterative_indices(self.indices, i, self.n_recurse_neighbors, self.max_neighs)

            if self.t0 is not None:
                t0, t1 = self.t0[i], self.t1[i]
                if t0 >= 0 and t1 > 0:
                    t1_idx = np.where(self.t0 == t1)[0]
                    if len(t1_idx) > len(neighs_idx):
                        t1_idx = np.random.choice(t1_idx, len(neighs_idx), replace=False)
                    if len(t1_idx) > 0:
                        neighs_idx = np.unique(np.concatenate([neighs_idx, t1_idx]))

            if self.V[i].max() != 0 or self.V[i].min() != 0:
                dX = self.X[neighs_idx] - self.X[i, None]  # 60% of runtime
                if self.sqrt_transform: dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                val = cosine_correlation(dX, self.V[i])  # 40% of runtime
            else:
                val = np.zeros(len(neighs_idx))
            vals.extend(val)
            rows.extend(np.ones(len(neighs_idx)) * i)
            cols.extend(neighs_idx)
            if self.report: progress.update()
        if self.report: progress.finish()

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 1e-10  # actually zero; just to store these entries in sparse matrix.

        self.graph, self.graph_neg = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs), split_negative=True)


def velocity_graph(data, vkey='velocity', xkey='Ms', tkey=None, basis=None, n_neighbors=None, n_recurse_neighbors=None,
                   random_neighbors_at_max=None, sqrt_transform=False, approx=False, copy=False):
    """Computes velocity graph based on cosine similarities.

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

    vgraph = VelocityGraph(adata, vkey=vkey, xkey=xkey, tkey=tkey, basis=basis, n_neighbors=n_neighbors, approx=approx,
                           n_recurse_neighbors=n_recurse_neighbors, random_neighbors_at_max=random_neighbors_at_max,
                           sqrt_transform=sqrt_transform, report=True)

    logg.info('computing velocity graph', r=True)
    vgraph.compute_cosines()

    adata.uns[vkey+'_graph'] = vgraph.graph
    adata.uns[vkey+'_graph_neg'] = vgraph.graph_neg

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + vkey + '_graph\', sparse matrix with cosine correlations (adata.uns)')

    return adata if copy else None
