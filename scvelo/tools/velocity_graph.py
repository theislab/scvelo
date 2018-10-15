from ..logging import logg, settings
from .utils import cosine_correlation, get_indices, get_iterative_indices
from .velocity import velocity
from .rank_velocity_genes import rank_velocity_genes
from scanpy.api import Neighbors
from scipy.sparse import coo_matrix
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
    def __init__(self, adata, xkey='Ms', vkey='velocity', n_neighbors=None, n_recurse_neighbors=None,
                 n_top_genes=None, sqrt_trafo=False, use_copy=False):
        subset = rank_velocity_genes(adata, n_genes=n_top_genes, min_r2=.1, min_dispersion=0, min_counts=10) \
            if (n_top_genes is not None and n_top_genes < adata.var['velocity_genes'].sum()) \
            else adata.var['velocity_genes'] if 'velocity_genes' in adata.var.keys() \
            else np.ones(adata.n_vars, dtype=bool)

        self.X = np.array(adata[:, subset].layers[xkey].copy(), dtype=np.float32)
        self.V = np.array(adata[:, subset].layers[vkey].copy(), dtype=np.float32)

        self.sqrt_trafo = sqrt_trafo
        if sqrt_trafo: self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= self.V.mean(1)[:, None]

        self.n_recurse_neighbors = 1 if n_neighbors is not None \
            else 2 if n_recurse_neighbors is None else n_recurse_neighbors

        if n_neighbors is None and 'neighbors' in adata.uns.keys():
            self.indices = get_indices(dist=adata.uns['neighbors']['distances'])[0]
        else:
            keys = [key for key in ['X_pca', 'X_tsne', 'X_umap'] if key in adata.obsm.keys()]
            neighs = Neighbors(adata)
            neighs.compute_neighbors(n_neighbors=n_neighbors, use_rep=keys[-1], n_pcs=10)
            self.indices = get_indices(dist=neighs.distances)[0]

        self.graph = adata.uns[vkey + '_graph'] if vkey + '_graph' in adata.uns.keys() else []
        self.graph_neg = adata.uns[vkey + '_graph_neg'] if vkey + '_graph_neg' in adata.uns.keys() else []

    def compute_cosines(self):
        vals, rows, cols, n_obs = [], [], [], self.X.shape[0]
        for i in range(n_obs):
            neighs_idx = get_iterative_indices(self.indices, i, self.n_recurse_neighbors)
            if self.V[i].max() != 0 or self.V[i].min() != 0:
                dX = self.X[neighs_idx] - self.X[i, None]  # 60% of runtime
                if self.sqrt_trafo: dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                val = cosine_correlation(dX, self.V[i])  # 40% of runtime
            else:
                val = np.zeros(len(neighs_idx))
            vals.extend(val)
            rows.extend(np.ones(len(neighs_idx)) * i)
            cols.extend(neighs_idx)

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 1e-10  # actually zero; just to store these entries in sparse matrix.

        self.graph, self.graph_neg = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs))


def velocity_graph(adata, vkey='velocity', n_neighbors=None, n_recurse_neighbors=None, sqrt_transform=False,
                   n_top_genes=None, copy=False, use_copy=False):
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

    vgraph = VelocityGraph(adata, vkey=vkey, n_neighbors=n_neighbors, n_recurse_neighbors=n_recurse_neighbors,
                           n_top_genes=n_top_genes, sqrt_trafo=sqrt_transform, use_copy=use_copy)
    vgraph.compute_cosines()

    adata.uns[vkey+'_graph'] = vgraph.graph
    adata.uns[vkey+'_graph_neg'] = vgraph.graph_neg

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.uns`\n'
        '    \'' + vkey + '_graph\', sparse matrix with cosine correlations')

    return adata if copy else None
