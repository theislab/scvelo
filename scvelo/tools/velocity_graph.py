from .. import settings
from .. import logging as logg
from ..preprocessing.neighbors import pca, neighbors, neighbors_to_be_recomputed
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
    def __init__(self, adata, vkey='velocity', xkey='Ms', tkey=None, basis=None, n_neighbors=None, sqrt_transform=None,
                 n_recurse_neighbors=None, random_neighbors_at_max=None, gene_subset=None, approx=None, report=False,
                 mode_neighbors='distances'):

        subset = np.ones(adata.n_vars, bool)
        if gene_subset is not None:
            subset &= adata.var_names.isin(gene_subset) if len(adata.var_names.isin(gene_subset)) > 0 else gene_subset
        elif vkey + '_genes' in adata.var.keys():
            subset &= np.array(adata.var[vkey + '_genes'].values, dtype=bool)

        xkey = xkey if xkey in adata.layers.keys() else 'spliced'

        X = np.array(adata.layers[xkey].A[:, subset] if issparse(adata.layers[xkey]) else adata.layers[xkey][:, subset])
        V = np.array(adata.layers[vkey].A[:, subset] if issparse(adata.layers[vkey]) else adata.layers[vkey][:, subset])

        nans = np.isnan(np.sum(V, axis=0))
        if np.any(nans):
            X = X[:, ~nans]
            V = V[:, ~nans]

        if approx is True and X.shape[1] > 100:
            X_pca, PCs, _, _ = pca(X,  n_comps=30, svd_solver='arpack', return_info=True)
            self.X = np.array(X_pca, dtype=np.float32)
            self.V = (V - V.mean(0)).dot(PCs.T)
            self.V[V.sum(1) == 0] = 0
        else:
            self.X = np.array(X, dtype=np.float32)
            self.V = np.array(V, dtype=np.float32)

        self.sqrt_transform = (adata.uns[vkey + '_settings']['mode'] is 'stochastic') if sqrt_transform is None else sqrt_transform
        if self.sqrt_transform: self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= np.nanmean(self.V, axis=1)[:, None]

        self.n_recurse_neighbors = 1 if n_neighbors is not None \
            else 2 if n_recurse_neighbors is None else n_recurse_neighbors

        if 'neighbors' not in adata.uns.keys(): neighbors(adata)
        if np.min((adata.uns['neighbors']['distances'] > 0).sum(1).A1) == 0:
            raise ValueError('Your neighbor graph seems to be corrupted. Consider recomputing via pp.neighbors.')
        if n_neighbors is None or n_neighbors <= adata.uns['neighbors']['params']['n_neighbors']:
            self.indices = get_indices(dist=adata.uns['neighbors']['distances'], n_neighbors=n_neighbors,
                                       mode_neighbors=mode_neighbors)[0]
        else:
            if basis is None: basis = [key for key in ['X_pca', 'X_tsne', 'X_umap'] if key in adata.obsm.keys()][-1]
            elif 'X_' + basis in adata.obsm.keys(): basis = 'X_' + basis

            if isinstance(approx, str) and approx in adata.obsm.keys():
                from sklearn.neighbors import NearestNeighbors
                neighs = NearestNeighbors(n_neighbors=n_neighbors + 1)
                neighs.fit(adata.obsm[approx])
                self.indices = neighs.kneighbors_graph(mode='connectivity').indices.reshape((-1, n_neighbors + 1))
            else:
                from .. import Neighbors
                neighs = Neighbors(adata)
                neighs.compute_neighbors(n_neighbors=n_neighbors, use_rep=basis, n_pcs=10)
                self.indices = get_indices(dist=neighs.distances, mode_neighbors=mode_neighbors)[0]

        self.max_neighs = random_neighbors_at_max

        self.graph = adata.uns[vkey + '_graph'] if vkey + '_graph' in adata.uns.keys() else []
        self.graph_neg = adata.uns[vkey + '_graph_neg'] if vkey + '_graph_neg' in adata.uns.keys() else []

        if tkey in adata.obs.keys():
            self.t0 = adata.obs[tkey].copy()
            init = min(self.t0) if isinstance(min(self.t0), int) else 0
            self.t0.cat.categories = np.arange(init, len(self.t0.cat.categories))
            self.t1 = self.t0.copy()
            self.t1.cat.categories = self.t0.cat.categories + 1
        else: self.t0 = None

        self.report = report
        self.self_prob = None

    def compute_cosines(self):
        vals, rows, cols, n_obs = [], [], [], self.X.shape[0]
        progress = logg.ProgressReporter(n_obs)

        for i in range(n_obs):
            if self.V[i].max() != 0 or self.V[i].min() != 0:
                neighs_idx = get_iterative_indices(self.indices, i, self.n_recurse_neighbors, self.max_neighs)

                if self.t0 is not None:
                    t0, t1 = self.t0[i], self.t1[i]
                    if t0 >= 0 and t1 > 0:
                        t1_idx = np.where(self.t0 == t1)[0]
                        if len(t1_idx) > len(neighs_idx):
                            t1_idx = np.random.choice(t1_idx, len(neighs_idx), replace=False)
                        if len(t1_idx) > 0:
                            neighs_idx = np.unique(np.concatenate([neighs_idx, t1_idx]))

                dX = self.X[neighs_idx] - self.X[i, None]  # 60% of runtime
                if self.sqrt_transform: dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                val = cosine_correlation(dX, self.V[i])  # 40% of runtime

                vals.extend(val)
                rows.extend(np.ones(len(neighs_idx)) * i)
                cols.extend(neighs_idx)
                if self.report: progress.update()
        if self.report: progress.finish()

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 0

        self.graph, self.graph_neg = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs), split_negative=True)

        confidence = self.graph.max(1).A.flatten()
        self.self_prob = np.clip(np.percentile(confidence, 98) - confidence, 0, 1)


def velocity_graph(data, vkey='velocity', xkey='Ms', tkey=None, basis=None, n_neighbors=None, n_recurse_neighbors=None,
                   random_neighbors_at_max=None, sqrt_transform=None, variance_stabilization=None, gene_subset=None,
                   approx=None, mode_neighbors='distances', copy=False):
    """Computes velocity graph based on cosine similarities.

    The cosine similarities are computed between velocities and potential cell state transitions, i.e. it measures how
    well a corresponding change in gene expression :math:`\\delta_{ij} = x_j - x_i` matches the predicted change
    according to the velocity vector :math:`\\nu_i`,

    .. math::
        \\pi_{ij} = \\cos\\angle(\\delta_{ij}, \\nu_i)
        = \\frac{\\delta_{ij}^T \\nu_i}{\\left\\lVert\\delta_{ij}\\right\\rVert \\left\\lVert \\nu_i \\right\\rVert}.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    xkey: `str` (default: `'Ms'`)
        Layer key to extract count data from.
    tkey: `str` (default: `None`)
        Observation key to extract time data from.
    basis: `str` (default: `None`)
        Basis / Embedding to use.
    n_neighbors: `int` or `None` (default: None)
        Use fixed number of neighbors or do recursive neighbor search (if `None`).
    n_recurse_neighbors: `int` (default: 2)
        Number of recursions to be done for neighbors search.
    random_neighbors_at_max: `int` or `None` (default: `None`)
        If number of iterative neighbors for an individual cell is higher than this threshold,
        a random selection of such are chosen as reference neighbors.
    sqrt_transform: `bool` (default: `False`)
        Whether to variance-transform the cell states changes and velocities before computing cosine similarities.
    gene_subset: `list` of `str`, subset of adata.var_names or `None`(default: `None`)
        Subset of genes to compute velocity graph on exclusively.
    approx: `bool` or `None` (default: `None`)
        If True, first 30 pc's are used instead of the full count matrix
    mode_neighbors: 'str' (default: `'distances'`)
        Determines the type of KNN graph used. Options are 'distances' or 'connectivities'. The latter yields a
        symmetric KNN graph while the former does not.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_graph: `.uns`
        sparse matrix with transition probabilities
    """
    adata = data.copy() if copy else data
    if neighbors_to_be_recomputed(adata): neighbors(adata)
    if vkey not in adata.layers.keys(): velocity(adata, vkey=vkey)
    if sqrt_transform is None: sqrt_transform = variance_stabilization

    vgraph = VelocityGraph(adata, vkey=vkey, xkey=xkey, tkey=tkey, basis=basis, n_neighbors=n_neighbors, approx=approx,
                           n_recurse_neighbors=n_recurse_neighbors, random_neighbors_at_max=random_neighbors_at_max,
                           sqrt_transform=sqrt_transform, gene_subset=gene_subset, report=True,
                           mode_neighbors=mode_neighbors)

    if isinstance(basis, str):
        logg.warn(
            'The velocity graph is computed on ' + basis + 'embedding coordinates. Consider computing \n'
            '         the graph in an unbiased manner on full expression space by not specifying basis.\n')

    logg.info('computing velocity graph', r=True)
    vgraph.compute_cosines()

    adata.uns[vkey+'_graph'] = vgraph.graph
    adata.uns[vkey+'_graph_neg'] = vgraph.graph_neg

    adata.obs[vkey+'_self_transition'] = vgraph.self_prob

    if vkey + '_settings' in adata.uns.keys() and 'embeddings' in adata.uns[vkey + '_settings']:
        del adata.uns[vkey + '_settings']['embeddings']

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + vkey + '_graph\', sparse matrix with cosine correlations (adata.uns)')

    return adata if copy else None
