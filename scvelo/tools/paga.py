# This is adapted from https://github.com/theislab/paga
from .. import settings
from .. import logging as logg
from .utils import strings_to_categoricals, most_common_in_list
from .velocity_graph import vals_to_csr
from .velocity_pseudotime import velocity_pseudotime
from .rank_velocity_genes import velocity_clusters
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from pandas.api.types import is_categorical
from scanpy.tools._paga import PAGA


def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shap[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warn(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g


def get_sparse_from_igraph(graph, weight_attr=None):
    from scipy.sparse import csr_matrix
    edges = graph.get_edgelist()
    if weight_attr is None:
        weights = [1] * len(edges)
    else:
        weights = graph.es[weight_attr]
    if not graph.is_directed():
        edges.extend([(v, u) for u, v in edges])
        weights.extend(weights)
    shape = graph.vcount()
    shape = (shape, shape)
    if len(edges) > 0:
        return csr_matrix((weights, zip(*edges)), shape=shape)
    else:
        return csr_matrix(shape)


def set_row_csr(csr, rows, value=0):
    """Set all nonzero elements to the given value. Useful to set to 0 mostly.
    """
    for row in rows:
        csr.data[csr.indptr[row]:csr.indptr[row+1]] = value
    if value == 0:
        csr.eliminate_zeros()


class PAGA_tree(PAGA):
    def __init__(self, adata, groups=None, vkey=None, use_time_prior=None, use_root_prior=None,
                 minimum_spanning_tree=None):
        super().__init__(adata=adata, groups=groups, model='v1.2')
        if groups is None:
            groups = 'clusters' if 'clusters' in adata.obs.keys() \
                else 'louvain' if 'louvain' in adata.obs.keys() else None
        self.groups = groups
        self.vkey = vkey
        self.use_time_prior = use_time_prior
        self.use_root_prior = use_root_prior
        self.minimum_spanning_tree = minimum_spanning_tree

    # overwrite to use flexible vkey
    def compute_transitions(self):
        vkey = self.vkey + '_graph'
        if vkey not in self._adata.uns:
            if 'velocyto_transitions' in self._adata.uns:
                self._adata.uns[vkey] = self._adata.uns['velocyto_transitions']
                logg.warn("The key 'velocyto_transitions' has been changed to 'velocity_graph'.")
            else:
                raise ValueError(
                    'The passed AnnData needs to have an `uns` annotation '
                    "with key 'velocity_graph' - a sparse matrix from RNA velocity."
                )
        if self._adata.uns[vkey].shape != (self._adata.n_obs, self._adata.n_obs):
            raise ValueError(
                f"The passed 'velocity_graph' have shape {self._adata.uns[vkey].shape} "
                f"but shoud have shape {(self._adata.n_obs, self._adata.n_obs)}"
            )
        import igraph
        root = None
        clusters = self._adata.obs[self.groups]
        cats = clusters.cat.categories
        vgraph = self._adata.uns[vkey] > .1

        if isinstance(self.use_time_prior, str) and self.use_time_prior in self._adata.obs.keys():
            vpt = self._adata.obs[self.use_time_prior].values
            rows, cols, vals = [], [], []
            for i in range(vgraph.shape[0]):
                indices = vgraph[i].indices
                idx_bool = vpt[i] < vpt[indices]
                cols.extend(indices[idx_bool])
                vals.extend(vgraph[i].data[idx_bool])
                rows.extend([i] * np.sum(idx_bool))
            vgraph = vals_to_csr(vals, rows, cols, shape=vgraph.shape)

            if 'final_cells' in self._adata.obs.keys() and is_categorical(self._adata.obs['final_cells']):
                #set_row_csr(vgraph, rows=np.where([isinstance(c, str) for c in self._adata.obs['final_cells']])[0])
                final_cells = self._adata.obs['final_cells'].cat.categories
                if isinstance(final_cells[0], str):
                    set_row_csr(vgraph, rows=np.where(clusters.values.isin(final_cells))[0])
            if 'end_points' in self._adata.obs.keys() and not is_categorical(self._adata.obs['end_points']):
                set_row_csr(vgraph, rows=np.where(self._adata.obs['end_points'] > .7)[0])

        if self.use_root_prior and 'root_cells' in self._adata.obs.keys():
            if is_categorical(self._adata.obs['root_cells']):
                # set_row_csr(vgraph, rows=np.where([isinstance(c, str) for c in self._adata.obs['root_cells']])[0])
                root_cells = self._adata.obs['root_cells'].cat.categories
                if isinstance(root_cells[0], str):
                    root = most_common_in_list(self._adata.obs['root_cells'])
                    vgraph[:, np.where(clusters.values == root)[0]] = 0
                    vgraph.eliminate_zeros()
            else:
                vgraph[:, np.where(self._adata.obs['root_cells'] > .7)[0]] = 0
                vgraph.eliminate_zeros()

        membership = self._adata.obs[self._groups_key].cat.codes.values
        g = get_igraph_from_adjacency(vgraph, directed=True)
        vc = igraph.VertexClustering(g, membership=membership)
        cg_full = vc.cluster_graph(combine_edges='sum')
        transitions = get_sparse_from_igraph(cg_full, weight_attr='weight')
        transitions = transitions - transitions.T
        transitions_conf = transitions.copy()
        transitions = transitions.tocoo()
        total_n = self._neighbors.n_neighbors * np.array(vc.sizes())
        for i, j, v in zip(transitions.row, transitions.col, transitions.data):
            reference = np.sqrt(total_n[i] * total_n[j])
            transitions_conf[i, j] = 0 if v < 0 else v / reference
        transitions_conf.eliminate_zeros()

        # remove non-confident direct paths if more confident indirect path is found.
        T = transitions_conf.A
        threshold = max(np.nanmin(np.nanmax(T / (T > 0), axis=0)) - 1e-6, .01)
        T *= (T > threshold)
        for i in range(len(T)):
            idx = T[i] > 0
            if np.any(idx):
                indirect = np.clip(T[idx], None, T[i][idx][:, None]).max(0)
                T[i, T[i] < indirect] = 0

        if self.minimum_spanning_tree:
            T_tmp = T.copy()
            T_num = T > 0
            T_sum = np.sum(T_num, 0)
            T_max = np.max(T_tmp)
            for i in range(len(T_tmp)):
                if T_sum[i] == 1:
                    T_tmp[np.where(T_num[:, i])[0][0], i] = T_max
            from scipy.sparse.csgraph import minimum_spanning_tree
            T_tmp = np.abs(minimum_spanning_tree(-T_tmp).A) > 0
            T = T_tmp * T

        transitions_conf = csr_matrix(T)
        self.transitions_confidence = transitions_conf.T

        # set threshold for minimal spanning tree.
        df = pd.DataFrame(T, index=cats, columns=cats)
        if root is not None:
            df.pop(root)
        self.threshold = max(np.nanmin(np.nanmax(df.values / (df.values > 0), axis=0)) - 1e-6, .01)


def paga(adata, groups=None, vkey='velocity', use_time_prior=True, use_root_prior=None,
         minimum_spanning_tree=True, copy=False):
    """PAGA graph with velocity-directed edges.

    Mapping out the coarse-grained connectivity structures of complex manifolds [Wolf19]_.
    By quantifying the connectivity of partitions (groups, clusters) of the
    single-cell graph, partition-based graph abstraction (PAGA) generates a much
    simpler abstracted graph (*PAGA graph*) of partitions, in which edge weights
    represent confidence in the presence of connections.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        An annotated data matrix.
    groups : key for categorical in `adata.obs`, optional (default: 'louvain')
        You can pass your predefined groups by choosing any categorical
        annotation of observations (`adata.obs`).
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    use_time_prior : `str` or bool, optional (default: True)
        Obs key for pseudo-time values. If True, 'velocity_pseudotime' is used if available.
    use_root_prior : `str` or bool, optional (default: True)
        Obs key for root/final states. If True, 'root_cells' and 'final_cells' are used if available.
    minimum_spanning_tree : bool, optional (default: True)
        Whether to prune the tree such that a path from A-to-B is removed if another more confident path exists.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy. Otherwise, perform computation inplace and return `None`.
    Returns
    -------
    **connectivities** : :class:`numpy.ndarray` (adata.uns['connectivities'])
        The full adjacency matrix of the abstracted graph, weights correspond to
        confidence in the connectivities of partitions.
    **connectivities_tree** : :class:`scipy.sparse.csr_matrix` (adata.uns['connectivities_tree'])
        The adjacency matrix of the tree-like subgraph that best explains the topology.
    **transitions_confidence** : :class:`scipy.sparse.csr_matrix` (adata.uns['transitions_confidence'])
        The adjacency matrix of the abstracted directed graph, weights correspond to
        confidence in the transitions between partitions.
    """
    if 'neighbors' not in adata.uns:
        raise ValueError('You need to run `pp.neighbors` first to compute a neighborhood graph.')
    adata = adata.copy() if copy else adata
    strings_to_categoricals(adata)

    if groups is None:
        groups = 'clusters' if 'clusters' in adata.obs.keys() else 'louvain' if 'louvain' in adata.obs.keys() else None
    elif groups == 'velocity_clusters' and 'velocity_clusters' not in adata.obs.keys():
        velocity_clusters(adata)
    if use_time_prior and not isinstance(use_time_prior, str):
        use_time_prior = 'velocity_pseudotime'
        if use_time_prior not in adata.obs.keys():
            velocity_pseudotime(adata, vkey=vkey)

    logg.info('running PAGA', r=True)
    paga = PAGA_tree(adata, groups, vkey=vkey, use_time_prior=use_time_prior, use_root_prior=use_root_prior,
                     minimum_spanning_tree=minimum_spanning_tree)

    if 'paga' not in adata.uns:
        adata.uns['paga'] = {}

    paga.compute_connectivities()
    adata.uns['paga']['connectivities'] = paga.connectivities
    adata.uns['paga']['connectivities_tree'] = paga.connectivities_tree
    adata.uns[groups + '_sizes'] = np.array(paga.ns)

    paga.compute_transitions()
    adata.uns['paga']['transitions_confidence'] = paga.transitions_confidence
    adata.uns['paga']['threshold'] = paga.threshold
    adata.uns['paga']['groups'] = groups

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n' +
              "    'paga/transitions_confidence', connectivities adjacency (adata.uns)\n"
              "    'paga/connectivities', connectivities adjacency (adata.uns)\n"
              "    'paga/connectivities_tree', connectivities subtree (adata.uns)")

    return adata if copy else None
