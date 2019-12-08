# This is adapted from https://github.com/theislab/paga
from .. import settings
from .. import logging as logg
from .utils import strings_to_categoricals
import numpy as np

try: from scanpy.tools.paga import PAGA
except ImportError:
    try: from scanpy.tools._paga import PAGA
    except ImportError: pass


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


class PAGA2(PAGA):
    def __init__(self, adata, groups, model='v1.2', vkey='velocity'):
        super().__init__(adata=adata, groups=groups, model=model)
        self.vkey = vkey

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
        # restore this at some point
        # if 'expected_n_edges_random' not in self._adata.uns['paga']:
        #     raise ValueError(
        #         'Before running PAGA with `use_rna_velocity=True`, run it with `False`.')
        import igraph
        g = get_igraph_from_adjacency(
            self._adata.uns[vkey].astype('bool'), directed=True)
        vc = igraph.VertexClustering(
            g, membership=self._adata.obs[self._groups_key].cat.codes.values)
        # set combine_edges to False if you want self loops
        cg_full = vc.cluster_graph(combine_edges='sum')
        transitions = get_sparse_from_igraph(cg_full, weight_attr='weight')
        transitions = transitions - transitions.T
        transitions_conf = transitions.copy()
        transitions = transitions.tocoo()
        total_n = self._neighbors.n_neighbors * np.array(vc.sizes())
        # total_n_sum = sum(total_n)
        # expected_n_edges_random = self._adata.uns['paga']['expected_n_edges_random']
        for i, j, v in zip(transitions.row, transitions.col, transitions.data):
            # if expected_n_edges_random[i, j] != 0:
            #     # factor 0.5 because of asymmetry
            #     reference = 0.5 * expected_n_edges_random[i, j]
            # else:
            #     # approximate
            #     reference = self._neighbors.n_neighbors * total_n[i] * total_n[j] / total_n_sum
            reference = np.sqrt(total_n[i] * total_n[j])
            transitions_conf[i, j] = 0 if v < 0 else v / reference
        transitions_conf.eliminate_zeros()
        # transpose in order to match convention of stochastic matrices
        # entry ij means transition from j to i
        self.transitions_confidence = transitions_conf.T


def paga(
        adata,
        vkey='velocity',
        groups=None,
        model='v1.2',
        copy=False):
    """Mapping out the coarse-grained connectivity structures of complex manifolds [Wolf19]_.
    By quantifying the connectivity of partitions (groups, clusters) of the
    single-cell graph, partition-based graph abstraction (PAGA) generates a much
    simpler abstracted graph (*PAGA graph*) of partitions, in which edge weights
    represent confidence in the presence of connections. By tresholding this
    confidence in :func:`~scanpy.pl.paga`, a much simpler representation of the
    manifold data is obtained, which is nonetheless faithful to the topology of
    the manifold.
    The confidence should be interpreted as the ratio of the actual versus the
    expected value of connetions under the null model of randomly connecting
    partitions. We do not provide a p-value as this null model does not
    precisely capture what one would consider "connected" in real data, hence it
    strongly overestimates the expected value. See an extensive discussion of
    this in [Wolf19]_.
    .. note::
        Note that you can use the result of :func:`~scanpy.pl.paga` in
        :func:`~scanpy.tl.umap` and :func:`~scanpy.tl.draw_graph` via
        `init_pos='paga'` to get single-cell embeddings that are typically more
        faithful to the global topology.
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        An annotated data matrix.
    groups : key for categorical in `adata.obs`, optional (default: 'louvain')
        You can pass your predefined groups by choosing any categorical
        annotation of observations (`adata.obs`).
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    model : {'v1.2', 'v1.0'}, optional (default: 'v1.2')
        The PAGA connectivity model.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy. Otherwise, perform
        computation inplace and return `None`.
    Returns
    -------
    **connectivities** : :class:`numpy.ndarray` (adata.uns['connectivities'])
        The full adjacency matrix of the abstracted graph, weights correspond to
        confidence in the connectivities of partitions.
    **connectivities_tree** : :class:`scipy.sparse.csr_matrix` (adata.uns['connectivities_tree'])
        The adjacency matrix of the tree-like subgraph that best explains
        the topology.
    Notes
    -----
    Together with a random walk-based distance measure
    (e.g. :func:`scanpy.tl.dpt`) this generates a partial coordinatization of
    data useful for exploring and explaining its variation.
    See Also
    --------
    pl.paga
    pl.paga_path
    pl.paga_compare
    """
    if groups is None:
        groups = 'clusters' if 'clusters' in adata.obs.keys() else 'louvain' if 'louvain' in adata.obs.keys() else 'grey'
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    adata = adata.copy() if copy else adata
    strings_to_categoricals(adata)
    start = logg.info('running PAGA')
    paga = PAGA2(adata, groups, model=model, vkey=vkey)
    # only add if not present
    if 'paga' not in adata.uns:
        adata.uns['paga'] = {}
    paga.compute_connectivities()
    adata.uns['paga']['connectivities'] = paga.connectivities
    adata.uns['paga']['connectivities_tree'] = paga.connectivities_tree
    # adata.uns['paga']['expected_n_edges_random'] = paga.expected_n_edges_random
    adata.uns[groups + '_sizes'] = np.array(paga.ns)
    paga.compute_transitions()
    adata.uns['paga']['transitions_confidence'] = paga.transitions_confidence
    # adata.uns['paga']['transitions_ttest'] = paga.transitions_ttest
    adata.uns['paga']['groups'] = groups
    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint('added\n' +
              "    'paga/transitions_confidence', connectivities adjacency (adata.uns)\n"
              "    'paga/connectivities', connectivities adjacency (adata.uns)\n"
              "    'paga/connectivities_tree', connectivities subtree (adata.uns)")

    return adata if copy else None
