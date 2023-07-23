import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from scanpy.tools._paga import PAGA

# This is adapted from https://github.com/theislab/paga
from scvelo import logging as logg
from scvelo import settings
from .rank_velocity_genes import velocity_clusters
from .utils import strings_to_categoricals
from .velocity_graph import vals_to_csr
from .velocity_pseudotime import velocity_pseudotime


# TODO: Finish docstrings
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
    g.es["weight"] = weights
    if g.vcount() != adjacency.shape[0]:
        logg.warn(
            f"The constructed graph has only {g.vcount()} nodes. "
            "Your adjacency matrix contained redundant nodes."
        )
    return g


# TODO: Add docstrings
def get_sparse_from_igraph(graph, weight_attr=None):
    """TODO."""
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


# TODO: Finish docstrings
def set_row_csr(csr, rows, value=0):
    """Set all nonzero elements to the given value. Useful to set to 0 mostly."""
    for row in rows:
        start = csr.indptr[row]
        end = csr.indptr[row + 1]
        csr.data[start:end] = value
    if value == 0:
        csr.eliminate_zeros()


# TODO: Add docstrings
class PAGA_tree(PAGA):
    """TODO."""

    def __init__(
        self,
        adata,
        groups=None,
        vkey=None,
        use_time_prior=None,
        root_key=None,
        end_key=None,
        threshold_root_end_prior=None,
        minimum_spanning_tree=None,
    ):
        super().__init__(adata=adata, groups=groups, model="v1.2")
        self.groups = groups
        self.vkey = vkey
        self.use_time_prior = use_time_prior
        self.root_key = root_key
        self.end_key = end_key
        self.threshold_root_end_prior = threshold_root_end_prior
        if self.threshold_root_end_prior is None:
            self.threshold_root_end_prior = 0.9
        self.minimum_spanning_tree = minimum_spanning_tree

    # TODO: Add docstrings
    def compute_transitions(self):
        """TODO."""
        try:
            import igraph
        except ImportError:
            raise ImportError("To run paga, you need to install `pip install igraph`")
        vkey = f"{self.vkey}_graph"
        if vkey not in self._adata.uns:
            raise ValueError(
                "The passed AnnData needs to have an `uns` annotation "
                "with key 'velocity_graph' - a sparse matrix from RNA velocity."
            )
        if self._adata.uns[vkey].shape != (self._adata.n_obs, self._adata.n_obs):
            raise ValueError(
                f"The passed 'velocity_graph' has shape {self._adata.uns[vkey].shape} "
                f"but shoud have shape {(self._adata.n_obs, self._adata.n_obs)}"
            )

        clusters = self._adata.obs[self.groups]
        cats = clusters.cat.categories
        vgraph = self._adata.uns[vkey] > 0.1
        time_prior = self.use_time_prior

        if isinstance(time_prior, str) and time_prior in self._adata.obs.keys():
            vpt = self._adata.obs[time_prior].values
            vpt_mean = self._adata.obs.groupby(self.groups)[time_prior].mean()
            vpt_means = np.array([vpt_mean[cat] for cat in clusters])
            rows, cols, vals = [], [], []
            for i in range(vgraph.shape[0]):
                indices = vgraph[i].indices
                idx_bool = vpt[i] < vpt[indices]
                idx_bool &= vpt_means[indices] > vpt_means[i] - 0.1
                cols.extend(indices[idx_bool])
                vals.extend(vgraph[i].data[idx_bool])
                rows.extend([i] * np.sum(idx_bool))
            vgraph = vals_to_csr(vals, rows, cols, shape=vgraph.shape)

        lb = self.threshold_root_end_prior  # cells to be consider as terminal states
        if isinstance(self.end_key, str) and self.end_key in self._adata.obs.keys():
            set_row_csr(vgraph, rows=np.where(self._adata.obs[self.end_key] > lb)[0])
        if isinstance(self.root_key, str) and self.root_key in self._adata.obs.keys():
            vgraph[:, np.where(self._adata.obs[self.root_key] > lb)[0]] = 0
            vgraph.eliminate_zeros()

        membership = self._adata.obs[self.groups].cat.codes.values
        g = get_igraph_from_adjacency(vgraph, directed=True)
        vc = igraph.VertexClustering(g, membership=membership)
        cg_full = vc.cluster_graph(combine_edges="sum")
        transitions = get_sparse_from_igraph(cg_full, weight_attr="weight")
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
        threshold = max(np.nanmin(np.nanmax(T / (T > 0), axis=0)) - 1e-6, 0.01)
        T *= T > threshold
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
        self.threshold = np.nanmin(np.nanmax(df.values / (df.values > 0), axis=0))
        self.threshold = max(self.threshold - 1e-6, 0.01)


def paga(
    adata,
    groups=None,
    vkey="velocity",
    use_time_prior=True,
    root_key=None,
    end_key=None,
    threshold_root_end_prior=None,
    minimum_spanning_tree=True,
    copy=False,
):
    """PAGA graph with velocity-directed edges.

    Mapping out the coarse-grained connectivity structures of complex manifolds
    :cite:p:`Wolf19`. By quantifying the connectivity of partitions (groups, clusters) of the
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
        Obs key for pseudo-time values.
        If True, 'velocity_pseudotime' is used if available.
    root_key : `str` or bool, optional (default: None)
        Obs key for root states.
    end_key : `str` or bool, optional (default: None)
        Obs key for end states.
    threshold_root_end_prior : `float` (default: 0.9)
        Threshold for root and final states priors, to be in the range of [0,1].
        Values above the threshold will be considered as terminal and included as prior.
    minimum_spanning_tree : bool, optional (default: True)
        Whether to prune the tree such that a path from A-to-B
        is removed if another more confident path exists.
    copy : `bool`, optional (default: `False`)
        Copy `adata` before computation and return a copy.
        Otherwise, perform computation inplace and return `None`.

    Returns
    -------
    connectivities: `.uns`
        The full adjacency matrix of the abstracted graph, weights correspond to
        confidence in the connectivities of partitions.
    connectivities_tree: `.uns`
        The adjacency matrix of the tree-like subgraph that best explains the topology.
    transitions_confidence: `.uns`
        The adjacency matrix of the abstracted directed graph, weights correspond to
        confidence in the transitions between partitions.
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "You need to run `pp.neighbors` first to compute a neighborhood graph."
        )

    adata = adata.copy() if copy else adata
    strings_to_categoricals(adata)

    if groups is None:
        groups = (
            "clusters"
            if "clusters" in adata.obs.keys()
            else "louvain"
            if "louvain" in adata.obs.keys()
            else None
        )
    elif groups == "velocity_clusters" and "velocity_clusters" not in adata.obs.keys():
        velocity_clusters(adata)
    if use_time_prior and not isinstance(use_time_prior, str):
        use_time_prior = "velocity_pseudotime"
        if use_time_prior not in adata.obs.keys():
            velocity_pseudotime(adata, vkey=vkey, root_key=root_key, end_key=end_key)

    priors = [p for p in [use_time_prior, root_key, end_key] if p in adata.obs.keys()]
    logg.info(
        "running PAGA",
        f"using priors: {priors}" if len(priors) > 0 else "",
        r=True,
    )
    paga = PAGA_tree(
        adata,
        groups,
        vkey=vkey,
        use_time_prior=use_time_prior,
        root_key=root_key,
        end_key=end_key,
        threshold_root_end_prior=threshold_root_end_prior,
        minimum_spanning_tree=minimum_spanning_tree,
    )

    if "paga" not in adata.uns:
        adata.uns["paga"] = {}

    paga.compute_connectivities()
    adata.uns["paga"]["connectivities"] = paga.connectivities
    adata.uns["paga"]["connectivities_tree"] = paga.connectivities_tree
    adata.uns[f"{groups}_sizes"] = np.array(paga.ns)

    paga.compute_transitions()
    adata.uns["paga"]["transitions_confidence"] = paga.transitions_confidence
    adata.uns["paga"]["threshold"] = paga.threshold
    adata.uns["paga"]["groups"] = groups

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added\n" + "    'paga/connectivities', connectivities adjacency (adata.uns)\n"
        "    'paga/connectivities_tree', connectivities subtree (adata.uns)\n"
        "    'paga/transitions_confidence', velocity transitions (adata.uns)"
    )

    return adata if copy else None
