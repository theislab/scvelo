import os

import numpy as np
from scipy.sparse import coo_matrix, issparse

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import get_n_jobs, l2_norm, parallelize
from scvelo.preprocessing.moments import get_moments
from scvelo.preprocessing.neighbors import (
    get_n_neighs,
    get_neighs,
    neighbors,
    pca,
    verify_neighbors,
)
from .utils import cosine_correlation, get_indices, get_iterative_indices
from .velocity import velocity


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
    def __init__(
        self,
        adata,
        vkey="velocity",
        xkey="Ms",
        tkey=None,
        basis=None,
        n_neighbors=None,
        sqrt_transform=None,
        n_recurse_neighbors=None,
        random_neighbors_at_max=None,
        gene_subset=None,
        approx=None,
        report=False,
        compute_uncertainties=None,
        mode_neighbors="distances",
    ):

        subset = np.ones(adata.n_vars, bool)
        if gene_subset is not None:
            var_names_subset = adata.var_names.isin(gene_subset)
            subset &= var_names_subset if len(var_names_subset) > 0 else gene_subset
        elif f"{vkey}_genes" in adata.var.keys():
            subset &= np.array(adata.var[f"{vkey}_genes"].values, dtype=bool)

        xkey = xkey if xkey in adata.layers.keys() else "spliced"

        X = np.array(
            adata.layers[xkey].A[:, subset]
            if issparse(adata.layers[xkey])
            else adata.layers[xkey][:, subset]
        )
        V = np.array(
            adata.layers[vkey].A[:, subset]
            if issparse(adata.layers[vkey])
            else adata.layers[vkey][:, subset]
        )

        nans = np.isnan(np.sum(V, axis=0))
        if np.any(nans):
            X = X[:, ~nans]
            V = V[:, ~nans]

        if approx is True and X.shape[1] > 100:
            X_pca, PCs, _, _ = pca(X, n_comps=30, svd_solver="arpack", return_info=True)
            self.X = np.array(X_pca, dtype=np.float32)
            self.V = (V - V.mean(0)).dot(PCs.T)
            self.V[V.sum(1) == 0] = 0
        else:
            self.X = np.array(X, dtype=np.float32)
            self.V = np.array(V, dtype=np.float32)
        self.V_raw = np.array(self.V)

        self.sqrt_transform = sqrt_transform
        uns_key = f"{vkey}_params"
        if self.sqrt_transform is None:
            if uns_key in adata.uns.keys() and "mode" in adata.uns[uns_key]:
                self.sqrt_transform = adata.uns[uns_key]["mode"] == "stochastic"
        if self.sqrt_transform:
            self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= np.nanmean(self.V, axis=1)[:, None]

        self.n_recurse_neighbors = n_recurse_neighbors
        if self.n_recurse_neighbors is None:
            if n_neighbors is not None or mode_neighbors == "connectivities":
                self.n_recurse_neighbors = 1
            else:
                self.n_recurse_neighbors = 2

        if "neighbors" not in adata.uns.keys():
            neighbors(adata)
        if np.min((get_neighs(adata, "distances") > 0).sum(1).A1) == 0:
            raise ValueError(
                "Your neighbor graph seems to be corrupted. "
                "Consider recomputing via pp.neighbors."
            )
        if n_neighbors is None or n_neighbors <= get_n_neighs(adata):
            self.indices = get_indices(
                dist=get_neighs(adata, "distances"),
                n_neighbors=n_neighbors,
                mode_neighbors=mode_neighbors,
            )[0]
        else:
            if basis is None:
                basis_keys = ["X_pca", "X_tsne", "X_umap"]
                basis = [key for key in basis_keys if key in adata.obsm.keys()][-1]
            elif f"X_{basis}" in adata.obsm.keys():
                basis = f"X_{basis}"

            if isinstance(approx, str) and approx in adata.obsm.keys():
                from sklearn.neighbors import NearestNeighbors

                neighs = NearestNeighbors(n_neighbors=n_neighbors + 1)
                neighs.fit(adata.obsm[approx])
                self.indices = neighs.kneighbors_graph(
                    mode="connectivity"
                ).indices.reshape((-1, n_neighbors + 1))
            else:
                from scvelo import Neighbors

                neighs = Neighbors(adata)
                neighs.compute_neighbors(
                    n_neighbors=n_neighbors, use_rep=basis, n_pcs=10
                )
                self.indices = get_indices(
                    dist=neighs.distances, mode_neighbors=mode_neighbors
                )[0]

        self.max_neighs = random_neighbors_at_max

        gkey, gkey_ = f"{vkey}_graph", f"{vkey}_graph_neg"
        self.graph = adata.uns[gkey] if gkey in adata.uns.keys() else []
        self.graph_neg = adata.uns[gkey_] if gkey_ in adata.uns.keys() else []

        if tkey in adata.obs.keys():
            self.t0 = adata.obs[tkey].copy()
            init = min(self.t0) if isinstance(min(self.t0), int) else 0
            self.t0.cat.categories = np.arange(init, len(self.t0.cat.categories))
            self.t1 = self.t0.copy()
            self.t1.cat.categories = self.t0.cat.categories + 1
        else:
            self.t0 = None

        self.compute_uncertainties = compute_uncertainties
        self.uncertainties = None
        self.self_prob = None
        self.report = report
        self.adata = adata

    def compute_cosines(self, n_jobs=None, backend="loky"):
        n_jobs = get_n_jobs(n_jobs=n_jobs)

        n_obs = self.X.shape[0]

        # TODO: Use batches and vectorize calculation of dX in self._calculate_cosines
        res = parallelize(
            self._compute_cosines,
            range(self.X.shape[0]),
            n_jobs=n_jobs,
            unit="cells",
            backend=backend,
        )()
        uncertainties, vals, rows, cols = map(_flatten, zip(*res))

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 0

        self.graph, self.graph_neg = vals_to_csr(
            vals, rows, cols, shape=(n_obs, n_obs), split_negative=True
        )
        if self.compute_uncertainties:
            uncertainties = np.hstack(uncertainties)
            uncertainties[np.isnan(uncertainties)] = 0
            self.uncertainties = vals_to_csr(
                uncertainties, rows, cols, shape=(n_obs, n_obs), split_negative=False
            )
            self.uncertainties.eliminate_zeros()

        confidence = self.graph.max(1).A.flatten()
        self.self_prob = np.clip(np.percentile(confidence, 98) - confidence, 0, 1)

    def _compute_cosines(self, obs_idx, queue):
        vals, rows, cols, uncertainties = [], [], [], []
        if self.compute_uncertainties:
            moments = get_moments(self.adata, np.sign(self.V_raw), second_order=True)

        for obs_id in obs_idx:
            if self.V[obs_id].max() != 0 or self.V[obs_id].min() != 0:
                neighs_idx = get_iterative_indices(
                    self.indices, obs_id, self.n_recurse_neighbors, self.max_neighs
                )

                if self.t0 is not None:
                    t0, t1 = self.t0[obs_id], self.t1[obs_id]
                    if t0 >= 0 and t1 > 0:
                        t1_idx = np.where(self.t0 == t1)[0]
                        if len(t1_idx) > len(neighs_idx):
                            t1_idx = np.random.choice(
                                t1_idx, len(neighs_idx), replace=False
                            )
                        if len(t1_idx) > 0:
                            neighs_idx = np.unique(np.concatenate([neighs_idx, t1_idx]))

                dX = self.X[neighs_idx] - self.X[obs_id, None]  # 60% of runtime
                if self.sqrt_transform:
                    dX = np.sqrt(np.abs(dX)) * np.sign(dX)
                val = cosine_correlation(dX, self.V[obs_id])  # 40% of runtime

                if self.compute_uncertainties:
                    dX /= l2_norm(dX)[:, None]
                    uncertainties.extend(
                        np.nansum(dX ** 2 * moments[obs_id][None, :], 1)
                    )

                vals.extend(val)
                rows.extend(np.ones(len(neighs_idx)) * obs_id)
                cols.extend(neighs_idx)

            if queue is not None:
                queue.put(1)

        if queue is not None:
            queue.put(None)

        return uncertainties, vals, rows, cols


def _flatten(iterable):
    return [i for it in iterable for i in it]


def velocity_graph(
    data,
    vkey="velocity",
    xkey="Ms",
    tkey=None,
    basis=None,
    n_neighbors=None,
    n_recurse_neighbors=None,
    random_neighbors_at_max=None,
    sqrt_transform=None,
    variance_stabilization=None,
    gene_subset=None,
    compute_uncertainties=None,
    approx=None,
    mode_neighbors="distances",
    copy=False,
    n_jobs=None,
    backend="loky",
):
    """Computes velocity graph based on cosine similarities.

    The cosine similarities are computed between velocities and potential cell state
    transitions, i.e. it measures how well a corresponding change in gene expression
    :math:`\\delta_{ij} = x_j - x_i` matches the predicted change according to the
    velocity vector :math:`\\nu_i`,

    .. math::
        \\pi_{ij} = \\cos\\angle(\\delta_{ij}, \\nu_i)
        = \\frac{\\delta_{ij}^T \\nu_i}{\\left\\lVert\\delta_{ij}\\right\\rVert
        \\left\\lVert \\nu_i \\right\\rVert}.

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
    n_recurse_neighbors: `int` (default: `None`)
        Number of recursions for neighbors search. Defaults to
        2 if mode_neighbors is 'distances', and 1 if mode_neighbors is 'connectivities'.
    random_neighbors_at_max: `int` or `None` (default: `None`)
        If number of iterative neighbors for an individual cell is higher than this
        threshold, a random selection of such are chosen as reference neighbors.
    sqrt_transform: `bool` (default: `False`)
        Whether to variance-transform the cell states changes
        and velocities before computing cosine similarities.
    gene_subset: `list` of `str`, subset of adata.var_names or `None`(default: `None`)
        Subset of genes to compute velocity graph on exclusively.
    compute_uncertainties: `bool` (default: `None`)
        Whether to compute uncertainties along with cosine correlation.
    approx: `bool` or `None` (default: `None`)
        If True, first 30 pc's are used instead of the full count matrix
    mode_neighbors: 'str' (default: `'distances'`)
        Determines the type of KNN graph used. Options are 'distances' or
        'connectivities'. The latter yields a symmetric graph.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.
    n_jobs: `int` or `None` (default: `None`)
        Number of parallel jobs.
    backend: `str` (default: "loky")
        Backend used for multiprocessing. See :class:`joblib.Parallel` for valid
        options.

    Returns
    -------
    velocity_graph: `.uns`
        sparse matrix with correlations of cell state transitions with velocities
    """

    adata = data.copy() if copy else data
    verify_neighbors(adata)
    if vkey not in adata.layers.keys():
        velocity(adata, vkey=vkey)
    if sqrt_transform is None:
        sqrt_transform = variance_stabilization

    vgraph = VelocityGraph(
        adata,
        vkey=vkey,
        xkey=xkey,
        tkey=tkey,
        basis=basis,
        n_neighbors=n_neighbors,
        approx=approx,
        n_recurse_neighbors=n_recurse_neighbors,
        random_neighbors_at_max=random_neighbors_at_max,
        sqrt_transform=sqrt_transform,
        gene_subset=gene_subset,
        compute_uncertainties=compute_uncertainties,
        report=True,
        mode_neighbors=mode_neighbors,
    )

    if isinstance(basis, str):
        logg.warn(
            f"The velocity graph is computed on {basis} embedding coordinates.\n"
            f"        Consider computing the graph in an unbiased manner \n"
            f"        on full expression space by not specifying basis.\n"
        )

    n_jobs = get_n_jobs(n_jobs=n_jobs)
    logg.info(
        f"computing velocity graph (using {n_jobs}/{os.cpu_count()} cores)", r=True
    )
    vgraph.compute_cosines(n_jobs=n_jobs, backend=backend)

    adata.uns[f"{vkey}_graph"] = vgraph.graph
    adata.uns[f"{vkey}_graph_neg"] = vgraph.graph_neg

    if vgraph.uncertainties is not None:
        adata.uns[f"{vkey}_graph_uncertainties"] = vgraph.uncertainties

    adata.obs[f"{vkey}_self_transition"] = vgraph.self_prob

    if f"{vkey}_params" in adata.uns.keys():
        if "embeddings" in adata.uns[f"{vkey}_params"]:
            del adata.uns[f"{vkey}_params"]["embeddings"]
    else:
        adata.uns[f"{vkey}_params"] = {}
    adata.uns[f"{vkey}_params"]["mode_neighbors"] = mode_neighbors
    adata.uns[f"{vkey}_params"]["n_recurse_neighbors"] = vgraph.n_recurse_neighbors

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        f"    '{vkey}_graph', sparse matrix with cosine correlations (adata.uns)"
    )

    return adata if copy else None
