import warnings
from collections import Counter
from typing import Dict, Literal, Optional

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, issparse

from anndata import AnnData
from scanpy import Neighbors
from scanpy.preprocessing import pca

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import get_initial_size


# TODO: Add docstrings
def _get_hnsw_neighbors(
    adata: AnnData,
    use_rep: str,
    n_pcs: int,
    n_neighbors: int,
    num_threads: int,
    **kwargs,
):
    X = adata.X if use_rep == "X" else adata.obsm[use_rep]
    neighbors = FastNeighbors(n_neighbors=n_neighbors, num_threads=num_threads)
    neighbors.fit(X if n_pcs is None else X[:, :n_pcs], **kwargs)

    return neighbors


# TODO: Add docstrings
def _get_scanpy_neighbors(adata: AnnData, **kwargs):
    logg.switch_verbosity("off", module="scanpy")
    with warnings.catch_warnings():  # ignore numba warning (umap/issues/252)
        warnings.simplefilter("ignore")
        neighbors = Neighbors(adata)
        neighbors.compute_neighbors(write_knn_indices=True, **kwargs)
    logg.switch_verbosity("on", module="scanpy")

    return neighbors


# TODO: Add docstrings
def _get_sklearn_neighbors(
    adata: AnnData, use_rep: str, n_pcs: Optional[int], n_neighbors: int, **kwargs
):
    from sklearn.neighbors import NearestNeighbors

    # TODO: Use `scv.core.get_modality`
    X = adata.X if use_rep == "X" else adata.obsm[use_rep]
    neighbors = NearestNeighbors(n_neighbors=n_neighbors - 1, **kwargs)
    neighbors.fit(X if n_pcs is None else X[:, :n_pcs])
    knn_distances, neighbors.knn_indices = neighbors.kneighbors()
    knn_distances, neighbors.knn_indices = set_diagonal(
        knn_distances, neighbors.knn_indices
    )
    neighbors.distances, neighbors.connectivities = compute_connectivities_umap(
        neighbors.knn_indices, knn_distances, X.shape[0], n_neighbors=n_neighbors
    )

    return neighbors


# TODO: Add docstrings
def _get_rep(adata: AnnData, use_rep: str, n_pcs: int):
    if use_rep is None:
        rep = "X" if adata.n_vars < 50 or n_pcs == 0 else "X_pca"
    elif use_rep not in adata.obsm.keys() and f"X_{use_rep}" in adata.obsm.keys():
        rep = f"X_{use_rep}"
    else:
        rep = use_rep

    if (rep == "X") and (n_pcs is not None):
        logg.warn(
            f"Unexpected pair of parameters: `use_rep='X'` but `n_pcs={n_pcs}`. "
            f"This will only consider the frist {n_pcs} variables when calculating the "
            "neighbor graph. To use all of `X`, pass `n_pcs=None`."
        )

    return rep


# TODO: Add docstrings
def _set_neighbors_data(
    adata: AnnData,
    neighbors,
    n_neighbors: int,
    method: str,
    metric: str,
    n_pcs: int,
    use_rep: str,
):
    adata.uns["neighbors"] = {}
    # TODO: Remove except statement. `n_obs` x `n_obs` arrays need to be written to
    # `AnnData.obsp`.
    try:
        adata.obsp["distances"] = neighbors.distances
        adata.obsp["connectivities"] = neighbors.connectivities
        adata.uns["neighbors"]["connectivities_key"] = "connectivities"
        adata.uns["neighbors"]["distances_key"] = "distances"
    except ValueError as e:
        logg.warning(f"Could not write neighbors to `AnnData.obsp`: {e}")
        adata.uns["neighbors"]["distances"] = neighbors.distances
        adata.uns["neighbors"]["connectivities"] = neighbors.connectivities

    if hasattr(neighbors, "knn_indices"):
        adata.uns["neighbors"]["indices"] = neighbors.knn_indices
    adata.uns["neighbors"]["params"] = {
        "n_neighbors": n_neighbors,
        "method": method,
        "metric": metric,
        "n_pcs": n_pcs,
        "use_rep": use_rep,
    }


# TODO: Add docstrings
def _set_pca(adata, n_pcs: Optional[int], use_highly_variable: bool):
    if (
        "X_pca" not in adata.obsm.keys()
        or n_pcs is not None
        and n_pcs > adata.obsm["X_pca"].shape[1]
    ):
        if use_highly_variable and "highly_variable" in adata.var.keys():
            n_vars = np.sum(adata.var["highly_variable"])
        else:
            n_vars = adata.n_vars

        n_comps = min(30 if n_pcs is None else n_pcs, n_vars - 1, adata.n_obs - 1)
        use_highly_variable &= "highly_variable" in adata.var.keys()
        pca(
            adata,
            n_comps=n_comps,
            use_highly_variable=use_highly_variable,
            svd_solver="arpack",
        )
    elif n_pcs is None and adata.obsm["X_pca"].shape[1] < 10:
        logg.warn(
            f"Neighbors are computed on {adata.obsm['X_pca'].shape[1]} "
            "principal components only."
        )


def neighbors(
    adata: AnnData,
    n_neighbors: int = 30,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    use_highly_variable: bool = True,
    knn: bool = True,
    random_state: int = 0,
    method: Literal["umap", "sklearn", "hnsw", "gauss", "rapids"] = "umap",
    metric: str = "euclidean",
    metric_kwds: Optional[Dict] = None,
    num_threads: int = -1,
    copy: bool = False,
):
    """Compute a neighborhood graph of observations.

    The neighbor graph methods (umap, hnsw, sklearn) only differ in runtime and
    yield the same result as Scanpy :cite:p:`Wolf18`. Connectivities are computed with
    adaptive kernel width as proposed in Haghverdi et al. 2016 (doi:10.1038/nmeth.3971).

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
        Use the indicated representation. If `None`, the representation is chosen
        automatically: for .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.
    use_highly_variable: `bool` (default: True)
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    random_state
        A numpy random seed.
    method : {{'umap', 'hnsw', 'sklearn'}}  (default: `'umap'`)
        Method to compute neighbors, only differs in runtime.
        The 'hnsw' method is most efficient and requires to `pip install hnswlib`.
        Connectivities are computed with adaptive kernel.
    metric
        A known metric’s name or a callable that returns a distance.
    metric_kwds
        Options for the metric.
    num_threads
        Number of threads to be used (for runtime).
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    connectivities : `.obsp`
        Sparse weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    distances : `.obsp`
        Sparse matrix of distances for each pair of neighbors.
    """
    adata = adata.copy() if copy else adata

    use_rep = _get_rep(adata=adata, use_rep=use_rep, n_pcs=n_pcs)

    if use_rep == "X_pca":
        _set_pca(adata=adata, n_pcs=n_pcs, use_highly_variable=use_highly_variable)

        n_duplicate_cells = len(get_duplicate_cells(adata))
        if n_duplicate_cells > 0:
            logg.warn(
                f"You seem to have {n_duplicate_cells} duplicate cells in your data.",
                "Consider removing these via pp.remove_duplicate_cells.",
            )

    if metric_kwds is None:
        metric_kwds = {}

    logg.info("computing neighbors", r=True)

    if method == "sklearn":
        neighbors = _get_sklearn_neighbors(
            adata=adata,
            use_rep=use_rep,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            metric=metric,
            metric_params=metric_kwds,
            n_jobs=num_threads,
        )
    elif method == "hnsw":
        neighbors = _get_hnsw_neighbors(
            adata=adata,
            use_rep=use_rep,
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            num_threads=num_threads,
            metric=metric,
            random_state=random_state,
            **metric_kwds,
        )
    elif method in ["umap", "gauss", "rapids"]:
        neighbors = _get_scanpy_neighbors(
            adata=adata,
            n_neighbors=n_neighbors,
            knn=knn,
            n_pcs=n_pcs,
            method=method,
            use_rep=use_rep,
            random_state=random_state,
            metric=metric,
            metric_kwds=metric_kwds,
        )
    else:
        raise ValueError(
            f"Provided `method={method}`. Admissible values are `'umap'`, `'sklearn'`, "
            "`'hnsw'`, `'gauss'`, and `'rapids'`."
        )

    _set_neighbors_data(
        adata=adata,
        neighbors=neighbors,
        n_neighbors=n_neighbors,
        method=method,
        metric=metric,
        n_pcs=n_pcs,
        use_rep=use_rep,
    )

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        "    'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)"
    )

    return adata if copy else None


# TODO: Add docstrings
class FastNeighbors:
    """TODO."""

    def __init__(self, n_neighbors=30, num_threads=-1):
        self.n_neighbors = n_neighbors
        self.num_threads = num_threads
        self.knn_indices, self.knn_distances = None, None
        self.distances, self.connectivities = None, None

    # TODO: Add docstrings
    def fit(self, X, metric="l2", M=16, ef=100, ef_construction=100, random_state=0):
        """TODO."""
        try:
            import hnswlib
        except ImportError:
            print(
                "In order to use fast approx neighbor search, "
                "you need to `pip install hnswlib`\n"
            )

        ef_c, ef = max(ef_construction, self.n_neighbors), max(self.n_neighbors, ef)
        metric = "l2" if metric == "euclidean" else metric

        X = X.A if issparse(X) else X
        ns, dim = X.shape

        knn = hnswlib.Index(space=metric, dim=dim)
        knn.init_index(
            max_elements=ns, ef_construction=ef_c, M=M, random_seed=random_state
        )
        knn.add_items(X)
        knn.set_ef(ef)

        knn_indices, knn_distances = knn.knn_query(
            X, k=self.n_neighbors, num_threads=self.num_threads
        )

        n_neighbors = self.n_neighbors

        if metric == "l2":
            knn_distances = np.sqrt(knn_distances)

        self.distances, self.connectivities = compute_connectivities_umap(
            knn_indices, knn_distances, ns, n_neighbors
        )
        self.knn_indices = knn_indices


# TODO: Add docstrings
def set_diagonal(knn_distances, knn_indices, remove_diag=False):
    """TODO."""
    if remove_diag and knn_distances[0, 0] == 0:
        knn_distances = knn_distances[:, 1:]
        knn_indices = knn_indices[:, 1:].astype(int)
    elif knn_distances[0, 0] != 0:
        knn_distances = np.hstack(
            [np.zeros(len(knn_distances))[:, None], knn_distances]
        )
        knn_indices = np.array(
            np.hstack([np.arange(len(knn_indices), dtype=int)[:, None], knn_indices]),
            dtype=int,
        )
    return knn_distances, knn_indices


# TODO: Add docstrings
def select_distances(dist, n_neighbors=None):
    """TODO."""
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if issparse(D) else (D > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()
    return D


# TODO: Add docstrings
def select_connectivities(connectivities, n_neighbors=None):
    """TODO."""
    C = connectivities.copy()
    n_counts = (C > 0).sum(1).A1 if issparse(C) else (C > 0).sum(1)
    n_neighbors = (
        n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    )
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = C.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[::-1][n_neighbors:]
        dat[rm_idx] = 0
    C.eliminate_zeros()
    return C


# TODO: Add docstrings
def get_neighs(adata, mode="distances"):
    """TODO."""
    if hasattr(adata, "obsp") and mode in adata.obsp.keys():
        return adata.obsp[mode]
    elif "neighbors" in adata.uns.keys() and mode in adata.uns["neighbors"]:
        return adata.uns["neighbors"][mode]
    else:
        raise ValueError("The selected mode is not valid.")


# TODO: Add docstrings
def get_n_neighs(adata):
    """TODO."""
    return adata.uns.get("neighbors", {}).get("params", {}).get("n_neighbors", 0)


# TODO: Add docstrings
def verify_neighbors(adata):
    """TODO."""
    valid = "neighbors" in adata.uns.keys() and "params" in adata.uns["neighbors"]
    if valid:
        n_neighs = (get_neighs(adata, "distances") > 0).sum(1)
        # test whether the graph is corrupted
        valid = n_neighs.min() * 2 > n_neighs.max()
    if not valid:
        logg.warn(
            "The neighbor graph has an unexpected format "
            "(e.g. computed outside scvelo) \n"
            "or is corrupted (e.g. due to subsetting). "
            "Consider recomputing with `pp.neighbors`."
        )


# TODO: Add docstrings
def neighbors_to_be_recomputed(adata, n_neighbors=None):  # deprecated
    """TODO."""
    # check whether neighbors graph is disrupted or has insufficient number of neighbors
    invalid_neighs = (
        "neighbors" not in adata.uns.keys()
        or "params" not in adata.uns["neighbors"]
        or (n_neighbors is not None and n_neighbors > get_n_neighs(adata))
    )
    if invalid_neighs:
        return True
    else:
        n_neighs = (get_neighs(adata, "distances") > 0).sum(1)
        return n_neighs.max() * 0.1 > n_neighs.min()


# TODO: Add docstrings
def get_connectivities(
    adata, mode="connectivities", n_neighbors=None, recurse_neighbors=False
):
    """TODO."""
    if "neighbors" in adata.uns.keys():
        C = get_neighs(adata, mode)
        if n_neighbors is not None and n_neighbors < get_n_neighs(adata):
            if mode == "connectivities":
                C = select_connectivities(C, n_neighbors)
            else:
                C = select_distances(C, n_neighbors)
        connectivities = C > 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            connectivities.setdiag(1)
            if recurse_neighbors:
                connectivities += connectivities.dot(connectivities * 0.5)
                connectivities.data = np.clip(connectivities.data, 0, 1)
            connectivities = connectivities.multiply(1.0 / connectivities.sum(1))
        return connectivities.tocsr().astype(np.float32)
    else:
        return None


# TODO: Add docstrings
def get_csr_from_indices(knn_indices, knn_dists, n_obs, n_neighbors):
    """TODO."""
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # we didn't get the full knn for i
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


# TODO: Finish docstrings
def compute_connectivities_umap(
    knn_indices,
    knn_dists,
    n_obs,
    n_neighbors,
    set_op_mix_ratio=1.0,
    local_connectivity=1.0,
):
    """Computes fuzzy simplical set associated with data.

    This is from umap.fuzzy_simplicial_set :cite:p:`McInnes18`.
    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    from umap.umap_ import fuzzy_simplicial_set

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(
        X,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if isinstance(connectivities, tuple):  # umap returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    distances = get_csr_from_indices(knn_indices, knn_dists, n_obs, n_neighbors)

    return distances, connectivities.tocsr()


# TODO: Add docstrings
def get_duplicate_cells(data):
    """TODO."""
    if isinstance(data, AnnData):
        X = data.X
        lst = list(np.sum(np.abs(data.obsm["X_pca"]), 1) + get_initial_size(data))
    else:
        X = data
        lst = list(np.sum(X, 1).A1 if issparse(X) else np.sum(X, 1))

    idx_dup = []
    if len(set(lst)) < len(lst):
        vals = [val for val, count in Counter(lst).items() if count > 1]
        idx_dup = np.where(pd.Series(lst).isin(vals))[0]

        X_new = np.array(X[idx_dup].A if issparse(X) else X[idx_dup])
        sorted_idx = np.lexsort(X_new.T)
        sorted_data = X_new[sorted_idx, :]

        row_mask = np.invert(np.append([True], np.any(np.diff(sorted_data, axis=0), 1)))
        idx = sorted_idx[row_mask]
        idx_dup = np.array(idx_dup)[idx]
    return idx_dup


# TODO: Add docstrings
def remove_duplicate_cells(adata):
    """TODO."""
    if "X_pca" not in adata.obsm.keys():
        pca(adata)
    idx_duplicates = get_duplicate_cells(adata)
    if len(idx_duplicates) > 0:
        mask = np.ones(adata.n_obs, bool)
        mask[idx_duplicates] = 0
        logg.info("Removed", len(idx_duplicates), "duplicate cells.")
        adata._inplace_subset_obs(mask)
        if "neighbors" in adata.uns.keys():
            neighbors(adata)
