import numpy as np
from scipy.sparse import csr_matrix, issparse

from scvelo import logging as logg
from scvelo import settings
from .neighbors import get_connectivities, get_n_neighs, neighbors, verify_neighbors
from .utils import normalize_per_cell, not_yet_normalized


def moments(
    data,
    n_neighbors=30,
    n_pcs=None,
    mode="connectivities",
    method="umap",
    use_rep=None,
    use_highly_variable=True,
    copy=False,
):
    """Computes moments for velocity estimation.

    First-/second-order moments are computed for each cell across its nearest neighbors,
    where the neighbor graph is obtained from euclidean distances in PCA space.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    n_neighbors: `int` (default: 30)
        Number of neighbors to use.
    n_pcs: `int` (default: None)
        Number of principal components to use.
        If not specified, the full space is used of a pre-computed PCA,
        or 30 components are used when PCA is computed internally.
    mode: `'connectivities'` or `'distances'`  (default: `'connectivities'`)
        Distance metric to use for moment computation.
    method : {{'umap', 'hnsw', 'sklearn', `None`}}  (default: `'umap'`)
        Method to compute neighbors, only differs in runtime.
        Connectivities are computed with adaptive kernel width as proposed in
        Haghverdi et al. 2016 (https://doi.org/10.1038/nmeth.3971).
    use_rep : `None`, `'X'` or any key for `.obsm` (default: None)
        Use the indicated representation. If `None`, the representation is chosen
        automatically: for .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.
    use_highly_variable: `bool` (default: True)
        Whether to use highly variable genes only, stored in .var['highly_variable'].
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Ms: `.layers`
        dense matrix with first order moments of spliced counts.
    Mu: `.layers`
        dense matrix with first order moments of unspliced counts.
    """

    adata = data.copy() if copy else data

    layers = [layer for layer in {"spliced", "unspliced"} if layer in adata.layers]
    if any([not_yet_normalized(adata.layers[layer]) for layer in layers]):
        normalize_per_cell(adata)

    if n_neighbors is not None and n_neighbors > get_n_neighs(adata):
        neighbors(
            adata,
            n_neighbors=n_neighbors,
            use_rep=use_rep,
            use_highly_variable=use_highly_variable,
            n_pcs=n_pcs,
            method=method,
        )
    verify_neighbors(adata)

    if "spliced" not in adata.layers.keys() or "unspliced" not in adata.layers.keys():
        logg.warn("Skipping moments, because un/spliced counts were not found.")
    else:
        logg.info(f"computing moments based on {mode}", r=True)
        connectivities = get_connectivities(
            adata, mode, n_neighbors=n_neighbors, recurse_neighbors=False
        )

        adata.layers["Ms"] = (
            csr_matrix.dot(connectivities, csr_matrix(adata.layers["spliced"]))
            .astype(np.float32)
            .A
        )
        adata.layers["Mu"] = (
            csr_matrix.dot(connectivities, csr_matrix(adata.layers["unspliced"]))
            .astype(np.float32)
            .A
        )
        # if renormalize: normalize_per_cell(adata, layers={'Ms', 'Mu'}, enforce=True)

        logg.info(
            "    finished", time=True, end=" " if settings.verbosity > 2 else "\n"
        )
        logg.hint(
            "added \n"
            "    'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)"
        )
    return adata if copy else None


def second_order_moments(adata, adjusted=False):
    """Computes second order moments for stochastic velocity estimation.

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    Returns
    -------
    Mss: Second order moments for spliced abundances
    Mus: Second order moments for spliced with unspliced abundances
    """

    if "neighbors" not in adata.uns:
        raise ValueError(
            "You need to run `pp.neighbors` first to compute a neighborhood graph."
        )

    connectivities = get_connectivities(adata)
    s, u = csr_matrix(adata.layers["spliced"]), csr_matrix(adata.layers["unspliced"])
    if s.shape[0] == 1:
        s, u = s.T, u.T
    Mss = csr_matrix.dot(connectivities, s.multiply(s)).astype(np.float32).A
    Mus = csr_matrix.dot(connectivities, s.multiply(u)).astype(np.float32).A
    if adjusted:
        Mss = 2 * Mss - adata.layers["Ms"].reshape(Mss.shape)
        Mus = 2 * Mus - adata.layers["Mu"].reshape(Mus.shape)
    return Mss, Mus


def second_order_moments_u(adata):
    """Computes second order moments for stochastic velocity estimation.

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    Returns
    -------
    Muu: Second order moments for unspliced abundances
    """

    if "neighbors" not in adata.uns:
        raise ValueError(
            "You need to run `pp.neighbors` first to compute a neighborhood graph."
        )

    connectivities = get_connectivities(adata)
    u = csr_matrix(adata.layers["unspliced"])
    Muu = csr_matrix.dot(connectivities, u.multiply(u)).astype(np.float32).A

    return Muu


def magic_impute(adata, knn=5, t=2, verbose=0, **kwargs):
    logg.info(
        "To be used carefully. Magic has not yet been tested for this application."
    )
    import magic

    magic_operator = magic.MAGIC(verbose=verbose, knn=knn, t=t, **kwargs)
    adata.layers["Ms"] = magic_operator.fit_transform(adata.layers["spliced"])
    adata.layers["Mu"] = magic_operator.transform(adata.layers["unspliced"])


def get_moments(
    adata, layer=None, second_order=None, centered=True, mode="connectivities"
):
    """Computes moments for a specified layer.

    First and second order moments.
    If centered, that corresponds to means and variances across nearest neighbors.

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.
    layer: `str` (default: `None`)
        Key of layer with abundances to consider for moment computation.
    second_order: `bool` (default: `None`)
        Whether to compute second order moments from abundances.
    centered: `bool` (default: `True`)
        Whether to compute centered (=variance) or uncentered second order moments.
    mode: `'connectivities'` or `'distances'`  (default: `'connectivities'`)
        Distance metric to use for moment computation.

    Returns
    -------
    Mx: first or second order moments
    """

    if "neighbors" not in adata.uns:
        raise ValueError(
            "You need to run `pp.neighbors` first to compute a neighborhood graph."
        )
    connectivities = get_connectivities(adata, mode=mode)
    X = (
        adata.X
        if layer is None
        else adata.layers[layer]
        if isinstance(layer, str)
        else layer
    )
    X = (
        csr_matrix(X)
        if isinstance(layer, str) and layer in {"spliced", "unspliced"}
        else np.array(X)
        if not issparse(X)
        else X
    )
    if not issparse(X):
        X = X[:, ~np.isnan(X.sum(0))]
    if second_order:
        X2 = X.multiply(X) if issparse(X) else X ** 2
        Mx = (
            csr_matrix.dot(connectivities, X2)
            if second_order
            else csr_matrix.dot(connectivities, X)
        )
        if centered:
            mu = csr_matrix.dot(connectivities, X)
            mu2 = mu.multiply(mu) if issparse(mu) else mu ** 2
            Mx = Mx - mu2
    else:
        Mx = csr_matrix.dot(connectivities, X)
    if issparse(X):
        Mx = Mx.astype(np.float32).A
    return Mx
