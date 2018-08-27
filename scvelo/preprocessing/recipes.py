from .moments import moments


def recipe_velocity(adata, min_counts=10, n_top_genes=None, n_pcs=30, n_neighbors=30, log=True, copy=False):
    """Filtering, normalization, neighbors graph and moments.

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    The recipe runs the following steps

    .. code:: python

        sc.pp.filter_genes(adata, min_counts=10)
        sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')
        sc.pp.filter_genes_dispersion(adata, n_top_genes=3000)
        sc.pp.normalize_per_cell(adata, layers='all)
        if log: sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=30)
        sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca')
        scv.pp.moments(adata)


    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    n_top_genes: `int` (default: 1000)
        Number of genes to keep.
    log: `bool` (default: `True`)
        Take logarithm.
    n_pcs: `int` (default: 30)
        Number of principal components to use.
    n_neighbors: `int` (default: 30)
        Number of neighbors to use.
    copy: `bool` (default: `False`)
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    from scanpy.api.pp import filter_genes, filter_genes_dispersion, normalize_per_cell, log1p, pca, neighbors

    filter_genes(adata, min_counts=min_counts)
    normalize_per_cell(adata, key_n_counts='n_counts_all')

    filter_genes_dispersion(adata, n_top_genes=n_top_genes)
    normalize_per_cell(adata, layers='all')
    if log: log1p(adata)
    pca(adata, n_comps=n_pcs)
    neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')

    moments(adata)
    return adata if copy else None
