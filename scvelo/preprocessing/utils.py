import numpy as np
from scipy.sparse import issparse


def show_proportions(adata):
    """Fraction of spliced/unspliced/ambiguous abundances

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.

    Returns
    -------
    Prints the fractions of abundances.
    """
    layers_keys = [key for key in ['spliced', 'unspliced', 'ambiguous'] if key in adata.layers.keys()]
    tot_mol_cell_layers = [adata.layers[key].sum(1) for key in layers_keys]

    mean_abundances = np.round(
        [np.mean(tot_mol_cell / np.sum(tot_mol_cell_layers, 0)) for tot_mol_cell in tot_mol_cell_layers], 2)

    print('Abundance of ' + str(layers_keys) + ': ' + str(mean_abundances))


def cleanup(data, clean='layers', keep={'spliced', 'unspliced'}, copy=False):
    """Deletes attributes not needed.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    cleanup: `str` or list of `str` (default: `layers`)
        Which attributes to consider for freeing memory.
    keep: `str` or list of `str` (default: `['spliced', unspliced']`)
        Which attributes to keep.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with selection of attributes kept.
    """
    adata = data.copy() if copy else data

    if any(['obs' in clean, 'all' in clean]):
        for key in list(adata.obs.keys()):
            if key not in keep: del adata.obs[key]

    if any(['var' in clean, 'all' in clean]):
        for key in list(adata.var.keys()):
            if key not in keep: del adata.var[key]

    if any(['uns' in clean, 'all' in clean]):
        for key in list(adata.uns.keys()):
            if key not in keep: del adata.uns[key]

    if any(['layers' in clean, 'all' in clean]):
        for key in list(adata.layers.keys()):  # remove layers that are not needed
            if key not in keep: del adata.layers[key]

    return adata if copy else None


def set_initial_size(adata, layers={'spliced', 'unspliced'}):
    if all([layer in adata.layers.keys() for layer in layers]):
        for layer in layers:
            X = adata.layers[layer]
            adata.obs['initial_size_' + layer] = X.sum(1).A1.copy() if issparse(X) else X.sum(1).copy()


def filter_genes(data, min_counts=3, min_counts_u=3, min_cells=None, min_cells_u=None, copy=False):
    """Filtering, normalization and log transform

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    Runs the following steps

    .. code:: python

        sc.pp.filter_genes(adata, min_counts=10)
        sc.pp.normalize_per_cell(adata)
        sc.pp.filter_genes_dispersion(adata, n_top_genes=10000)
        sc.pp.normalize_per_cell(adata)
        if log: sc.pp.log1p(adata)


    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    min_counts: `int` (default: 10)
        Minimum number of gene counts per cell.
    n_top_genes: `int` (default: 10000)
        Number of genes to keep.
    log: `bool` (default: `True`)
        Take logarithm.
    copy: `bool` (default: `False`)
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    adata = data.copy() if copy else data
    from scanpy.api.pp import filter_genes

    def filter_genes_u(adata, min_counts_u=None, min_cells_u=None):
        counts = adata.layers['unspliced'] if min_counts_u is not None else adata.layers['unspliced'] > 0
        counts = counts.sum(0).A1 if issparse(counts) else counts.sum(0)
        adata._inplace_subset_var(counts >= (min_counts_u if min_counts_u is not None else min_cells_u))

    # compute initial counts per cell first
    if all([key in adata.layers.keys() for key in {'spliced', 'unspliced'}]): set_initial_size(adata)

    if min_counts is not None: filter_genes(adata, min_counts=min_counts)
    if min_cells is not None: filter_genes(adata, min_cells=min_cells)

    if 'unspliced' in adata.layers.keys():
        if min_counts_u is not None: filter_genes_u(adata, min_counts_u=min_counts_u)
        if min_cells_u is not None: filter_genes_u(adata, min_cells_u=min_cells_u)

    return adata if copy else None


def filter_genes_dispersion(data, n_top_genes=None, flavor='seurat', copy=False):
    adata = data.copy() if copy else data
    if all([key in adata.layers.keys() for key in {'spliced', 'unspliced'}]): set_initial_size(adata)

    if n_top_genes is not None and n_top_genes < adata.shape[1]:
        #normalize_per_cell(adata)
        if flavor == 'svr':
            mu = adata.X.mean(0).A1 if issparse(adata.X) else adata.X.mean(0)
            sigma = np.sqrt(adata.X.multiply(adata.X).mean(0).A1 - mu**2) if issparse(adata.X) else adata.X.std(0)
            log_mu = np.log2(mu)
            log_cv = np.log2(sigma / mu)

            from sklearn.svm import SVR
            clf = SVR(gamma=150. / len(mu))
            clf.fit(log_mu[:, None], log_cv)
            score = log_cv - clf.predict(log_mu[:, None])
            nth_score = np.sort(score)[::-1][n_top_genes]
            adata._inplace_subset_var(score >= nth_score)
        else:
            from scanpy.api.pp import filter_genes_dispersion
            filter_genes_dispersion(adata, n_top_genes=n_top_genes, flavor=flavor)

    return adata if copy else None


def normalize_layers(data, layers={'spliced', 'unspliced'}, copy=False):
    """Normalize by total counts to median
    """
    adata = data.copy() if copy else data
    from scanpy.api.pp import normalize_per_cell

    def get_size_and_bool(adata, layer):
        if layer not in {'spliced', 'unspliced'}:
            return None, True
        else:
            size = adata.obs['initial_size_' + layer].copy() if 'initial_size_' + layer in adata.obs.keys() else None
            X = adata.layers[layer]
        return size, np.allclose((X.data[:10] if issparse(X) else X[0]) % 1, 0, atol=1e-3)

    for layer in layers:
        size, not_yet_normalized = get_size_and_bool(adata, layer)
        if not_yet_normalized:
            adata.layers[layer] = normalize_per_cell(adata.layers[layer], None, size, copy=True)
    return adata if copy else None


def normalize_per_cell(data, log=True, copy=False):
    adata = data.copy() if copy else data
    from scanpy.api.pp import normalize_per_cell, log1p
    size = adata.obs['initial_size_spliced'].copy() if 'initial_size_spliced' in adata.obs.keys() else None
    normalize_per_cell(adata, None, size)
    normalize_layers(adata)
    if log: log1p(adata)
    return adata if copy else None


def filter_and_normalize(data, min_counts=3, min_counts_u=3, min_cells=None, min_cells_u=None, n_top_genes=None,
                         log=True, flavor='seurat', copy=False):
    """Runs pp.filter_genes(), pp.filter_genes_dispersion() and pp.normalize_per_cell()
    """
    adata = data.copy() if copy else data
    filter_genes(adata, min_counts=min_counts, min_counts_u=min_counts_u, min_cells=min_cells, min_cells_u=min_cells_u)
    filter_genes_dispersion(adata, n_top_genes=n_top_genes, flavor=flavor)
    normalize_per_cell(adata, log=log)
    return adata if copy else None


def recipe_velocity(adata, min_counts=3, min_counts_u=3, n_top_genes=None, n_pcs=30, n_neighbors=30, log=True, copy=False):
    """Runs pp.filter_and_normalize() and pp.moments()
    """
    from .moments import moments
    filter_and_normalize(adata, min_counts=min_counts, min_counts_u=min_counts_u, n_top_genes=n_top_genes, log=log)
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata if copy else None
