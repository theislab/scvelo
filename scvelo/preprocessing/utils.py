import numpy as np


def show_proportions(adata, copy=False):
    """Fraction of spliced/unspliced/ambiguous abundances

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Prints the fractions of abundances.
    """
    layers_keys = [key for key in ['spliced', 'unspliced', 'ambiguous'] if key in adata.layers.keys()]
    tot_mol_cell_layers = [adata.layers[key].sum(1) for key in layers_keys]

    mean_abundances = np.round(
        [np.mean(tot_mol_cell / np.sum(tot_mol_cell_layers, 0)) for tot_mol_cell in tot_mol_cell_layers], 2)

    print('Abundance of ' + str(layers_keys) + ': ' + str(mean_abundances))

    return adata if copy else None


def cleanup(adata, clean='layers', keep={'spliced', 'unspliced'}, copy=False):
    """Deletes attributes not needed.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
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


def filter_and_normalize(adata, min_counts=10, n_top_genes=None, log=True, copy=False):
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
    adata: :class:`~anndata.AnnData`
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
    from scanpy.api.pp import filter_genes, filter_genes_dispersion, normalize_per_cell, log1p
    filter_genes(adata, min_counts=min_counts)
    if n_top_genes is not None and n_top_genes < adata.shape[1]:
        normalize_per_cell(adata)
        filter_genes_dispersion(adata, n_top_genes=n_top_genes)
    normalize_per_cell(adata)
    if log: log1p(adata)
    return adata if copy else None


def recipe_velocity(adata, min_counts=10, n_top_genes=None, n_pcs=30, n_neighbors=30, log=True, copy=False):
    """Runs pp.filter_and_normalize() and pp.moments()
    """
    from .moments import moments
    filter_and_normalize(adata, min_counts=min_counts, n_top_genes=n_top_genes, log=log)
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata if copy else None
