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

    print('abundance of ' + str(layers_keys) + ': ' + str(mean_abundances))

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

    if any(['obsm' in clean, 'all' in clean]):
        for key in list(adata.obsm.keys()):
            if key not in keep: del adata.obsm[key]

    if any(['varm' in clean, 'all' in clean]):
        for key in list(adata.varm.keys()):
            if key not in keep: del adata.varm[key]

    if any(['uns' in clean, 'all' in clean]):
        for key in list(adata.uns.keys()):
            if key not in keep: del adata.uns[key]

    if any(['layers' in clean, 'all' in clean]):
        for key in list(adata.layers.keys()):  # remove layers that are not needed
            if key not in keep: del adata.layers[key]

    return adata if copy else None
