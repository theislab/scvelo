import numpy as np


def print_fraction_of_abundances(adata, cleanup_layers=True, copy=False):
    """Computes the fraction of abundances of spliced, unspliced and ambiguous counts

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix

    cleanup_layers: `bool` (default: `True`)
        Delete all layers except 'spliced' and 'unspliced'

    copy: `bool`, (default: `False`)
        Return a copy instead of writing to adata

    Returns
    -------
    Prints the fractions of abundances.
    Returns or updates `adata` with all layers freed except `spliced` and `unspliced` (`cleanup_layers == True`).
    """
    layers_keys = [key for key in ['spliced', 'unspliced', 'ambiguous'] if key in adata.layers.keys()]
    tot_mol_cell_layers = [adata.layers[key].sum(1) for key in layers_keys]

    mean_abundances = np.round(
        [np.mean(tot_mol_cell / np.sum(tot_mol_cell_layers, 0)) for tot_mol_cell in tot_mol_cell_layers], 2)

    print('abundance of ' + str(layers_keys) + ': ' + str(mean_abundances))

    if cleanup_layers:
        for key in list(adata.layers.keys()):  # remove layers that are not needed
            if key not in ['spliced', 'unspliced']: del adata.layers[key]

    return adata if copy else None
