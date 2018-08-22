from .moments import moments
from .fractions import print_fraction_of_abundances


def read_loom_layers(file_name, backup_url=None):
    """Read loom file and return AnnData object

    Arguments
    ---------
    loom_filepath: str
        directory of loom file

    backup_url: str

    all_loom_layers: bool, default=False
        whether to load loom layers 'loom_obs' and 'loom_var' into AnnData object

    Returns
    -------
    Updates the attributes .X, .var, .obsm['Xs'], .obsm['Xs']
    """
    import os
    import loompy
    from scanpy.api import AnnData
    from collections import OrderedDict

    if not os.path.exists(file_name):
        try:
            from urllib.request import urlretrieve
            urlretrieve(backup_url, file_name)
        except OSError:
            print("OS error: {0}".format(OSError))

    with loompy.connect(file_name, 'r') as lc:
        X = lc.layer['spliced'].sparse().T.tocsr()

        layers = OrderedDict()
        layers['spliced'] = lc.layer["spliced"].sparse().T.tocsr()
        layers['unspliced'] = lc.layer["unspliced"].sparse().T.tocsr()

        obs = dict(lc.col_attrs)
        obs['obs_names'] = obs.pop('CellID')

        var = dict(lc.row_attrs)
        var['var_names'] = var.pop('Gene')

        adata = AnnData(X, obs=obs, var=var, layers=layers)

    return adata


def recipe_velocity(adata, min_counts=10, n_top_genes=3000, n_pcs=30, n_neighbors=15, copy=False):
    from scanpy.api.pp import \
        filter_genes, filter_genes_dispersion, normalize_per_cell, pca, neighbors
    filter_genes(adata, min_counts=min_counts)
    if n_top_genes < adata.n_vars: filter_genes_dispersion(adata, n_top_genes=n_top_genes)
    normalize_per_cell(adata, layers='all')
    pca(adata, n_comps=n_pcs)
    neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')

    moments(adata)
    return adata if copy else None
