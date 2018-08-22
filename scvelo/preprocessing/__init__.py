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
        X = lc.layer['spliced'][:, :].T

        layers = OrderedDict()
        layers['spliced'] = lc.layer["spliced"][:, :].T
        layers['unspliced'] = lc.layer["unspliced"][:, :].T

        obs = dict(lc.col_attrs)
        obs['obs_names'] = obs.pop('CellID')

        var = dict(lc.row_attrs)
        var['var_names'] = var.pop('Gene')

        adata = AnnData(X, obs=obs, var=var, layers=layers)

    return adata