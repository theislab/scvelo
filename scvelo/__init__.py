"""scvelo - stochastic single cell RNA velocity"""

from .get_version import get_version
__version__ = get_version(__file__)
del get_version

from .read_load import read, load, clean_obs_names, merge
from .logging import settings

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
from . import datasets
from . import logging


def run(data, basis=None, mode='deterministic', min_counts=10, n_pcs=30, n_neighbors=30, copy=False):
    from time import time
    start = time()

    adata = data.copy() if copy else data
    pp.filter_and_normalize(adata, min_counts=min_counts)
    pp.moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    tl.velocity(adata, mode=mode)
    tl.velocity_graph(adata)
    if basis is not None: tl.velocity_embedding(adata, basis=basis)

    print('/n Total time (seconds): ' + str(round(time() - start, 2)))

    return adata if copy else None
