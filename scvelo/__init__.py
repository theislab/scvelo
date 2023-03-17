"""scvelo - RNA velocity generalized through dynamical modeling."""
from anndata import AnnData
from scanpy import read, read_loom

from scvelo import datasets, logging
from scvelo import plotting as pl
from scvelo import preprocessing as pp
from scvelo import settings
from scvelo import tools as tl
from scvelo import utils
from scvelo.core import get_df
from scvelo.plotting.gridspec import GridSpec
from scvelo.preprocessing.neighbors import Neighbors
from scvelo.read_load import load, read_csv
from scvelo.settings import set_figure_params
from scvelo.tools.run import run_all, test
from scvelo.tools.utils import round
from scvelo.tools.velocity import Velocity
from scvelo.tools.velocity_graph import VelocityGraph

import sys  # isort:skip

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["tl", "pp", "pl"]})

try:
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    from importlib.metadata import version

    __version__ = version(__name__)
    del version


__all__ = [
    "AnnData",
    "datasets",
    "get_df",
    "GridSpec",
    "load",
    "logging",
    "Neighbors",
    "pl",
    "pp",
    "read",
    "read_csv",
    "read_loom",
    "round",
    "run_all",
    "set_figure_params",
    "settings",
    "test",
    "tl",
    "utils",
    "Velocity",
    "VelocityGraph",
]
