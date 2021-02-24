"""scvelo - RNA velocity generalized through dynamical modeling"""

try:
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    try:
        from importlib_metadata import version  # Python < 3.8
    except Exception:
        from importlib.metadata import version  # Python = 3.8
    __version__ = version(__name__)
    del version

from scanpy import read, read_loom
from .read_load import AnnData, load, read_csv, get_df, DataFrame
from .preprocessing.neighbors import Neighbors
from .tools.run import run_all, test
from .tools.utils import round
from .tools.velocity import Velocity
from .tools.velocity_graph import VelocityGraph
from .plotting.gridspec import GridSpec
from .settings import set_figure_params

from . import pp
from . import tl
from . import pl
from . import utils
from . import datasets
from . import logging
from . import settings


__all__ = [
    "AnnData",
    "DataFrame",
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
