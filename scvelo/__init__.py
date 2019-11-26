"""scvelo - RNA velocity using dynamical modeling"""

try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    try: from importlib_metadata import version  # Python < 3.8
    except: from importlib.metadata import version  # Python = 3.8
    __version__ = version(__name__)
    del version

from .read_load import AnnData, read, read_loom, load, read_csv
from .preprocessing.neighbors import Neighbors
from .tools.run import run_all, test
from .tools.velocity import Velocity
from .tools.velocity_graph import VelocityGraph
from .settings import set_figure_params

from . import pp
from . import tl
from . import pl
from . import datasets
from . import utils
from . import logging
from . import settings
