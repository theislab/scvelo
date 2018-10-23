"""scvelo - stochastic single cell RNA velocity"""

from .get_version import get_version
__version__ = get_version(__file__)
del get_version

from .read_load import AnnData, read, read_loom, load
from .logging import settings
from .tools.run import run_all
from .tools.velocity import Velocity
from .tools.velocity_graph import VelocityGraph

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
from . import datasets
from . import logging
from . import utils