"""scvelo - stochastic single cell RNA velocity"""

from .get_version import get_version
__version__ = get_version(__file__)
del get_version

from .read_load import read, load
from .logging import settings
settings.verbosity = 3   # global verbosity level: show errors(0), warnings(1), info(2) and hints(3)

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
from . import datasets
from . import logging
