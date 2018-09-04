"""scvelo - stochastic single cell RNA velocity"""

from get_version import get_version
__version__ = get_version(__file__)
del get_version

from . import preprocessing as pp
from . import tools as tl
from . import plotting as pl
from . import datasets
from . import logging
