from scanpy import logging as logg
from scanpy import settings
settings.verbosity = 3   # global verbosity level: show errors(0), warnings(1), info(2) and hints(3)
settings._frameon = False

from datetime import datetime
from time import time


def get_passed_time():
    now = time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def print_version():
    from . import __version__
    logg._write_log('Running scvelo', __version__, 'on {}.'.format(datetime.now().strftime("%Y-%m-%d %H:%M")))


def print_versions():
    for mod in ['scvelo', 'scanpy', 'anndata', 'loompy', 'numpy', 'scipy', 'matplotlib', 'sklearn', 'pandas']:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try: print('{}=={}'.format(mod_install, __import__(mod_name).__version__), end='  ')
        except (ImportError, AttributeError): pass
    print()


from anndata.logging import print_memory_usage
print_version_and_date = print_version
