from scanpy import settings
from scanpy import logging as logg
from scanpy.logging import print_memory_usage

from datetime import datetime
from time import time


def get_passed_time():
    now = time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def get_date_string():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


def print_version_and_date():
    from . import __version__
    logg._write_log('Running scvelo', __version__, 'on {}.'.format(get_date_string()))


def print_versions():
    print_version_and_date()
    print('Dependencies:', end='  ')
    for mod in ['scanpy', 'anndata', 'numpy', 'pandas']:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try: print('{}=={}'.format(mod_install, __import__(mod_name).__version__), end='  ')
        except (ImportError, AttributeError): pass
    print()
