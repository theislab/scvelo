from scanpy import settings
from scanpy import logging as logg

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
