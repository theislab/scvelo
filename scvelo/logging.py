from scanpy import logging as logg
from scanpy import settings
settings.verbosity = 3   # global verbosity level: show errors(0), warnings(1), info(2) and hints(3)
settings._frameon = False

import sys
from datetime import datetime
from time import time


def get_passed_time():
    now = time()
    elapsed = now - settings._previous_time
    settings._previous_time = now

    from functools import reduce
    elapsed = "%d:%02d:%02d.%02d" % reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(elapsed*100,), 100, 60, 60])
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


class ProgressReporter:
    def __init__(self, total, interval=3):
        self.count = 0
        self.total = total
        self.timestamp = time()
        self.interval = interval

    def update(self):
        self.count += 1
        if time() - self.timestamp > self.interval or self.count == self.total:
            self.timestamp = time()
            percent = int(self.count * 100 / self.total)
            sys.stdout.write('\r' + '... %d%%' % percent)
            sys.stdout.flush()

    def finish(self):
        sys.stdout.write('\r')
        sys.stdout.flush()


from anndata.logging import print_memory_usage
print_version_and_date = print_version
