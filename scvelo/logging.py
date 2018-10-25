"""Logging and Profiling
"""

from . import settings
from datetime import datetime
from time import time as get_time


_VERBOSITY_LEVELS_FROM_STRINGS = {'error': 0, 'warn': 1, 'info': 2, 'hint': 3}


def info(*args, **kwargs):
    return msg(*args, v='info', **kwargs)


def error(*args, **kwargs):
    args = ('Error:',) + args
    return msg(*args, v='error', **kwargs)


def warn(*args, **kwargs):
    args = ('WARNING:',) + args
    return msg(*args, v='warn', **kwargs)


def hint(*args, **kwargs):
    return msg(*args, v='hint', **kwargs)


def _settings_verbosity_greater_or_equal_than(v):
    if isinstance(settings.verbosity, str):
        settings_v = _VERBOSITY_LEVELS_FROM_STRINGS[settings.verbosity]
    else:
        settings_v = settings.verbosity
    return settings_v >= v


def msg(*msg, v=4, time=False, memory=False, reset=False, end='\n',
        no_indent=False, t=None, m=None, r=None):
    """Write message to logging output.
    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.
    v : {'error', 'warn', 'info', 'hint'} or int, (default: 4)
        0/'error', 1/'warn', 2/'info', 3/'hint', 4, 5, 6...
    time, t : bool, optional (default: False)
        Print timing information; restart the clock.
    memory, m : bool, optional (default: Faulse)
        Print memory information.
    reset, r : bool, optional (default: False)
        Reset timing and memory measurement. Is automatically reset
        when passing one of ``time`` or ``memory``.
    end : str (default: '\n')
        Same meaning as in builtin ``print()`` function.
    no_indent : bool (default: False)
        Do not indent for ``v >= 4``.
    """
    # variable shortcuts
    if t is not None: time = t
    if m is not None: memory = m
    if r is not None: reset = r
    if isinstance(v, str):
        v = _VERBOSITY_LEVELS_FROM_STRINGS[v]
    if v == 3:  # insert "--> " before hints
        msg = ('-->',) + msg
    if v >= 4 and not no_indent:
        msg = ('   ',) + msg
    if _settings_verbosity_greater_or_equal_than(v):
        if not time and not memory and len(msg) > 0:
            _write_log(*msg, end=end)
        if reset:
            try:
                settings._previous_memory_usage, _ = get_memory_usage()
            except:
                pass
            settings._previous_time = get_time()
        if time:
            elapsed = get_passed_time()
            msg = msg + ('({})'.format(_sec_to_str(elapsed)),)
            _write_log(*msg, end=end)
        if memory:
            _write_log(get_memory_usage(), end=end)

m = msg


def _write_log(*msg, end='\n'):
    """Write message to log output, ignoring the verbosity level.
    This is the most basic function.
    Parameters
    ----------
    *msg :
        One or more arguments to be formatted as string. Same behavior as print
        function.
    """
    from .settings import logfile
    if logfile == '':
        print(*msg, end=end)
    else:
        out = ''
        for s in msg:
            out += str(s) + ' '
        with open(logfile, 'a') as f:
            f.write(out + end)


def _sec_to_str(t):
    """Format time in seconds.
    Parameters
    ----------
    t : int
        Time in seconds.
    """
    from functools import reduce
    return "%d:%02d:%02d.%02d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*100,), 100, 60, 60])


def get_passed_time():
    now = get_time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def print_passed_time():
    now = get_time()
    elapsed = now - settings._previous_time
    settings._previous_time = now

    from functools import reduce
    elapsed = "%d:%02d:%02d.%02d" % reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(elapsed*100,), 100, 60, 60])
    return elapsed


def print_version():
    from . import __version__
    _write_log('Running scvelo', __version__, 'on {}.'.format(datetime.now().strftime("%Y-%m-%d %H:%M")))


def print_versions():
    for mod in ['scvelo', 'scanpy', 'anndata', 'loompy', 'numpy', 'scipy', 'matplotlib', 'sklearn', 'pandas']:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try: print('{}=={}'.format(mod_install, __import__(mod_name).__version__), end='  ')
        except (ImportError, AttributeError): pass
    print()


def get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


from anndata.logging import print_memory_usage
from anndata.logging import get_memory_usage


from sys import stdout
class ProgressReporter:
    def __init__(self, total, interval=3):
        self.count = 0
        self.total = total
        self.timestamp = get_time()
        self.interval = interval

    def update(self):
        self.count += 1
        if settings.verbosity > 1 and (get_time() - self.timestamp > self.interval or self.count == self.total):
            self.timestamp = get_time()
            percent = int(self.count * 100 / self.total)
            stdout.write('\r' + '... %d%%' % percent)
            stdout.flush()

    def finish(self):
        if settings.verbosity > 1:
            stdout.write('\r')
            stdout.flush()
