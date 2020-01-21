"""Logging and Profiling
"""

from . import settings
from sys import stdout
from datetime import datetime
from time import time as get_time
from platform import python_version
from anndata.logging import get_memory_usage
from anndata.logging import print_memory_usage


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


def _sec_to_str(t, show_microseconds=False):
    """Format time in seconds.
    Parameters
    ----------
    t : int
        Time in seconds.
    """
    from functools import reduce
    t_str = "%d:%02d:%02d.%02d" % reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t*100,), 100, 60, 60])
    return t_str if show_microseconds else t_str[:-3]


def get_passed_time():
    now = get_time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


def print_passed_time():
    return _sec_to_str(get_passed_time())


def timeout(func, args=(), kwargs={}, timeout_duration=2, default=None):
    """This will spwan a thread and run the given function using the args, kwargs and
    return the given default value if the timeout_duration is exceeded
    """
    import threading

    class InterruptableThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.result = default

        def run(self):
            try: self.result = func(*args, **kwargs)
            except: pass

    it = InterruptableThread()
    it.start()
    it.join(timeout_duration)
    return it.result


def get_latest_pypi_version():
    from subprocess import check_output, CalledProcessError
    try:  # needs to work offline as well
        result = check_output(["pip", "search", "scvelo"])
        return str(result.split()[-1])[2:-1]
    except CalledProcessError:
        return '0.0.0'


def check_if_latest_version():
    from . import __version__
    latest_version = timeout(get_latest_pypi_version, timeout_duration=2, default='0.0.0')
    if __version__.rsplit('.dev')[0] < latest_version.rsplit('.dev')[0]:
        warn('There is a newer scvelo version available on PyPI:\n',
             'Your version: \t\t', __version__, '\n',
             'Latest version: \t', latest_version)


def print_version():
    from . import __version__
    _write_log('Running scvelo', __version__, '(python ' + python_version() + ')', 'on {}.'.format(get_date_string()))
    check_if_latest_version()


def print_versions():
    for mod in ['scvelo', 'scanpy', 'anndata', 'loompy', 'numpy', 'scipy', 'matplotlib', 'sklearn', 'pandas']:
        mod_name = mod[0] if isinstance(mod, tuple) else mod
        mod_install = mod[1] if isinstance(mod, tuple) else mod
        try: print('{}=={}'.format(mod_install, __import__(mod_name).__version__), end='  ')
        except (ImportError, AttributeError): pass
    check_if_latest_version()


def get_date_string():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


def switch_verbosity(mode='on', module=None):
    if module is None: from . import settings
    elif module is 'scanpy': from scanpy import settings
    else: exec('from ' + module + ' import settings')

    if mode is 'on' and hasattr(settings, 'tmp_verbosity'):
        settings.verbosity = settings.tmp_verbosity
        del settings.tmp_verbosity

    elif mode is 'off':
        settings.tmp_verbosity = settings.verbosity
        settings.verbosity = 0

    elif not isinstance(mode, str):
        settings.tmp_verbosity = settings.verbosity
        settings.verbosity = mode


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



