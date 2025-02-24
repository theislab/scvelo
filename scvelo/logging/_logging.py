"""Logging and Profiling."""
from collections.abc import Iterable
from datetime import datetime
from platform import python_version
from sys import stdout
from time import time as get_time

from anndata.logging import get_memory_usage

from scvelo import settings

_VERBOSITY_LEVELS_FROM_STRINGS = {"error": 0, "warn": 1, "info": 2, "hint": 3}


# TODO: Add docstrings
def info(*args, **kwargs):
    """TODO."""
    return msg(*args, v="info", **kwargs)


# TODO: Add docstrings
def error(*args, **kwargs):
    """TODO."""
    args = ("Error:",) + args
    return msg(*args, v="error", **kwargs)


# TODO: Add docstrings
def warn(*args, **kwargs):
    """TODO."""
    args = ("WARNING:",) + args
    return msg(*args, v="warn", **kwargs)


# TODO: Add docstrings
def hint(*args, **kwargs):
    """TODO."""
    return msg(*args, v="hint", **kwargs)


# TODO: Add docstrings
def _settings_verbosity_greater_or_equal_than(v):
    """TODO."""
    if isinstance(settings.verbosity, str):
        settings_v = _VERBOSITY_LEVELS_FROM_STRINGS[settings.verbosity]
    else:
        settings_v = settings.verbosity
    return settings_v >= v


def msg(
    *msg,
    v=None,
    time=False,
    memory=False,
    reset=False,
    end="\n",
    no_indent=False,
    t=None,
    m=None,
    r=None,
):
    r"""Write message to logging output.

    Log output defaults to standard output but can be set to a file
    by setting `sc.settings.log_file = 'mylogfile.txt'`.

    Parameters
    ----------
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
    if t is not None:
        time = t
    if m is not None:
        memory = m
    if r is not None:
        reset = r
    if v is None:
        v = 4
    if isinstance(v, str):
        v = _VERBOSITY_LEVELS_FROM_STRINGS[v]
    if v == 3:  # insert "--> " before hints
        msg = ("-->",) + msg
    if v >= 4 and not no_indent:
        msg = ("   ",) + msg
    if _settings_verbosity_greater_or_equal_than(v):
        if not time and not memory and len(msg) > 0:
            _write_log(*msg, end=end)
        if reset:
            try:
                settings._previous_memory_usage, _ = get_memory_usage()
            except ImportError as e:
                ImportError(e)
            settings._previous_time = get_time()
        if time:
            elapsed = get_passed_time()
            msg = msg + (f"({_sec_to_str(elapsed)})",)
            _write_log(*msg, end=end)
        if memory:
            _write_log(get_memory_usage(), end=end)


# TODO: Add docstrings
def _write_log(*msg, end="\n"):
    """Write message to log output, ignoring the verbosity level.

    This is the most basic function.

    Parameters
    ----------
    msg
        One or more arguments to be formatted as string. Same behavior as print function.
    """
    from scvelo.settings import logfile

    if logfile == "":
        print(*msg, end=end)
    else:
        out = ""
        for s in msg:
            out += f"{s} "
        with open(logfile, "a") as f:
            f.write(out + end)


# TODO: Add docstrings
def _sec_to_str(t, show_microseconds=False):
    """Format time in seconds.

    Parameters
    ----------
    t : int
        Time in seconds.
    """
    from functools import reduce

    t_str = "%d:%02d:%02d.%02d" % reduce(
        lambda ll, b: divmod(ll[0], b) + ll[1:], [(t * 100,), 100, 60, 60]
    )
    return t_str if show_microseconds else t_str[:-3]


# TODO: Add docstrings
def get_passed_time():
    """TODO."""
    now = get_time()
    elapsed = now - settings._previous_time
    settings._previous_time = now
    return elapsed


# TODO: Finish docstrings
def timeout(func, args=(), timeout_duration=2, default=None, **kwargs):
    """Spwans thread and runs the given function using the args, kwargs, and return default value on timeout."""
    import threading

    class InterruptableThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.result = default

        def run(self):
            self.result = func(*args, **kwargs)

    it = InterruptableThread()
    it.start()
    it.join(timeout_duration)
    return it.result


# TODO: Add docstrings
def print_version():
    """TODO."""
    from scvelo import __version__

    _write_log(
        f"Running scvelo {__version__} "
        f"(Python {python_version()}) on {datetime.now().strftime('%Y-%m-%d %H:%M')}.",
    )


_DEPENDENCIES_NUMERICS = [
    "anndata",
    "loompy",
    "numba",
    "numpy",
    "pandas",
    "scanpy",
    "scipy",
    ("sklearn", "scikit-learn"),
]


_DEPENDENCIES_PLOTTING = ["matplotlib"]


# Adapted from https://github.com/theislab/cellrank/blob/main/src/cellrank/logging/_logging.py#L98-L128
def _versions_dependencies(dependencies: Iterable[str]):
    # this is not the same as the requirements!
    for mod in dependencies:
        mod_name, dist_name = mod if isinstance(mod, tuple) else (mod, mod)
        try:
            imp = __import__(mod_name)
            yield dist_name, imp.__version__
        except (ImportError, AttributeError):
            pass


def print_versions():
    """Print versions of relevant packages."""
    from scvelo import settings

    modules = ["scvelo"] + _DEPENDENCIES_NUMERICS + _DEPENDENCIES_PLOTTING
    print(
        " ".join(
            f"{module}=={version}"
            for module, version in _versions_dependencies(modules)
        ),
        file=settings.logfile,
    )


# TODO: Add docstrings
def switch_verbosity(mode="on", module=None):
    """TODO."""
    if module is None:
        from . import settings
    elif module == "scanpy":
        from scanpy import settings
    else:
        exec(f"from {module} import settings")

    if mode == "on" and hasattr(settings, "tmp_verbosity"):
        settings.verbosity = settings.tmp_verbosity
        del settings.tmp_verbosity

    elif mode == "off":
        settings.tmp_verbosity = settings.verbosity
        settings.verbosity = 0

    elif not isinstance(mode, str):
        settings.tmp_verbosity = settings.verbosity
        settings.verbosity = mode


# TODO: Add docstrings
class ProgressReporter:
    """TODO."""

    def __init__(self, total, interval=3):
        self.count = 0
        self.total = total
        self.timestamp = get_time()
        self.interval = interval

    # TODO: Add docstrings
    def update(self):
        """TODO."""
        self.count += 1
        if settings.verbosity > 1 and (
            get_time() - self.timestamp > self.interval or self.count == self.total
        ):
            self.timestamp = get_time()
            percent = int(self.count * 100 / self.total)
            stdout.write(f"\r... {percent}%")
            stdout.flush()

    # TODO: Add docstrings
    def finish(self):
        """TODO."""
        if settings.verbosity > 1:
            stdout.write("\r")
            stdout.flush()
