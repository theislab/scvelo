import warnings
from functools import wraps
from typing import Mapping


# Modified from https://github.com/scverse/scanpy/blob/master/scanpy/_utils/__init__.py
def deprecated_arg_names(arg_mapping: Mapping[str, str]):
    """Marks a functions keyword arguments as deprecated.

    This decorator marks a functions keyword arguments as deprecated. It will
    result in a warning being emitted when the deprecated keyword argument is
    used, and the function being called with the new argument.

    Parameters
    ----------
    arg_mapping
        Mapping from deprecated argument name to current argument name.
    """

    def decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            warnings.simplefilter("always", DeprecationWarning)  # turn off filter
            for old, new in arg_mapping.items():
                if old in kwargs:
                    warnings.warn(
                        f"Keyword argument '{old}' has been "
                        f"deprecated in favour of '{new}'. "
                        f"'{old}' will be removed in a future version.",
                        category=DeprecationWarning,
                        stacklevel=2,
                    )
                    val = kwargs.pop(old)
                    if old == "copy":
                        kwargs[new] = not val
                    else:
                        kwargs[new] = val
            # reset filter
            warnings.simplefilter("default", DeprecationWarning)
            return func(*args, **kwargs)

        return func_wrapper

    return decorator
