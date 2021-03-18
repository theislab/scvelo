from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import ndarray

from scvelo.core import l2_norm


# TODO: Extend test to generate sparse inputs as well
@given(
    a=arrays(
        float,
        shape=st.integers(min_value=1, max_value=100),
        elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
    ),
    axis=st.integers(min_value=0, max_value=1),
)
def test_l2_norm(a: ndarray, axis: int):
    if a.ndim == 1:
        np.allclose(np.linalg.norm(a), l2_norm(a, axis=axis))
    else:
        np.allclose(np.linalg.norm(a, axis=axis), l2_norm(a, axis=axis))
