from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix

from scvelo.core import l2_norm


# TODO: Extend test to generate sparse inputs as well
class TestL2Norm:
    @given(
        a=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_1d_array(self, a: ndarray, axis: int):
        np.allclose(np.linalg.norm(a), l2_norm(a, axis=axis))

    @given(
        a=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_2d_array(self, a: ndarray, axis: int):
        np.allclose(np.linalg.norm(a, axis=axis), l2_norm(a, axis=axis))

    @given(
        a=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_sparse_input(self, a: ndarray, axis: int):
        np.allclose(np.linalg.norm(a, axis=axis), l2_norm(csr_matrix(a), axis=axis))
