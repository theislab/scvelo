from typing import List

from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import ndarray
from numpy.testing import assert_almost_equal, assert_array_equal

from scvelo.core import clipped_log, invert, prod_sum, sum


class TestClippedLog:
    @given(
        a=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        bounds=st.lists(
            st.floats(
                min_value=0, max_value=100, allow_infinity=False, allow_nan=False
            ),
            min_size=2,
            max_size=2,
            unique=True,
        ),
        eps=st.floats(
            min_value=1e-6, max_value=1, allow_infinity=False, allow_nan=False
        ),
    )
    def test_flat_arrays(self, a: ndarray, bounds: List[float], eps: float):
        lb = min(bounds)
        ub = max(bounds) + 2 * eps

        a_logged = clipped_log(a, lb=lb, ub=ub, eps=eps)

        assert a_logged.shape == a.shape
        if (a <= lb).any():
            assert_almost_equal(np.abs(a_logged - np.log(lb + eps)).min(), 0)
        else:
            assert (a_logged >= np.log(lb + eps)).all()
        if (a >= ub).any():
            assert_almost_equal(np.abs(a_logged - np.log(ub - eps)).min(), 0)
        else:
            assert (a_logged <= np.log(ub - eps)).all()

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
        bounds=st.lists(
            st.floats(
                min_value=0, max_value=100, allow_infinity=False, allow_nan=False
            ),
            min_size=2,
            max_size=2,
            unique=True,
        ),
        eps=st.floats(
            min_value=1e-6, max_value=1, allow_infinity=False, allow_nan=False
        ),
    )
    def test_2d_arrays(self, a: ndarray, bounds: List[float], eps: float):
        lb = min(bounds)
        ub = max(bounds) + 2 * eps

        a_logged = clipped_log(a, lb=lb, ub=ub, eps=eps)

        assert a_logged.shape == a.shape
        if (a <= lb).any():
            assert_almost_equal(np.abs(a_logged - np.log(lb + eps)).min(), 0)
        else:
            assert (a_logged >= np.log(lb + eps)).all()
        if (a >= ub).any():
            assert_almost_equal(np.abs(a_logged - np.log(ub - eps)).min(), 0)
        else:
            assert (a_logged <= np.log(ub - eps)).all()


class TestInvert:
    @given(
        a=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        )
    )
    def test_flat_arrays(self, a: ndarray):
        a_inv = invert(a)

        if a[a != 0].size == 0:
            assert a_inv[a != 0].size == 0
        else:
            assert_array_equal(a_inv[a != 0], 1 / a[a != 0])

        if 0 in a:
            assert np.isnan(a_inv[a == 0]).all()
        else:
            assert set(a_inv[a == 0]) == set()

    @given(
        a=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        )
    )
    def test_2d_arrays(self, a: ndarray):
        a_inv = invert(a)

        if a[a != 0].size == 0:
            assert a_inv[a != 0].size == 0
        else:
            assert_array_equal(a_inv[a != 0], 1 / a[a != 0])

        if 0 in a:
            assert np.isnan(a_inv[a == 0]).all()
        else:
            assert set(a_inv[a == 0]) == set()


# TODO: Extend test to generate sparse inputs as well
# TODO: Make test to generate two different arrays a1, a2
# TODO: Check why tests fail with assert_almost_equal
class TestProdSum:
    @given(
        a=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_flat_array(self, a: ndarray, axis: int):
        assert np.allclose((a * a).sum(axis=0), prod_sum(a, a, axis=axis))

    @given(
        a=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_2d_array(self, a: ndarray, axis: int):
        assert np.allclose((a * a).sum(axis=axis), prod_sum(a, a, axis=axis))


# TODO: Extend test to generate sparse inputs as well
class TestSum:
    @given(
        a=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        ),
    )
    def test_flat_arrays(self, a: ndarray):
        a_summed = sum(a=a, axis=0)

        assert_array_equal(a_summed, a.sum(axis=0))

    @given(
        a=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(max_value=1e3, allow_infinity=False, allow_nan=False),
        ),
        axis=st.integers(min_value=0, max_value=1),
    )
    def test_2d_arrays(self, a: ndarray, axis: int):
        a_summed = sum(a=a, axis=axis)

        if a.ndim == 1:
            axis = 0

        assert_array_equal(a_summed, a.sum(axis=axis))
