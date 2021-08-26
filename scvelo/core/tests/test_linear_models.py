import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from numpy import ndarray
from numpy.testing import assert_almost_equal, assert_array_equal

from scvelo.core import LinearRegression


class TestLinearRegression:
    @given(
        x=arrays(
            float,
            shape=st.integers(min_value=1, max_value=100),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        coef=st.floats(
            min_value=-1000, max_value=1000, allow_infinity=False, allow_nan=False
        ),
    )
    def test_perfect_fit(self, x: ndarray, coef: float):
        lr = LinearRegression()
        lr.fit(x, x * coef)

        assert lr.intercept_ == 0
        if set(x) != {0}:  # fit is only unique if x is non-trivial
            assert_almost_equal(lr.coef_, coef)

    @given(
        x=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        coef=arrays(
            float,
            shape=100,
            elements=st.floats(
                min_value=-1000, max_value=1000, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    # TODO: Extend test to use `percentile`. Zero columns (after trimming) make the
    # previous implementation of the unit test fail
    # TODO: Check why test fails if number of columns is increased to e.g. 1000 (500)
    def test_perfect_fit_2d(self, x: ndarray, coef: ndarray):
        coef = coef[: x.shape[1]]
        lr = LinearRegression()
        lr.fit(x, x * coef)

        assert lr.coef_.shape == (x.shape[1],)
        assert lr.intercept_.shape == (x.shape[1],)
        assert_array_equal(lr.intercept_, np.zeros(x.shape[1]))
        if set(x.flatten()) != {0}:  # fit is only unique if x is non-trivial
            assert_almost_equal(lr.coef_, coef)

    # TODO: Use hypothesis
    # TODO: Integrate into `test_perfect_fit_2d`
    @pytest.mark.parametrize(
        "x, coef, intercept",
        [
            (np.array([[0], [1], [2], [3]]), 0, 1),
            (np.array([[0], [1], [2], [3]]), 2, 1),
            (np.array([[0], [1], [2], [3]]), 2, -1),
        ],
    )
    def test_perfect_fit_with_intercept(
        self, x: ndarray, coef: float, intercept: float
    ):
        lr = LinearRegression(fit_intercept=True, positive_intercept=False)
        lr.fit(x, x * coef + intercept)

        assert lr.coef_.shape == (x.shape[1],)
        assert lr.intercept_.shape == (x.shape[1],)
        assert_array_equal(lr.intercept_, intercept)
        assert_array_equal(lr.coef_, coef)
