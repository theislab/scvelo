from typing import Callable

import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, spmatrix

from scvelo.preprocessing.utils import csr_vcorrcoef, get_mean_var


class TestCsrVcorrcoef:
    @pytest.mark.parametrize(
        "X",
        (
            np.zeros(3),
            np.array([1, 0, -4]),
            np.array([-0.3, 0.5, 0.93]),
            np.zeros(shape=(3, 3)),
            np.eye(3),
            np.array([[1, 2, 3], [1, -1, 1]]),
            np.array([[0.1, -0.3, 7.5], [8.3, 0.4, -0.9]]),
        ),
    )
    @pytest.mark.parametrize(
        "y",
        (
            np.zeros(3),
            np.ones(3),
            np.array([1, 0, 0]),
            np.array([1, 2, 3]),
            np.array([-0.24, 0.7, 0.4]),
        ),
    )
    def test_dense_arrays(self, X: np.ndarray, y: np.ndarray):
        pearsonr = csr_vcorrcoef(X=X, y=y)

        if X.ndim == 1:
            np.testing.assert_almost_equal(np.corrcoef(X, y)[0, 1], pearsonr)
        else:
            assert all(
                np.allclose(np.corrcoef(sample, y)[0, 1], corr, equal_nan=True)
                for corr, sample in zip(pearsonr, X)
            )

    @pytest.mark.parametrize(
        "X",
        (
            csr_matrix(np.zeros(3)),
            csr_matrix(np.array([1, 0, -4])),
            csr_matrix(np.array([-0.3, 0.5, 0.93])),
            csr_matrix(np.zeros(shape=(3, 3))),
            csr_matrix(np.eye(3)),
            csr_matrix(np.array([[1, 2, 3], [1, -1, 1]])),
            csr_matrix(np.array([[0.1, -0.3, 7.5], [8.3, 0.4, -0.9]])),
        ),
    )
    @pytest.mark.parametrize(
        "y",
        (
            np.zeros(3),
            np.ones(3),
            np.array([1, 0, 0]),
            np.array([1, 2, 3]),
            np.array([-0.24, 0.7, 0.4]),
        ),
    )
    def test_sparse_arrays(self, X: spmatrix, y: np.ndarray):
        pearsonr = csr_vcorrcoef(X=X, y=y)

        X_dense = X.A.squeeze()

        if X_dense.ndim == 1:
            np.testing.assert_almost_equal(np.corrcoef(X_dense, y)[0, 1], pearsonr)
        else:
            assert all(
                np.allclose(np.corrcoef(sample, y)[0, 1], corr, equal_nan=True)
                for corr, sample in zip(pearsonr, X_dense)
            )


class TestGetMeanVar:
    @given(
        X=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=-1e3, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    @pytest.mark.parametrize("sparse", (True, False))
    @pytest.mark.parametrize("flat_array", (True, False))
    def test_mean_var_w_random_arrays(self, X, sparse, flat_array):
        if flat_array:
            X = X.flatten()
        numpy_mean = X.mean(axis=0)
        numpy_var = np.nan_to_num(X.var(axis=0, ddof=1))

        if sparse and not flat_array:
            X = csr_matrix(X)
        mean, var = get_mean_var(X=X)

        if flat_array:
            assert np.isscalar(mean)
            assert np.isscalar(var)
        else:
            assert mean.shape == (X.shape[1],)
            assert var.shape == (X.shape[1],)
        np.testing.assert_almost_equal(numpy_mean, mean)
        np.testing.assert_almost_equal(numpy_var, var)

    @pytest.mark.parametrize(
        "X, analytic_mean, analytic_variance",
        (
            (np.array([1, 2, 3, np.nan]), 2, 1),
            (np.array([1, 2, 2, np.nan, np.inf]), 5 / 3, 1 / 3),
            (np.array([1, 2, 2, 2, np.nan, np.inf, -np.inf]), 7 / 4, 0.25),
            (
                np.array([[1, 2, 3, np.nan], [1, 2, np.inf, 3]]),
                np.array([1, 2, 3, 3]),
                np.array([0, 0, 0, 0]),
            ),
            (
                np.array([[1, 2, 3, np.nan], [1, 2, np.inf, 3], [2, -np.inf, 1, 3]]),
                np.array([4 / 3, 2, 2, 3]),
                np.array([1 / 3, 0, 2, 0]),
            ),
        ),
    )
    def test_with_nan_and_inf_entries(
        self, X: np.ndarray, analytic_mean: np.ndarray, analytic_variance: np.ndarray
    ):
        nan_or_inf_entries = np.isnan(X) | np.isinf(X) | np.isneginf(X)
        mean, var = get_mean_var(X=X)

        assert (X[nan_or_inf_entries] == 0).all()

        np.testing.assert_almost_equal(mean, analytic_mean)
        np.testing.assert_almost_equal(var, analytic_variance)

    @pytest.mark.parametrize(
        "X, analytic_mean, analytic_variance",
        (
            (
                np.array([[1, 2, 3, np.nan], [1, 2, np.inf, 3]]),
                np.array([1, 2, 3, 3]),
                np.array([0, 0, 0, 0]),
            ),
            (
                np.array([[1, 2, 3, np.nan], [1, 2, np.inf, 3], [2, -np.inf, 1, 3]]),
                np.array([4 / 3, 2, 2, 3]),
                np.array([1 / 3, 0, 2, 0]),
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (csr_matrix, csc_matrix))
    def test_sparse_with_nan_and_inf_entries(
        self,
        X: spmatrix,
        analytic_mean: np.ndarray,
        analytic_variance: np.ndarray,
        sparse_format: Callable,
    ):
        X = sparse_format(X)
        nan_or_inf_entries = np.isnan(X.data) | np.isinf(X.data) | np.isneginf(X.data)
        mean, var = get_mean_var(X=X)

        assert (X.data[nan_or_inf_entries] == 0).all()

        np.testing.assert_almost_equal(mean, analytic_mean)
        np.testing.assert_almost_equal(var, analytic_variance)

    @pytest.mark.parametrize(
        "X",
        (
            np.arange(0, 101),
            np.array([np.arange(0, 101), np.arange(0, 101)]),
            np.array([np.arange(0, 101), np.arange(0, 101)[::-1]]),
        ),
    )
    @pytest.mark.parametrize(
        "percentile, lower_percentile, upper_percentile",
        ((10, 10, 100), (60, 0, 60), ([10, 90], 10, 90), (np.array([10, 90]), 10, 90)),
    )
    def test_percentile_dense_input(
        self, X, percentile, lower_percentile, upper_percentile
    ):
        clipped_array = X.copy()
        clipped_array[clipped_array <= lower_percentile] = lower_percentile
        clipped_array[clipped_array >= upper_percentile] = upper_percentile

        numpy_mean = clipped_array.mean(axis=0)
        numpy_var = np.var(clipped_array, axis=0, ddof=1)

        mean, var = get_mean_var(X=X, perc=percentile)

        np.testing.assert_almost_equal(mean, numpy_mean)
        np.testing.assert_almost_equal(var, numpy_var)

    @pytest.mark.parametrize(
        "X",
        (
            csr_matrix(np.array([np.arange(0, 101), np.arange(0, 101)])),
            csr_matrix(np.array([np.arange(0, 101), np.arange(0, 101)][::-1])),
        ),
    )
    @pytest.mark.parametrize(
        "percentile, lower_percentile, upper_percentile",
        (
            (10, 10.9, 100),
            (60, 1, 60.4),
            ([10, 90], 10.9, 90.1),
            (np.array([10, 90]), 10.9, 90.1),
        ),
    )
    def test_percentile_sparse_input(
        self, X, percentile, lower_percentile, upper_percentile
    ):
        clipped_array = X.A.copy().astype(float)
        clipped_array[
            (clipped_array <= lower_percentile) & (clipped_array != 0)
        ] = lower_percentile
        clipped_array[
            (clipped_array >= upper_percentile) & (clipped_array != 0)
        ] = upper_percentile

        numpy_mean = clipped_array.mean(axis=0)
        numpy_var = np.var(clipped_array, axis=0, ddof=1)

        mean, var = get_mean_var(X=X.copy(), perc=percentile)

        np.testing.assert_almost_equal(mean, numpy_mean)
        np.testing.assert_almost_equal(var, numpy_var)

    @pytest.mark.parametrize(
        "X, analytic_mean, analytic_var",
        (
            (np.array([1, 0, 3]), 2, 2),
            (
                np.array([[1, 2, 3], [0, 5, 0], [0, -1, 1]]),
                np.array([1, 2, 2]),
                np.array([0, 9, 2]),
            ),
            (
                np.array([[1, 2, 3], [np.nan, 5, np.inf], [0, -1, 1]]),
                np.array([1, 2, 2]),
                np.array([0, 9, 2]),
            ),
            (
                csr_matrix(np.array([[1, 2, 3], [0, 5, 0], [-np.inf, -1, 1]])),
                np.array([1, 2, 2]),
                np.array([0, 9, 2]),
            ),
        ),
    )
    def test_ignore_zeros(self, X, analytic_mean, analytic_var):
        mean, var = get_mean_var(X=X, ignore_zeros=True)

        np.testing.assert_almost_equal(mean, analytic_mean)
        np.testing.assert_almost_equal(var, analytic_var)
