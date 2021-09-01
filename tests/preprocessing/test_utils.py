from typing import Callable

import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, spmatrix

from anndata import AnnData

from scvelo.preprocessing.utils import (
    _filter,
    check_if_valid_dtype,
    csr_vcorrcoef,
    get_mean_var,
    log1p,
)
from tests.core import get_adata


class TestCheckIfValidDtype:
    @given(adata=get_adata(max_obs=5, max_vars=5))
    def test_check_x(self, adata: AnnData):
        if "X" in adata.layers:
            adata.layers["_X"] = adata.layers.pop("X")
        # anndata stores array as float in `AnnData.X`
        adata.X = adata.X.astype(int)

        check_if_valid_dtype(adata, layer="X")

        assert adata.X.dtype == "float32"
        for layer in adata.layers:
            assert np.issubdtype(adata.layers[layer].dtype, np.integer)

    @given(data=st.data(), adata=get_adata(max_obs=5, max_vars=5))
    def test_check_layers(self, data, adata: AnnData):
        if "X" in adata.layers:
            adata.layers["_X"] = adata.layers.pop("X")
        layer_to_convert = data.draw(st.sampled_from([*adata.layers]))
        # anndata stores array as float in `AnnData.X`
        adata.X = adata.X.astype(int)

        check_if_valid_dtype(adata, layer=layer_to_convert)

        assert np.issubdtype(adata.X.dtype, np.integer)
        for layer in adata.layers:
            if layer == layer_to_convert:
                assert adata.layers[layer].dtype == "float32"
            else:
                assert np.issubdtype(adata.layers[layer].dtype, np.integer)


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


class TestFilter:
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
        min_max_counts=st.lists(
            st.floats(allow_infinity=True, allow_nan=False),
            min_size=2,
            max_size=2,
            unique=True,
        ),
    )
    def test_filter_based_on_counts(self, X, min_max_counts):
        min_counts = min(min_max_counts)
        max_counts = min(min_max_counts)

        filter_mask, counts_per_var = _filter(
            X=X, min_counts=min_counts, max_counts=max_counts
        )

        assert filter_mask.dtype == np.dtype("bool")
        assert filter_mask.shape == (X.shape[1],)
        assert isinstance(counts_per_var, np.ndarray)
        assert counts_per_var.shape == (X.shape[1],)
        assert counts_per_var.dtype == X.dtype

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
        min_max_cells=st.lists(
            st.floats(allow_infinity=True, allow_nan=False),
            min_size=2,
            max_size=2,
            unique=True,
        ),
    )
    def test_filter_based_on_cells(self, X, min_max_cells):
        min_cells = min(min_max_cells)
        max_cells = min(min_max_cells)

        filter_mask, counts_per_var = _filter(
            X=X, min_cells=min_cells, max_cells=max_cells
        )

        assert filter_mask.dtype == np.dtype("bool")
        assert filter_mask.shape == (X.shape[1],)
        assert isinstance(counts_per_var, np.ndarray)
        assert counts_per_var.shape == (X.shape[1],)
        assert np.issubdtype(counts_per_var.dtype, np.integer)

    @pytest.mark.parametrize(
        "min_counts, filtered_out_low",
        (
            (None, np.ones(6, dtype=bool)),
            (0, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([True, True, True, False, True, True])),
            (3, np.array([True, True, True, False, False, True])),
            (4, np.array([True, True, True, False, False, False])),
        ),
    )
    @pytest.mark.parametrize(
        "max_counts, filtered_out_upper",
        (
            (None, np.ones(6, dtype=bool)),
            (9, np.array([False, False, False, True, True, True])),
            (10, np.ones(6, dtype=bool)),
            (11, np.ones(6, dtype=bool)),
        ),
    )
    def test_filter_mask(
        self, min_counts, filtered_out_low, max_counts, filtered_out_upper
    ):
        filter_mask, _ = _filter(
            X=np.diag([10, 10, 10, 1, 2, 3]),
            min_counts=min_counts,
            max_counts=max_counts,
        )

        np.testing.assert_equal(filter_mask, filtered_out_low & filtered_out_upper)


class TestLog1p:
    @given(adata=get_adata(max_obs=5, max_vars=5), copy=st.booleans())
    def test_dense_adata(self, adata: AnnData, copy: bool):
        original_X = adata.X.copy()
        returned_adata = log1p(data=adata, copy=copy)

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        np.testing.assert_almost_equal(adata.X, np.log1p(original_X))

    @given(
        adata=get_adata(max_obs=5, max_vars=5, sparse_entries=True), copy=st.booleans()
    )
    def test_sparse_adata(self, adata: AnnData, copy: bool):
        original_X = adata.X.copy()
        returned_adata = log1p(data=adata, copy=copy)

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        np.testing.assert_almost_equal(adata.X.data, np.log1p(original_X.data))

    @given(
        data=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        copy=st.booleans(),
    )
    def test_array(self, data: np.ndarray, copy: bool):
        original_array = data.copy()
        returned_data = log1p(data=data, copy=copy)

        if copy:
            assert isinstance(returned_data, np.ndarray)
            data = returned_data
        else:
            assert returned_data is None

        np.testing.assert_almost_equal(data, np.log1p(original_array))
