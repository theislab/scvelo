from typing import Callable, Dict, List, Optional, Union

import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix, issparse, spmatrix

from anndata import AnnData

from scvelo.preprocessing.utils import (
    _filter,
    check_if_valid_dtype,
    counts_per_cell_quantile,
    csr_vcorrcoef,
    filter_and_normalize,
    filter_genes,
    filter_genes_dispersion,
    get_mean_var,
    log1p,
    materialize_as_ndarray,
    normalize_per_cell,
    not_yet_normalized,
    recipe_velocity,
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


class TestCountsPerCellQuantile:
    @given(
        data=st.data(),
        X=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        max_proportion_per_cell=st.floats(
            min_value=0, max_value=1, allow_infinity=False, allow_nan=False
        ),
        provide_counts_per_cell=st.booleans(),
        sparse=st.sampled_from([False, csr_matrix, csc_matrix]),
    )
    def test_with_random_input(
        self,
        data,
        X: np.ndarray,
        max_proportion_per_cell,
        provide_counts_per_cell,
        sparse: Union[bool, Callable],
    ):
        n_obs = X.shape[0]
        if provide_counts_per_cell:
            counts_per_cell = data.draw(
                arrays(
                    int,
                    shape=st.integers(min_value=X.shape[0], max_value=X.shape[0]),
                    elements=st.integers(min_value=0, max_value=1e3),
                )
            )
        else:
            counts_per_cell = None
        if sparse:
            X = sparse(X)
        counts = counts_per_cell_quantile(
            X=X,
            max_proportion_per_cell=max_proportion_per_cell,
            counts_per_cell=counts_per_cell,
        )

        if sparse:
            assert issparse(X)
        else:
            assert isinstance(X, np.ndarray)
        assert counts.shape == (n_obs,)

    @pytest.mark.parametrize(
        "X",
        (
            np.ones(shape=(3, 100)),
            np.hstack((np.ones((5, 100)), np.zeros((5, 5)))),
            np.array([[1, 4, 95], [4, 7, 89]]),
        ),
    )
    @pytest.mark.parametrize(
        "max_proportion_per_cell", (0.0, 0.005, 0.1, 0.3, 0.5, 0.9, 1.0)
    )
    @pytest.mark.parametrize("sparse", (False, csr_matrix, csc_matrix))
    def test_max_proportion_per_cell(self, X, max_proportion_per_cell, sparse):
        ground_truth = X[:, (X <= max_proportion_per_cell * 100).all(axis=0)].sum(
            axis=1
        )

        if sparse:
            X = sparse(X)

        counts = counts_per_cell_quantile(
            X=X, max_proportion_per_cell=max_proportion_per_cell
        )

        np.testing.assert_almost_equal(counts, ground_truth)

    @given(
        X=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    @pytest.mark.parametrize("sparse", (False, csr_matrix, csc_matrix))
    def test_include_all_variables(self, X: np.ndarray, sparse: Union[bool, Callable]):
        counts_with_all_vars = X.sum(axis=1)

        if sparse:
            X = sparse(X)
        counts = counts_per_cell_quantile(X=X, max_proportion_per_cell=1)

        np.testing.assert_almost_equal(counts_with_all_vars, counts)

    @given(
        X=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
    )
    @pytest.mark.parametrize("sparse", (False, csr_matrix, csc_matrix))
    def test_exclude_all_variables(self, X: np.ndarray, sparse: Union[bool, Callable]):
        if sparse:
            X = sparse(X)
        counts = counts_per_cell_quantile(X=X, max_proportion_per_cell=0)

        np.testing.assert_almost_equal(np.zeros(X.shape[0]), counts)

    @given(
        X=arrays(
            float,
            shape=st.tuples(
                st.integers(min_value=1, max_value=100),
                st.integers(min_value=1, max_value=100),
            ),
            elements=st.floats(
                min_value=0, max_value=1e3, allow_infinity=False, allow_nan=False
            ),
        ),
        max_proportion_per_cell=st.floats(min_value=0, max_value=1),
        sparse=st.sampled_from([False, csr_matrix, csc_matrix]),
    )
    def test_counts_per_cell_specified(
        self,
        X: np.ndarray,
        max_proportion_per_cell: float,
        sparse: Union[bool, Callable],
    ):
        ground_truth = X[:, (X <= max_proportion_per_cell * 100).all(axis=0)].sum(
            axis=1
        )

        if sparse:
            X = sparse(X)
        counts = counts_per_cell_quantile(
            X=X,
            max_proportion_per_cell=max_proportion_per_cell,
            counts_per_cell=100 * np.ones(X.shape[0]),
        )

        np.testing.assert_almost_equal(ground_truth, counts)


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


class TestFilterAndNormalize:
    @pytest.mark.parametrize("sparse_X", (False, csr_matrix, csc_matrix))
    @pytest.mark.parametrize("sparse_layers", (False, csr_matrix, csc_matrix))
    @pytest.mark.parametrize("log", (True, False))
    def test_X_looks_processed(
        self, capfd, sparse_X: bool, sparse_layers: bool, log: bool
    ):
        original_X = sparse_X(np.eye(10)) if sparse_X else np.eye(10)
        if sparse_layers:
            layers_value = sparse_layers(np.triu(np.ones(10), k=0))
        else:
            layers_value = np.triu(np.ones(10), k=0)
        adata = AnnData(
            X=original_X.copy(),
            layers={"unspliced": layers_value, "spliced": layers_value},
        )

        filter_and_normalize(adata, log=log)

        assert type(adata.X) is type(original_X)
        if not log:
            if issparse(original_X):
                assert (adata.X != original_X).sum() == 0
            else:
                assert (adata.X == original_X).all()

        expected_log = "Normalized count data: X, spliced, unspliced.\n"
        if log:
            expected_log += "Logarithmized X.\n"

        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize("sparse_X", (False, csr_matrix, csc_matrix))
    @pytest.mark.parametrize("sparse_layers", (False, csr_matrix, csc_matrix))
    @pytest.mark.parametrize("log", (True, False))
    def test_X_log_advised(self, capfd, sparse_X: bool, sparse_layers: bool, log: bool):
        original_X = sparse_X(np.eye(10)) if sparse_X else np.eye(10)
        if sparse_layers:
            layers_value = sparse_layers(np.eye(10))
        else:
            layers_value = np.eye(10)
        adata = AnnData(
            X=original_X.copy(),
            layers={"unspliced": layers_value, "spliced": layers_value},
        )

        filter_and_normalize(adata, log=log)

        assert type(adata.X) is type(original_X)
        if not log and issparse(original_X):
            assert (adata.X != original_X).sum() == 0
        elif not log:
            assert (adata.X == original_X).all()

        expected_log = "Normalized count data: X, spliced, unspliced.\n"
        if log:
            expected_log += "Logarithmized X.\n"

        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_pancreas_50obs(self, capfd, pancreas_50obs: AnnData):
        filter_and_normalize(
            pancreas_50obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
            n_top_genes=5,
        )

        assert pancreas_50obs.shape == (50, 5)
        assert pancreas_50obs.var_names.equals(
            pd.Index(["Ppy", "Isl1", "Ppp1r1a", "Lrpprc", "Nnat"])
        )
        expected_log = (
            "Filtered out 27530 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            "Extracted 5 highly variable genes.\n"
            "Logarithmized X.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_pancreas_100obs(self, capfd, pancreas_100obs: AnnData):
        filter_and_normalize(
            pancreas_100obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
            n_top_genes=5,
        )

        assert pancreas_100obs.shape == (100, 5)
        assert pancreas_100obs.var_names.equals(
            pd.Index(["Ppy", "Gch1", "Ppp1r1a", "Sst", "Maged2"])
        )
        expected_log = (
            "Filtered out 27029 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            "Extracted 5 highly variable genes.\n"
            "Logarithmized X.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_counts_after_normalization_pancreas_50obs(self, pancreas_50obs: AnnData):
        filter_and_normalize(
            pancreas_50obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
            log=False,
            use_initial_size=False,
        )
        np.testing.assert_almost_equal(pancreas_50obs.X.sum(axis=1), 1)
        np.testing.assert_almost_equal(
            pancreas_50obs.layers["unspliced"].sum(axis=1), 1, decimal=6
        )
        np.testing.assert_almost_equal(pancreas_50obs.layers["spliced"].sum(axis=1), 1)

    @pytest.mark.parametrize("log", (True, False))
    def test_without_spliced(self, capfd, log: bool):
        adata = AnnData(np.triu(np.ones(10), k=0))
        filter_and_normalize(adata, counts_per_cell_after=1, log=log)

        expected_log = (
            "WARNING: Could not find spliced / unspliced counts.\n"
            "Normalized count data: X.\n"
        )
        if log:
            expected_log += "Logarithmized X.\n"

        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "layers",
        (
            {},
            {"unspliced": np.eye(10)},
            {"unspliced": np.eye(10), "random_layer": np.eye(10)},
            {"spliced": np.eye(10)},
        ),
    )
    def test_unspliced_spliced_not_found(self, capfd, layers: Dict):
        adata = AnnData(np.eye(10), layers=layers)
        filter_and_normalize(adata)

        expected_log_start = "WARNING: Could not find spliced / unspliced counts.\n"
        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log_start)

    @pytest.mark.parametrize("layers_to_normalize", (None, "all"))
    def test_enforce_normalization(self, capfd, layers_to_normalize: Optional[str]):
        adata = AnnData(
            0.1 * np.eye(10),
            layers={
                "unspliced": 0.1 * np.eye(10),
                "spliced": 0.1 * np.eye(10),
                "random_layer": 0.1 * np.eye(10),
            },
        )
        filter_and_normalize(adata, layers_normalize=layers_to_normalize)

        if layers_to_normalize is not None:
            expected_log = (
                "Normalized count data: X, unspliced, spliced, random_layer.\n"
            )
        else:
            expected_log = (
                "WARNING: Did not normalize X as it looks processed already. To "
                "enforce normalization, set `enforce=True`.\n"
                "WARNING: Did not normalize spliced as it looks processed already. To "
                "enforce normalization, set `enforce=True`.\n"
                "WARNING: Did not normalize unspliced as it looks processed already. "
                "To enforce normalization, set `enforce=True`.\n"
            )
        expected_log += "Logarithmized X.\n"
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log


class TestFilterGenes:
    def get_vars_to_keep(self, data, adata: AnnData, as_string: bool):
        """Draw subset from variable names."""

        if as_string:
            return data.draw(st.sampled_from(adata.var_names.to_list()))
        else:
            return data.draw(
                st.lists(
                    st.sampled_from(adata.var_names.to_list()),
                    min_size=1,
                    max_size=len(adata.var_names),
                    unique=True,
                )
            )

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced"],
        ),
        copy=st.booleans(),
        min_counts=st.integers(min_value=0, max_value=100),
        max_counts=st.integers(min_value=0, max_value=100),
        min_counts_u=st.integers(min_value=0, max_value=100),
        max_counts_u=st.integers(min_value=0, max_value=100),
        min_cells=st.integers(min_value=0, max_value=100),
        max_cells=st.integers(min_value=0, max_value=100),
        min_cells_u=st.integers(min_value=0, max_value=100),
        max_cells_u=st.integers(min_value=0, max_value=100),
        pass_var_to_keep_as_string=st.booleans(),
    )
    def test_retain_genes_and_copy(
        self,
        data,
        adata: AnnData,
        copy: bool,
        min_counts: int,
        max_counts: int,
        min_counts_u: int,
        max_counts_u: int,
        min_cells: int,
        max_cells: int,
        min_cells_u: int,
        max_cells_u: int,
        pass_var_to_keep_as_string: bool,
    ):
        vars_to_keep = self.get_vars_to_keep(
            data=data,
            adata=adata,
            as_string=pass_var_to_keep_as_string,
        )

        returned_adata = filter_genes(
            adata,
            retain_genes=vars_to_keep,
            copy=copy,
            min_counts=min_counts,
            max_counts=max_counts,
            min_counts_u=min_counts_u,
            max_counts_u=max_counts_u,
            min_cells=min_cells,
            max_cells=max_cells,
            min_cells_u=min_cells_u,
            max_cells_u=max_cells_u,
        )

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        if isinstance(vars_to_keep, str):
            assert pd.Index([vars_to_keep]).isin(adata.var_names).all()
        else:
            assert pd.Index(vars_to_keep).isin(adata.var_names).all()

    @given(
        data=st.data(),
        adata=get_adata(
            max_obs=5,
            max_vars=5,
            layer_keys=["unspliced", "spliced"],
        ),
        copy=st.booleans(),
        min_shared_counts=st.integers(min_value=0, max_value=100),
        min_shared_cells=st.integers(min_value=0, max_value=100),
        pass_var_to_keep_as_string=st.booleans(),
    )
    def test_retain_genes_and_copy_w_shared_counts(
        self,
        data,
        adata: AnnData,
        copy: bool,
        min_shared_counts: int,
        min_shared_cells: int,
        pass_var_to_keep_as_string: bool,
    ):
        vars_to_keep = self.get_vars_to_keep(
            data=data,
            adata=adata,
            as_string=pass_var_to_keep_as_string,
        )

        returned_adata = filter_genes(
            adata,
            retain_genes=vars_to_keep,
            copy=copy,
            min_shared_counts=min_shared_counts,
            min_shared_cells=min_shared_cells,
        )

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        if isinstance(vars_to_keep, str):
            assert pd.Index([vars_to_keep]).isin(adata.var_names).all()
        else:
            assert pd.Index(vars_to_keep).isin(adata.var_names).all()

    @pytest.mark.parametrize(
        "X, unspliced, spliced",
        (
            (
                np.diag([10, 10, 10, 1, 2, 3]),
                np.diag([10, 10, 10, 1, 2, 3]),
                np.diag([10, 10, 10, 1, 2, 3]),
            ),
            (
                csr_matrix(np.diag([10, 10, 10, 1, 2, 3])),
                csr_matrix(np.diag([10, 10, 10, 1, 2, 3])),
                csr_matrix(np.diag([10, 10, 10, 1, 2, 3])),
            ),
        ),
    )
    @pytest.mark.parametrize(
        "min_counts, filtered_out_low_spliced",
        (
            (None, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([True, True, True, False, True, True])),
            (3, np.array([True, True, True, False, False, True])),
            (4, np.array([True, True, True, False, False, False])),
        ),
    )
    @pytest.mark.parametrize(
        "max_counts, filtered_out_upper_spliced",
        (
            (None, np.ones(6, dtype=bool)),
            (9, np.array([False, False, False, True, True, True])),
            (10, np.ones(6, dtype=bool)),
            (11, np.ones(6, dtype=bool)),
        ),
    )
    @pytest.mark.parametrize(
        "min_counts_u, filtered_out_low_unspliced",
        (
            (None, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([True, True, True, False, True, True])),
            (3, np.array([True, True, True, False, False, True])),
            (4, np.array([True, True, True, False, False, False])),
        ),
    )
    @pytest.mark.parametrize(
        "max_counts_u, filtered_out_upper_unspliced",
        (
            (None, np.ones(6, dtype=bool)),
            (9, np.array([False, False, False, True, True, True])),
            (10, np.ones(6, dtype=bool)),
            (11, np.ones(6, dtype=bool)),
        ),
    )
    def test_min_max_counts_filter(
        self,
        X: Union[np.ndarray, spmatrix],
        unspliced: Union[np.ndarray, spmatrix],
        spliced: Union[np.ndarray, spmatrix],
        min_counts: Optional[int],
        max_counts: Optional[int],
        min_counts_u: Optional[int],
        max_counts_u: Optional[int],
        filtered_out_low_spliced: bool,
        filtered_out_upper_spliced: bool,
        filtered_out_low_unspliced: bool,
        filtered_out_upper_unspliced: bool,
    ):
        var_names = pd.Index([f"var_{var_id}" for var_id in range(X.shape[1])])
        adata = AnnData(
            X=X,
            layers={"unspliced": unspliced, "spliced": spliced},
            var=pd.DataFrame(index=var_names),
        )

        filter_genes(
            adata,
            min_counts=min_counts,
            max_counts=max_counts,
            min_counts_u=min_counts_u,
            max_counts_u=max_counts_u,
        )
        assert pd.Index.equals(
            adata.var_names,
            var_names[
                filtered_out_low_spliced
                & filtered_out_upper_spliced
                & filtered_out_low_unspliced
                & filtered_out_upper_unspliced
            ],
        )

    @pytest.mark.parametrize(
        "X, unspliced, spliced",
        (
            (
                np.triu(np.ones(6), k=0),
                np.triu(np.ones(6), k=0),
                np.triu(np.ones(6), k=0),
            ),
            (
                csr_matrix(np.triu(np.ones(6), k=0)),
                csr_matrix(np.triu(np.ones(6), k=0)),
                csr_matrix(np.triu(np.ones(6), k=0)),
            ),
        ),
    )
    @pytest.mark.parametrize(
        "min_cells, filtered_out_low_spliced",
        (
            (None, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([False, True, True, True, True, True])),
            (3, np.array([False, False, True, True, True, True])),
            (4, np.array([False, False, False, True, True, True])),
        ),
    )
    @pytest.mark.parametrize(
        "max_cells, filtered_out_upper_spliced",
        (
            (None, np.ones(6, dtype=bool)),
            (5, np.array([True, True, True, True, True, False])),
            (6, np.ones(6, dtype=bool)),
            (7, np.ones(6, dtype=bool)),
        ),
    )
    @pytest.mark.parametrize(
        "min_cells_u, filtered_out_low_unspliced",
        (
            (None, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([False, True, True, True, True, True])),
            (3, np.array([False, False, True, True, True, True])),
            (4, np.array([False, False, False, True, True, True])),
        ),
    )
    @pytest.mark.parametrize(
        "max_cells_u, filtered_out_upper_unspliced",
        (
            (None, np.ones(6, dtype=bool)),
            (5, np.array([True, True, True, True, True, False])),
            (6, np.ones(6, dtype=bool)),
            (7, np.ones(6, dtype=bool)),
        ),
    )
    def test_min_max_cells_filter(
        self,
        X: Union[np.ndarray, spmatrix],
        unspliced: Union[np.ndarray, spmatrix],
        spliced: Union[np.ndarray, spmatrix],
        min_cells,
        max_cells,
        min_cells_u,
        max_cells_u,
        filtered_out_low_spliced,
        filtered_out_upper_spliced,
        filtered_out_low_unspliced,
        filtered_out_upper_unspliced,
    ):
        var_names = pd.Index([f"var_{var_id}" for var_id in range(X.shape[1])])
        adata = AnnData(
            X=X,
            layers={"unspliced": unspliced, "spliced": spliced},
            var=pd.DataFrame(index=var_names),
        )

        filter_genes(
            adata,
            min_cells=min_cells,
            max_counts=max_cells,
            min_cells_u=min_cells_u,
            max_cells_u=max_cells_u,
        )

        assert pd.Index.equals(
            adata.var_names,
            var_names[
                filtered_out_low_spliced
                & filtered_out_upper_spliced
                & filtered_out_low_unspliced
                & filtered_out_upper_unspliced
            ],
        )

    @pytest.mark.parametrize(
        "X, unspliced, spliced",
        (
            (
                np.triu(np.ones(6), k=0),
                np.triu(np.ones(6), k=0)
                + np.diag([1, 2, -1, 1, 4, 2])
                + np.diag([0, 1, 0, 2, 1], k=-1),
                np.triu(np.ones(6), k=0),
            ),
            (
                csr_matrix(np.triu(np.ones(6), k=0)),
                csr_matrix(
                    np.triu(np.ones(6), k=0)
                    + np.diag([1, 2, -1, 1, 4, 2])
                    + np.diag([0, 1, 0, 2, 1], k=-1)
                ),
                csr_matrix(np.triu(np.ones(6), k=0)),
            ),
        ),
    )
    @pytest.mark.parametrize(
        "min_shared_counts, filter_mask_shared_counts",
        (
            (None, np.ones(6, dtype=bool)),
            (0, np.ones(6, dtype=bool)),
            (3, np.array([True, True, True, True, True, True])),
            (4, np.array([False, True, True, True, True, True])),
            (6, np.array([False, True, False, True, True, True])),
        ),
    )
    @pytest.mark.parametrize(
        "min_shared_cells, filter_mask_shared_cells",
        (
            (None, np.ones(6, dtype=bool)),
            (1, np.ones(6, dtype=bool)),
            (2, np.array([False, True, True, True, True, True])),
            (4, np.array([False, False, False, True, True, True])),
        ),
    )
    def test_shared_counts_cells_filter(
        self,
        X: Union[np.ndarray, spmatrix],
        unspliced: Union[np.ndarray, spmatrix],
        spliced: Union[np.ndarray, spmatrix],
        min_shared_counts,
        filter_mask_shared_counts,
        min_shared_cells,
        filter_mask_shared_cells,
    ):
        var_names = pd.Index([f"var_{var_id}" for var_id in range(X.shape[1])])
        adata = AnnData(
            X=X,
            layers={"unspliced": unspliced, "spliced": spliced},
            var=pd.DataFrame(index=var_names),
        )

        filter_genes(
            adata,
            min_shared_counts=min_shared_counts,
            min_shared_cells=min_shared_cells,
        )

        assert pd.Index.equals(
            adata.var_names,
            var_names[filter_mask_shared_counts & filter_mask_shared_cells],
        )

    @pytest.mark.parametrize(
        "min_counts, min_counts_u, n_vars_filtered_out_spliced, "
        "n_vars_filtered_out_unspliced",
        ((0, 1, 0, 0), (2, 0, 1, 0), (0, 3, 0, 2), (4, 3, 3, 0)),
    )
    def test_min_counts_logging(
        self,
        capfd,
        min_counts: int,
        min_counts_u: int,
        n_vars_filtered_out_spliced: int,
        n_vars_filtered_out_unspliced: int,
    ):
        expected_log = ""
        if n_vars_filtered_out_spliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_spliced} genes that are detected "
                f"{min_counts} counts (spliced).\n"
            )
        if n_vars_filtered_out_unspliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_unspliced} genes that are detected "
                f"{min_counts_u} counts (unspliced).\n"
            )

        adata = AnnData(
            X=np.triu(np.ones(6), k=0),
            layers={
                "unspliced": np.triu(np.ones(6), k=0),
                "spliced": np.triu(np.ones(6), k=0),
            },
        )

        filter_genes(adata, min_counts=min_counts, min_counts_u=min_counts_u)
        actual_log, _ = capfd.readouterr()

        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "max_counts, max_counts_u, n_vars_filtered_out_spliced, "
        "n_vars_filtered_out_unspliced",
        ((6, 7, 0, 0), (5, 6, 1, 0), (6, 4, 0, 2), (4, 3, 2, 1)),
    )
    def test_max_counts_logging(
        self,
        capfd,
        max_counts: int,
        max_counts_u: int,
        n_vars_filtered_out_spliced: int,
        n_vars_filtered_out_unspliced: int,
    ):
        expected_log = ""
        if n_vars_filtered_out_spliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_spliced} genes that are detected "
                f"{max_counts} counts (spliced).\n"
            )
        if n_vars_filtered_out_unspliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_unspliced} genes that are detected "
                f"{max_counts_u} counts (unspliced).\n"
            )

        adata = AnnData(
            X=np.triu(np.ones(6), k=0),
            layers={
                "unspliced": np.triu(np.ones(6), k=0),
                "spliced": np.triu(np.ones(6), k=0),
            },
        )

        filter_genes(adata, max_counts=max_counts, max_counts_u=max_counts_u)
        actual_log, _ = capfd.readouterr()

        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "min_cells, min_cells_u, n_vars_filtered_out_spliced, "
        "n_vars_filtered_out_unspliced",
        ((0, 1, 0, 0), (2, 0, 1, 0), (0, 3, 0, 2), (4, 3, 3, 0)),
    )
    def test_min_cells_logging(
        self,
        capfd,
        min_cells: int,
        min_cells_u: int,
        n_vars_filtered_out_spliced: int,
        n_vars_filtered_out_unspliced: int,
    ):
        expected_log = ""
        if n_vars_filtered_out_spliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_spliced} genes that are detected "
                f"in less than {min_cells} cells (spliced).\n"
            )
        if n_vars_filtered_out_unspliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_unspliced} genes that are detected "
                f"in less than {min_cells_u} cells (unspliced).\n"
            )

        adata = AnnData(
            X=np.triu(np.ones(6), k=0),
            layers={
                "unspliced": np.triu(np.ones(6), k=0),
                "spliced": np.triu(np.ones(6), k=0),
            },
        )

        filter_genes(adata, min_cells=min_cells, min_cells_u=min_cells_u)
        actual_log, _ = capfd.readouterr()

        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "max_cells, max_cells_u, n_vars_filtered_out_spliced, "
        "n_vars_filtered_out_unspliced",
        ((6, 7, 0, 0), (5, 6, 1, 0), (6, 4, 0, 2), (4, 3, 2, 1)),
    )
    def test_max_cells_logging(
        self,
        capfd,
        max_cells: int,
        max_cells_u: int,
        n_vars_filtered_out_spliced: int,
        n_vars_filtered_out_unspliced: int,
    ):
        expected_log = ""
        if n_vars_filtered_out_spliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_spliced} genes that are detected "
                f"in more than {max_cells} cells (spliced).\n"
            )
        if n_vars_filtered_out_unspliced > 0:
            expected_log += (
                f"Filtered out {n_vars_filtered_out_unspliced} genes that are detected "
                f"in more than {max_cells_u} cells (unspliced).\n"
            )

        adata = AnnData(
            X=np.triu(np.ones(6), k=0),
            layers={
                "unspliced": np.triu(np.ones(6), k=0),
                "spliced": np.triu(np.ones(6), k=0),
            },
        )

        filter_genes(adata, max_cells=max_cells, max_cells_u=max_cells_u)
        actual_log, _ = capfd.readouterr()

        assert actual_log == expected_log


class TestFilterGenesDispersion:
    @given(
        adata=get_adata(
            max_obs=5,
            max_vars=5,
        ),
        subset=st.booleans(),
        copy=st.booleans(),
    )
    def test_subsetting_and_copy(self, adata: AnnData, subset: bool, copy: bool):
        original_n_obs = adata.n_obs

        returned_adata = filter_genes_dispersion(data=adata, subset=subset, copy=copy)

        if copy:
            assert isinstance(returned_adata, AnnData)
            adata = returned_adata
        else:
            assert returned_adata is None

        if subset:
            assert original_n_obs <= adata.n_obs
        else:
            assert original_n_obs == adata.n_obs

    def test_wrong_flavor(self):
        with pytest.raises(
            ValueError, match=r'`flavor` needs to be "seurat" or "cell_ranger"'
        ):
            filter_genes_dispersion(data=AnnData(np.eye(2)), flavor="random_flavor")

    @pytest.mark.parametrize(
        "flavor, n_top_genes, vars_after_filter",
        (
            ("svr", 5, pd.Index(["Ppy", "Ppp1r1a", "Gcg", "Nnat", "Ins2"])),
            (
                "svr",
                10,
                pd.Index(
                    [
                        "Ppy",
                        "Ppp1r1a",
                        "Lrpprc",
                        "Ttr",
                        "Rbp4",
                        "Gcg",
                        "Chgb",
                        "Nnat",
                        "Iapp",
                        "Ins2",
                    ]
                ),
            ),
            ("seurat", 5, pd.Index(["Ppy", "Ppp1r1a", "Nnat", "Cyr61", "Tyms"])),
            (
                "seurat",
                10,
                pd.Index(
                    [
                        "Dtymk",
                        "Col18a1",
                        "Ppy",
                        "Ppp1r1a",
                        "Lrpprc",
                        "Gcg",
                        "Nnat",
                        "Cyr61",
                        "Tyms",
                        "Hmgb2",
                    ]
                ),
            ),
            ("cell_ranger", 5, pd.Index(["Ppy", "Ppp1r1a", "Rbp4", "Nnat", "Ins2"])),
            (
                "cell_ranger",
                10,
                pd.Index(
                    [
                        "Ppy",
                        "Cdc14b",
                        "Ghr",
                        "Ppp1r1a",
                        "Ttr",
                        "Rbp4",
                        "Gcg",
                        "Nnat",
                        "Spp1",
                        "Ins2",
                    ]
                ),
            ),
        ),
    )
    def test_n_top_genes_pancreas_50obs(
        self,
        capfd,
        pancreas_50obs: AnnData,
        flavor: str,
        n_top_genes: int,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(pancreas_50obs, min_shared_counts=20, copy=True)
        filter_genes_dispersion(adata, flavor=flavor, n_top_genes=n_top_genes)

        assert adata.shape == (50, n_top_genes)
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 27530 genes that are detected 20 counts (shared).\n"
            f"Extracted {n_top_genes} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, n_top_genes, vars_after_filter",
        (
            ("svr", 5, pd.Index(["Ppy", "Sst", "Gcg", "Ghrl", "Ins2"])),
            (
                "svr",
                10,
                pd.Index(
                    [
                        "Ppy",
                        "Ppp1r1a",
                        "Sst",
                        "Gcg",
                        "Chgb",
                        "Nnat",
                        "Ghrl",
                        "Iapp",
                        "Ins2",
                        "Cck",
                    ]
                ),
            ),
            ("seurat", 5, pd.Index(["Cdk1", "Top2a", "Ppy", "Ppp1r1a", "Stmn1"])),
            (
                "seurat",
                10,
                pd.Index(
                    [
                        "Cdk1",
                        "Top2a",
                        "Ppy",
                        "Ppp1r1a",
                        "Sst",
                        "Papss2",
                        "Igfbpl1",
                        "Stmn1",
                        "Tyms",
                        "Sytl4",
                    ]
                ),
            ),
            ("cell_ranger", 5, pd.Index(["Ppy", "Sst", "Gcg", "Ghrl", "Iapp"])),
            (
                "cell_ranger",
                10,
                pd.Index(
                    [
                        "Cdk1",
                        "Ppy",
                        "Pyy",
                        "Ppp1r1a",
                        "Sst",
                        "Gcg",
                        "Igfbpl1",
                        "Ghrl",
                        "Iapp",
                        "Sytl4",
                    ]
                ),
            ),
        ),
    )
    def test_n_top_genes_pancreas_100obs(
        self,
        capfd,
        pancreas_100obs: AnnData,
        flavor: str,
        n_top_genes: int,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(pancreas_100obs, min_shared_counts=20, copy=True)
        filter_genes_dispersion(adata, flavor=flavor, n_top_genes=n_top_genes)

        assert adata.shape == (100, n_top_genes)
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 27029 genes that are detected 20 counts (shared).\n"
            f"Extracted {n_top_genes} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, n_top_genes, vars_after_filter",
        (
            ("svr", 5, pd.Index(["Cxcl14", "Atp9b", "Nlgn1", "Sirt2", "Slc17a7"])),
            (
                "svr",
                10,
                pd.Index(
                    [
                        "Adarb2",
                        "Cxcl14",
                        "Atp9b",
                        "Gad2",
                        "Nlgn1",
                        "Dab1",
                        "Tmsb10",
                        "Sirt2",
                        "Slc17a7",
                        "Spock3",
                    ]
                ),
            ),
            ("seurat", 5, pd.Index(["Adarb2", "Cxcl14", "Atp9b", "Nlgn1", "Sirt2"])),
            (
                "seurat",
                10,
                pd.Index(
                    [
                        "Syt1",
                        "Prkca",
                        "Adarb2",
                        "Cxcl14",
                        "Atp9b",
                        "Apba1",
                        "Nlgn1",
                        "Igfbpl1",
                        "Sirt2",
                        "Slc17a7",
                    ]
                ),
            ),
            (
                "cell_ranger",
                5,
                pd.Index(["Gria1", "Adarb2", "Atp9b", "Nlgn1", "Sirt2"]),
            ),
            (
                "cell_ranger",
                10,
                pd.Index(
                    [
                        "Gria1",
                        "Adarb2",
                        "Cxcl14",
                        "Atp9b",
                        "Rtn3",
                        "Nlgn1",
                        "Tmsb10",
                        "Sirt2",
                        "Spock3",
                        "Rps25",
                    ]
                ),
            ),
        ),
    )
    def test_n_top_genes_dentategyrus_50obs(
        self,
        capfd,
        dentategyrus_50obs: AnnData,
        flavor: str,
        n_top_genes: int,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(dentategyrus_50obs, min_shared_counts=20, copy=True)
        filter_genes_dispersion(adata, flavor=flavor, n_top_genes=n_top_genes)

        assert adata.shape == (50, n_top_genes)
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 13719 genes that are detected 20 counts (shared).\n"
            f"Extracted {n_top_genes} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, n_top_genes, vars_after_filter",
        (
            ("svr", 5, pd.Index(["Cxcl14", "Atp9b", "Cst3", "Npy", "Sirt2"])),
            (
                "svr",
                10,
                pd.Index(
                    [
                        "Cxcl14",
                        "Atp9b",
                        "Gad2",
                        "Slc1a2",
                        "Cst3",
                        "Npy",
                        "Sirt2",
                        "Hs3st4",
                        "Spock3",
                        "Sgcz",
                    ]
                ),
            ),
            ("seurat", 5, pd.Index(["Adarb2", "Cxcl14", "Atp9b", "Npy", "Cpne4"])),
            (
                "seurat",
                10,
                pd.Index(
                    [
                        "Adarb2",
                        "Cxcl14",
                        "Atp9b",
                        "Gad2",
                        "Slc1a2",
                        "Nlgn1",
                        "Npy",
                        "Sirt2",
                        "Hs3st4",
                        "Cpne4",
                    ]
                ),
            ),
            (
                "cell_ranger",
                5,
                pd.Index(["Adarb2", "Atp9b", "Gad2", "Sirt2", "Spock3"]),
            ),
            (
                "cell_ranger",
                10,
                pd.Index(
                    [
                        "Adarb2",
                        "Cxcl14",
                        "Atp9b",
                        "Gad2",
                        "Cst3",
                        "Nlgn1",
                        "Tmsb10",
                        "Sirt2",
                        "Hs3st4",
                        "Spock3",
                    ]
                ),
            ),
        ),
    )
    def test_n_top_genes_dentategyrus_100obs(
        self,
        capfd,
        dentategyrus_100obs: AnnData,
        flavor: str,
        n_top_genes: int,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(dentategyrus_100obs, min_shared_counts=20, copy=True)
        filter_genes_dispersion(adata, flavor=flavor, n_top_genes=n_top_genes)

        assert adata.shape == (100, n_top_genes)
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 13615 genes that are detected 20 counts (shared).\n"
            f"Extracted {n_top_genes} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, min_disp, max_disp, min_mean, max_mean, vars_after_filter",
        (
            ("seurat", None, None, None, None, pd.Index([])),
            (
                "seurat",
                0.0125,
                np.inf,
                0.004,
                3,
                pd.Index(["Pyy", "Gcg", "Gnas", "Iapp", "Rps9"]),
            ),
            ("cell_ranger", None, None, None, None, pd.Index([])),
            (
                "cell_ranger",
                0.0125,
                np.inf,
                0.0055,
                3,
                pd.Index(["Pyy", "Malat1", "Iapp", "Rpl13a"]),
            ),
        ),
    )
    def test_min_max_disp_min_max_mean_pancreas_50obs(
        self,
        capfd,
        pancreas_50obs: AnnData,
        min_disp: float,
        max_disp: float,
        min_mean: float,
        max_mean: float,
        flavor: str,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(pancreas_50obs, min_shared_counts=20, copy=True)
        normalize_per_cell(adata, counts_per_cell_after=1)
        log1p(adata)

        filter_genes_dispersion(
            adata,
            flavor=flavor,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )

        assert adata.shape == (50, len(vars_after_filter))
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 27530 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            f"Extracted {len(vars_after_filter)} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, min_disp, max_disp, min_mean, max_mean, vars_after_filter",
        (
            ("seurat", None, None, None, None, pd.Index([])),
            (
                "seurat",
                0.0125,
                np.inf,
                0.005,
                3,
                pd.Index(["Pyy", "Malat1", "Iapp", "Rpl13a"]),
            ),
            ("cell_ranger", None, None, None, None, pd.Index([])),
            (
                "cell_ranger",
                0.0125,
                np.inf,
                0.0055,
                3,
                pd.Index(["Pyy", "Malat1", "Rpl13a"]),
            ),
        ),
    )
    def test_min_max_disp_min_max_mean_pancreas_100obs(
        self,
        capfd,
        pancreas_100obs: AnnData,
        min_disp: float,
        max_disp: float,
        min_mean: float,
        max_mean: float,
        flavor: str,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(pancreas_100obs, min_shared_counts=20, copy=True)
        normalize_per_cell(adata, counts_per_cell_after=1)
        log1p(adata)

        filter_genes_dispersion(
            adata,
            flavor=flavor,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )

        assert adata.shape == (100, len(vars_after_filter))
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 27029 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            f"Extracted {len(vars_after_filter)} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, min_disp, max_disp, min_mean, max_mean, vars_after_filter",
        (
            ("seurat", None, None, None, None, pd.Index([])),
            ("seurat", 0.0125, np.inf, 0.005, 3, pd.Index(["Ubb", "Fth1"])),
            ("cell_ranger", None, None, None, None, pd.Index([])),
            (
                "cell_ranger",
                0.0125,
                np.inf,
                0.005,
                3,
                pd.Index(["Ubb", "Rpl3", "Fth1"]),
            ),
        ),
    )
    def test_min_max_disp_min_max_mean_dentategyrus_50obs(
        self,
        capfd,
        dentategyrus_50obs: AnnData,
        flavor: str,
        min_disp: float,
        max_disp: float,
        min_mean: float,
        max_mean: float,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(dentategyrus_50obs, min_shared_counts=20, copy=True)
        normalize_per_cell(adata, counts_per_cell_after=1)
        log1p(adata)

        filter_genes_dispersion(
            adata,
            flavor=flavor,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )

        assert adata.shape == (50, len(vars_after_filter))
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 13719 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            f"Extracted {len(vars_after_filter)} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize(
        "flavor, min_disp, max_disp, min_mean, max_mean, vars_after_filter",
        (
            ("seurat", None, None, None, None, pd.Index([])),
            ("seurat", 0.0125, np.inf, 0.005, 3, pd.Index(["Ubb", "Rpl3", "Fth1"])),
            ("cell_ranger", None, None, None, None, pd.Index([])),
            ("cell_ranger", 0.0125, np.inf, 0.005, 3, pd.Index(["Ubb", "Fth1"])),
        ),
    )
    def test_min_max_disp_min_max_mean_dentategyrus_100obs(
        self,
        capfd,
        dentategyrus_100obs: AnnData,
        flavor: str,
        min_disp: float,
        max_disp: float,
        min_mean: float,
        max_mean: float,
        vars_after_filter: pd.Index,
    ):
        adata = filter_genes(dentategyrus_100obs, min_shared_counts=20, copy=True)
        normalize_per_cell(adata, counts_per_cell_after=1)
        log1p(adata)

        filter_genes_dispersion(
            adata,
            flavor=flavor,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )

        assert adata.shape == (100, len(vars_after_filter))
        assert adata.var_names.equals(vars_after_filter)
        if flavor == "svr":
            assert adata.var.columns.equals(pd.Index(["highly_variable"]))
        elif flavor == "seurat":
            assert adata.var.columns.equals(
                pd.Index(
                    ["means", "dispersions", "dispersions_norm", "highly_variable"]
                )
            )

        expected_log = (
            "Filtered out 13615 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            f"Extracted {len(vars_after_filter)} highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize("retain_genes", ("Tram1", ["Tram1", "Ankrd44", "Cryba2"]))
    def test_retain_genes(
        self, capfd, pancreas_50obs: AnnData, retain_genes: Union[str, List[str]]
    ):
        adata = filter_genes(pancreas_50obs, min_shared_counts=20, copy=True)
        normalize_per_cell(adata, counts_per_cell_after=1)
        log1p(adata)

        filter_genes_dispersion(
            adata, flavor="seurat", retain_genes=retain_genes, n_top_genes=2
        )

        if isinstance(retain_genes, str):
            retain_genes = [retain_genes]

        assert adata.shape == (50, len(retain_genes + ["Ppy", "Lrpprc"]))
        assert adata.var_names.equals(pd.Index(retain_genes + ["Ppy", "Lrpprc"]))

        expected_log = (
            "Filtered out 27530 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            f"Extracted {len(retain_genes + ['Ppy', 'Lrpprc'])} highly variable "
            "genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_more_n_top_genes_than_vars(self, capfd, pancreas_50obs: AnnData):
        filter_genes_dispersion(pancreas_50obs, n_top_genes=100000)

        expected_log = (
            "Skip filtering by dispersion since number of variables are less than "
            "`n_top_genes`.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_passing_n_top_genes_and_mean_disp(self, capfd, pancreas_50obs: AnnData):
        filter_genes_dispersion(pancreas_50obs, n_top_genes=100, min_mean=0, max_mean=1)

        expected_log = (
            "If you pass `n_top_genes`, all cutoffs are ignored.\n"
            "Extracted 100 highly variable genes.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log


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


class TestMaterializeAsNdarray:
    @pytest.mark.parametrize(
        "key", (np.array([0, 1, 2]), np.eye(2), np.ones(shape=(2, 3, 4)))
    )
    def test_array(self, key):
        arr = materialize_as_ndarray(key)

        assert isinstance(arr, np.ndarray)

    @pytest.mark.parametrize(
        "key",
        (
            [np.array([0, 1, 2])],
            [np.eye(2), np.ones(shape=(2, 3, 4))],
            (np.eye(2), np.ones(shape=(2, 3, 4))),
            [0, 1, 2],
            (0, 1, 2),
            [["a", "b"], [1.3, 2.7], np.zeros(shape=(2, 3))],
            (["a", "b"], [1.3, 2.7], np.zeros(shape=(2, 3))),
            (("a", "b"), [1.3, 2.7], np.zeros(shape=(2, 3))),
        ),
    )
    def test_list(self, key):
        arr = materialize_as_ndarray(key)

        assert all(isinstance(entry, np.ndarray) for entry in arr)


class TestNormalizePerCell:
    @given(data=st.data(), adata=get_adata(max_obs=5, max_vars=5))
    def test_target_sum_dense(self, data, adata):
        counts_per_cell_after = data.draw(
            arrays(
                float,
                shape=st.integers(min_value=adata.n_obs, max_value=adata.n_obs),
                elements=st.floats(
                    min_value=1e-3, max_value=1e3, allow_infinity=False, allow_nan=False
                ),
            )
        )
        zero_obs = (adata.X == 0).all(axis=1)

        normalize_per_cell(adata, counts_per_cell_after=counts_per_cell_after)

        assert adata.X[zero_obs, :].sum() == 0
        if np.any(~zero_obs):
            np.testing.assert_almost_equal(
                adata.X[~zero_obs, :].sum(axis=1),
                counts_per_cell_after[~zero_obs],
                decimal=4,
            )

    @given(data=st.data(), adata=get_adata(max_obs=5, max_vars=5, sparse_entries=True))
    def test_target_sum_sparse(self, data, adata):
        counts_per_cell_after = data.draw(
            arrays(
                float,
                shape=st.integers(min_value=adata.n_obs, max_value=adata.n_obs),
                elements=st.floats(
                    min_value=1e-3, max_value=1e3, allow_infinity=False, allow_nan=False
                ),
            )
        )
        zero_obs = adata.X.getnnz(axis=1) == 0

        normalize_per_cell(adata, counts_per_cell_after=counts_per_cell_after)

        assert adata.X[zero_obs, :].sum() == 0
        if np.any(~zero_obs):
            np.testing.assert_almost_equal(
                adata.X[~zero_obs, :].sum(axis=1).A1,
                counts_per_cell_after[~zero_obs],
                decimal=4,
            )

    @pytest.mark.parametrize(
        "X, add_column",
        (
            (np.eye(3), False),
            (np.ones((3, 10000)), True),
            (1e-2 * np.ones((3, 10000)), False),
        ),
    )
    def test_adding_gene_count_corr(self, X, add_column):
        adata = AnnData(X=X)
        normalize_per_cell(adata)

        if add_column:
            assert "gene_count_corr" in adata.var.columns
        else:
            assert "gene_count_corr" not in adata.var.columns

    @pytest.mark.parametrize(
        "X", (np.eye(3), csr_matrix(np.eye(3)), csc_matrix(np.eye(3)))
    )
    @pytest.mark.parametrize(
        "obs, counts_per_cell, normed_counts",
        (
            (None, None, np.eye(3)),
            (None, "cell_size", np.eye(3)),
            (pd.DataFrame({"cell_size": 2 * np.ones(3)}), "cell_size", np.eye(3)),
            (
                pd.DataFrame(
                    {"cell_size": 2 * np.ones(3), "initial_size": np.array([2, 1, 2])}
                ),
                "initial_size",
                np.diag([1, 2, 1]),
            ),
        ),
    )
    def test_counts_per_cell_size(self, X, obs, counts_per_cell, normed_counts):
        adata = AnnData(X=X, obs=obs)

        normalize_per_cell(adata, counts_per_cell=counts_per_cell)
        if issparse(adata.X):
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, normed_counts)
        else:
            np.testing.assert_almost_equal(adata.X, normed_counts)

    @pytest.mark.parametrize(
        "X", (np.eye(3), csr_matrix(np.eye(3)), csc_matrix(np.eye(3)))
    )
    @pytest.mark.parametrize(
        "obs, use_initial_size, normed_counts",
        (
            (None, True, np.eye(3)),
            (None, False, np.eye(3)),
            (
                pd.DataFrame({"initial_size": np.array([2, 1, 2])}),
                True,
                np.diag([1, 2, 1]),
            ),
            (pd.DataFrame({"initial_size": np.array([2, 1, 2])}), False, np.eye(3)),
        ),
    )
    def test_use_initial_size(self, X, obs, use_initial_size, normed_counts):
        adata = AnnData(X=X, obs=obs)

        normalize_per_cell(adata, use_initial_size=use_initial_size)
        if issparse(adata.X):
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, normed_counts)
        else:
            np.testing.assert_almost_equal(adata.X, normed_counts)

    @pytest.mark.parametrize(
        "X", (np.eye(2), csr_matrix(np.eye(2)), csc_matrix(np.eye(2)))
    )
    @pytest.mark.parametrize(
        "layers",
        (
            {"unspliced": np.triu(np.ones((2, 2)), k=0)},
            {
                "unspliced": csr_matrix(np.triu(np.ones((2, 2)), k=0)),
                "spliced": csc_matrix(np.triu(np.ones((2, 2)), k=0)),
            },
            {
                "random_name": csr_matrix(np.triu(np.ones((2, 2)), k=0)),
                "spliced": csr_matrix(np.triu(np.ones((2, 2)), k=0)),
            },
        ),
    )
    @pytest.mark.parametrize(
        "layers_to_normalize",
        (
            None,
            "all",
            "unspliced",
            "random_name",
            ["spliced"],
            ["unspliced", "random_name"],
            ["random_name", "spliced"],
            ["random_name", "non_existing_layer"],
        ),
    )
    def test_layers(self, X, layers, layers_to_normalize):
        adata = AnnData(X=X, layers=layers)

        normalize_per_cell(
            data=adata,
            layers=layers_to_normalize,
            counts_per_cell_after=np.array([1, 0.5]),
        )

        if issparse(X):
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, np.diag([1, 0.5]))
        else:
            np.testing.assert_almost_equal(adata.X, np.diag([1, 0.5]))

        if layers_to_normalize is None:
            layers_to_normalize = ["unspliced", "spliced"]
        elif layers_to_normalize == "all":
            layers_to_normalize = [*adata.layers]
        elif isinstance(layers_to_normalize, str):
            layers_to_normalize = [layers_to_normalize]

        normalized_layer = np.array([[0.5, 0.5], [0, 0.5]])
        for layer in adata.layers:
            if layer in layers_to_normalize:
                if issparse(layers[layer]):
                    assert issparse(adata.layers[layer])
                    np.testing.assert_almost_equal(
                        adata.layers[layer].A, normalized_layer
                    )
                else:
                    np.testing.assert_almost_equal(
                        adata.layers[layer], normalized_layer
                    )
            else:
                if issparse(layers[layer]):
                    assert issparse(adata.layers[layer])
                    np.testing.assert_almost_equal(
                        adata.layers[layer].A, layers[layer].A
                    )
                else:
                    np.testing.assert_almost_equal(adata.layers[layer], layers[layer])

    @pytest.mark.parametrize(
        "X, max_proportion_per_cell, normalization_constant",
        (
            (
                np.array([[1, 4, 95], [4, 7, 89]]),
                0.9,
                np.array([5 / 8, 11 / 8])[:, None],
            ),
            (np.array([[1, 4, 95], [4, 7, 89]]), 1, 1),
        ),
    )
    def test_max_proportion_per_cell(
        self, X, max_proportion_per_cell, normalization_constant
    ):
        adata = AnnData(X)
        normalize_per_cell(adata, max_proportion_per_cell=max_proportion_per_cell)

        np.testing.assert_almost_equal(adata.X, X / normalization_constant, decimal=6)

    @pytest.mark.parametrize(
        "X, X_normalized", ((np.eye(3), False), (0.01 * np.eye(3), True))
    )
    def test_logging(self, capfd, X, X_normalized):
        if X_normalized:
            expected_log = (
                "WARNING: Did not normalize X as it looks processed already. To "
                "enforce normalization, set `enforce=True`.\n"
            )
        else:
            expected_log = "Normalized count data: X.\n"

        adata = AnnData(X=X)
        normalize_per_cell(adata)

        actual_log, _ = capfd.readouterr()

        assert actual_log == expected_log


class TestNotYetNormalized:
    @pytest.mark.parametrize(
        "X, normalized",
        (
            (
                (np.eye(3), True),
                (1.001 * np.eye(3), True),
                (csr_matrix(np.eye(3)), True),
                (1.01 * np.eye(3), False),
                (0.1 * np.eye(3), False),
                (csr_matrix(0.1 * np.eye(3)), False),
            )
        ),
    )
    def test_not_yet_normalized(self, X: Union[np.ndarray, spmatrix], normalized: bool):
        normalize_check = not_yet_normalized(X)

        assert normalized == normalize_check


class TestRecipeVelocity:
    def test_pancreas50obs(self, capfd, pancreas_50obs):
        recipe_velocity(pancreas_50obs)

        assert pancreas_50obs.shape == (50, 3571)
        assert pancreas_50obs.obs.columns.equals(
            pd.Index(
                [
                    "initial_size_unspliced",
                    "initial_size_spliced",
                    "initial_size",
                    "n_counts",
                ]
            )
        )
        assert pancreas_50obs.var.columns.equals(pd.Index(["gene_count_corr"]))
        assert [*pancreas_50obs.uns] == ["log1p", "pca", "neighbors"]
        assert [*pancreas_50obs.obsm] == ["X_pca"]
        assert [*pancreas_50obs.varm] == ["PCs"]
        assert [*pancreas_50obs.layers] == ["spliced", "unspliced", "Ms", "Mu"]
        assert [*pancreas_50obs.obsp] == ["distances", "connectivities"]

        expected_log = (
            "Filtered out 19269 genes that are detected 3 counts (spliced).\n"
            "Filtered out 5158 genes that are detected 3 counts (unspliced).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            "Logarithmized X.\n"
            "computing neighbors\n"
            "    finished ("
        )

        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'distances' and 'connectivities', weighted adjacency matrices "
            "(adata.obsp)\n"
            "computing moments based on connectivities\n"
            "    finished ("
        )
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n    'Ms' and 'Mu', moments of un/spliced abundances "
            "(adata.layers)\n"
        )
        assert actual_log.startswith(expected_log)

    def test_dentategyrus50obs(self, capfd, dentategyrus_50obs):
        recipe_velocity(dentategyrus_50obs)

        assert dentategyrus_50obs.shape == (50, 1150)
        assert dentategyrus_50obs.obs.columns.equals(
            pd.Index(
                [
                    "initial_size_unspliced",
                    "initial_size_spliced",
                    "initial_size",
                    "n_counts",
                ]
            )
        )
        assert dentategyrus_50obs.var.columns.equals(pd.Index([]))
        assert [*dentategyrus_50obs.uns] == ["log1p", "pca", "neighbors"]
        assert [*dentategyrus_50obs.obsm] == ["X_pca"]
        assert [*dentategyrus_50obs.varm] == ["PCs"]
        assert [*dentategyrus_50obs.layers] == [
            "ambiguous",
            "spliced",
            "unspliced",
            "Ms",
            "Mu",
        ]
        assert [*dentategyrus_50obs.obsp] == ["distances", "connectivities"]

        expected_log = (
            "Filtered out 7068 genes that are detected 3 counts (spliced).\n"
            "Filtered out 5695 genes that are detected 3 counts (unspliced).\n"
            "Normalized count data: X, spliced, unspliced.\n"
            "Logarithmized X.\n"
            "computing neighbors\n"
            "    finished ("
        )

        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'distances' and 'connectivities', weighted adjacency matrices "
            "(adata.obsp)\n"
            "computing moments based on connectivities\n"
            "    finished ("
        )
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n    'Ms' and 'Mu', moments of un/spliced abundances "
            "(adata.layers)\n"
        )
        assert actual_log.startswith(expected_log)
