from typing import Callable, Optional, Union

import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix, spmatrix

from anndata import AnnData

from scvelo.preprocessing.utils import (
    _filter,
    check_if_valid_dtype,
    csr_vcorrcoef,
    filter_genes,
    get_mean_var,
    log1p,
    materialize_as_ndarray,
    not_yet_normalized,
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
