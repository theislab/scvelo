from typing import Callable, Dict, Optional, Union

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
    normalize_per_cell,
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

        X_dense = X.toarray().squeeze()

        if X_dense.ndim == 1:
            np.testing.assert_almost_equal(np.corrcoef(X_dense, y)[0, 1], pearsonr)
        else:
            assert all(
                np.allclose(np.corrcoef(sample, y)[0, 1], corr, equal_nan=True)
                for corr, sample in zip(pearsonr, X_dense)
            )


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
    def test_X_looks_processed(self, capfd, sparse_X: bool, sparse_layers: bool):
        original_X = sparse_X(np.eye(10)) if sparse_X else np.eye(10)
        if sparse_layers:
            layers_value = sparse_layers(np.triu(np.ones(10), k=0))
        else:
            layers_value = np.triu(np.ones(10), k=0)
        adata = AnnData(
            X=original_X.copy(),
            layers={"unspliced": layers_value, "spliced": layers_value},
        )

        filter_and_normalize(adata)

        assert type(adata.X) is type(original_X)

        expected_log = "Normalized count data: X, spliced, unspliced.\n"

        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize("sparse_X", (False, csr_matrix, csc_matrix))
    @pytest.mark.parametrize("sparse_layers", (False, csr_matrix, csc_matrix))
    def test_X_log_advised(self, capfd, sparse_X: bool, sparse_layers: bool):
        original_X = sparse_X(np.eye(10)) if sparse_X else np.eye(10)
        if sparse_layers:
            layers_value = sparse_layers(np.eye(10))
        else:
            layers_value = np.eye(10)
        adata = AnnData(
            X=original_X.copy(),
            layers={"unspliced": layers_value, "spliced": layers_value},
        )

        filter_and_normalize(adata)

        assert type(adata.X) is type(original_X)

        expected_log = "Normalized count data: X, spliced, unspliced.\n"

        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_pancreas_50obs(self, capfd, pancreas_50obs: AnnData):
        filter_and_normalize(
            pancreas_50obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
        )

        expected_log = (
            "Filtered out 27530 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_pancreas_100obs(self, capfd, pancreas_100obs: AnnData):
        filter_and_normalize(
            pancreas_100obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
        )

        expected_log = (
            "Filtered out 27029 genes that are detected 20 counts (shared).\n"
            "Normalized count data: X, spliced, unspliced.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    def test_counts_after_normalization_pancreas_50obs(self, pancreas_50obs: AnnData):
        filter_and_normalize(
            pancreas_50obs,
            min_shared_counts=20,
            counts_per_cell_after=1,
            use_initial_size=False,
        )
        np.testing.assert_almost_equal(pancreas_50obs.X.sum(axis=1), 1)
        np.testing.assert_almost_equal(
            pancreas_50obs.layers["unspliced"].sum(axis=1), 1, decimal=6
        )
        np.testing.assert_almost_equal(pancreas_50obs.layers["spliced"].sum(axis=1), 1)

    def test_without_spliced(self, capfd):
        adata = AnnData(np.triu(np.ones(10), k=0))
        filter_and_normalize(adata, counts_per_cell_after=1)

        expected_log = (
            "WARNING: Could not find spliced / unspliced counts.\n"
            "Normalized count data: X.\n"
        )

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
            np.testing.assert_almost_equal(adata.X.toarray(), normed_counts)
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
            np.testing.assert_almost_equal(adata.X.toarray(), normed_counts)
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
            np.testing.assert_almost_equal(adata.X.toarray(), np.diag([1, 0.5]))
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
                        adata.layers[layer].toarray(), normalized_layer
                    )
                else:
                    np.testing.assert_almost_equal(
                        adata.layers[layer], normalized_layer
                    )
            else:
                if issparse(layers[layer]):
                    assert issparse(adata.layers[layer])
                    np.testing.assert_almost_equal(
                        adata.layers[layer].toarray(), layers[layer].toarray()
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
