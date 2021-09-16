from typing import Callable, Optional

import pytest

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, issparse

from anndata import AnnData

from scvelo.preprocessing.neighbors import get_duplicate_cells, get_n_neighs, get_neighs


class TestGetDuplicateCells:
    @pytest.mark.parametrize(
        "X, true_duplicate_row_idx",
        (
            (np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]), np.array([2])),
            (np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]), np.array([2, 3])),
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
                np.array([3]),
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (None, csr_matrix, csc_matrix))
    def test_array(
        self,
        X: np.ndarray,
        true_duplicate_row_idx: np.ndarray,
        sparse_format: Optional[Callable],
    ):
        if sparse_format:
            X = sparse_format(X)
        returned_duplicate_idx = get_duplicate_cells(data=X)

        np.testing.assert_almost_equal(returned_duplicate_idx, true_duplicate_row_idx)

    @pytest.mark.parametrize(
        "X, X_pca, true_duplicate_row_idx",
        (
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [2, 7], [1, 0]]),
                np.array([2]),
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [2, 7], [0, 0]]),
                np.array([]),
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 0], [0, 1]]),
                np.array([2, 3]),
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 0], [1, 1]]),
                np.array([2]),
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 1], [0, 1]]),
                np.array([3]),
            ),
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
                np.array([[0], [1], [0.1], [0]]),
                np.array([3]),
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (None, csr_matrix, csc_matrix))
    def test_anndata(
        self,
        X: np.ndarray,
        X_pca: np.ndarray,
        true_duplicate_row_idx: np.ndarray,
        sparse_format: Optional[Callable],
    ):
        if sparse_format:
            X = sparse_format(X)

        adata = AnnData(X=X, obsm={"X_pca": X_pca})
        returned_duplicate_idx = get_duplicate_cells(data=adata)

        np.testing.assert_almost_equal(returned_duplicate_idx, true_duplicate_row_idx)


class TestGetNeighs:
    @pytest.mark.parametrize(
        "adata, ground_truth",
        (
            (
                AnnData(
                    np.eye(5),
                    obsp={"distances": np.eye(5), "connectivities": np.eye(5)},
                    uns={"neighbors": {"distances": 2}},
                ),
                np.eye(5),
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={
                        "distances": csr_matrix(np.eye(5)),
                        "connectivities": np.eye(5),
                    },
                    uns={"neighbors": {"distances": 2}},
                ),
                csr_matrix(np.eye(5)),
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={
                        "connectivities": np.eye(5),
                    },
                    uns={"neighbors": {"distances": csr_matrix(np.eye(5))}},
                ),
                csr_matrix(np.eye(5)),
            ),
            (
                AnnData(
                    np.eye(5),
                    uns={"neighbors": {"distances": csr_matrix(np.eye(5))}},
                ),
                csr_matrix(np.eye(5)),
            ),
        ),
    )
    def test_with_default_values(self, adata: AnnData, ground_truth):
        returned_value = get_neighs(adata=adata)

        if issparse(ground_truth):
            assert (returned_value != ground_truth).sum() == 0
        else:
            np.testing.assert_almost_equal(returned_value, ground_truth)

    @pytest.mark.parametrize(
        "adata, mode",
        (
            (AnnData(np.eye(5), obsp={"distances": np.eye(5)}), "distances"),
            (
                AnnData(
                    np.eye(5),
                    obsp={"distances": np.eye(5)},
                    uns={"neigbors": {"distances": 2}},
                ),
                "distances",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={"distances": np.eye(5), "random_obsp": np.ones(shape=(5, 5))},
                ),
                "random_obsp",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={
                        "distances": csr_matrix(np.eye(5)),
                        "random_obsp": csr_matrix(np.triu(np.ones(shape=(5, 5)), k=0)),
                    },
                ),
                "random_obsp",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={
                        "distances": csr_matrix(np.eye(5)),
                        "random_obsp": csr_matrix(np.triu(np.ones(shape=(5, 5)), k=0)),
                    },
                ),
                "distances",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={
                        "distances": csr_matrix(np.eye(5)),
                        "random_obsp": csc_matrix(np.triu(np.ones(shape=(5, 5)), k=0)),
                    },
                ),
                "random_obsp",
            ),
        ),
    )
    def test_get_from_obsp(self, adata: AnnData, mode: str):
        returned_value = get_neighs(adata=adata, mode=mode)

        if issparse(returned_value):
            assert (returned_value != adata.obsp[mode]).sum() == 0
        else:
            np.testing.assert_almost_equal(returned_value, adata.obsp[mode])

    @pytest.mark.parametrize(
        "adata, mode",
        (
            (
                AnnData(np.eye(5), uns={"neighbors": {"distances": np.eye(5)}}),
                "distances",
            ),
            (
                AnnData(
                    np.eye(5), uns={"neighbors": {"distances": csr_matrix(np.eye(5))}}
                ),
                "distances",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={"distances": csr_matrix(np.eye(5))},
                    uns={"neighbors": {"n_neighbors": 1}},
                ),
                "n_neighbors",
            ),
            (
                AnnData(
                    np.eye(5),
                    obsp={"distances": csr_matrix(np.eye(5))},
                    uns={"neighbors": {"n_neighbors": 1, "random_entry": np.eye(2)}},
                ),
                "n_neighbors",
            ),
        ),
    )
    def test_get_from_uns(self, adata: AnnData, mode: str):
        returned_value = get_neighs(adata=adata, mode=mode)

        if issparse(returned_value):
            assert (returned_value != adata.uns["neighbors"][mode]).sum() == 0
        else:
            np.testing.assert_almost_equal(returned_value, adata.uns["neighbors"][mode])

    def test_selected_mode_not_available(self):
        adata = AnnData(np.eye(2))

        with pytest.raises(ValueError, match=r"The selected mode is not valid."):
            _ = get_neighs(adata=adata, mode="distances")


class TestGetNNeighs:
    @pytest.mark.parametrize(
        "adata, expected_return_value",
        (
            (AnnData(np.eye(2)), 0),
            (AnnData(np.eye(2), uns={"neighbors": {}}), 0),
            (AnnData(np.eye(2), uns={"neighbors": {"random_key": 0}}), 0),
            (AnnData(np.eye(2), uns={"neighbors": {"params": {}}}), 0),
            (AnnData(np.eye(2), uns={"neighbors": {"params": {"random_key": 2}}}), 0),
            (AnnData(np.eye(2), uns={"neighbors": {"params": {"n_neighbors": 5}}}), 5),
        ),
    )
    def test_get_n_neighs(self, adata: AnnData, expected_return_value: int):
        n_neigbors = get_n_neighs(adata=adata)

        assert n_neigbors == expected_return_value
