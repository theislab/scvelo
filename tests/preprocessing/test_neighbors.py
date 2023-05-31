from typing import Callable, Dict, Optional

import hypothesis.strategies as st
import pytest
from hypothesis import given

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix, issparse, load_npz, spmatrix
from sklearn.neighbors import NearestNeighbors

from anndata import AnnData

from scvelo.preprocessing.neighbors import (
    _get_hnsw_neighbors,
    _get_rep,
    _get_scanpy_neighbors,
    _get_sklearn_neighbors,
    _set_pca,
    compute_connectivities_umap,
    get_connectivities,
    get_csr_from_indices,
    get_duplicate_cells,
    get_n_neighs,
    get_neighs,
    neighbors,
    neighbors_to_be_recomputed,
    remove_duplicate_cells,
    select_connectivities,
    select_distances,
    set_diagonal,
    verify_neighbors,
)
from tests.core import get_adata


class TestComputeConnectivitiesUmap:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_real_data(self, adata, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        knn_indices = adata.uns["neighbors"]["indices"]

        knn_distances = []
        for row_distance, row_index in zip(adata.obsp["distances"], knn_indices):
            knn_distances.append(row_distance.A[0, row_index])
        knn_distances = np.array(knn_distances)

        distance_matrix, connectivity_matrix = compute_connectivities_umap(
            knn_indices=knn_indices,
            knn_dists=knn_distances,
            n_obs=n_obs,
            n_neighbors=adata.uns["neighbors"]["params"]["n_neighbors"],
        )

        assert isinstance(distance_matrix, csr_matrix)
        np.testing.assert_almost_equal(
            distance_matrix.A, adata.obsp["distances"].A, decimal=4
        )

        assert isinstance(connectivity_matrix, csr_matrix)
        np.testing.assert_almost_equal(
            connectivity_matrix.A, adata.obsp["connectivities"].A, decimal=4
        )


class TestGetConnectivities:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("raw", [True, False])
    @pytest.mark.parametrize("uns", [{}, {"random": 0}])
    @pytest.mark.parametrize("n_neighbors", [None, 15, 30])
    @pytest.mark.parametrize("recurse_neighbors", [True, False])
    def test_neighbors_not_present(
        self,
        adata,
        dataset: str,
        n_obs: int,
        raw: bool,
        uns: Dict,
        n_neighbors: Optional[int],
        recurse_neighbors: bool,
    ):
        if raw:
            adata = adata(dataset=dataset, n_obs=n_obs, raw=True, preprocessed=False)
        else:
            adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.uns = uns

        returned_val = get_connectivities(
            adata=adata, n_neighbors=n_neighbors, recurse_neighbors=recurse_neighbors
        )
        assert returned_val is None

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_neighbors", [5, 15, 29])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    def test_connectivities_with_trimmed_neighbors(
        self, adata, dataset: str, n_obs: int, n_neighbors: Optional[int], mode: str
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        connectivities = get_connectivities(
            adata=adata, mode=mode, n_neighbors=n_neighbors
        )

        assert issparse(connectivities)
        assert (connectivities.getnnz(axis=1) == (n_neighbors + 1)).all()
        assert (connectivities.data == 1 / (n_neighbors + 1)).all()

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_neighbors", [None, 30, 45])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    def test_connectivities_with_original_neighbors(
        self, adata, dataset: str, n_obs: int, n_neighbors: Optional[int], mode: str
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        connectivities = get_connectivities(
            adata=adata, mode=mode, n_neighbors=n_neighbors
        )

        assert issparse(connectivities)
        assert (
            connectivities.getnnz(axis=1) == adata.obsp[mode].getnnz(axis=1) + 1
        ).all()

        for row in connectivities:
            np.testing.assert_equal(row.data, 1 / row.getnnz())

    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize(
        "adjacency_matrix, recursed_neighbors_matrix",
        [
            (
                np.array([[1, 0, 1], [0, 1, 1], [0, 1, 1]]),
                csr_matrix(
                    np.array([[0.4, 0.2, 0.4], [0, 0.5, 0.5], [0, 0.5, 0.5]]),
                ).astype(np.float32),
            ),
            (
                np.array([[1, 0, 1, 0], [1, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]]),
                csr_matrix(
                    np.array(
                        [
                            [0.4, 0.2, 0.4, 0],
                            [0.4, 0.4, 0.2, 0],
                            [0.2, 0.4, 0.4, 0],
                            [0.2, 0.4, 0, 0.4],
                        ]
                    )
                ).astype(np.float32),
            ),
        ],
    )
    def test_recursed_neighbors(
        self,
        mode: str,
        adjacency_matrix: np.ndarray,
        recursed_neighbors_matrix: spmatrix,
    ):
        adata = AnnData(
            np.eye(*adjacency_matrix.shape),
            obsp={mode: csr_matrix(adjacency_matrix)},
            uns={
                "neighbors": {"params": {"n_neighbors": adjacency_matrix[0, :].sum()}}
            },
        )

        connectivities = get_connectivities(
            adata=adata,
            mode=mode,
            recurse_neighbors=True,
        )

        assert issparse(connectivities)
        assert (connectivities != recursed_neighbors_matrix).getnnz() == 0

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_neighbors", [None, 15])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    def test_recursed_neighbors_real_data(
        self, adata, dataset: str, n_obs: int, n_neighbors: Optional[int], mode: str
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        connectivities = get_connectivities(
            adata=adata, mode=mode, n_neighbors=n_neighbors
        )
        ground_truth = load_npz(
            "tests/_data/test_neighbors/get_connectivities/"
            f"dataset={dataset}-n_obs={n_obs}-n_neighbors={n_neighbors}-mode={mode}.npz"
        )

        assert issparse(connectivities)
        np.testing.assert_almost_equal(connectivities.A, ground_truth.A)


class TestGetCsrFromIndices:
    @pytest.mark.parametrize(
        "knn_indices, knn_dists, n_obs, n_neighbors, ground_truth",
        (
            [
                (
                    np.array([[0, 1, 2], [1, 3, 0], [2, 1, 3], [3, 0, 1]]),
                    np.array(
                        [[0, 0.1, 0.2], [0, 0.5, 1.7], [0, 0.01, 0.02], [0, 0.5, 1]]
                    ),
                    4,
                    3,
                    csr_matrix(
                        [
                            [0.0, 0.1, 0.2, 0.0],
                            [1.7, 0.0, 0.0, 0.5],
                            [0.0, 0.01, 0.0, 0.02],
                            [0.5, 1.0, 0.0, 0.0],
                        ]
                    ),
                ),
                (
                    np.array([[0, 1, 2], [1, 3, 0], [2, 1, 3], [3, 0, 1]]),
                    np.array(
                        [[0, 0.1, 0.2], [0, 0.5, 1.7], [0, 0.01, 0.3], [0, 0.5, 1]]
                    ),
                    4,
                    2,
                    csr_matrix(
                        [
                            [0.0, 0.1, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.5],
                            [0.0, 0.01, 0.0, 0.0],
                            [0.5, 0.0, 0.0, 0.0],
                        ]
                    ),
                ),
                (
                    np.array([[0, 1, -1], [1, -1, 0], [2, 1, 3], [3, 0, 1]]),
                    np.array(
                        [[0, 0.1, 0.2], [0, 0.5, 1.7], [0, 0.01, 0.3], [0, 0.5, 1]]
                    ),
                    4,
                    3,
                    csr_matrix(
                        [
                            [0.0, 0.1, 0.0, 0.0],
                            [1.7, 0.0, 0.0, 0.0],
                            [0.0, 0.01, 0.0, 0.3],
                            [0.5, 1.0, 0.0, 0.0],
                        ]
                    ),
                ),
            ]
        ),
    )
    def test_manual_data(
        self,
        knn_indices: np.ndarray,
        knn_dists: np.ndarray,
        n_obs: int,
        n_neighbors: int,
        ground_truth: spmatrix,
    ):
        returned_matrix = get_csr_from_indices(
            knn_indices=knn_indices,
            knn_dists=knn_dists,
            n_obs=n_obs,
            n_neighbors=n_neighbors,
        )

        assert isinstance(returned_matrix, csr_matrix)
        assert (returned_matrix != ground_truth).getnnz() == 0

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_real_data(self, adata, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        knn_indices = adata.uns["neighbors"]["indices"]

        knn_distances = []
        for row_distance, row_index in zip(adata.obsp["distances"], knn_indices):
            knn_distances.append(row_distance.A[0, row_index])
        knn_distances = np.array(knn_distances)

        returned_matrix = get_csr_from_indices(
            knn_indices=knn_indices,
            knn_dists=knn_distances,
            n_obs=n_obs,
            n_neighbors=adata.uns["neighbors"]["params"]["n_neighbors"],
        )

        assert isinstance(returned_matrix, csr_matrix)
        assert (returned_matrix != adata.obsp["distances"]).getnnz() == 0


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


class TestGetHnswNeighbors:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 15])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X_pca(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_hnsw_neighbors(
            adata=adata,
            use_rep="X_pca",
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            num_threads=-1,
        )
        if n_pcs is None:
            n_pcs = adata.obsm["X_pca"].shape[1]

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_hnsw_neighbors/dataset={dataset}-"
                f"n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_hnsw_neighbors/dataset={dataset}-"
                f"n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [15, 30])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_hnsw_neighbors(
            adata=adata,
            use_rep="X",
            n_pcs=n_pcs,
            n_neighbors=n_neighbors,
            num_threads=-1,
        )

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_hnsw_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_hnsw_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))


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

        if isinstance(returned_value, (int, float)):
            assert returned_value == 1
        else:
            assert (returned_value != adata.obsp[mode]).sum() == 0

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


class TestGetRep:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize(
        "use_rep, expected_rep",
        [(None, "X_pca"), ("X", "X"), ("pca", "X_pca"), ("X_pca", "X_pca")],
    )
    @pytest.mark.parametrize("n_pcs", [None, 10, 30, 100])
    def test_pca_not_yet_calculated(
        self,
        adata,
        dataset: str,
        n_obs: int,
        use_rep: Optional[str],
        expected_rep: Optional[str],
        n_pcs: Optional[int],
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        returned_rep = _get_rep(adata=adata, use_rep=use_rep, n_pcs=n_pcs)
        assert returned_rep == expected_rep

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 10, 30, 100])
    @pytest.mark.parametrize("n_vars", [5, 10, 49, 50])
    def test_small_n_vars(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_vars: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata = adata[:, adata.var_names[:n_vars]]

        returned_rep = _get_rep(adata=adata, use_rep=None, n_pcs=n_pcs)
        if n_vars < 50:
            assert returned_rep == "X"
        else:
            assert returned_rep == "X_pca"

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_zero_n_pcs(
        self,
        adata,
        dataset: str,
        n_obs: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        returned_rep = _get_rep(adata=adata, use_rep=None, n_pcs=0)
        assert returned_rep == "X"

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [5, 10, 30])
    def test_X_and_n_pcs_specified(
        self,
        adata,
        capfd,
        dataset: str,
        n_obs: int,
        n_pcs: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        returned_rep = _get_rep(adata=adata, use_rep="X", n_pcs=n_pcs)
        assert returned_rep == "X"

        expected_log = (
            "WARNING: Unexpected pair of parameters: `use_rep='X'` but "
            f"`n_pcs={n_pcs}`. This will only consider the frist {n_pcs} variables "
            f"when calculating the neighbor graph. To use all of `X`, pass "
            "`n_pcs=None`.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log


class TestGetScanpyNeighbors:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 15])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X_pca(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_scanpy_neighbors(
            adata=adata, use_rep="X_pca", n_pcs=n_pcs, n_neighbors=n_neighbors
        )
        if n_pcs is None:
            n_pcs = adata.obsm["X_pca"].shape[1]

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_scanpy_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_scanpy_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [15, 30])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_scanpy_neighbors(
            adata=adata, use_rep="X", n_pcs=n_pcs, n_neighbors=n_neighbors
        )

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_scanpy_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_scanpy_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))


class TestGetSklearnNeighbors:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 15])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X_pca(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_sklearn_neighbors(
            adata=adata, use_rep="X_pca", n_pcs=n_pcs, n_neighbors=n_neighbors
        )
        if n_pcs is None:
            n_pcs = adata.obsm["X_pca"].shape[1]

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_sklearn_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_sklearn_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X_pca'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert isinstance(neighbors, NearestNeighbors)

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [15, 30])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    def test_neighbors_with_X(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        n_neighbors: int,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        neighbors = _get_sklearn_neighbors(
            adata=adata, use_rep="X", n_pcs=n_pcs, n_neighbors=n_neighbors
        )

        ground_truth_distances = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_sklearn_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_distances.npz"
            ),
        )
        ground_truth_connectivities = load_npz(
            file=(
                f"tests/_data/test_neighbors/_get_sklearn_neighbors/dataset={dataset}"
                f"-n_obs={n_obs}-rep='X'-n_pcs={n_pcs}-n_neighbors={n_neighbors}"
                "_connectivites.npz"
            ),
        )

        assert isinstance(neighbors, NearestNeighbors)

        assert hasattr(neighbors, "distances")
        assert issparse(neighbors.distances)
        assert (neighbors.distances.getnnz(axis=1) == n_neighbors - 1).all()
        np.testing.assert_almost_equal(
            neighbors.distances.A, ground_truth_distances.A, decimal=4
        )

        assert hasattr(neighbors, "connectivities")
        assert issparse(neighbors.connectivities)
        assert (neighbors.connectivities.getnnz(axis=1) >= n_neighbors - 1).all()
        assert (neighbors.connectivities != neighbors.connectivities.T).getnnz() == 0
        np.testing.assert_almost_equal(
            neighbors.connectivities.A, ground_truth_connectivities.A, decimal=4
        )

        assert hasattr(neighbors, "knn_indices")
        assert neighbors.knn_indices.shape == (adata.n_obs, n_neighbors)
        np.testing.assert_equal(neighbors.knn_indices[:, 0], np.arange(adata.n_obs))


class TestNeighbors:
    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        method=st.sampled_from(["random", "Scanpy", "UMAP"]),
    )
    def test_provide_not_supported_method(self, adata: AnnData, method: str):
        expected_value_error = (
            f"Provided `method={method}`. Admissible values are `'umap'`, `'sklearn'`, "
            "`'hnsw'`, `'gauss'`, and `'rapids'`."
        )
        with pytest.raises(ValueError, match=rf"{expected_value_error}"):
            neighbors(adata=adata, use_rep="random", method=method)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [15, 30])
    @pytest.mark.parametrize("method", ["hnsw", "sklearn", "umap"])
    @pytest.mark.parametrize("n_neighbors", [15, 30])
    @pytest.mark.parametrize("copy", [True, False])
    def test_output(
        self,
        adata,
        capfd,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        method: str,
        n_neighbors: int,
        copy: bool,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        del adata.layers["Mu"]
        del adata.layers["Ms"]
        adata.obsm = {}
        adata.obsp = {}
        adata.uns = {}
        adata.varm = {}

        returned_adata = neighbors(
            adata=adata, method=method, n_pcs=n_pcs, n_neighbors=n_neighbors, copy=copy
        )

        # Check returned value
        if copy:
            assert isinstance(returned_adata, AnnData)
        else:
            assert returned_adata is None
            returned_adata = adata.copy()

        assert returned_adata.obs_names.equals(adata.obs_names)
        assert returned_adata.var_names.equals(adata.var_names)

        # Check PCA
        assert set(returned_adata.obsm) == {"X_pca"}
        assert returned_adata.obsm["X_pca"].shape == (adata.n_obs, n_pcs)

        assert "pca" in returned_adata.uns
        assert set(returned_adata.uns["pca"]) == {
            "params",
            "variance",
            "variance_ratio",
        }
        assert isinstance(returned_adata.uns["pca"]["params"], Dict)
        assert set(returned_adata.uns["pca"]["params"]) == {
            "use_highly_variable",
            "zero_center",
        }
        if "highly_variable" not in adata.var:
            assert returned_adata.uns["pca"]["params"]["use_highly_variable"] is False
        else:
            assert returned_adata.uns["pca"]["params"]["use_highly_variable"] is True

        # Check data related to neighbor graph
        assert set(returned_adata.obsp) == {"distances", "connectivities"}
        assert issparse(returned_adata.obsp["connectivities"])
        assert issparse(returned_adata.obsp["distances"])
        np.testing.assert_equal(
            returned_adata.obsp["distances"].getnnz(axis=1), n_neighbors - 1
        )

        assert "neighbors" in returned_adata.uns
        assert returned_adata.uns["neighbors"]["connectivities_key"] == "connectivities"
        assert returned_adata.uns["neighbors"]["distances_key"] == "distances"
        assert "indices" in returned_adata.uns["neighbors"]
        assert returned_adata.uns["neighbors"]["indices"].shape == (
            returned_adata.n_obs,
            n_neighbors,
        )
        np.testing.assert_equal(
            returned_adata.uns["neighbors"]["indices"][:, 0],
            np.arange(returned_adata.n_obs),
        )
        assert set(returned_adata.uns["neighbors"]["params"]) == {
            "n_neighbors",
            "method",
            "metric",
            "n_pcs",
            "use_rep",
        }
        assert returned_adata.uns["neighbors"]["params"]["n_neighbors"] == n_neighbors
        assert returned_adata.uns["neighbors"]["params"]["method"] == method
        assert returned_adata.uns["neighbors"]["params"]["metric"] == "euclidean"
        assert returned_adata.uns["neighbors"]["params"]["n_pcs"] == n_pcs
        assert returned_adata.uns["neighbors"]["params"]["use_rep"] == "X_pca"

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_log(self, adata, capfd, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        del adata.layers["Mu"]
        del adata.layers["Ms"]
        adata.obsm = {}
        adata.obsp = {}
        adata.uns = {}
        adata.varm = {}

        neighbors(adata=adata)

        expected_log = "computing neighbors\n    finished ("

        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'distances' and 'connectivities', weighted adjacency matrices "
            "(adata.obsp)\n"
        )
        assert actual_log == expected_log

    # FIXME: See https://github.com/theislab/scvelo/issues/922
    """
    # TODO: Make test more sophisticated to test multiple data matrices
    # TODO: Use additional representations besides `X`
    def test_duplicate_cells(self, capfd):
        adata = AnnData(
            X=np.array(
                [
                    [1.3, 0.2, -0.7],
                    [0.5, 1, -10],
                    [1.31, 0.21, -0.71],
                    [1.3, 0.2, -0.7],
                ]
            ),
            layers={
                "unspliced": np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
                "spliced": np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
            },
        )

        neighbors(adata=adata, use_rep="X_pca")

        expected_log = (
            "WARNING: You seem to have 1 duplicate cells in your "
            "data. Consider removing these via pp.remove_duplicate_cells.\n"
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
        )
        assert actual_log == expected_log
        """


class TestNeighborsToBeRecomputed:
    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        n_neighbors=st.integers(),
        uns=st.sampled_from([{}, {"neighbors": []}, {"neighbors": {"param": []}}]),
    )
    def test_incomplete_uns(self, adata: AnnData, n_neighbors: int, uns: Dict):
        adata.uns = uns
        assert neighbors_to_be_recomputed(adata=adata, n_neighbors=n_neighbors)

    @given(
        adata=get_adata(max_obs=5, max_vars=5),
        n_neighbors_original=st.integers(),
        n_additional_neighbors=st.integers(min_value=1),
    )
    def test_more_neighbors_than_originally_used(
        self, adata: AnnData, n_neighbors_original: int, n_additional_neighbors: int
    ):
        adata.uns = {"neighbors": {"params": {"n_neighbors": n_neighbors_original}}}

        assert neighbors_to_be_recomputed(
            adata=adata, n_neighbors=n_neighbors_original + n_additional_neighbors
        )

    def test_unexpected_distance_matrix(self):
        distance_matrix = np.eye(11)
        distance_matrix[0, :] = 1
        adata = AnnData(np.eye(11), obsp={"distances": csr_matrix(distance_matrix)})

        assert neighbors_to_be_recomputed(adata=adata)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_with_real_data(self, adata, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        assert not neighbors_to_be_recomputed(adata=adata)


class TestRemoveDuplicateCells:
    @pytest.mark.parametrize(
        "X, X_pca, X_without_duplicates, X_pca_without_duplicates, n_duplicates",
        (
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [2, 7], [1, 0]]),
                np.array([[1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [2, 7]]),
                1,
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [2, 7], [0, 0]]),
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [2, 7], [0, 0]]),
                0,
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 0], [0, 1]]),
                np.array([[1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1]]),
                2,
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 0], [1, 1]]),
                np.array([[1, 0, 0], [0, 1, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 1]]),
                1,
            ),
            (
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0], [0, 1, 0]]),
                np.array([[1, 0], [0, 1], [1, 1], [0, 1]]),
                np.array([[1, 0, 0], [0, 1, 0], [1, 0, 0]]),
                np.array([[1, 0], [0, 1], [1, 1]]),
                1,
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
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                np.array([[0], [1], [0.1]]),
                1,
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (None, csr_matrix, csc_matrix))
    def test_with_pca_present(
        self,
        capfd,
        X,
        X_pca,
        X_without_duplicates,
        X_pca_without_duplicates,
        n_duplicates,
        sparse_format,
    ):
        if sparse_format:
            X = sparse_format(X)
        adata = AnnData(X=X, obsm={"X_pca": X_pca})
        remove_duplicate_cells(adata=adata)

        if sparse_format:
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, X_without_duplicates)
        else:
            np.testing.assert_almost_equal(adata.X, X_without_duplicates)
        np.testing.assert_almost_equal(adata.obsm["X_pca"], X_pca_without_duplicates)

        if n_duplicates > 0:
            expected_log = f"Removed {n_duplicates} duplicate cells.\n"
        else:
            expected_log = ""
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    # FIXME: See https://github.com/theislab/scvelo/issues/922
    """
    @pytest.mark.parametrize(
        "X, X_without_duplicates, n_duplicates",
        (
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                1,
            ),
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                0,
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (None, csr_matrix, csc_matrix))
    def test_without_pca_present(
        self, capfd, X, X_without_duplicates, n_duplicates, sparse_format
    ):
        if sparse_format:
            X = sparse_format(X)
        adata = AnnData(X=X)
        remove_duplicate_cells(adata=adata)

        assert "X_pca" in adata.obsm
        assert "pca" in adata.uns
        assert "PCs" in adata.varm

        if sparse_format:
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, X_without_duplicates)
        else:
            np.testing.assert_almost_equal(adata.X, X_without_duplicates)

        if n_duplicates > 0:
            expected_log = f"Removed {n_duplicates} duplicate cells.\n"
        else:
            expected_log = ""
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log
    """

    # FIXME: See https://github.com/theislab/scvelo/issues/922
    """
    @pytest.mark.parametrize(
        "X, X_without_duplicates, n_duplicates",
        (
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                        [1.3, 0.2, -0.7],
                    ]
                ),
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                1,
            ),
            (
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                np.array(
                    [
                        [1.3, 0.2, -0.7],
                        [0.5, 1, -10],
                        [1.31, 0.21, -0.71],
                    ]
                ),
                0,
            ),
        ),
    )
    @pytest.mark.parametrize("sparse_format", (None, csr_matrix, csc_matrix))
    def test_neighbors_recalculated(
        self, capfd, X, X_without_duplicates, n_duplicates, sparse_format
    ):
        if sparse_format:
            X = sparse_format(X)
        adata = AnnData(
            X=X,
            uns={"neighbors": {}},
            obsp={
                "distances": np.eye(X.shape[0]),
                "connectivities": np.eye(X.shape[0]),
            },
        )
        remove_duplicate_cells(adata=adata)
        actual_log, _ = capfd.readouterr()

        assert "X_pca" in adata.obsm
        assert "pca" in adata.uns
        assert "PCs" in adata.varm

        if sparse_format:
            assert issparse(adata.X)
            np.testing.assert_almost_equal(adata.X.A, X_without_duplicates)
        else:
            np.testing.assert_almost_equal(adata.X, X_without_duplicates)

        if n_duplicates > 0:
            assert issparse(adata.obsp["distances"])
            assert issparse(adata.obsp["connectivities"])
            assert adata.uns["neighbors"]["connectivities_key"] == "connectivities"
            assert adata.uns["neighbors"]["distances_key"] == "distances"
            assert isinstance(adata.uns["neighbors"]["indices"], np.ndarray)
            assert adata.uns["neighbors"]["params"] == {
                "n_neighbors": 30,
                "method": "umap",
                "metric": "euclidean",
                "n_pcs": None,
                "use_rep": "X",
            }
            expected_log = (
                f"Removed {n_duplicates} duplicate cells.\n"
                "computing neighbors\n"
                "    finished ("
            )
            # `[7:]` removes execution time
            actual_log = actual_log.split(expected_log)[1][7:]
            expected_log = (
                ") --> added \n"
                "    'distances' and 'connectivities', weighted adjacency matrices "
                "(adata.obsp)\n"
            )
            assert actual_log == expected_log
        else:
            np.testing.assert_almost_equal(adata.obsp["distances"], np.eye(adata.n_obs))
            np.testing.assert_almost_equal(
                adata.obsp["connectivities"], np.eye(adata.n_obs)
            )
            assert adata.uns["neighbors"] == {}
            expected_log = ""
            assert actual_log == expected_log
    """


class TestSelectConnectivities:
    @pytest.mark.parametrize(
        "connectivities, ground_truth_result",
        (
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                csr_matrix(
                    np.array([[0, 0, 0, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                csr_matrix(
                    np.array([[0, 0, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 0, 1]])),
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 0, 1]])),
            ),
        ),
    )
    def test_default(self, connectivities, ground_truth_result):
        adjusted_connectivities = select_connectivities(connectivities=connectivities)

        assert issparse(adjusted_connectivities)
        assert (ground_truth_result != adjusted_connectivities).getnnz() == 0

    @pytest.mark.parametrize(
        "connectivities, n_neighbors, ground_truth_result",
        (
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                0,
                csr_matrix(np.zeros((4, 4))),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 0, 0, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                2,
                csr_matrix(
                    np.array([[0, 0, 0, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 0, 0, 3], [2, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 2, 0, 0], [1, 0, 2, 0]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 0, 0, 3], [2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                2,
                csr_matrix(
                    np.array([[0, 0, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 1, 1]])),
                1,
                csr_matrix(np.array([[0, 0, 2], [2, 0, 0], [0, 1, 0], [0, 0, 1]])),
            ),
        ),
    )
    def test_n_neighbors(self, connectivities, n_neighbors, ground_truth_result):
        adjusted_connectivities = select_connectivities(
            connectivities=connectivities,
            n_neighbors=n_neighbors,
        )

        assert issparse(adjusted_connectivities)
        assert (adjusted_connectivities.getnnz(axis=1) <= n_neighbors).all()
        assert (adjusted_connectivities != ground_truth_result).getnnz() == 0

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_neighbors", [0, 1, 5, 10, 30, None])
    def test_real_data(self, adata, dataset, n_obs, n_neighbors):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        adjusted_connectivities = select_connectivities(
            adata.obsp["connectivities"],
            n_neighbors=n_neighbors,
        )

        assert issparse(adjusted_connectivities)
        if (
            n_neighbors is not None
            and n_neighbors <= adata.obsp["connectivities"].getnnz(axis=1).min()
        ):
            assert (adjusted_connectivities.getnnz(axis=1) == n_neighbors).all()
            assert all(
                [
                    all(
                        adjusted_row.data
                        >= np.sort(original_row.data)[-n_neighbors - 1]
                    )
                    for adjusted_row, original_row in zip(
                        adjusted_connectivities, adata.obsp["connectivities"]
                    )
                ]
            )
        else:
            n_neighbors = adata.obsp["connectivities"].getnnz(axis=1).min()
            assert (adjusted_connectivities.getnnz(axis=1) == n_neighbors).all()
            assert all(
                [
                    all(adjusted_row.data >= np.sort(original_row.data)[-n_neighbors])
                    for adjusted_row, original_row in zip(
                        adjusted_connectivities, adata.obsp["connectivities"]
                    )
                ]
            )


class TestSelectDistances:
    @pytest.mark.parametrize(
        "distances, ground_truth_result",
        (
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                csr_matrix(
                    np.array([[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                csr_matrix(
                    np.array([[0, 1, 2, 0], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 0, 1]])),
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 0, 1]])),
            ),
        ),
    )
    def test_default(self, distances, ground_truth_result):
        adjusted_distances = select_distances(dist=distances)

        assert issparse(adjusted_distances)
        assert (ground_truth_result != adjusted_distances).getnnz() == 0

    @pytest.mark.parametrize(
        "distances, n_neighbors, ground_truth_result",
        (
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                0,
                csr_matrix(np.zeros((4, 4))),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
                2,
                csr_matrix(
                    np.array([[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [1, 0, 0, 0]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 2, 0, 0], [1, 0, 2, 0]])
                ),
                1,
                csr_matrix(
                    np.array([[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [1, 0, 0, 0]])
                ),
            ),
            (
                csr_matrix(
                    np.array([[0, 1, 2, 3], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
                2,
                csr_matrix(
                    np.array([[0, 1, 2, 0], [2, 0, 0, 1], [1, 1, 0, 0], [1, 0, 0, 1]])
                ),
            ),
            (
                csr_matrix(np.array([[0, 1, 2], [2, 0, 1], [1, 1, 0], [1, 1, 1]])),
                1,
                csr_matrix(np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0]])),
            ),
        ),
    )
    def test_n_neighbors(self, distances, n_neighbors, ground_truth_result):
        adjusted_dinstances = select_distances(dist=distances, n_neighbors=n_neighbors)

        assert issparse(adjusted_dinstances)
        assert (adjusted_dinstances.getnnz(axis=1) <= n_neighbors).all()
        assert (adjusted_dinstances != ground_truth_result).getnnz() == 0

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_neighbors", [0, 1, 5, 10, 30, None])
    def test_real_data(self, adata, dataset, n_obs, n_neighbors):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        adjusted_distances = select_distances(
            adata.obsp["distances"],
            n_neighbors=n_neighbors,
        )

        assert issparse(adjusted_distances)
        if (
            n_neighbors is not None
            and n_neighbors <= adata.obsp["distances"].getnnz(axis=1).min()
        ):
            assert (adjusted_distances.getnnz(axis=1) == n_neighbors).all()
            assert all(
                [
                    all(adjusted_row.data <= np.sort(original_row.data)[n_neighbors])
                    for adjusted_row, original_row in zip(
                        adjusted_distances, adata.obsp["distances"]
                    )
                ]
            )
        else:
            n_neighbors = adata.obsp["distances"].getnnz(axis=1).min()
            assert (adjusted_distances.getnnz(axis=1) == n_neighbors).all()
            assert all(
                [
                    all(
                        adjusted_row.data <= np.sort(original_row.data)[n_neighbors - 1]
                    )
                    for adjusted_row, original_row in zip(
                        adjusted_distances, adata.obsp["distances"]
                    )
                ]
            )


class TestSetDiagonal:
    @pytest.mark.parametrize(
        "knn_distances",
        [
            np.array([[0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 4, 2]]),
            np.array(
                [
                    [0, 0.21, 2.4, 0.4],
                    [0, 0.327, 0.3, 0.22],
                    [0, 0.3, 0.5, 1.7],
                ]
            ),
        ],
    )
    @pytest.mark.parametrize(
        "knn_indices",
        [
            np.array([[0, 1, 2, 4], [1, 4, 5, 2], [2, 7, 3, 1]]),
            np.array([[0, 2, 1, 4], [1, 4, 2, 3], [2, 3, 4, 1]]),
        ],
    )
    @pytest.mark.parametrize("remove_diag", [True, False])
    def test_remove_diag(self, knn_distances, knn_indices, remove_diag):
        knn_distances_, knn_indices_ = set_diagonal(
            knn_distances=knn_distances,
            knn_indices=knn_indices,
            remove_diag=remove_diag,
        )

        if remove_diag:
            assert knn_distances_.shape == (3, 3)
            np.testing.assert_equal(knn_distances_, knn_distances[:, 1:])
            np.testing.assert_equal(knn_indices_, knn_indices[:, 1:])
            assert knn_indices_.shape == (3, 3)
        else:
            np.testing.assert_equal(knn_distances_, knn_distances)
            np.testing.assert_equal(knn_indices_, knn_indices)

    @pytest.mark.parametrize(
        "knn_distances",
        [
            np.array([[1, 2, 3], [2, 3, 1], [3, 4, 2]]),
            np.array(
                [
                    [0.21, 2.4, 0.4],
                    [0.327, 0.3, 0.22],
                    [0.3, 0.5, 1.7],
                ]
            ),
        ],
    )
    @pytest.mark.parametrize(
        "knn_indices",
        [
            np.array([[1, 2, 4], [4, 5, 2], [7, 3, 1]]),
            np.array([[2, 1, 4], [4, 2, 3], [3, 4, 1]]),
        ],
    )
    @pytest.mark.parametrize("remove_diag", [True, False])
    def test_set_diagonal(self, knn_distances, knn_indices, remove_diag):
        knn_distances_, knn_indices_ = set_diagonal(
            knn_distances=knn_distances,
            knn_indices=knn_indices,
            remove_diag=remove_diag,
        )

        assert knn_distances_.shape == (3, 4)
        np.testing.assert_equal(knn_distances_[:, 0], np.zeros(3))
        np.testing.assert_equal(knn_distances_[:, 1:], knn_distances)

        assert knn_indices_.shape == (3, 4)
        np.testing.assert_equal(knn_indices_[:, 0], np.arange(3))
        np.testing.assert_equal(knn_indices_[:, 1:], knn_indices)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("remove_diag", [True, False])
    def test_real_data(self, adata, dataset, n_obs, remove_diag):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]

        knn_distances = adata.obsp["distances"][
            np.repeat(np.arange(adata.n_obs), n_neighbors).reshape(adata.n_obs, -1),
            adata.uns["neighbors"]["indices"],
        ].A
        knn_indices = adata.uns["neighbors"]["indices"]

        knn_distances_, knn_indices_ = set_diagonal(
            knn_distances=knn_distances,
            knn_indices=knn_indices,
            remove_diag=remove_diag,
        )

        if remove_diag:
            np.testing.assert_equal(knn_distances_, knn_distances[:, 1:])
            np.testing.assert_equal(knn_indices_, knn_indices[:, 1:])
        else:
            np.testing.assert_equal(knn_distances_, knn_distances)
            np.testing.assert_equal(knn_indices_, knn_indices)


class TestSetPCA:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 5, 10, 30])
    @pytest.mark.parametrize("use_highly_variable", [True, False])
    def test_pca_not_yet_calculated(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        use_highly_variable: bool,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        original_adata = adata.copy()

        del adata.layers["Mu"]
        del adata.layers["Ms"]
        adata.obsm = {}
        adata.obsp = {}
        adata.uns = {}
        adata.varm = {}
        cleaned_adata = adata.copy()

        returned_val = _set_pca(
            adata=adata, n_pcs=n_pcs, use_highly_variable=use_highly_variable
        )

        assert returned_val is None

        assert adata.obs_names.identical(cleaned_adata.obs_names)
        assert adata.var_names.identical(cleaned_adata.var_names)
        assert set(adata.layers) == set(cleaned_adata.layers)
        assert (adata.obs.columns == cleaned_adata.obs.columns).all()
        pd.testing.assert_frame_equal(adata.obs, cleaned_adata.obs)
        pd.testing.assert_frame_equal(adata.var, cleaned_adata.var)
        np.testing.assert_almost_equal(adata.X.A, cleaned_adata.X.A)
        for layer in adata.layers:
            np.testing.assert_almost_equal(
                adata.layers[layer].A,
                cleaned_adata.layers[layer].A,
            )

        if n_pcs is None:
            n_pcs = 30

        assert set(adata.obsm) == {"X_pca"}
        assert adata.obsm["X_pca"].shape == (adata.n_obs, n_pcs)
        np.testing.assert_almost_equal(
            adata.obsm["X_pca"][:, :n_pcs],
            original_adata.obsm["X_pca"][:, :n_pcs],
            decimal=2,
        )

        assert "pca" in adata.uns
        assert set(adata.uns["pca"]) == {"params", "variance", "variance_ratio"}
        assert isinstance(adata.uns["pca"]["params"], Dict)
        assert set(adata.uns["pca"]["params"]) == {"use_highly_variable", "zero_center"}
        if "highly_variable" not in adata.var:
            assert adata.uns["pca"]["params"]["use_highly_variable"] is False
        else:
            assert (
                adata.uns["pca"]["params"]["use_highly_variable"] is use_highly_variable
            )

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [1, 5, 9])
    @pytest.mark.parametrize("use_highly_variable", [True, False])
    @pytest.mark.parametrize("pass_n_pcs", [True, False])
    def test_small_n_pcs(
        self,
        adata,
        capfd,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        use_highly_variable: bool,
        pass_n_pcs: bool,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :n_pcs]

        if pass_n_pcs:
            returned_val = _set_pca(
                adata=adata, n_pcs=n_pcs, use_highly_variable=use_highly_variable
            )
        else:
            returned_val = _set_pca(
                adata=adata, n_pcs=None, use_highly_variable=use_highly_variable
            )

        assert returned_val is None

        if pass_n_pcs:
            expected_log = ""
        else:
            expected_log = (
                f"WARNING: Neighbors are computed on {n_pcs} principal components "
                "only.\n"
            )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 50])
    @pytest.mark.parametrize("use_highly_variable", [True, False])
    def test_n_pcs_exceed_n_vars(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        use_highly_variable: bool,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        del adata.layers["Mu"]
        del adata.layers["Ms"]
        adata.obsm = {}
        adata.obsp = {}
        adata.uns = {}
        adata.varm = {}
        adata = adata[:, adata.var_names[:15]]

        returned_val = _set_pca(
            adata=adata, n_pcs=n_pcs, use_highly_variable=use_highly_variable
        )
        assert returned_val is None

        assert set(adata.obsm) == {"X_pca"}
        assert adata.obsm["X_pca"].shape == (adata.n_obs, adata.n_vars - 1)

        assert "pca" in adata.uns
        assert set(adata.uns["pca"]) == {"params", "variance", "variance_ratio"}
        assert isinstance(adata.uns["pca"]["params"], Dict)
        assert set(adata.uns["pca"]["params"]) == {"use_highly_variable", "zero_center"}
        if "highly_variable" not in adata.var:
            assert adata.uns["pca"]["params"]["use_highly_variable"] is False
        else:
            assert (
                adata.uns["pca"]["params"]["use_highly_variable"] is use_highly_variable
            )

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("n_pcs", [None, 50])
    @pytest.mark.parametrize("use_highly_variable", [True, False])
    def test_n_pcs_exceed_n_obs(
        self,
        adata,
        dataset: str,
        n_obs: int,
        n_pcs: Optional[int],
        use_highly_variable: bool,
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        del adata.layers["Mu"]
        del adata.layers["Ms"]
        adata.obsm = {}
        adata.obsp = {}
        adata.uns = {}
        adata.varm = {}
        adata = adata[adata.obs_names[:15], :]

        returned_val = _set_pca(
            adata=adata, n_pcs=n_pcs, use_highly_variable=use_highly_variable
        )

        assert returned_val is None

        assert set(adata.obsm) == {"X_pca"}
        assert adata.obsm["X_pca"].shape == (adata.n_obs, adata.n_obs - 1)

        assert "pca" in adata.uns
        assert set(adata.uns["pca"]) == {"params", "variance", "variance_ratio"}
        assert isinstance(adata.uns["pca"]["params"], Dict)
        assert set(adata.uns["pca"]["params"]) == {"use_highly_variable", "zero_center"}
        if "highly_variable" not in adata.var:
            assert adata.uns["pca"]["params"]["use_highly_variable"] is False
        else:
            assert (
                adata.uns["pca"]["params"]["use_highly_variable"] is use_highly_variable
            )


class TestVerifyNeighbors:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_invalid_graph(self, capfd, adata, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.obsp["distances"][0, :] = 0
        adata.obsp["distances"].eliminate_zeros()

        returned_val = verify_neighbors(adata=adata)
        assert returned_val is None

        actual_warning, _ = capfd.readouterr()
        expected_warning = (
            "WARNING: The neighbor graph has an unexpected format "
            "(e.g. computed outside scvelo) \n"
            "or is corrupted (e.g. due to subsetting). "
            "Consider recomputing with `pp.neighbors`.\n"
        )
        assert actual_warning == expected_warning

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("raw", [True, False])
    @pytest.mark.parametrize(
        "uns", [{}, {"neighbors": []}, {"neighbors": {"param": []}}, {"random": 0}]
    )
    def test_neighbors_or_params_not_present(
        self, capfd, adata, dataset, n_obs, raw, uns
    ):
        if raw:
            adata = adata(dataset=dataset, n_obs=n_obs, raw=True, preprocessed=False)
        else:
            adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.uns = uns

        returned_val = verify_neighbors(adata=adata)
        assert returned_val is None

        actual_warning, _ = capfd.readouterr()
        expected_warning = (
            "WARNING: The neighbor graph has an unexpected format "
            "(e.g. computed outside scvelo) \n"
            "or is corrupted (e.g. due to subsetting). "
            "Consider recomputing with `pp.neighbors`.\n"
        )
        assert actual_warning == expected_warning

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_valid_graph(self, capfd, adata, dataset, n_obs):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        returned_val = verify_neighbors(adata=adata)
        assert returned_val is None

        actual_warning, _ = capfd.readouterr()
        expected_warning = ""
        assert actual_warning == expected_warning
