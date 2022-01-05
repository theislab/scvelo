import pytest

import numpy as np
from scipy.sparse import csr_matrix

from anndata import AnnData

from scvelo.preprocessing.moments import get_moments, second_order_moments_u


class TestGetMoments:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_neighbor_graph_not_present(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        del adata.uns["neighbors"]

        with pytest.raises(
            ValueError,
            match=(
                "You need to run `pp.neighbors` first to compute a neighborhood graph."
            ),
        ):
            _ = get_moments(adata=adata)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("layer", [None, "unspliced", "spliced"])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize("dense", [True, False])
    def test_first_moments(
        self, adata, dataset: str, n_obs: int, layer: bool, mode: str, dense: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        if dense:
            if layer is None:
                adata.X = adata.X.A
            else:
                adata.layers[layer] = adata.layers[layer].A

        first_order_moment = get_moments(adata=adata, layer=layer, mode=mode)
        assert isinstance(first_order_moment, np.ndarray)

        ground_truth = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer={layer}-mode={mode}_first_moment.npy"
            )
        )
        np.testing.assert_almost_equal(first_order_moment, ground_truth)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("layer", [None, "unspliced", "spliced"])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize("dense", [True, False])
    def test_second_moments(
        self, adata, dataset: str, n_obs: int, layer: bool, mode: str, dense: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        if dense:
            if layer is None:
                adata.X = adata.X.A
            else:
                adata.layers[layer] = adata.layers[layer].A

        second_order_moment = get_moments(
            adata=adata, layer=layer, mode=mode, second_order=True, centered=False
        )
        assert isinstance(second_order_moment, np.ndarray)

        ground_truth = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer={layer}-mode={mode}_second_moment.npy"
            )
        )
        np.testing.assert_almost_equal(second_order_moment, ground_truth)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("layer", [None, "unspliced", "spliced"])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize("dense", [True, False])
    def test_passing_array_for_layer(
        self, adata, dataset: str, n_obs: int, layer: bool, mode: str, dense: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        if dense:
            if layer is None:
                adata.X = adata.X.A
            else:
                adata.layers[layer] = adata.layers[layer].A

        if layer is None:
            first_order_moment = get_moments(adata=adata, layer=adata.X, mode=mode)
        else:
            first_order_moment = get_moments(
                adata=adata, layer=adata.layers[layer], mode=mode
            )

        assert isinstance(first_order_moment, np.ndarray)

        ground_truth = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer={layer}-mode={mode}_first_moment.npy"
            )
        )
        np.testing.assert_almost_equal(first_order_moment, ground_truth)

    @pytest.mark.parametrize("sparse", [True, False])
    def test_analytic_example(self, sparse: bool):
        adata = AnnData(
            X=np.array([[1, 2, 0], [2, 3, 1], [1, 0.5, 2]]),
            obsp={
                "connectivities": csr_matrix(
                    np.array([[0, 0.5, 0.1], [0.5, 0, 0], [0.5, 0, 0]])
                )
            },
            uns={"neighbors": []},
        )
        if sparse:
            adata.X = csr_matrix(adata.X)

        first_order_moment = get_moments(adata=adata)
        first_order_moment_ground_truth = np.array(
            [[4 / 3, 5.5 / 3, 1], [1.5, 2.5, 0.5], [1, 1.25, 1]]
        )
        np.testing.assert_almost_equal(
            first_order_moment, first_order_moment_ground_truth
        )

        second_order_moment_uncentered = get_moments(
            adata=adata, second_order=True, centered=False
        )
        second_order_moment_uncentered_ground_truth = np.array(
            [[2, 13.25 / 3, 5 / 3], [2.5, 6.5, 0.5], [1, 2.125, 2]]
        )
        np.testing.assert_almost_equal(
            second_order_moment_uncentered,
            second_order_moment_uncentered_ground_truth,
            decimal=5,
        )

        second_order_moment_centered = get_moments(adata=adata, second_order=True)
        np.testing.assert_almost_equal(
            second_order_moment_centered,
            second_order_moment_uncentered - first_order_moment_ground_truth ** 2,
        )


class TestSecondOrderMomentsU:
    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_neighbor_graph_not_present(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        del adata.uns["neighbors"]

        with pytest.raises(
            ValueError,
            match=(
                "You need to run `pp.neighbors` first to compute a neighborhood graph."
            ),
        ):
            _ = second_order_moments_u(adata=adata)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_output(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        second_order_moment = second_order_moments_u(adata=adata)
        assert isinstance(second_order_moment, np.ndarray)

        ground_truth = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=unspliced-mode=connectivities_second_moment.npy"
            )
        )
        np.testing.assert_almost_equal(second_order_moment, ground_truth)
