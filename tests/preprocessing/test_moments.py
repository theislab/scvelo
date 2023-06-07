from typing import List

import pytest

import numpy as np
from scipy.sparse import csr_matrix, issparse

from anndata import AnnData

from scvelo.preprocessing.moments import (
    get_moments,
    moments,
    second_order_moments,
    second_order_moments_u,
)


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
        np.testing.assert_allclose(
            first_order_moment,
            first_order_moment_ground_truth,
            rtol=1e-6,
            atol=1e-6,
        )

        second_order_moment_uncentered = get_moments(
            adata=adata, second_order=True, centered=False
        )
        second_order_moment_uncentered_ground_truth = np.array(
            [[2, 13.25 / 3, 5 / 3], [2.5, 6.5, 0.5], [1, 2.125, 2]]
        )
        np.testing.assert_allclose(
            second_order_moment_uncentered,
            second_order_moment_uncentered_ground_truth,
            rtol=1e-6,
            atol=1e-6,
        )

        second_order_moment_centered = get_moments(adata=adata, second_order=True)
        np.testing.assert_allclose(
            second_order_moment_centered,
            second_order_moment_uncentered - first_order_moment_ground_truth**2,
            rtol=1e-6,
            atol=1e-6,
        )


class TestMoments:
    def _compare_adatas(self, adata_1, adata_2):
        # Check `.layers`
        assert set(adata_1.layers) == set(adata_2.layers).union(["Mu", "Ms"])
        for layer in set(adata_1.layers).difference(["Mu", "Ms"]):
            assert (adata_1.layers[layer] != adata_2.layers[layer]).getnnz() == 0

        # Check `.obsm` is unchanged
        assert set(adata_1.obsm) == {"X_pca"}
        np.testing.assert_equal(adata_1.obsm["X_pca"], adata_2.obsm["X_pca"])

        # Check `.obsp` is unchanged
        assert set(adata_1.obsp) == {"distances", "connectivities"}
        assert issparse(adata_1.obsp["connectivities"])
        np.testing.assert_almost_equal(
            adata_1.obsp["connectivities"].A,
            adata_2.obsp["connectivities"].A,
            decimal=4,
        )
        assert issparse(adata_1.obsp["distances"])
        np.testing.assert_almost_equal(
            adata_1.obsp["distances"].A, adata_2.obsp["distances"].A, decimal=4
        )

        # Check `.uns` is unchanged
        assert set(adata_1.uns["pca"]) == {"params", "variance", "variance_ratio"}
        assert adata_1.uns["pca"]["params"] == adata_2.uns["pca"]["params"]
        np.testing.assert_equal(
            adata_1.uns["pca"]["variance"], adata_2.uns["pca"]["variance"]
        )
        np.testing.assert_equal(
            adata_1.uns["pca"]["variance_ratio"],
            adata_2.uns["pca"]["variance_ratio"],
        )

        assert set(adata_1.uns["neighbors"]) == {
            "connectivities_key",
            "distances_key",
            "indices",
            "params",
        }
        assert (
            adata_1.uns["neighbors"]["connectivities_key"]
            == adata_2.uns["neighbors"]["connectivities_key"]
        )
        assert (
            adata_1.uns["neighbors"]["distances_key"]
            == adata_2.uns["neighbors"]["distances_key"]
        )
        np.testing.assert_equal(
            adata_1.uns["neighbors"]["indices"],
            adata_2.uns["neighbors"]["indices"],
        )
        assert adata_1.uns["neighbors"]["params"] == adata_2.uns["neighbors"]["params"]

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize(
        "layers_to_remove", [["unspliced"], ["spliced"], ["unspliced", "spliced"]]
    )
    def test_skip_moment_calculation(
        self, adata, capfd, dataset: str, n_obs: int, layers_to_remove: List[str]
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        for layer in layers_to_remove:
            del adata.layers[layer]

        moments(data=adata)

        expected_log = (
            "WARNING: Skipping moments, because un/spliced counts were not found.\n"
        )
        actual_log, _ = capfd.readouterr()
        assert actual_log == expected_log

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize("copy", [True, False])
    def test_moment_calculation(
        self, adata, dataset: str, n_obs: int, mode: str, copy: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        original_adata = adata.copy()
        del adata.layers["Mu"]
        del adata.layers["Ms"]

        returned_adata = moments(data=adata, mode=mode, copy=copy)

        if copy:
            assert isinstance(returned_adata, AnnData)
        else:
            assert returned_adata is None
            returned_adata = adata.copy()

        self._compare_adatas(returned_adata, original_adata)

        # Check calculated moments
        ground_truth_unspliced = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=unspliced-mode={mode}_first_moment.npy"
            )
        )
        np.testing.assert_almost_equal(
            returned_adata.layers["Mu"], ground_truth_unspliced
        )

        ground_truth_spliced = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=spliced-mode={mode}_first_moment.npy"
            )
        )
        np.testing.assert_almost_equal(
            returned_adata.layers["Ms"], ground_truth_spliced
        )

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("mode", ["connectivities", "distances"])
    @pytest.mark.parametrize("copy", [True, False])
    def test_log(self, adata, capfd, dataset: str, n_obs: int, mode: str, copy: bool):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        del adata.layers["Mu"]
        del adata.layers["Ms"]

        _ = moments(data=adata, mode=mode, copy=copy)

        expected_log = f"computing moments based on {mode}\n    finished ("

        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)\n"
        )
        assert actual_log == expected_log

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    @pytest.mark.parametrize("copy", [True, False])
    def test_neighbors_and_moments_calculation(
        self, adata, capfd, dataset: str, n_obs: int, copy: bool
    ):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.uns["neighbors"]["params"]["n_pcs"] = None
        original_adata = adata.copy()
        del adata.layers["Mu"]
        del adata.layers["Ms"]
        del adata.uns["neighbors"]
        adata.obsp = {}

        returned_adata = moments(data=adata, copy=copy)
        if copy:
            assert isinstance(returned_adata, AnnData)
        else:
            assert returned_adata is None
            returned_adata = adata.copy()

        self._compare_adatas(returned_adata, original_adata)

        np.testing.assert_almost_equal(
            returned_adata.layers["Mu"], original_adata.layers["Mu"]
        )
        np.testing.assert_almost_equal(
            returned_adata.layers["Ms"], original_adata.layers["Ms"]
        )

        expected_log = "computing neighbors\n    finished ("

        actual_log, _ = capfd.readouterr()
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'distances' and 'connectivities', weighted adjacency matrices "
            "(adata.obsp)\n"
            "computing moments based on connectivities\n    finished ("
        )
        assert actual_log.startswith(expected_log)

        # `[7:]` removes execution time
        actual_log = actual_log.split(expected_log)[1][7:]
        expected_log = (
            ") --> added \n"
            "    'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)\n"
        )
        assert actual_log == expected_log

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_raw_input(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=True, preprocessed=False)

        moments(data=adata)

        ground_truth_unspliced = np.load(
            file=(
                f"tests/_data/test_moments/moments/dataset={dataset}-n_obs={n_obs}"
                f"first_moment_unspliced.npy"
            )
        )
        np.testing.assert_almost_equal(adata.layers["Mu"], ground_truth_unspliced)

        ground_truth_spliced = np.load(
            file=(
                f"tests/_data/test_moments/moments/dataset={dataset}-n_obs={n_obs}"
                f"first_moment_spliced.npy"
            )
        )
        np.testing.assert_almost_equal(adata.layers["Ms"], ground_truth_spliced)


class TestSecondOrderMoments:
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
            _ = second_order_moments(adata=adata)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_output(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)

        second_order_moment_spliced, second_order_moment_mixed = second_order_moments(
            adata=adata
        )
        assert isinstance(second_order_moment_spliced, np.ndarray)
        assert isinstance(second_order_moment_mixed, np.ndarray)

        ground_truth_spliced = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=spliced-mode=connectivities_second_moment.npy"
            )
        )
        np.testing.assert_almost_equal(
            second_order_moment_spliced, ground_truth_spliced
        )

        ground_truth_mixed = np.load(
            file=(
                f"tests/_data/test_moments/second_order_moments/dataset={dataset}"
                f"-n_obs={n_obs}-mode=connectivities_second_moment_mixed.npy"
            )
        )
        np.testing.assert_almost_equal(second_order_moment_mixed, ground_truth_mixed)

    @pytest.mark.parametrize("dataset", ["pancreas", "dentategyrus"])
    @pytest.mark.parametrize("n_obs", [50, 100])
    def test_adjusted(self, adata, dataset: str, n_obs: int):
        adata = adata(dataset=dataset, n_obs=n_obs, raw=False, preprocessed=True)
        adata.layers["Mu"] = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=unspliced-mode=connectivities_first_moment.npy"
            )
        )
        adata.layers["Ms"] = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=spliced-mode=connectivities_first_moment.npy"
            )
        )

        second_order_moment_spliced, second_order_moment_mixed = second_order_moments(
            adata=adata, adjusted=True
        )
        assert isinstance(second_order_moment_spliced, np.ndarray)
        assert isinstance(second_order_moment_mixed, np.ndarray)

        second_order_spliced = np.load(
            file=(
                f"tests/_data/test_moments/get_moments/dataset={dataset}-n_obs={n_obs}"
                f"-layer=spliced-mode=connectivities_second_moment.npy"
            )
        )
        np.testing.assert_almost_equal(
            second_order_moment_spliced, 2 * second_order_spliced - adata.layers["Ms"]
        )

        second_order_mixed = np.load(
            file=(
                f"tests/_data/test_moments/second_order_moments/dataset={dataset}"
                f"-n_obs={n_obs}-mode=connectivities_second_moment_mixed.npy"
            )
        )
        np.testing.assert_almost_equal(
            second_order_moment_mixed, 2 * second_order_mixed - adata.layers["Mu"]
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
