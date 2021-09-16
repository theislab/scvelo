import pytest

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, issparse

from anndata import AnnData

from scvelo.preprocessing.neighbors import get_n_neighs, get_neighs


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
