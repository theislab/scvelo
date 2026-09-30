import pytest

import numpy as np
from scipy.sparse import csr_matrix

from anndata import AnnData

import scvelo as scv


@pytest.fixture
def velocity_graphs():
    positive = csr_matrix(
        [[0, 0.8, 0.5, 0], [0, 0, 0.7, 0.4], [0.3, 0, 0, 0.9], [0.6, 0.4, 0, 0]]
    )
    negative = csr_matrix(
        [[0, 0, 0, -0.4], [-0.8, 0, 0, 0], [0, -0.2, 0, 0], [0, 0, -0.7, 0]]
    )
    return positive, negative


def _adata(positive, negative, obsp=False):
    adata = AnnData(np.ones((4, 3)))
    adata.obs_names = ["a", "b", "c", "d"]
    adata.uns["velocity_graph"] = positive.copy()
    adata.uns["velocity_graph_neg"] = negative.copy()
    if obsp:
        adata.obsp["velocity_graph"] = positive.copy()
        adata.obsp["velocity_graph_neg"] = negative.copy()
    return adata


@pytest.mark.parametrize("use_negative_cosines", [False, True])
@pytest.mark.parametrize("self_transitions", [False, True])
def test_transition_matrix_obsp_tracks_reordered_cells(
    velocity_graphs, use_negative_cosines, self_transitions
):
    adata = _adata(*velocity_graphs, obsp=True)
    kwargs = {
        "use_negative_cosines": use_negative_cosines,
        "self_transitions": self_transitions,
    }
    reference = scv.tl.transition_matrix(adata, **kwargs).toarray()
    permutation = np.array([2, 0, 3, 1])
    reordered = adata[permutation].copy()

    observed = scv.tl.transition_matrix(reordered, **kwargs).toarray()

    np.testing.assert_allclose(observed, reference[np.ix_(permutation, permutation)])


@pytest.mark.parametrize("use_negative_cosines", [False, True])
def test_transition_matrix_obsp_negative_does_not_require_uns_copy(
    velocity_graphs, use_negative_cosines
):
    adata = _adata(*velocity_graphs, obsp=True)
    expected = scv.tl.transition_matrix(
        adata, self_transitions=False, use_negative_cosines=use_negative_cosines
    ).toarray()
    del adata.uns["velocity_graph_neg"]

    observed = scv.tl.transition_matrix(
        adata, self_transitions=False, use_negative_cosines=use_negative_cosines
    ).toarray()

    np.testing.assert_allclose(observed, expected)


@pytest.mark.parametrize("use_negative_cosines", [False, True])
def test_transition_matrix_uns_only_matches_aligned_obsp(
    velocity_graphs, use_negative_cosines
):
    uns_only = _adata(*velocity_graphs)
    aligned_obsp = _adata(*velocity_graphs, obsp=True)
    kwargs = {
        "self_transitions": False,
        "use_negative_cosines": use_negative_cosines,
    }

    np.testing.assert_allclose(
        scv.tl.transition_matrix(uns_only, **kwargs).toarray(),
        scv.tl.transition_matrix(aligned_obsp, **kwargs).toarray(),
    )


@pytest.mark.parametrize("use_negative_cosines", [False, True])
def test_transition_matrix_explicit_vgraph_ignores_stored_negative_graphs(
    velocity_graphs, use_negative_cosines
):
    positive, negative = velocity_graphs
    adata = _adata(positive, negative, obsp=True)
    expected = scv.tl.transition_matrix(
        adata,
        vgraph=positive,
        self_transitions=False,
        use_negative_cosines=use_negative_cosines,
    ).toarray()
    adata.uns["velocity_graph_neg"] = negative * 0.5
    adata.obsp["velocity_graph_neg"] = negative * 0.25

    observed = scv.tl.transition_matrix(
        adata,
        vgraph=positive,
        self_transitions=False,
        use_negative_cosines=use_negative_cosines,
    ).toarray()

    np.testing.assert_allclose(observed, expected)
