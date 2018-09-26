import scvelo as scv
import numpy as np


def test_einsum():
    from scvelo.tools.utils import prod_sum_obs, prod_sum_var, norm
    Ms, Mu = np.random.rand(5, 4), np.random.rand(5, 4)
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


# def test_velocity_graph():
#     adata = scv.datasets.toy_data(n_obs=500)
#     scv.pp.recipe_velocity(adata, n_top_genes=300)
#     scv.tl.velocity(adata)
#
#     scv.tl.velocity_graph(adata)
#     graph = adata.uns['velocity_graph'] + adata.uns['velocity_graph_neg']
#     graph.setdiag(0)
#     graph.eliminate_zeros()
#
#     T = scv.tl.transition_matrix(adata)
#     T.setdiag(0)
#     T.eliminate_zeros()
#
#     assert np.allclose((T > 0).toarray(), (graph > 0).toarray())
