from scvelo.tools.utils import prod_sum_obs, prod_sum_var, norm
import numpy as np
import scvelo as scv


def test_einsum():
    adata = scv.datasets.toy_data(n_obs=100, n_vars=100)
    adata = scv.pp.recipe_velocity(adata)
    Ms, Mu = adata.layers['Ms'], adata.layers['Mu']
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_velocity_graph():
    adata = scv.datasets.toy_data(n_obs=1000, n_vars=100)
    scv.tl.velocity_graph(adata)

    graph1 = adata.uns['velocity_graph'].copy()

    scv.tl.velocity_graph(adata, n_jobs=2)
    graph2 = adata.uns['velocity_graph'].copy()

    assert np.allclose((scv.tl.transition_matrix(adata) > 0).toarray(), (graph1 > 0).toarray())
    assert np.allclose(graph1.toarray(), graph2.toarray())
