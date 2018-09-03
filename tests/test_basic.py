from scvelo.tools.utils import prod_sum_obs, prod_sum_var, norm
import numpy as np
import scvelo as scv


def test_einsum():
    adata = scv.datasets.toy_data(n_obs=100)
    scv.pp.recipe_velocity(adata, n_top_genes=300)
    Ms, Mu = adata.layers['Ms'], adata.layers['Mu']
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_velocity_graph():
    adata = scv.datasets.toy_data(n_obs=1000)
    scv.pp.recipe_velocity(adata, n_top_genes=300)
    scv.tl.velocity(adata)

    scv.tl.velocity_graph(adata)
    graph = adata.uns['velocity_graph'].copy()

    assert np.allclose((scv.tl.transition_matrix(adata) > 0).toarray(), (graph > 0).toarray())