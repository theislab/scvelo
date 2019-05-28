import scvelo as scv
import numpy as np


def test_einsum():
    from scvelo.tools.utils import prod_sum_obs, prod_sum_var, norm
    Ms, Mu = np.random.rand(5, 4), np.random.rand(5, 4)
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_pipeline():
    adata = scv.datasets.simulation(n_obs=100)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata)
