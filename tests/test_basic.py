import scvelo as scv
import numpy as np


def test_einsum():
    from scvelo.tools.utils import prod_sum_obs, prod_sum_var, norm
    Ms, Mu = np.random.rand(5, 4), np.random.rand(5, 4)
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_neighbors():
    adata = scv.datasets.simulation(random_seed=0, n_vars=100)
    scv.pp.filter_and_normalize(adata)

    scv.pp.pca(adata)
    scv.pp.neighbors(adata)
    adata_ = scv.pp.neighbors(adata, method='sklearn', copy=True)
    assert np.all(np.round(adata.obsp['distances'][0].data, 2) == np.round(adata_.obsp['distances'][0].data, 2))


def test_dynamical_model():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.recover_dynamics(adata, var_names=adata.var_names[0])
    assert np.round(adata[:, adata.var_names[0]].var['fit_alpha'][0], 4) == 4.7409


def test_pipeline():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)

    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata)
    adata.var.velocity_genes = True

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata)

    scv.tl.velocity_confidence(adata)
    scv.tl.latent_time(adata)

    with scv.GridSpec() as pl:
        pl.velocity_graph(adata)
        pl.velocity_embedding(adata, arrow_length=3, arrow_size=3, c='latent_time')
        pl.velocity_embedding_grid(adata, scale=.5, density=.5, c='latent_time', cmap='gnuplot')
        pl.velocity_embedding_stream(adata, c=adata.var_names[0], layer='velocity')
        pl.scatter(adata, basis=adata.var_names[0], c='velocity', use_raw=True)
        pl.hist([adata.obs.initial_size_spliced, adata.obs.initial_size_unspliced])
