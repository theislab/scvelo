import numpy as np

import scvelo as scv


def test_einsum():
    from scvelo.tools.utils import norm, prod_sum_obs, prod_sum_var

    Ms, Mu = np.random.rand(5, 4), np.random.rand(5, 4)
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_neighbors():
    adata = scv.datasets.simulation(random_seed=0, n_vars=100)
    scv.pp.filter_and_normalize(adata)
    scv.pp.pca(adata)
    scv.pp.neighbors(adata)
    adata_ = scv.pp.neighbors(adata, method="sklearn", copy=True)
    dists = np.round(adata.obsp["distances"][0].data, 2)
    dists_ = np.round(adata_.obsp["distances"][0].data, 2)
    assert np.all(dists == dists_)


def test_dynamical_model():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.recover_dynamics(adata, var_names=adata.var_names[0])
    assert np.round(adata[:, adata.var_names[0]].var["fit_alpha"][0], 4) == 4.7409


def test_pipeline():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)

    scv.pp.filter_and_normalize(adata, n_top_genes=5)
    scv.pp.pca(adata)
    scv.pp.moments(adata)

    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity(adata, vkey="dynamical_velocity", mode="dynamical")
    adata.var.velocity_genes = True

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata)

    scv.tl.velocity_confidence(adata)
    scv.tl.latent_time(adata)

    with scv.GridSpec() as pl:
        pl.velocity_graph(adata)
        pl.velocity_embedding(adata, arrow_length=3, arrow_size=3, c="latent_time")
        pl.velocity_embedding_grid(adata, scale=0.5, c="latent_time", cmap="gnuplot")
        pl.velocity_embedding_stream(adata, c=adata.var_names[0], layer="velocity")
        pl.scatter(adata, basis=adata.var_names[0], c="velocity", use_raw=True)
        pl.hist([adata.obs.initial_size_spliced, adata.obs.initial_size_unspliced])

    Ms, Mu = adata.layers["Ms"][0], adata.layers["Mu"][0]
    Vs, Vd = adata.layers["velocity"][0], adata.layers["dynamical_velocity"][0]
    Vgraph = adata.uns["velocity_graph"].data[:5]
    pars = adata[:, 0].var[["fit_alpha", "fit_gamma"]].values

    assert np.allclose(Ms, [0.8269, 1.0772, 0.9396, 1.0846, 1.0398], rtol=1e-2)
    assert np.allclose(Mu, [3.8412, 3.1976, 3.5523, 3.3433, 3.8006], rtol=1e-2)
    assert np.allclose(adata.X[0], [0.0, 0.0, 0.0, 0.4981, 0.0], rtol=1e-2)
    # assert np.allclose(Vpca, [0.0163, 0.0185, 0.0472, 0.0025], rtol=1e-2)
    assert np.allclose(Vd, [1.7582, 2.0214, 1.73, 0.6615, 1.5118], rtol=1e-2)
    assert np.allclose(Vs, [3.2961, 2.5254, 2.9926, 2.634, 3.1352], rtol=1e-2)
    assert np.allclose(Vgraph, [0.915, 0.5997, 0.8494, 0.1615, 0.7708], rtol=1e-2)
    assert np.allclose(pars, [4.9257, 0.3239], rtol=1e-2)


def test_highly_variable_subset():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    bdata = adata.copy()

    scv.pp.filter_and_normalize(adata, n_top_genes=5, subset_highly_variable=True)
    scv.pp.filter_and_normalize(bdata, n_top_genes=5, subset_highly_variable=False)

    scv.pp.pca(adata)
    scv.pp.pca(bdata)

    scv.pp.moments(adata, use_rep="pca")
    scv.pp.moments(bdata, use_rep="pca")

    scv.tl.velocity_graph(adata)
    scv.tl.velocity_graph(bdata)

    bdata._inplace_subset_var(bdata.var["highly_variable"])

    assert np.allclose(adata.layers["Ms"][0], bdata.layers["Ms"][0])
    assert np.allclose(adata.layers["velocity"][0], bdata.layers["velocity"][0])
    assert np.allclose(
        adata.uns["velocity_graph"].data[:5], bdata.uns["velocity_graph"].data[:5]
    )
