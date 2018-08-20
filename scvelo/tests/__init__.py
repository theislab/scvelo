from ..tools.velocity_graph import *
from ..tools.solver import *
from ..tools.velocity_embedding import transition_matrix

from anndata import AnnData
from scanpy.api.pp import normalize_per_cell, pca, neighbors
from time import time
from scipy import sparse, stats


def test_data(n_obs=10000, n_vars=100):
    adata = AnnData(sparse.random(n_obs, n_vars, data_rvs=stats.poisson(1).rvs, density=.6, format='csr'))
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = \
        .5 * adata.X + .3 * sparse.random(n_obs, n_vars, density=.6, data_rvs=stats.poisson(1).rvs, format='csr')

    normalize_per_cell(adata, layers='all')
    pca(adata)
    neighbors(adata, use_rep='X_pca')
    moments(adata)
    return adata


def test_einsum():
    adata = test_data()
    Ms, Mu = adata.layers['Ms'], adata.layers['Mu']
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))

    start = time()

    prod_sum_obs(Ms, Mu)
    stop1 = time()

    np.sum(Ms * Mu, 0)
    stop2 = time()

    print('einsum product is ', round((stop2 - start) / (stop1 - start), 2), ' times faster than numpy.')

    start = time()

    norm(Ms)
    stop1 = time()

    np.linalg.norm(Ms, axis=1)
    stop2 = time()

    print('einsum norm is ', round((stop2 - start) / (stop1 - start), 2), ' times faster than numpy.')


# def test_solvers():
#     adata = test_data()
#     Ms, Mu = adata.layers['Ms'], adata.layers['Mu']
#
#     assert np.allclose(solve_cov(Ms, Mu), solve_inv(Ms, Mu))


def test_velocity_graph():
    adata = test_data()
    velocity(adata)
    velocity_graph(adata)

    graph1 = adata.uns['velocity_graph'].copy()

    velocity_graph(adata, n_jobs=2)
    graph2 = adata.uns['velocity_graph'].copy()

    assert np.allclose((transition_matrix(adata) > 0).toarray(), (graph1 > 0).toarray())
    assert np.allclose(graph1.toarray(), graph2.toarray())