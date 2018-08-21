import scvelo as scv
import scanpy.api as sc
from scvelo.tools.utils import *
from scipy import sparse, stats


def test_data(n_obs=1000, n_vars=100):
    adata = sc.AnnData(sparse.random(n_obs, n_vars, data_rvs=stats.poisson(1).rvs, density=.6, format='csr'))
    adata.layers['spliced'] = adata.X
    adata.layers['unspliced'] = \
       .5 * adata.X + .3 * sparse.random(n_obs, n_vars, density=.6, data_rvs=stats.poisson(1).rvs, format='csr')

    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.filter_genes_dispersion(adata)
    sc.pp.normalize_per_cell(adata, layers='all')
    sc.pp.pca(adata, n_comps=10)
    sc.pp.neighbors(adata, use_rep='X_pca')
    scv.pp.moments(adata)
    return adata


def test_einsum():
    adata = test_data()
    Ms, Mu = adata.layers['Ms'], adata.layers['Mu']
    assert np.allclose(prod_sum_obs(Ms, Mu), np.sum(Ms * Mu, 0))
    assert np.allclose(prod_sum_var(Ms, Mu), np.sum(Ms * Mu, 1))
    assert np.allclose(norm(Ms), np.linalg.norm(Ms, axis=1))


def test_velocity_graph():
    adata = test_data()
    scv.tl.velocity_graph(adata)

    graph1 = adata.uns['velocity_graph'].copy()

    scv.tl.velocity_graph(adata, n_jobs=2)
    graph2 = adata.uns['velocity_graph'].copy()

    assert np.allclose((scv.tl.transition_matrix(adata) > 0).toarray(), (graph1 > 0).toarray())
    assert np.allclose(graph1.toarray(), graph2.toarray())