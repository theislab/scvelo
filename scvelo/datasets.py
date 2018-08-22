"""Builtin Datasets.
"""

import scanpy.api as sc
import numpy as np


def dentategyrus():
    """Dentate Gyrus dataset from Hochgerner et al. (2018)
    Dentate gyrus is part of the hippocampus involved in learning, episodic memory formation and spatial coding.
    It is measured using 10X Genomics Chromium and described in Hochgerner et al. (2018).
    The data consists of 25,919 genes across 3,396 cells and provides several interesting characteristics.
    """
    filename = 'data/DentateGyrus/10X43_1.loom'
    url = 'http://pklab.med.harvard.edu/velocyto/DG1/10X43_1.loom'
    adata = sc.read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
    return adata


def toy_data(n_obs):
    # from scipy import sparse, stats
    # adata = sc.AnnData(sparse.random(n_obs, n_vars, data_rvs=stats.poisson(1).rvs, density=.6, format='csr'))
    # adata.layers['spliced'] = adata.X
    # adata.layers['unspliced'] = \
    #    .5 * adata.X + .3 * sparse.random(n_obs, n_vars, density=.6, data_rvs=stats.poisson(1).rvs, format='csr')
    adata = dentategyrus()
    indices = np.random.choice(adata.n_obs, n_obs)
    adata = adata[indices]
    adata.obs_names = range(adata.n_obs)
    return adata
