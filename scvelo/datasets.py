"""Builtin Datasets.
"""
from .read_load import read, load
from .preprocessing.utils import cleanup
import numpy as np
import pandas as pd


def toy_data(n_obs):
    """
    Randomly samples from the Dentate Gyrus dataset.

    Arguments
    ---------
    n_obs: `int`
        Size of the sampled dataset

    Returns
    -------
    Returns `adata` object
    """

    """Random samples from Dentate Gyrus.
    """
    adata = dentategyrus()
    indices = np.random.choice(adata.n_obs, n_obs)
    adata = adata[indices]
    adata.obs_names = np.array(range(adata.n_obs), dtype='str')
    adata.var_names_make_unique()
    return adata


def dentategyrus():
    """Dentate Gyrus dataset from Hochgerner et al. (2018).

    Dentate gyrus is part of the hippocampus involved in learning, episodic memory formation and spatial coding.
    It is measured using 10X Genomics Chromium and described in Hochgerner et al. (2018).
    The data consists of 25,919 genes across 3,396 cells and provides several interesting characteristics.

    Returns
    -------
    Returns `adata` object
    """
    filename = 'data/DentateGyrus/10X43_1.loom'
    url = 'http://pklab.med.harvard.edu/velocyto/DG1/10X43_1.loom'
    adata = read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
    cleanup(adata, clean='all', keep={'spliced', 'unspliced', 'ambiguous'})

    url_louvain = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/DG_clusters.npy'
    url_umap = 'https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/DG_umap.npy'

    adata.obs['clusters'] = load('./data/DentateGyrus/DG_clusters.npy', url_louvain)
    adata.obsm['X_umap'] = load('./data/DentateGyrus/DG_umap.npy', url_umap)

    adata.obs['clusters'] = pd.Categorical(adata.obs['clusters'])

    return adata


def forebrain():
    """Developing human forebrain.
    Forebrain tissue of a week 10 embryo, focusing on the glutamatergic neuronal lineage.

    Returns
    -------
    Returns `adata` object
    """
    filename = 'data/ForebrainGlut/hgForebrainGlut.loom'
    url = 'http://pklab.med.harvard.edu/velocyto/hgForebrainGlut/hgForebrainGlut.loom'
    adata = read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
    adata.var_names_make_unique()
    return adata
