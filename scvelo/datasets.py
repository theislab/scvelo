"""Builtin Datasets.
"""

from scanpy.api import read
import numpy as np


def dentategyrus():
    """Dentate Gyrus dataset from Hochgerner et al. (2018).
    Dentate gyrus is part of the hippocampus involved in learning, episodic memory formation and spatial coding.
    It is measured using 10X Genomics Chromium and described in Hochgerner et al. (2018).
    The data consists of 25,919 genes across 3,396 cells and provides several interesting characteristics.
    """
    filename = 'data/DentateGyrus/10X43_1.loom'
    url = 'http://pklab.med.harvard.edu/velocyto/DG1/10X43_1.loom'
    adata = read(filename, backup_url=url, cleanup=True, sparse=True, cache=True)
    return adata


def toy_data(n_obs):
    """Random samples from Dentate Gyrus.
    """
    adata = dentategyrus()
    indices = np.random.choice(adata.n_obs, n_obs)
    adata = adata[indices]
    adata.obs_names = np.array(range(adata.n_obs), dtype='str')
    return adata
