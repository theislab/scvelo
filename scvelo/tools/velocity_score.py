from ..logging import logg
from .utils import prod_sum_var, norm
from ..preprocessing.moments import moments
from .transition_matrix import transition_matrix

from scanpy.api.pp import neighbors
import numpy as np


def score_transition(adata, vkey='velocity', scale=10, copy=False):
    # compute velocity score as in correlation of mean transition with velocity
    if vkey not in adata.layers.keys():
        raise ValueError(
            'You need to run `tl.velocity` first.')

    dX = transition_matrix(adata, vkey, scale).multiply(adata.layers['Ms']) - adata.layers['Ms']
    dX -= dX.mean(1)[:, None]

    V = adata.layers[vkey].copy()
    V -= V.mean(1)[:, None]

    adata.obs[vkey+'_score_transition'] = prod_sum_var(dX, V) / (norm(dX) * norm(V))

    logg.hint(
        'added to `.obs`\n'
        '    ' + vkey + '_score_transition')

    return adata if copy else None


def score_smoothness(adata, vkey='velocity', copy=False):
    if vkey not in adata.layers.keys():
        raise ValueError(
            'You need to run `tl.velocity` first.')

    X, V = adata.layers['Ms'].copy(), adata.layers[vkey].copy()

    n_neighbors = adata.uns['neighbors']['params']['n_neighbors'] - 1
    indices = adata.uns['neighbors']['distances'].indices.reshape((-1, n_neighbors))

    V -= V.mean(1)[:, None]
    V_norm = norm(V)
    R = np.zeros(adata.n_obs)

    for i in range(adata.n_obs):
        Vi_neighs = V[indices[i]]
        Vi_neighs -= Vi_neighs.mean(1)[:, None]
        R[i] = np.mean(np.einsum('ij, j', Vi_neighs, V[i]) / (norm(Vi_neighs) * V_norm[i])[None, :])

    adata.obs[vkey + '_score_smoothness'] = R

    logg.hint(
        'added to `.obs`\n'
        '    ' + vkey + '_score_smoothness')

    return adata if copy else None


def random_subsample(adata, frac=.5):
    subset = np.random.choice([True, False], size=adata.n_obs, p=[frac, 1-frac]).sum()
    adata.obs['subset'] = subset

    adata_subset = adata[subset].copy()
    neighbors(adata_subset)
    moments(adata_subset)

    return adata_subset


def score_robustness(adata, adata_subset, vkey='velocity', copy=False):
    V = adata[adata.obs['subset']].layers[vkey]
    V_subset = adata_subset.layers[vkey]

    adata_subset.obs[vkey + '_score_robustness'] = prod_sum_var(V, V_subset) / (norm(V) * norm(V_subset))

    return adata_subset if copy else None
