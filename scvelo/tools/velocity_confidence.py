from ..logging import logg
from .utils import prod_sum_var, norm, get_indices
from ..preprocessing.moments import moments
from .transition_matrix import transition_matrix

from scanpy.api.pp import neighbors
import numpy as np


def velocity_confidence(adata, vkey='velocity', copy=False):
    if vkey not in adata.layers.keys():
        raise ValueError(
            'You need to run `tl.velocity` first.')

    idx = adata.var['velocity_genes']
    X, V = adata.layers['Ms'][:, idx].copy(), adata.layers[vkey][:, idx].copy()
    indices = get_indices(dist=adata.uns['neighbors']['distances'])[0]

    V -= V.mean(1)[:, None]
    V_norm = norm(V)
    R = np.zeros(adata.n_obs)

    for i in range(adata.n_obs):
        Vi_neighs = V[indices[i]]
        Vi_neighs -= Vi_neighs.mean(1)[:, None]
        R[i] = np.mean(np.einsum('ij, j', Vi_neighs, V[i]) / (norm(Vi_neighs) * V_norm[i])[None, :])

    adata.obs[vkey + '_length'] = V_norm.round(2)
    adata.obs[vkey + '_confidence'] = R.round(2)

    logg.hint(
        'added to `.obs`\n'
        '    ' + vkey + '_confidence')

    return adata if copy else None


def velocity_confidence_transition(adata, vkey='velocity', scale=10, copy=False):
    if vkey not in adata.layers.keys():
        raise ValueError(
            'You need to run `tl.velocity` first.')

    idx = adata.var['velocity_genes']
    T = transition_matrix(adata, vkey=vkey, scale=scale)
    dX = T.dot(adata.layers['Ms'][:, idx]) - adata.layers['Ms'][:, idx]
    dX -= dX.mean(1)[:, None]

    V = adata.layers[vkey][:, idx].copy()
    V -= V.mean(1)[:, None]

    adata.obs[vkey + '_confidence_transition'] = prod_sum_var(dX, V) / (norm(dX) * norm(V))

    logg.hint(
        'added to `.obs`\n'
        '    ' + vkey + '_confidence_transition')

    return adata if copy else None


def random_subsample(adata, frac=.5):
    subset = np.random.choice([True, False], size=adata.n_obs, p=[frac, 1-frac]).sum()
    adata.obs['subset'] = subset

    adata_subset = adata[subset].copy()
    neighbors(adata_subset)
    moments(adata_subset)

    return adata_subset


def score_robustness(adata, adata_subset=None, vkey='velocity', copy=False):
    if adata_subset is None: adata_subset = random_subsample(adata)
    V = adata[adata.obs['subset']].layers[vkey]
    V_subset = adata_subset.layers[vkey]

    adata_subset.obs[vkey + '_score_robustness'] = prod_sum_var(V, V_subset) / (norm(V) * norm(V_subset))

    return adata_subset if copy else None
