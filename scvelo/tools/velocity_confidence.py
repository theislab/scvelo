from .. import logging as logg
from .utils import prod_sum_var, norm, get_indices
from ..preprocessing.moments import moments
from .transition_matrix import transition_matrix

from scanpy.api.pp import neighbors
import numpy as np


def velocity_confidence(data, vkey='velocity', copy=False):
    """Computes confidences of velocities.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_length: `.obs`
        Length of the velocity vectors for each individual cell
    velocity_confidence: `.obs`
        Confidence for each cell
    """
    adata = data.copy() if copy else data
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
    adata.obs[vkey + '_confidence'] = R

    logg.hint('added \'' + vkey + '_confidence\' (adata.obs)')

    if vkey + '_confidence_transition' not in adata.obs.keys():
        velocity_confidence_transition(adata, vkey)

    return adata if copy else None


def velocity_confidence_transition(data, vkey='velocity', scale=10, copy=False):
    """Computes confidences of velocity transitions.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    scale: `float` (default: 10)
        Scale parameter of gaussian kernel.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_confidence_transition: `.obs`
        Confidence of transition for each cell
    """
    adata = data.copy() if copy else data
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

    logg.hint('added \'' + vkey + '_confidence_transition\' (adata.obs)')

    return adata if copy else None


def random_subsample(adata, frac=.5):
    subset = np.random.choice([True, False], size=adata.n_obs, p=[frac, 1-frac]).sum()
    adata.obs['subset'] = subset

    adata_subset = adata[subset].copy()
    neighbors(adata_subset)
    moments(adata_subset)

    return adata_subset


def score_robustness(data, adata_subset=None, vkey='velocity', copy=False):
    adata = data.copy() if copy else data
    if adata_subset is None: adata_subset = random_subsample(adata)
    V = adata[adata.obs['subset']].layers[vkey]
    V_subset = adata_subset.layers[vkey]

    adata_subset.obs[vkey + '_score_robustness'] = prod_sum_var(V, V_subset) / (norm(V) * norm(V_subset))

    return adata_subset if copy else None
