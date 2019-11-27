from .. import logging as logg
from .utils import prod_sum_var, norm, get_indices, random_subsample
from .transition_matrix import transition_matrix

import numpy as np


def velocity_confidence(data, vkey='velocity', copy=False):
    """Computes confidences of velocities.

    .. code:: python

        scv.tl.velocity_confidence(adata)
        scv.pl.scatter(adata, color='velocity_confidence', perc=[2,98])

    .. image:: https://user-images.githubusercontent.com/31883718/69626334-b6df5200-1048-11ea-9171-495845c5bc7a.png
       :width: 600px


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

    V = np.array(adata.layers[vkey])
    if vkey + '_genes' in adata.var.keys():
        V = V[:, np.array(adata.var[vkey + '_genes'], dtype=bool)]

    nans = np.isnan(np.sum(V, axis=0))
    if np.any(nans):
        V = V[:, ~nans]

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

    if vkey + '_genes' in adata.var.keys():
        idx = np.array(adata.var[vkey + '_genes'], dtype=bool)
        X, V = adata.layers['Ms'][:, idx].copy(), adata.layers[vkey][:, idx].copy()
    else:
        X, V = adata.layers['Ms'].copy(), adata.layers[vkey].copy()

    nans = np.isnan(np.sum(V, axis=0))
    if np.any(nans):
        X = X[:, ~nans]
        V = V[:, ~nans]

    T = transition_matrix(adata, vkey=vkey, scale=scale)
    dX = T.dot(X) - X
    dX -= dX.mean(1)[:, None]
    V -= V.mean(1)[:, None]

    norms = norm(dX) * norm(V)
    norms += norms == 0

    adata.obs[vkey + '_confidence_transition'] = prod_sum_var(dX, V) / norms

    logg.hint('added \'' + vkey + '_confidence_transition\' (adata.obs)')

    return adata if copy else None


def score_robustness(data, adata_subset=None, fraction=.5, vkey='velocity', copy=False):
    adata = data.copy() if copy else data

    if adata_subset is None:
        from ..preprocessing.moments import moments
        from ..preprocessing.neighbors import neighbors
        from .velocity import velocity

        logg.switch_verbosity('off')
        adata_subset = adata.copy()
        subset = random_subsample(adata_subset, fraction=fraction, return_subset=True)
        neighbors(adata_subset)
        moments(adata_subset)
        velocity(adata_subset, vkey=vkey)
        logg.switch_verbosity('on')
    else:
        subset = adata.obs_names.isin(adata_subset.obs_names)

    V = adata[subset].layers[vkey]
    V_subset = adata_subset.layers[vkey]

    score = np.nan * (subset == False)
    score[subset] = prod_sum_var(V, V_subset) / (norm(V) * norm(V_subset))
    adata.obs[vkey + '_score_robustness'] = score

    return adata_subset if copy else None
