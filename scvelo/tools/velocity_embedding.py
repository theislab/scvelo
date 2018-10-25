from .. import settings
from .. import logging as logg
from .utils import norm
from .transition_matrix import transition_matrix

from scipy.sparse import issparse
import numpy as np
import warnings


def quiver_autoscale(X, Y, U, V):
    import matplotlib.pyplot as pl
    Q = pl.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=None)
    Q._init()
    pl.clf()
    return Q.scale


def velocity_embedding(data, basis=None, vkey='velocity', scale=10, self_transitions=True, use_negative_cosines=True,
                       retain_scale=False, copy=False):
    """Computes the single cell velocities in the embedding

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'tsne'`)
        Which embedding to use.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    scale: `int` (default: 10)
        Scale parameter of gaussian kernel for transition matrix.
    retain_scale: `bool` (default: `False`)
        Whether to retain scale from high dimensional space in embedding.

    Returns
    -------
    Returns or updates `adata` with the attributes
    velocity_umap: `.obsm`
        coordinates of velocity projection on embedding
    """
    adata = data.copy() if copy else data
    T = transition_matrix(adata, vkey=vkey, scale=scale, self_transitions=self_transitions, use_negative_cosines=use_negative_cosines)
    T.setdiag(0)
    T.eliminate_zeros()

    logg.info('computing velocity embedding', r=True)

    if basis is None:
        keys = [key for key in ['pca', 'tsne', 'umap'] if 'X_' + key in adata.obsm.keys()]
        if len(keys) > 0: basis = keys[-1]
        else: raise ValueError('No basis specified')

    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError('You need compute the embedding first.')

    X_emb = adata.obsm['X_' + basis][:, :2]
    V_emb = np.zeros((adata.n_obs, 2))

    TA = T.A if adata.n_obs < 8192 else None
    sparse = False if adata.n_obs < 8192 else True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(adata.n_obs):
            indices = T[i].indices
            dX = X_emb[indices] - X_emb[i, None]  # shape (n_neighbors, 2)
            if not retain_scale: dX /= norm(dX)[:, None]
            dX[np.isnan(dX)] = 0  # zero diff in a steady-state
            probs = T[i].data if sparse else TA[i, indices]
            V_emb[i] = probs.dot(dX) - probs.mean() * dX.sum(0)  # probs.sum() / len(indices)

    if retain_scale:
        delta = T.dot(adata.X) - adata.X
        if issparse(delta): delta = delta.A
        cos_proj = (adata.layers[vkey] * delta).sum(1) / norm(delta)
        V_emb *= np.clip(cos_proj[:, None] * 10, 0, 1)

    V_emb /= (3 * quiver_autoscale(X_emb[:, 0], X_emb[:, 1], V_emb[:, 0], V_emb[:, 1]))

    vkey += '_' + basis
    adata.obsm[vkey] = V_emb

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added\n'
        '    \'' + vkey + '\', embedded velocity vectors (adata.obsm)')

    return adata if copy else None

