from scipy.spatial.distance import pdist, squareform
from .utils import normalize_sparse
import numpy as np


def transition_matrix(adata, vkey='velocity', basis=None, direction="forward", diffuse_prob=0, scale=10):
    """Computes transition probabilities by applying Gaussian kernel to cosine similarities x scale

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    basis: `str` or `None` (default: `None`)
        Restrict transition to embedding if specified
    direction: `'forward'` or `'backward'` (default: `'foward'`)
        Whether to use the transition matrix to push forward or to pull backward
    diffuse_prob:
        Fraction to account for diffusion kernel (Brownian motion)
    scale: `float` (default: 10)
        Scale parameter of gaussian kernel.

    Returns
    -------
    Returns sparse matrix with transition probabilities.
    """
    if vkey+'_graph' not in adata.uns:
        raise ValueError(
            'You need to run `tl.velocity_graph` first to compute cosine correlations.')

    T = np.expm1(adata.uns[vkey + '_graph'] * scale)
    T = .5 * T + .5 * (adata.uns['neighbors']['distances'] > 0).multiply(T)

    if direction == 'backward': T = T.T
    T = normalize_sparse(T)

    if basis is not None and 'X_' + basis in adata.obsm.keys():
        dists_emb = (T > 0).multiply(squareform(pdist(adata.obsm['X_' + basis])))
        sigma = dists_emb.mean()
        kernel_locality = np.expm1(-.5 * dists_emb.multiply(dists_emb) / sigma**2)
        T = T.multiply(kernel_locality)
        T = (1-diffuse_prob) * T + diffuse_prob * np.expm1(-.5 * dists_emb.multiply(dists_emb) / (sigma/2)**2)
        T = normalize_sparse(T)

    return T
