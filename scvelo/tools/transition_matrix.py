from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
import numpy as np
import warnings


def transition_matrix(adata, vkey='velocity', basis=None, direction="forward", scale=10, p = 0.5):
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
    scale: `float` (default: 10)
        Scale parameter of gaussian kernel.
    p: 'float' (default: 0.5)
        Relative weight to be given to the cosine correlations. Must be in [0, 1].

    Returns
    -------
    Returns sparse matrix with transition probabilities.
    """
    if vkey+'_graph' not in adata.uns:
        raise ValueError(
            'You need to run `tl.velocity_graph` first to compute cosine correlations.')
        
    if p< 0 or p>1:
        raise ValueError(
            'You need to choose p in [0, 1].')
        
    if direction not in ['forward', 'backward']:
        raise ValueError(
            'You need to choose a valid direction')

    T = np.expm1(adata.uns[vkey + '_graph'] * scale)

    dists = adata.uns['neighbors']['distances']
    T = p * T + (1-p) * (dists > 0).multiply(T)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    # compute distances in the embedding to smooth the matrix
    if basis in adata.obsm.keys():
        if direction == 'backward': T = T.T
        X_emb = adata.obsm['X_' + basis]
        dists_emb = (T > 0).multiply(squareform(pdist(X_emb)))
        kernel = np.expm1(-dists_emb**2 / dists_emb.data.mean()**2 / 2)
        T = T.multiply(kernel)
        
        # must renormalise
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    return T
