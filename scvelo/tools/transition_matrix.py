from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
import numpy as np
import warnings


def transition_matrix(adata, vkey='velocity', basis=None, direction="forward", scale=10, locality = 2):
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
    
    #%% check inputs
    if vkey+'_graph' not in adata.uns:
        raise ValueError(
            'You need to run `tl.velocity_graph` first to compute cosine correlations.')
        
    if direction not in ['forward', 'backward']:
        raise ValueError(
            'You need to choose a valid direction')
        
    if basis != None and 'X_' + str(basis) not in adata.obsm.keys():
        raise ValueError(
                'You need to chose a valid basis')
    #%% contruct the desired transition matrix

    # get the veclocity based transition matrix
    T = np.expm1(adata.uns[vkey + '_graph'] * scale)

    # get the distance based transition matrix
    dists = adata.uns['neighbors']['distances']
    
    # combine the two
    T = 0.5 * T + .5 * (dists > 0).multiply(T)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T = T.multiply(csr_matrix(1. / T.sum(1)))
        
    if direction == 'backward': 
        T = T.T
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    # compute distances in the embedding to smooth the matrix
    if 'X_' + str(basis) in adata.obsm.keys():
        X_emb = adata.obsm['X_' + basis]
        dists_emb = (T > 0).multiply(squareform(pdist(X_emb)))
        kernel = np.expm1(-dists_emb**2 / dists_emb.data.mean()**2 / 2 * 1/locality)
        T = T.multiply(kernel)
        
        # must renormalise
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    return T
