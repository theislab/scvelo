from scipy.spatial.distance import pdist, squareform
from .utils import normalize_sparse
import numpy as np


def transition_matrix(adata, vkey='velocity', basis=None, backward=False, scale=10, weight_diffusion=0, scale_diffusion=1):
    """Computes transition probabilities by applying Gaussian kernel to cosine similarities x scale

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    basis: `str` or `None` (default: `None`)
        Restrict transition to embedding if specified
    backward: `bool` (default: `False`)
        Whether to use the transition matrix to push forward (`False`) or to pull backward (`True`)
    scale: `float` (default: 10)
        Scale parameter of gaussian kernel.
    weight_diffusion: `float` (default: 0)
        Relative weight to be given to diffusion kernel (Brownian motion)
    scale_diffusion: `float` (default: 1)
        Scale of diffusion kernel. 

    Returns
    -------
    Returns sparse matrix with transition probabilities.
    """
    if vkey+'_graph' not in adata.uns:
        raise ValueError('You need to run `tl.velocity_graph` first to compute cosine correlations.')
 
    T = np.expm1(adata.uns[vkey + '_graph'] * scale)

    # weight direct neighbors with 1 and indirect (recurse) neighbors with 0.5
    T = .5 * T + .5 * (adata.uns['neighbors']['distances'] > 0).multiply(T)

    if backward: T = T.T
    T = normalize_sparse(T)

    if 'X_' + str(basis) in adata.obsm.keys():
        dists_emb = (T > 0).multiply(squareform(pdist(adata.obsm['X_' + basis])))
        scale_diffusion *= dists_emb.data.mean()
        
        diffusion_kernel = dists_emb.copy()
        diffusion_kernel.data = np.exp(-.5 * dists_emb.data ** 2 / scale_diffusion ** 2)
        T = T.multiply(diffusion_kernel)  # combine velocity based kernel with diffusion based kernel

        if 0 < weight_diffusion < 1:  # add another diffusion kernel (Brownian motion - like)
            diffusion_kernel.data = np.exp(-.5 * dists_emb.data ** 2 / (scale_diffusion/2) ** 2)
            T = (1-weight_diffusion) * T + weight_diffusion * diffusion_kernel

        T = normalize_sparse(T)

    return T