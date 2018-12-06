from ..preprocessing.neighbors import get_connectivities
from .utils import normalize

from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
import numpy as np


def transition_matrix(adata, vkey='velocity', basis=None, backward=False, self_transitions=True, scale=10, perc=None,
                      use_negative_cosines=False, weight_diffusion=0, scale_diffusion=1, weight_indirect_neighbors=.5,
                      n_neighbors=None, vgraph=None):
    """Computes transition probabilities from velocity graph

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

    graph = csr_matrix(adata.uns[vkey + '_graph']).copy() if vgraph is None else vgraph.copy()

    if self_transitions:
        confidence = graph.max(1).A.flatten()
        ub = np.percentile(confidence, 98)
        self_prob = np.clip(ub - confidence, 0, 1)
        graph.setdiag(self_prob)

    T = np.expm1(graph * scale)  # equivalent to np.exp(adata.uns[vkey + '_graph'].A * scale) - 1
    if vkey + '_graph_neg' in adata.uns.keys():
        T = T - np.expm1(-adata.uns[vkey + '_graph_neg'] * scale) if use_negative_cosines \
            else T + (adata.uns[vkey + '_graph_neg'] < 0) * 1e-10

    # weight direct neighbors with 1 and indirect (recurse) neighbors with 0.5
    if 'neighbors' in adata.uns.keys():
        direct_neighbors = adata.uns['neighbors']['distances'] > 0
        direct_neighbors.setdiag(1)
        w = weight_indirect_neighbors
        T = w * T + (1-w) * direct_neighbors.multiply(T)

    if backward: T = T.T
    T = normalize(T)

    if n_neighbors is not None:
        T = T.multiply(get_connectivities(adata, mode='distances', n_neighbors=n_neighbors, recurse_neighbors=True))

    if perc is not None:
        threshold = np.percentile(T.data, perc)
        T.data[T.data < threshold] = 0
        T.eliminate_zeros()

    if 'X_' + str(basis) in adata.obsm.keys():
        dists_emb = (T > 0).multiply(squareform(pdist(adata.obsm['X_' + basis])))
        scale_diffusion *= dists_emb.data.mean()
        
        diffusion_kernel = dists_emb.copy()
        diffusion_kernel.data = np.exp(-.5 * dists_emb.data ** 2 / scale_diffusion ** 2)
        T = T.multiply(diffusion_kernel)  # combine velocity based kernel with diffusion based kernel

        if 0 < weight_diffusion < 1:  # add another diffusion kernel (Brownian motion - like)
            diffusion_kernel.data = np.exp(-.5 * dists_emb.data ** 2 / (scale_diffusion/2) ** 2)
            T = (1-weight_diffusion) * T + weight_diffusion * diffusion_kernel

        T = normalize(T)

    return T