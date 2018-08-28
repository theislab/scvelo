from scipy.sparse import csr_matrix
import numpy as np
import warnings


def transition_matrix(adata, vkey='velocity', scale=10):
    """Computes transition probabilities by applying Gaussian kernel to cosine similarities x scale

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    return T


def transition_matrix_markov(adata, direction="forward"):
    """Prepare a transition probability for Markov process
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    direction: `str` (default: `'forward'`
        Whether to diffuse forward of backwards
    Returns
    -------
    Returns sparse matrix with Markov transition probabilities.
    """
    T = transition_matrix(adata) if direction == 'forward' else transition_matrix(adata).T
    T.data[T.data < 1e-10] = 0

    dists = adata.uns['neighbors']['distances']
    T = .5 * T + .5 * (dists > 0).multiply(T)  # + .2 * np.expm1(dists * dists.mean() / 2)

    T = T.multiply(csr_matrix(1. / T.sum(1)))

    return T
