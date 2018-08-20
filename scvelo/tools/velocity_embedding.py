from .velocity_graph import *


def transition_matrix(adata, vkey='velocity', scale=10):
    """Computes transition probabilities by applying Gaussian kernel to cosine similarities x scale

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix

    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used

    scale: `float` (default: 10)
        scale parameter of gaussian kernel

    Returns
    -------
    sparse matrix with transition probabilities
    """
    if vkey+'_graph' not in adata.uns:
        raise ValueError(
            'You need to run `tl.velocity_graph` first to compute cosine correlations.')

    T = np.expm1(adata.uns[vkey + '_graph'] * scale)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T = T.multiply(csr_matrix(1. / T.sum(1)))

    return T


def velocity_embedding(adata, basis='tsne', vkey='velocity', scale=10, retain_scale=False, copy=False):
    """Computes the single cell velocities in the embedding

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix

    basis: `string` (default: `'tsne'`)
        Which embedding to use

    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used

    scale: `int` (default: 10)
        scale parameter of gaussian kernel for transition matrix

    retain_scale: `bool` (default: `False`)
        Whether to retain scale from high dimensional space in embedding

    Returns
    -------
    Returns or updates `adata` with the attributes
    `'velocity_' + basis`: coordinates of quivers on embedding (`.obsm`)
    """
    T = transition_matrix(adata, vkey, scale)

    logg.info('computing velocity embedding', r=True)

    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError(
            'You need to run `tl.{}` first to compute embedded velocity vectors.').format(basis)

    X_emb = adata.obsm['X_' + basis][:, :2]
    V_emb = np.zeros((adata.n_obs, 2))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if adata.n_obs < 8192:
            TA = T.A
            for i in range(adata.n_obs):
                indices = T[i].indices
                diff = X_emb[indices] - X_emb[i, None]
                diff /= norm(diff)[:, None]
                diff[np.isnan(diff)] = 0  # zero diff in a steady-state
                V_emb[i] = TA[i, indices].dot(diff) - diff.sum(0) / len(indices)
        else:
            for i in range(adata.n_obs):
                indices = T[i].indices
                diff = X_emb[indices] - X_emb[i, None]
                diff /= norm(diff)[:, None]
                diff[np.isnan(diff)] = 0  # zero diff in a steady-state
                V_emb[i] = T[i].data.dot(diff) - diff.sum(0) / len(indices)

    if retain_scale:
        delta = T.dot(adata.X) - adata.X
        cos_proj = (adata.obsm[vkey] * delta).sum(1) / norm(delta)
        V_emb *= np.clip(cos_proj, 0, 1)

    vkey += '_' + basis
    adata.obsm[vkey] = V_emb

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.obsm`\n'
        '    \'' + vkey + '\', embedded velocity vectors')

    return adata if copy else None

