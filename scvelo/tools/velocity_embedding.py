import warnings

import numpy as np
from scipy.sparse import issparse

from scvelo import logging as logg
from scvelo import settings
from scvelo.core import l2_norm
from .transition_matrix import transition_matrix


# TODO: Addd docstrings
def quiver_autoscale(X_emb, V_emb):
    """TODO."""
    import matplotlib.pyplot as pl

    scale_factor = np.abs(X_emb).max()  # just so that it handles very large values
    fig, ax = pl.subplots()
    Q = ax.quiver(
        X_emb[:, 0] / scale_factor,
        X_emb[:, 1] / scale_factor,
        V_emb[:, 0],
        V_emb[:, 1],
        angles="xy",
        scale_units="xy",
        scale=None,
    )
    Q._init()
    fig.clf()
    pl.close(fig)
    return Q.scale / scale_factor


def velocity_embedding(
    data,
    basis=None,
    vkey="velocity",
    scale=10,
    self_transitions=True,
    use_negative_cosines=True,
    direct_pca_projection=None,
    retain_scale=False,
    autoscale=True,
    all_comps=True,
    T=None,
    copy=False,
):
    r"""Projects the single cell velocities into any embedding.

    Given normalized difference of the embedding positions

    .. math::
        \tilde \delta_{ij} = \frac{x_j-x_i}{\left\lVert x_j-x_i \right\rVert},

    the projections are obtained as expected displacements with respect to the
    transition matrix :math:`\tilde \pi_{ij}` as

    .. math::
        \tilde \nu_i = E_{\tilde \pi_{i\cdot}} [\tilde \delta_{i \cdot}]
        = \sum_{j \neq i} \left( \tilde \pi_{ij} - \frac1n \right) \tilde
        \delta_{ij}.


    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'tsne'`)
        Which embedding to use.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    scale: `int` (default: 10)
        Scale parameter of gaussian kernel for transition matrix.
    self_transitions: `bool` (default: `True`)
        Whether to allow self transitions, based on the confidences of transitioning to
        neighboring cells.
    use_negative_cosines: `bool` (default: `True`)
        Whether to project cell-to-cell transitions with negative cosines into
        negative/opposite direction.
    direct_pca_projection: `bool` (default: `None`)
        Whether to directly project the velocities into PCA space,
        thus skipping the velocity graph.
    retain_scale: `bool` (default: `False`)
        Whether to retain scale from high dimensional space in embedding.
    autoscale: `bool` (default: `True`)
        Whether to scale the embedded velocities by a scalar multiplier,
        which simply ensures that the arrows in the embedding are properly scaled.
    all_comps: `bool` (default: `True`)
        Whether to compute the velocities on all embedding components.
    T: `csr_matrix` (default: `None`)
        Allows the user to directly pass a transition matrix.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    velocity_umap: `.obsm`
        coordinates of velocity projection on embedding (e.g., basis='umap')
    """
    adata = data.copy() if copy else data

    if basis is None:
        keys = [
            key for key in ["pca", "tsne", "umap"] if f"X_{key}" in adata.obsm.keys()
        ]
        if len(keys) > 0:
            basis = "pca" if direct_pca_projection else keys[-1]
        else:
            raise ValueError("No basis specified")

    if f"X_{basis}" not in adata.obsm_keys():
        raise ValueError("You need to compute the embedding first.")

    if direct_pca_projection and "pca" in basis:
        logg.warn(
            "Directly projecting velocities into PCA space is for exploratory analysis "
            "on principal components.\n"
            "         It does not reflect the actual velocity field from high "
            "dimensional gene expression space.\n"
            "         To visualize velocities, consider applying "
            "`direct_pca_projection=False`.\n"
        )

    logg.info("computing velocity embedding", r=True)

    V = np.array(adata.layers[vkey])
    vgenes = np.ones(adata.n_vars, dtype=bool)
    if f"{vkey}_genes" in adata.var.keys():
        vgenes &= np.array(adata.var[f"{vkey}_genes"], dtype=bool)
    vgenes &= ~np.isnan(V.sum(0))
    V = V[:, vgenes]

    if direct_pca_projection and "pca" in basis:
        PCs = adata.varm["PCs"] if all_comps else adata.varm["PCs"][:, :2]
        PCs = PCs[vgenes]

        X_emb = adata.obsm[f"X_{basis}"]
        V_emb = (V - V.mean(0)).dot(PCs)

    else:
        X_emb = (
            adata.obsm[f"X_{basis}"] if all_comps else adata.obsm[f"X_{basis}"][:, :2]
        )
        V_emb = np.zeros(X_emb.shape)

        T = (
            transition_matrix(
                adata,
                vkey=vkey,
                scale=scale,
                self_transitions=self_transitions,
                use_negative_cosines=use_negative_cosines,
            )
            if T is None
            else T
        )
        T.setdiag(0)
        T.eliminate_zeros()

        densify = adata.n_obs < 1e4
        TA = T.A if densify else None

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(adata.n_obs):
                indices = T[i].indices
                dX = X_emb[indices] - X_emb[i, None]  # shape (n_neighbors, 2)
                if not retain_scale:
                    dX /= l2_norm(dX)[:, None]
                dX[np.isnan(dX)] = 0  # zero diff in a steady-state
                probs = TA[i, indices] if densify else T[i].data
                V_emb[i] = probs.dot(dX) - probs.mean() * dX.sum(0)

        if retain_scale:
            X = (
                adata.layers["Ms"]
                if "Ms" in adata.layers.keys()
                else adata.layers["spliced"]
            )
            delta = T.dot(X[:, vgenes]) - X[:, vgenes]
            if issparse(delta):
                delta = delta.A
            cos_proj = (V * delta).sum(1) / l2_norm(delta)
            V_emb *= np.clip(cos_proj[:, None] * 10, 0, 1)

    if autoscale:
        V_emb /= 3 * quiver_autoscale(X_emb, V_emb)

    if f"{vkey}_params" in adata.uns.keys():
        adata.uns[f"{vkey}_params"]["embeddings"] = (
            []
            if "embeddings" not in adata.uns[f"{vkey}_params"]
            else list(adata.uns[f"{vkey}_params"]["embeddings"])
        )
        adata.uns[f"{vkey}_params"]["embeddings"].extend([basis])

    vkey += f"_{basis}"
    adata.obsm[vkey] = V_emb

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint("added\n" f"    '{vkey}', embedded velocity vectors (adata.obsm)")

    return adata if copy else None
