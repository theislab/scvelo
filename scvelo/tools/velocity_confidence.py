import numpy as np

from scvelo import logging as logg
from scvelo.core import l2_norm, prod_sum
from scvelo.preprocessing.neighbors import get_neighs
from .transition_matrix import transition_matrix
from .utils import get_indices, random_subsample


def velocity_confidence(data, vkey="velocity", copy=False):
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
    velocity_length: `.obs`
        Length of the velocity vectors for each individual cell
    velocity_confidence: `.obs`
        Confidence for each cell
    """  # noqa E501

    adata = data.copy() if copy else data
    if vkey not in adata.layers.keys():
        raise ValueError("You need to run `tl.velocity` first.")

    V = np.array(adata.layers[vkey])

    tmp_filter = np.invert(np.isnan(np.sum(V, axis=0)))
    if f"{vkey}_genes" in adata.var.keys():
        tmp_filter &= np.array(adata.var[f"{vkey}_genes"], dtype=bool)
    if "spearmans_score" in adata.var.keys():
        tmp_filter &= adata.var["spearmans_score"].values > 0.1

    V = V[:, tmp_filter]

    V -= V.mean(1)[:, None]
    V_norm = l2_norm(V, axis=1)
    R = np.zeros(adata.n_obs)

    indices = get_indices(dist=get_neighs(adata, "distances"))[0]
    for i in range(adata.n_obs):
        Vi_neighs = V[indices[i]]
        Vi_neighs -= Vi_neighs.mean(1)[:, None]
        R[i] = np.mean(
            np.einsum("ij, j", Vi_neighs, V[i])
            / (l2_norm(Vi_neighs, axis=1) * V_norm[i])[None, :]
        )

    adata.obs[f"{vkey}_length"] = V_norm.round(2)
    adata.obs[f"{vkey}_confidence"] = np.clip(R, 0, None)

    logg.hint(f"added '{vkey}_length' (adata.obs)")
    logg.hint(f"added '{vkey}_confidence' (adata.obs)")

    if f"{vkey}_confidence_transition" not in adata.obs.keys():
        velocity_confidence_transition(adata, vkey)

    return adata if copy else None


def velocity_confidence_transition(data, vkey="velocity", scale=10, copy=False):
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
    velocity_confidence_transition: `.obs`
        Confidence of transition for each cell
    """

    adata = data.copy() if copy else data
    if vkey not in adata.layers.keys():
        raise ValueError("You need to run `tl.velocity` first.")

    X = np.array(adata.layers["Ms"])
    V = np.array(adata.layers[vkey])

    tmp_filter = np.invert(np.isnan(np.sum(V, axis=0)))
    if f"{vkey}_genes" in adata.var.keys():
        tmp_filter &= np.array(adata.var[f"{vkey}_genes"], dtype=bool)
    if "spearmans_score" in adata.var.keys():
        tmp_filter &= adata.var["spearmans_score"].values > 0.1

    V = V[:, tmp_filter]
    X = X[:, tmp_filter]

    T = transition_matrix(adata, vkey=vkey, scale=scale)
    dX = T.dot(X) - X
    dX -= dX.mean(1)[:, None]
    V -= V.mean(1)[:, None]

    norms = l2_norm(dX, axis=1) * l2_norm(V, axis=1)
    norms += norms == 0

    adata.obs[f"{vkey}_confidence_transition"] = prod_sum(dX, V, axis=1) / norms

    logg.hint(f"added '{vkey}_confidence_transition' (adata.obs)")

    return adata if copy else None


def score_robustness(
    data, adata_subset=None, fraction=0.5, vkey="velocity", copy=False
):
    adata = data.copy() if copy else data

    if adata_subset is None:
        from scvelo.preprocessing.moments import moments
        from scvelo.preprocessing.neighbors import neighbors
        from .velocity import velocity

        logg.switch_verbosity("off")
        adata_subset = adata.copy()
        subset = random_subsample(adata_subset, fraction=fraction, return_subset=True)
        neighbors(adata_subset)
        moments(adata_subset)
        velocity(adata_subset, vkey=vkey)
        logg.switch_verbosity("on")
    else:
        subset = adata.obs_names.isin(adata_subset.obs_names)

    V = adata[subset].layers[vkey]
    V_subset = adata_subset.layers[vkey]

    score = np.nan * (subset is False)
    score[subset] = prod_sum(V, V_subset, axis=1) / (
        l2_norm(V, axis=1) * l2_norm(V_subset, axis=1)
    )
    adata.obs[f"{vkey}_score_robustness"] = score

    return adata_subset if copy else None
