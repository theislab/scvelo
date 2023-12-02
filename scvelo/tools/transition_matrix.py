import warnings

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
from scipy.spatial.distance import pdist, squareform

from scvelo.preprocessing.neighbors import get_connectivities, get_neighs
from .utils import normalize

warnings.simplefilter("ignore", SparseEfficiencyWarning)


def transition_matrix(
    adata,
    vkey="velocity",
    basis=None,
    backward=False,
    self_transitions=True,
    scale=10,
    perc=None,
    threshold=None,
    use_negative_cosines=False,
    weight_diffusion=0,
    scale_diffusion=1,
    weight_indirect_neighbors=None,
    n_neighbors=None,
    vgraph=None,
    basis_constraint=None,
):
    r"""Computes cell-to-cell transition probabilities.

    .. math::
        \tilde \pi_{ij} = \frac1{z_i} \exp( \pi_{ij} / \sigma),

    from the velocity graph :math:`\pi_{ij}`, with row-normalization :math:`z_i` and
    kernel width :math:`\sigma` (scale parameter :math:`\lambda = \sigma^{-1}`).

    Alternatively, use :func:`cellrank.tl.transition_matrix` to account for uncertainty
    in the velocity estimates.

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    basis: `str` or `None` (default: `None`)
        Restrict transition to embedding if specified
    backward: `bool` (default: `False`)
        Whether to use the transition matrix to
        push forward (`False`) or to pull backward (`True`)
    self_transitions: `bool` (default: `True`)
        Allow transitions from one node to itself.
    scale: `float` (default: 10)
        Scale parameter of gaussian kernel.
    perc: `float` between `0` and `100` or `None` (default: `None`)
        Determines threshold of transitions to include.
    use_negative_cosines: `bool` (default: `False`)
        If True, negatively similar transitions are taken into account.
    weight_diffusion: `float` (default: 0)
        Relative weight to be given to diffusion kernel (Brownian motion)
    scale_diffusion: `float` (default: 1)
        Scale of diffusion kernel.
    weight_indirect_neighbors: `float` between `0` and `1` or `None` (default: `None`)
        Weight to be assigned to indirect neighbors (i.e. neighbors of higher degrees).
    n_neighbors:`int` (default: None)
        Number of nearest neighbors to consider around each cell.
    vgraph: csr matrix or `None` (default: `None`)
        Velocity graph representation to use instead of adata.uns[f'{vkey}_graph'].

    Returns
    -------
    Returns sparse matrix with transition probabilities.
    """
    if f"{vkey}_graph" not in adata.uns:
        raise ValueError(
            "You need to run `tl.velocity_graph` first to compute cosine correlations."
        )

    graph_neg = None
    if vgraph is not None:
        graph = vgraph.copy()
    else:
        if hasattr(adata, "obsp") and f"{vkey}_graph" in adata.obsp.keys():
            graph = csr_matrix(adata.obsp[f"{vkey}_graph"]).copy()
            if f"{vkey}_graph_neg" in adata.obsp.keys():
                graph_neg = adata.obsp[f"{vkey}_graph_neg"]
        else:
            graph = csr_matrix(adata.uns[f"{vkey}_graph"]).copy()
            if f"{vkey}_graph_neg" in adata.uns.keys():
                graph_neg = adata.uns[f"{vkey}_graph_neg"]

    if basis_constraint is not None and f"X_{basis_constraint}" in adata.obsm.keys():
        from sklearn.neighbors import NearestNeighbors

        neighs = NearestNeighbors(n_neighbors=100)
        neighs.fit(adata.obsm[f"X_{basis_constraint}"])
        basis_graph = neighs.kneighbors_graph(mode="connectivity") > 0
        graph = graph.multiply(basis_graph)

    if self_transitions:
        confidence = graph.max(1).A.flatten()
        ub = np.percentile(confidence, 98)
        self_prob = np.clip(ub - confidence, 0, 1)
        graph.setdiag(self_prob)

    T = np.expm1(graph * scale)  # equivalent to np.exp(graph.A * scale) - 1
    if graph_neg is not None:
        graph_neg = adata.uns[f"{vkey}_graph_neg"]
        if use_negative_cosines:
            T -= np.expm1(-graph_neg * scale)
        else:
            T += np.expm1(graph_neg * scale)
            T.data += 1

    # weight direct and indirect (recursed) neighbors
    if weight_indirect_neighbors is not None and weight_indirect_neighbors < 1:
        direct_neighbors = get_neighs(adata, "distances") > 0
        direct_neighbors.setdiag(1)
        w = weight_indirect_neighbors
        T = w * T + (1 - w) * direct_neighbors.multiply(T)

    if n_neighbors is not None:
        T = T.multiply(
            get_connectivities(
                adata, mode="distances", n_neighbors=n_neighbors, recurse_neighbors=True
            )
        )

    if perc is not None or threshold is not None:
        if threshold is None:
            threshold = np.percentile(T.data, perc)
        T.data[T.data < threshold] = 0
        T.eliminate_zeros()

    if backward:
        T = T.T
    T = normalize(T)

    if f"X_{basis}" in adata.obsm.keys():
        dists_emb = (T > 0).multiply(squareform(pdist(adata.obsm[f"X_{basis}"])))
        scale_diffusion *= dists_emb.data.mean()

        diffusion_kernel = dists_emb.copy()
        diffusion_kernel.data = np.exp(
            -0.5 * dists_emb.data**2 / scale_diffusion**2
        )
        T = T.multiply(diffusion_kernel)  # combine velocity kernel & diffusion kernel

        if 0 < weight_diffusion < 1:  # add diffusion kernel (Brownian motion - like)
            diffusion_kernel.data = np.exp(
                -0.5 * dists_emb.data**2 / (scale_diffusion / 2) ** 2
            )
            T = (1 - weight_diffusion) * T + weight_diffusion * diffusion_kernel

        T = normalize(T)

    return T


def get_cell_transitions(
    adata,
    starting_cell=0,
    basis=None,
    n_steps=100,
    n_neighbors=30,
    backward=False,
    random_state=None,
    **kwargs,
):
    """Simulate cell transitions.

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    starting_cell: `int` (default: `0`)
        Index (`int`) or name (`obs_names`) of starting cell.
    n_steps: `int` (default: `100`)
        Number of transitions/steps to be simulated.
    backward: `bool` (default: `False`)
        Whether to use the transition matrix to push forward (`False`) or to pull backward (`True`)
    random_state: `int` or `None` (default: `None`)
        Set to `int` for reproducibility, otherwise `None` for a random seed.
    **kwargs:
        To be passed to tl.transition_matrix.

    Returns
    -------
    Returns embedding coordinates (if basis is specified),
    otherwise return indices of simulated cell transitions.
    """
    np.random.seed(random_state)
    if isinstance(starting_cell, str) and starting_cell in adata.obs_names:
        starting_cell = adata.obs_names.get_loc(starting_cell)
    X = [starting_cell]
    T = transition_matrix(
        adata,
        backward=backward,
        basis_constraint=basis,
        self_transitions=False,
        **kwargs,
    )
    for _ in range(n_steps):
        t = T[X[-1]]
        indices, p = t.indices, t.data
        if n_neighbors is not None and n_neighbors < len(p):
            idx = np.argsort(t.data)[::-1][:n_neighbors]
            indices, p = indices[idx], p[idx]
        if len(p) == 0:
            indices, p = [X[-1]], [1]
        p /= np.sum(p)
        ix = np.random.choice(indices, p=p)
        X.append(ix)
    X = pd.unique(X)
    if basis is not None and f"X_{basis}" in adata.obsm.keys():
        X = adata.obsm[f"X_{basis}"][X].T
    if backward:
        X = np.flip(X, axis=-1)
    return X
