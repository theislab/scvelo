import numpy as np
from scipy.sparse import csr_matrix, issparse, linalg

from scvelo import logging as logg
from scvelo import settings
from scvelo.preprocessing.moments import get_connectivities
from scvelo.preprocessing.neighbors import verify_neighbors
from .transition_matrix import transition_matrix
from .utils import get_plasticity_score, groups_to_bool, scale, strings_to_categoricals
from .velocity_graph import VelocityGraph


def cell_fate(
    data,
    groupby="clusters",
    disconnected_groups=None,
    self_transitions=False,
    n_neighbors=None,
    copy=False,
):
    """Computes individual cell endpoints.

    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    groupby: `str` (default: `'clusters'`)
        Key to which to assign the fates.
    disconnected_groups: list of `str` (default: `None`)
        Which groups to treat as disconnected for fate assignment.
    self_transitions: `bool` (default: `False`)
        Whether to include self-transitions.
    n_neighbors: `int` (default: `None`)
        Number of neighbors to restrict transitions to.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    cell_fate: `.obs`
        most likely cell fate for each individual cell
    cell_fate_confidence: `.obs`
        confidence of transitioning to the assigned fate
    """
    adata = data.copy() if copy else data
    logg.info("computing cell fates", r=True)

    n_neighbors = 10 if n_neighbors is None else n_neighbors
    _adata = adata.copy()
    vgraph = VelocityGraph(
        _adata, n_neighbors=n_neighbors, approx=True, n_recurse_neighbors=1
    )
    vgraph.compute_cosines()
    _adata.uns["velocity_graph"] = vgraph.graph
    _adata.uns["velocity_graph_neg"] = vgraph.graph_neg

    T = transition_matrix(_adata, self_transitions=self_transitions)
    fate = np.linalg.inv(np.eye(_adata.n_obs) - T)
    if issparse(T):
        fate = fate.A
    cell_fates = np.array(_adata.obs[groupby][fate.argmax(1)])
    if disconnected_groups is not None:
        idx = _adata.obs[groupby].isin(disconnected_groups)
        cell_fates[idx] = _adata.obs[groupby][idx]

    adata.obs["cell_fate"] = cell_fates
    adata.obs["cell_fate_confidence"] = fate.max(1) / fate.sum(1)
    strings_to_categoricals(adata)

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added\n"
        "    'cell_fate', most likely cell fate (adata.obs)\n"
        "    'cell_fate_confidence', confidence of fate transition (adata.obs)"
    )
    return adata if copy else None


def cell_origin(
    data,
    groupby="clusters",
    disconnected_groups=None,
    self_transitions=False,
    n_neighbors=None,
    copy=False,
):
    """Computes individual cell root points.

    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    groupby: `str` (default: `'clusters'`)
        Key to which to assign the fates.
    disconnected_groups: list of `str` (default: `None`)
        Which groups to treat as disconnected for fate assignment.
    n_neighbors: `int` (default: `None`)
        Number of neighbors to restrict transitions to.
    self_transitions: `bool` (default: `False`)
        Whether to include self-transitions.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    cell_origin: `.obs`
        most likely cell origin for each individual cell
    cell_origin_confidence: `.obs`
        confidence of coming from assigned origin
    """
    adata = data.copy() if copy else data
    logg.info("computing cell fates", r=True)

    n_neighbors = 10 if n_neighbors is None else n_neighbors
    _adata = adata.copy()
    vgraph = VelocityGraph(
        _adata, n_neighbors=n_neighbors, approx=True, n_recurse_neighbors=1
    )
    vgraph.compute_cosines()
    _adata.uns["velocity_graph"] = vgraph.graph
    _adata.uns["velocity_graph_neg"] = vgraph.graph_neg

    T = transition_matrix(_adata, self_transitions=self_transitions, backward=True)
    fate = np.linalg.inv(np.eye(_adata.n_obs) - T)
    if issparse(T):
        fate = fate.A
    cell_fates = np.array(_adata.obs[groupby][fate.argmax(1)])
    if disconnected_groups is not None:
        idx = _adata.obs[groupby].isin(disconnected_groups)
        cell_fates[idx] = _adata.obs[groupby][idx]

    adata.obs["cell_origin"] = cell_fates
    adata.obs["cell_origin_confidence"] = fate.max(1) / fate.sum(1)
    strings_to_categoricals(adata)

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added\n"
        "    'cell_origin', most likely cell origin (adata.obs)\n"
        "    'cell_origin_confidence', confidence of assigned origin (adata.obs)"
    )


# TODO: Add docstrings
def eigs(T, k=10, eps=1e-3, perc=None, random_state=None, v0=None):
    """TODO."""
    if random_state is not None:
        np.random.seed(random_state)
        v0 = np.random.rand(min(T.shape))
    # TODO: Find proper exception and handle appropriately
    try:
        # find k eigs with largest real part, and sort in descending order of eigenvals
        eigvals, eigvecs = linalg.eigs(T.T, k=k, which="LR", v0=v0)
        p = np.argsort(eigvals)[::-1]
        eigvals = eigvals.real[p]
        eigvecs = eigvecs.real[:, p]

        # select eigenvectors with eigenvalue of 1 - eps.
        idx = eigvals >= 1 - eps
        eigvals = eigvals[idx]
        eigvecs = np.absolute(eigvecs[:, idx])

        if perc is not None:
            lbs, ubs = np.percentile(eigvecs, perc, axis=0)
            eigvecs[eigvecs < lbs] = 0
            eigvecs = np.clip(eigvecs, 0, ubs)
            eigvecs /= eigvecs.max(0)

    except ValueError as e:
        logg.warn(f"Failed to fine `k=`{k} egenvalues with real part: {e}")
        eigvals, eigvecs = np.empty(0), np.zeros(shape=(T.shape[0], 0))

    return eigvals, eigvecs


# TODO: Add docstrings
def verify_roots(adata, roots, modality="Ms"):
    """TODO."""
    if "gene_count_corr" in adata.var.keys():
        p = get_plasticity_score(adata, modality)
        p_ub, root_ub = p > 0.5, roots > 0.9
        n_right_assignments = np.sum(root_ub * p_ub) / np.sum(p_ub)
        n_false_assignments = np.sum(root_ub * np.invert(p_ub)) / np.sum(
            np.invert(p_ub)
        )
        n_randn_assignments = np.mean(root_ub)
        if n_right_assignments > 3 * n_randn_assignments:  # mu + 2*mu (std=mu)
            roots *= p_ub
        elif (
            n_false_assignments > n_randn_assignments
            or n_right_assignments < n_randn_assignments
        ):
            logg.warn("Uncertain or fuzzy root cell identification. Please verify.")
    return roots


# TODO: Add docstrings
def write_to_obs(adata, key, vals, cell_subset=None):
    """TODO."""
    if cell_subset is None:
        adata.obs[key] = vals
    else:
        vals_all = (
            adata.obs[key].copy() if key in adata.obs.keys() else np.zeros(adata.n_obs)
        )
        vals_all[cell_subset] = vals
        adata.obs[key] = vals_all


def terminal_states(
    data,
    vkey="velocity",
    modality="Ms",
    groupby=None,
    groups=None,
    self_transitions=False,
    eps=1e-3,
    random_state=0,
    copy=False,
    **kwargs,
):
    r"""Computes terminal states (root and end points).

    The end points and root cells are obtained as stationary states of the
    velocity-inferred transition matrix and its transposed, respectively,
    which is given by left eigenvectors corresponding to an eigenvalue of 1, i.e.

    .. math::
        μ^{\textrm{end}}=μ^{\textrm{end}} \pi, \quad
        μ^{\textrm{root}}=μ^{\textrm{root}} \pi^{\small \textrm{T}}.

    .. code:: python

        scv.tl.terminal_states(adata)
        scv.pl.scatter(adata, color=["root_cells", "end_points"])

    .. image:: https://user-images.githubusercontent.com/31883718/69496183-bcfdf300-0ecf-11ea-9aae-685300a0b1ba.png

    Alternatively, we recommend to use :func:`cellrank.tl.terminal_states`
    providing an improved/generalized approach of identifying terminal states.

    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    modality: `str` (default: `'Ms'`)
        Layer used to calculate terminal states.
    groupby: `str`, `list` or `np.ndarray` (default: `None`)
        Key of observations grouping to consider. Only to be set, if each group is
        assumed to have a distinct lineage with an independent root and end point.
    groups: `str`, `list` or `np.ndarray` (default: `None`)
        Groups selected to find terminal states on. Must be an element of .obs[groupby].
        To be specified only for very distinct/disconnected clusters.
    self_transitions: `bool` (default: `False`)
        Allow transitions from one node to itself.
    eps: `float` (default: 1e-3)
        Tolerance for eigenvalue selection.
    random_state: `int` or None (default: 0)
        Seed used by the random number generator.
        If `None`, use the `RandomState` instance by `np.random`.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.
    **kwargs:
        Passed to scvelo.tl.transition_matrix(), e.g. basis, weight_diffusion.

    Returns
    -------
    root_cells: `.obs`
        sparse matrix with transition probabilities.
    end_points: `.obs`
        sparse matrix with transition probabilities.
    """
    adata = data.copy() if copy else data
    verify_neighbors(adata)

    logg.info("computing terminal states", r=True)

    strings_to_categoricals(adata)
    if groupby is not None:
        logg.warn(
            "Only set groupby, when you have evident distinct clusters/lineages,"
            " each with an own root and end point."
        )

    kwargs.update({"self_transitions": self_transitions})
    categories = [None]
    if groupby is not None and groups is None:
        categories = adata.obs[groupby].cat.categories
    for cat in categories:
        groups = cat if cat is not None else groups
        cell_subset = groups_to_bool(adata, groups=groups, groupby=groupby)
        _adata = adata if groups is None else adata[cell_subset]
        connectivities = get_connectivities(_adata, "distances")

        T = transition_matrix(_adata, vkey=vkey, backward=True, **kwargs)
        eigvecs_roots = eigs(T, eps=eps, perc=[2, 98], random_state=random_state)[1]
        roots = csr_matrix.dot(connectivities, eigvecs_roots).sum(1)
        roots = scale(np.clip(roots, 0, np.percentile(roots, 98)))
        roots = verify_roots(_adata, roots, modality)
        write_to_obs(adata, "root_cells", roots, cell_subset)

        T = transition_matrix(_adata, vkey=vkey, backward=False, **kwargs)
        eigvecs_ends = eigs(T, eps=eps, perc=[2, 98], random_state=random_state)[1]
        ends = csr_matrix.dot(connectivities, eigvecs_ends).sum(1)
        ends = scale(np.clip(ends, 0, np.percentile(ends, 98)))
        write_to_obs(adata, "end_points", ends, cell_subset)

        n_roots, n_ends = eigvecs_roots.shape[1], eigvecs_ends.shape[1]
        groups_str = f" ({groups})" if isinstance(groups, str) else ""
        roots_str = f"{n_roots} {'regions' if n_roots > 1 else 'region'}"
        ends_str = f"{n_ends} {'regions' if n_ends > 1 else 'region'}"

        logg.info(
            f"    identified {roots_str} of root cells "
            f"and {ends_str} of end points {groups_str}."
        )

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added\n"
        "    'root_cells', root cells of Markov diffusion process (adata.obs)\n"
        "    'end_points', end points of Markov diffusion process (adata.obs)"
    )
    return adata if copy else None
