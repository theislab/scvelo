import numpy as np
from scipy.sparse import issparse, linalg, spdiags

from scanpy.tools._dpt import DPT

from scvelo import logging as logg
from scvelo.preprocessing.moments import get_connectivities
from .terminal_states import terminal_states
from .utils import groups_to_bool, scale, strings_to_categoricals


def principal_curve(data, basis="pca", n_comps=4, clusters_list=None, copy=False):
    """Computes the principal curve.

    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'pca'`)
        Basis to use for computing the principal curve.
    n_comps: `int` (default: 4)
        Number of pricipal components to be used.
    copy: `bool`, (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    principal_curve: `.uns`
        dictionary containing `projections`, `ixsort` and `arclength`
    """
    adata = data.copy() if copy else data
    from rpy2 import robjects
    from rpy2.robjects.packages import importr

    if clusters_list is not None:
        cell_subset = np.array(
            [label in clusters_list for label in adata.obs["clusters"]]
        )
        X_emb = adata[cell_subset].obsm[f"X_{basis}"][:, :n_comps]
    else:
        cell_subset = None
        X_emb = adata.obsm[f"X_{basis}"][:, :n_comps]

    n_obs, n_dim = X_emb.shape

    # convert array to R matrix
    xvec = robjects.FloatVector(X_emb.T.reshape(X_emb.size))
    X_R = robjects.r.matrix(xvec, nrow=n_obs, ncol=n_dim)

    fit = importr("princurve").principal_curve(X_R)

    adata.uns["principal_curve"] = {}
    adata.uns["principal_curve"]["ixsort"] = ixsort = np.array(fit[1]) - 1
    adata.uns["principal_curve"]["projections"] = np.array(fit[0])[ixsort]
    adata.uns["principal_curve"]["arclength"] = np.array(fit[2])
    adata.uns["principal_curve"]["cell_subset"] = cell_subset

    return adata if copy else None


# TODO: Add docstrings
def velocity_map(adata=None, T=None, n_dcs=10, return_model=False):
    """TODO."""
    vpt = VPT(adata, n_dcs=n_dcs)
    if T is None:
        T = adata.uns["velocity_graph"] - adata.uns["velocity_graph_neg"]
        vpt._connectivities = T + T.T
    vpt.compute_transitions()
    vpt.compute_eigen(n_dcs)
    adata.obsm["X_vmap"] = vpt.eigen_basis
    return vpt if return_model else None


# TODO: Add docstrings
class VPT(DPT):
    """TODO."""

    # TODO: Add docstrings
    def set_iroot(self, root=None):
        """TODO."""
        if (
            isinstance(root, str)
            and root in self._adata.obs.keys()
            and self._adata.obs[root].max() != 0
        ):
            self.iroot = get_connectivities(self._adata).dot(self._adata.obs[root])
            self.iroot = scale(self.iroot).argmax()
        elif isinstance(root, str) and root in self._adata.obs_names:
            self.iroot = self._adata.obs_names.get_loc(root)
        elif isinstance(root, (int, np.integer)) and root < self._adata.n_obs:
            self.iroot = root
        else:
            self.iroot = None

    # TODO: Add docstrings
    def compute_transitions(self, density_normalize=True):
        """TODO."""
        T = self._connectivities
        if density_normalize:
            q = np.asarray(T.sum(axis=0))
            q += q == 0
            Q = (
                spdiags(1.0 / q, 0, T.shape[0], T.shape[0])
                if issparse(T)
                else np.diag(1.0 / q)
            )
            K = Q.dot(T).dot(Q)
        else:
            K = T
        z = np.sqrt(np.asarray(K.sum(axis=0)))
        Z = (
            spdiags(1.0 / z, 0, K.shape[0], K.shape[0])
            if issparse(K)
            else np.diag(1.0 / z)
        )
        self._transitions_sym = Z.dot(K).dot(Z)

    # TODO: Add docstrings
    def compute_eigen(self, n_comps=10, sym=None, sort="decrease"):
        """TODO."""
        if self._transitions_sym is None:
            raise ValueError("Run `.compute_transitions` first.")
        n_comps = min(self._transitions_sym.shape[0] - 1, n_comps)
        evals, evecs = linalg.eigsh(self._transitions_sym, k=n_comps, which="LM")
        self._eigen_values = evals[::-1]
        self._eigen_basis = evecs[:, ::-1]

    # TODO: Add docstrings
    def compute_pseudotime(self, inverse=False):
        """TODO."""
        if self.iroot is not None:
            self._set_pseudotime()
            self.pseudotime = 1 - self.pseudotime if inverse else self.pseudotime
            self.pseudotime[~np.isfinite(self.pseudotime)] = np.nan
        else:
            self.pseudotime = np.empty(self._adata.n_obs)
            self.pseudotime[:] = np.nan


def velocity_pseudotime(
    adata,
    vkey="velocity",
    groupby=None,
    groups=None,
    root_key=None,
    end_key=None,
    n_dcs=10,
    use_velocity_graph=True,
    save_diffmap=None,
    return_model=None,
    **kwargs,
):
    """Computes a pseudotime based on the velocity graph.

    Velocity pseudotime is a random-walk based distance measures on the velocity graph.
    After computing a distribution over root cells obtained from the velocity-inferred
    transition matrix, it measures the average number of steps it takes to reach a cell
    after start walking from one of the root cells. Contrarily to diffusion pseudotime,
    it implicitly infers the root cells and is based on the directed velocity graph
    instead of the similarity-based diffusion kernel.

    .. code:: python

        scv.tl.velocity_pseudotime(adata)
        scv.pl.scatter(adata, color="velocity_pseudotime", color_map="gnuplot")

    .. image:: https://user-images.githubusercontent.com/31883718/69545487-33fbc000-0f92-11ea-969b-194dc68400b0.png
       :width: 600px

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    groupby: `str`, `list` or `np.ndarray` (default: `None`)
        Key of observations grouping to consider.
    groups: `str`, `list` or `np.ndarray` (default: `None`)
        Groups selected to find terminal states on. Must be an element of
        adata.obs[groupby]. Only to be set, if each group is assumed to have a distinct
        lineage with an independent root and end point.
    root_key: `int` (default: `None`)
        Index of root cell to be used.
        Computed from velocity-inferred transition matrix if not specified.
    end_key: `int` (default: `None`)
        Index of end point to be used.
        Computed from velocity-inferred transition matrix if not specified.
    n_dcs: `int` (default: 10)
        The number of diffusion components to use.
    use_velocity_graph: `bool` (default: `True`)
        Whether to use the velocity graph.
        If False, it uses the similarity-based diffusion kernel.
    save_diffmap: `bool` (default: `None`)
        Whether to store diffmap coordinates.
    return_model: `bool` (default: `None`)
        Whether to return the vpt object for further inspection.
    **kwargs:
        Further arguments to pass to VPT (e.g. min_group_size, allow_kendall_tau_shift).

    Returns
    -------
    velocity_pseudotime: `.obs`
        Velocity pseudotime obtained from velocity graph.
    """
    strings_to_categoricals(adata)
    if root_key is None and "root_cells" in adata.obs.keys():
        root0 = adata.obs["root_cells"][0]
        if not np.isnan(root0) and not isinstance(root0, str):
            root_key = "root_cells"
    if end_key is None and "end_points" in adata.obs.keys():
        end0 = adata.obs["end_points"][0]
        if not np.isnan(end0) and not isinstance(end0, str):
            end_key = "end_points"

    groupby = (
        "cell_fate" if groupby is None and "cell_fate" in adata.obs.keys() else groupby
    )
    if groupby is not None:
        logg.warn(
            "Only set groupby, when you have evident distinct clusters/lineages,"
            " each with an own root and end point."
        )
    categories = (
        adata.obs[groupby].cat.categories
        if groupby is not None and groups is None
        else [None]
    )
    for cat in categories:
        groups = cat if cat is not None else groups
        if (
            root_key is None
            or root_key in adata.obs.keys()
            and np.max(adata.obs[root_key]) == np.min(adata.obs[root_key])
        ):
            terminal_states(adata, vkey=vkey, groupby=groupby, groups=groups)
            root_key, end_key = "root_cells", "end_points"
        cell_subset = groups_to_bool(adata, groups=groups, groupby=groupby)
        data = adata.copy() if cell_subset is None else adata[cell_subset].copy()
        if "allow_kendall_tau_shift" not in kwargs:
            kwargs["allow_kendall_tau_shift"] = True
        vpt = VPT(data, n_dcs=n_dcs, **kwargs)

        if use_velocity_graph:
            T = data.uns[f"{vkey}_graph"] - data.uns[f"{vkey}_graph_neg"]
            vpt._connectivities = T + T.T

        vpt.compute_transitions()
        vpt.compute_eigen(n_comps=n_dcs)

        vpt.set_iroot(root_key)
        vpt.compute_pseudotime()
        dpt_root = vpt.pseudotime

        if end_key is not None:
            vpt.set_iroot(end_key)
            vpt.compute_pseudotime(inverse=True)
            dpt_end = vpt.pseudotime

            # merge dpt_root and inverse dpt_end together
            vpt.pseudotime = np.nan_to_num(dpt_root) + np.nan_to_num(dpt_end)
            vpt.pseudotime[np.isfinite(dpt_root) & np.isfinite(dpt_end)] /= 2
            vpt.pseudotime = scale(vpt.pseudotime)
            vpt.pseudotime[np.isnan(dpt_root) & np.isnan(dpt_end)] = np.nan

        if "n_branchings" in kwargs and kwargs["n_branchings"] > 0:
            vpt.branchings_segments()
        else:
            vpt.indices = vpt.pseudotime.argsort()

        if f"{vkey}_pseudotime" not in adata.obs.keys():
            pseudotime = np.empty(adata.n_obs)
            pseudotime[:] = np.nan
        else:
            pseudotime = adata.obs[f"{vkey}_pseudotime"].values
        pseudotime[cell_subset] = vpt.pseudotime
        adata.obs[f"{vkey}_pseudotime"] = np.array(pseudotime, dtype=np.float64)

        if save_diffmap:
            diffmap = np.empty(shape=(adata.n_obs, n_dcs))
            diffmap[:] = np.nan
            diffmap[cell_subset] = vpt.eigen_basis
            adata.obsm[f"X_diffmap_{groups}"] = diffmap

    return vpt if return_model else None
