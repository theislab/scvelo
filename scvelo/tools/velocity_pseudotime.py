import numpy as np
from scanpy.tools.dpt import DPT
from scipy.sparse import issparse, spdiags, linalg

from .utils import groups_to_bool, scale, strings_to_categoricals
from ..preprocessing.moments import get_connectivities


def principal_curve(data, basis='pca', n_comps=4, clusters_list=None, copy=False):
    """Computes the principal curve
    
    Arguments
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
    Returns or updates `adata` with the attributes
    principal_curve: `.uns`
        dictionary containing `projections`, `ixsort` and `arclength`
    """
    adata = data.copy() if copy else data
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr

    if clusters_list is not None:
        cell_subset = np.array([label in clusters_list for label in adata.obs['clusters']])
        X_emb = adata[cell_subset].obsm['X_' + basis][:, :n_comps]
    else:
        cell_subset = None
        X_emb = adata.obsm['X_' + basis][:, :n_comps]

    n_obs, n_dim = X_emb.shape

    # convert array to R matrix
    xvec = robjects.FloatVector(X_emb.T.reshape((X_emb.size)))
    X_R = robjects.r.matrix(xvec, nrow=n_obs, ncol=n_dim)

    fit = importr("princurve").principal_curve(X_R)

    adata.uns['principal_curve'] = dict()
    adata.uns['principal_curve']['ixsort'] = ixsort = np.array(fit[1])-1
    adata.uns['principal_curve']['projections'] = np.array(fit[0])[ixsort]
    adata.uns['principal_curve']['arclength'] = np.array(fit[2])
    adata.uns['principal_curve']['cell_subset'] = cell_subset

    return adata if copy else None


def velocity_map(adata=None, T=None, n_dcs=10, return_model=False):
    vpt = VPT(adata, n_dcs=n_dcs)
    if T is None:
        T = adata.uns['velocity_graph'] - adata.uns['velocity_graph_neg']
        vpt._connectivities = T + T.T
    vpt.compute_transitions()
    vpt.compute_eigen(n_dcs)
    adata.obsm['X_vmap'] = vpt.eigen_basis
    return vpt if return_model else None


class VPT(DPT):
    def set_iroot(self, root=None):
        if isinstance(root, str) and root in self._adata.obs.keys() and self._adata.obs[root].max() != 0:
            self.iroot = scale(get_connectivities(self._adata).dot(self._adata.obs[root])).argmax()
        elif isinstance(root, str) and root in self._adata.obs_names:
            self.iroot = np.where(self._adata.obs_names == root)[0][0]
        elif isinstance(root, (int, np.integer)) and root < self._adata.n_obs:
            self.iroot = root
        else:
            self.iroot = None

    def compute_transitions(self, density_normalize=True):
        T = self._connectivities
        if density_normalize:
            q = np.asarray(T.sum(axis=0))
            Q = spdiags(1.0 / q, 0, T.shape[0], T.shape[0]) if issparse(T) else np.diag(1.0 / q)
            K = Q.dot(T).dot(Q)
        else:
            K = T
        z = np.sqrt(np.asarray(K.sum(axis=0)))
        Z = spdiags(1.0 / z, 0, K.shape[0], K.shape[0]) if issparse(K) else np.diag(1.0 / z)
        self._transitions_sym = Z.dot(K).dot(Z)

    def compute_eigen(self, n_comps=10, sym=None, sort='decrease'):
        if self._transitions_sym is None:
            raise ValueError('Run `.compute_transitions` first.')
        n_comps = min(self._transitions_sym.shape[0] - 1, n_comps)
        evals, evecs = linalg.eigsh(self._transitions_sym, k=n_comps, which='LM')
        self._eigen_values = evals[::-1]
        self._eigen_basis = evecs[:, ::-1]

    def compute_pseudotime(self, inverse=False):
        if self.iroot is not None:
            self._set_pseudotime()
            self.pseudotime = 1 - self.pseudotime if inverse else self.pseudotime
            self.pseudotime[~np.isfinite(self.pseudotime)] = np.nan
        else:
            self.pseudotime = np.empty(self._adata.n_obs)
            self.pseudotime[:] = np.nan


def velocity_pseudotime(adata, groupby=None, groups=None, root=None, end=None, n_dcs=10, n_branchings=0,
                        min_group_size=0.01, allow_kendall_tau_shift=True, use_velocity_field=True,
                        save_diffmap=False, return_model=False):
    strings_to_categoricals(adata)
    root = 'root_cells' if root is None and 'root_cells' in adata.obs.keys() else root
    end = 'end_points' if end is None and 'end_points' in adata.obs.keys() else end
    groupby = 'cell_fate' if groupby is None and 'cell_fate' in adata.obs.keys() else groupby
    categories = adata.obs[groupby].cat.categories if groupby is not None and groups is None else [None]
    for cat in categories:
        groups = cat if cat is not None else groups
        cell_subset = groups_to_bool(adata, groups=groups, groupby=groupby)
        data = adata.copy() if cell_subset is None else adata[cell_subset].copy()
        vpt = VPT(data, n_dcs=n_dcs, min_group_size=min_group_size,
                  n_branchings=n_branchings, allow_kendall_tau_shift=allow_kendall_tau_shift)

        if use_velocity_field:
            T = data.uns['velocity_graph'] - data.uns['velocity_graph_neg']
            vpt._connectivities = T + T.T

        vpt.compute_transitions()
        vpt.compute_eigen(n_comps=n_dcs)

        vpt.set_iroot(root)
        vpt.compute_pseudotime()
        dpt_root = vpt.pseudotime

        vpt.set_iroot(end)
        vpt.compute_pseudotime(inverse=True)
        dpt_end = vpt.pseudotime

        # merge dpt_root and inverse dpt_end together
        vpt.pseudotime = np.nan_to_num(dpt_root) + np.nan_to_num(dpt_end)
        vpt.pseudotime[np.isfinite(dpt_root) & np.isfinite(dpt_end)] /= 2
        vpt.pseudotime = scale(vpt.pseudotime)
        vpt.pseudotime[np.isnan(dpt_root) & np.isnan(dpt_end)] = np.nan

        if n_branchings > 0: vpt.branchings_segments()
        else: vpt.indices = vpt.pseudotime.argsort()

        if 'velocity_pseudotime' not in adata.obs.keys():
            pseudotime = np.empty(adata.n_obs)
            pseudotime[:] = np.nan
        else:
            pseudotime = adata.obs['velocity_pseudotime'].copy()
        pseudotime[cell_subset] = vpt.pseudotime
        adata.obs['velocity_pseudotime'] = pseudotime

        if save_diffmap:
            diffmap = np.empty(shape=(adata.n_obs, n_dcs))
            diffmap[:] = np.nan
            diffmap[cell_subset] = vpt.eigen_basis
            adata.obsm['X_diffmap_' + groups] = diffmap

    return vpt if return_model else None
