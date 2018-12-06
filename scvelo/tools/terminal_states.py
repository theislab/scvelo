from .. import settings
from .. import logging as logg
from ..preprocessing.moments import get_connectivities
from .velocity_graph import VelocityGraph
from .transition_matrix import transition_matrix
from .utils import scale, groups_to_bool, strings_to_categoricals

from scipy.sparse import linalg, csr_matrix, issparse
import numpy as np


def cell_fate(data, groupby='clusters', disconnected_groups=None, n_neighbors=None, copy=False):
    """Computes individual cell endpoints

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    groupby: `str` (default: `'clusters'`)
        Key to which to assign the fates.
    disconnected_groups: list of `str` (default: `None`)
        Which groups to treat as disconnected for fate assignment.
    n_neighbors: `int` (default: `None`)
        Number of neighbors to restrict transitions to.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Returns or updates `adata` with the attributes
    cell_fate: `.obs`
        most likely cell fate for each individual cell
    cell_fate_confidence: `.obs`
        confidence of transitioning to the assigned fate
    """
    adata = data.copy() if copy else data
    logg.info('computing cell fates', r=True)

    n_neighbors = 10 if n_neighbors is None else n_neighbors
    _adata = adata.copy()
    vgraph = VelocityGraph(_adata, n_neighbors=n_neighbors, approx=True, n_recurse_neighbors=1)
    vgraph.compute_cosines()
    _adata.uns['velocity_graph'] = vgraph.graph
    _adata.uns['velocity_graph_neg'] = vgraph.graph_neg

    T = transition_matrix(_adata)
    I = np.eye(_adata.n_obs)
    fate = np.linalg.inv(I - T)
    if issparse(T): fate = fate.A
    cell_fates = np.array(_adata.obs[groupby][fate.argmax(1)])
    if disconnected_groups is not None:
        idx = _adata.obs[groupby].isin(disconnected_groups)
        cell_fates[idx] = _adata.obs[groupby][idx]
    adata.obs['cell_fate'] = cell_fates
    adata.obs['cell_fate_confidence'] = fate.max(1) / fate.sum(1)
    strings_to_categoricals(adata)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added\n'
        '    \'cell_fate\', most likely cell fate (adata.obs)\n'
        '    \'cell_fate_confidence\', confidence of transitioning to the assigned fate (adata.obs)')
    return adata if copy else None


def cell_origin(data, groupby='clusters', distinct_clusters=None, self_transitions=False, n_neighbors=None):
    n_neighbors = 10 if n_neighbors is None else n_neighbors
    adata = data.copy()
    vgraph = VelocityGraph(adata, n_neighbors=n_neighbors, approx=True, n_recurse_neighbors=1)
    vgraph.compute_cosines()
    adata.uns['velocity_graph'] = vgraph.graph
    adata.uns['velocity_graph_neg'] = vgraph.graph_neg

    T = transition_matrix(adata, self_transitions=self_transitions, backward=True)
    I = np.eye(data.n_obs)
    fate = np.linalg.inv(I - T)
    if issparse(T): fate = fate.A
    cell_fates = np.array(data.obs[groupby][fate.argmax(1)])
    if distinct_clusters is not None:
        idx = data.obs[groupby].isin(distinct_clusters)
        cell_fates[idx] = data.obs[groupby][idx]
    data.obs['cell_origin'] = cell_fates
    data.obs['cell_origin_confidence'] = fate.max(1) / fate.sum(1)
    strings_to_categoricals(data)


def eigs(T, k=10, eps=1e-3, perc=None):
    try:
        eigvals, eigvecs = linalg.eigs(T.T, k=k, which='LR')  # find k eigs with largest real part

        p = np.argsort(eigvals)[::-1]                        # sort in descending order of eigenvalues
        eigvals = eigvals.real[p]
        eigvecs = eigvecs.real[:, p]

        idx = (eigvals >= 1 - eps)                           # select eigenvectors with eigenvalue of 1
        eigvals = eigvals[idx]
        eigvecs = np.absolute(eigvecs[:, idx])

        if perc is not None:
            lbs, ubs = np.percentile(eigvecs, perc, axis=0)
            eigvecs[eigvecs < lbs] = 0
            eigvecs = np.clip(eigvecs, 0, ubs)
            eigvecs /= eigvecs.max(0)

    except:
        eigvals, eigvecs = np.empty(0), np.zeros(shape=(T.shape[0], 0))

    return eigvals, eigvecs


def write_to_obs(adata, key, vals, cell_subset=None):
    if cell_subset is None:
        adata.obs[key] = vals
    else:
        vals_all = adata.obs[key].copy() if key in adata.obs.keys() else np.zeros(adata.n_obs)
        vals_all[cell_subset] = vals
        adata.obs[key] = vals_all


def terminal_states(data, vkey='velocity', groupby=None, groups=None, self_transitions=False, basis=None,
                    weight_diffusion=0, scale_diffusion=1, eps=1e-3, copy=False):
    """Computes terminal states (root and end points).

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    self_transitions: `bool` (default: `False`)
        Allow transitions from one node to itself.
    basis: `str` (default: `None`)
        Basis to use.
    weight_diffusion: `float` (default: 0)
        Relative weight to be given to diffusion kernel (Brownian motion)
    scale_diffusion: `float` (default: 1)
        Scale of diffusion kernel.
    eps: `float` (default: 1e-3)
        Tolerance for eigenvalue selection.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.

    Returns
    -------
    Returns or updates `data` with the attributes
    root: `.obs`
        sparse matrix with transition probabilities.
    end: `.obs`
        sparse matrix with transition probabilities.
    """
    adata = data.copy() if copy else data
    logg.info('computing terminal states', r=True)

    strings_to_categoricals(adata)

    groupby = 'cell_fate' if groupby is None and 'cell_fate' in adata.obs.keys() else groupby
    categories = adata.obs[groupby].cat.categories if groupby is not None and groups is None else [None]
    for cat in categories:
        groups = cat if cat is not None else groups
        cell_subset = groups_to_bool(adata, groups=groups, groupby=groupby)
        _adata = adata if groups is None else adata[cell_subset]
        connectivities = get_connectivities(_adata, 'distances')

        T = transition_matrix(_adata, vkey=vkey, basis=basis, weight_diffusion=weight_diffusion,
                              scale_diffusion=scale_diffusion, self_transitions=self_transitions, backward=True)
        eigvecs_roots = eigs(T, eps=eps, perc=[2, 98])[1]
        roots = csr_matrix.dot(connectivities, eigvecs_roots).sum(1)
        roots = scale(np.clip(roots, 0, np.percentile(roots, 98)))
        write_to_obs(adata, 'root_cells', roots, cell_subset)

        T = transition_matrix(_adata, vkey=vkey, basis=basis, weight_diffusion=weight_diffusion,
                              scale_diffusion=scale_diffusion, self_transitions=self_transitions, backward=False)
        eigvecs_ends = eigs(T, eps=eps, perc=[2, 98])[1]
        ends = csr_matrix.dot(connectivities, eigvecs_ends).sum(1)
        ends = scale(np.clip(ends, 0, np.percentile(ends, 98)))
        write_to_obs(adata, 'end_points', ends, cell_subset)

        n_roots, n_ends = eigvecs_roots.shape[1], eigvecs_ends.shape[1]
        groups_str = ' (' + groups + ')' if isinstance(groups, str) else ''
        logg.info('    identified ' + str(n_roots) + ' root cells and ' + str(n_ends) + ' end points' + groups_str)

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added\n'
        '    \'root_cells\', root cells of Markov diffusion process (adata.obs)\n'
        '    \'end_points\', end points of Markov diffusion process (adata.obs)')
    return adata if copy else None
