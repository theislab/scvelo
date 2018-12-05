from scipy.sparse import csr_matrix, issparse
import pandas as pd
import numpy as np
import warnings


def mean(x, axis=0):
    return x.mean(axis).A1 if issparse(x) else x.mean(axis)


def make_dense(X):
    return X.A if issparse(X) and X.ndim==2 else X.A1 if issparse(X) else X


def prod_sum_obs(A, B):
    """dot product and sum over axis 0 (obs) equivalent to np.sum(A * B, 0)
    """
    return A.multiply(B).sum(0).A1 if issparse(A) else np.einsum('ij, ij -> j', A, B)


def prod_sum_var(A, B):
    """dot product and sum over axis 1 (var) equivalent to np.sum(A * B, 1)
    """
    return A.multiply(B).sum(1).A1 if issparse(A) else np.einsum('ij, ij -> i', A, B)


def norm(A):
    """computes the L2-norm along axis 1 (e.g. genes or embedding dimensions) equivalent to np.linalg.norm(A, axis=1)
    """
    return np.sqrt(A.multiply(A).sum(1).A1) if issparse(A) else np.sqrt(np.einsum('ij, ij -> i', A, A))


def vector_norm(x):
    """computes the L2-norm along axis 1 (e.g. genes or embedding dimensions) equivalent to np.linalg.norm(A, axis=1)
    """
    return np.sqrt(np.einsum('i, i -> ', x, x))


def R_squared(residual, total):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = np.ones(residual.shape[1]) - prod_sum_obs(residual, residual) / prod_sum_obs(total, total)
    r2[np.isnan(r2)] = 0
    return r2


def cosine_correlation(dX, Vi):
    dX -= dX.mean(-1)[:, None]
    Vi_norm = vector_norm(Vi)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.zeros(dX.shape[0]) if Vi_norm == 0 else np.einsum('ij, j', dX, Vi) / (norm(dX) * Vi_norm)[None, :]
    return result


def normalize(X):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X = X.multiply(csr_matrix(1. / np.abs(X).sum(1))) if issparse(X) else X / X.sum(1)
    return X


def scale(X, min=0, max=1):
    idx = np.isfinite(X)
    if any(idx):
        X = X - X[idx].min() + min
        xmax = X[idx].max()
        X = X / xmax * max if xmax != 0 else X * max
    return X


def get_indices(dist, n_neighbors=None):
    D = dist.copy()
    n_counts = (D > 0).sum(1).A1 if issparse(D) else (D > 0).sum(1)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data

    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()
    indices = D.indices.reshape((-1, n_neighbors))
    return indices, D


def get_iterative_indices(indices, index, n_recurse_neighbors=2, max_neighs=None):
    def iterate_indices(indices, index, n_recurse_neighbors):
        return indices[iterate_indices(indices, index, n_recurse_neighbors - 1)] \
            if n_recurse_neighbors > 1 else indices[index]
    indices = np.unique(iterate_indices(indices, index, n_recurse_neighbors))
    if max_neighs is not None and len(indices) > max_neighs:
        indices = np.random.choice(indices, max_neighs, replace=False)
    return indices


def groups_to_bool(adata, groups, groupby=None):
    groups = [groups] if isinstance(groups, str) else groups
    if isinstance(groups, (list, tuple, np.ndarray, np.record)):
        groupby = groupby if groupby in adata.obs.keys() else 'clusters' if 'clusters' in adata.obs.keys() \
            else 'louvain' if 'louvain' in adata.obs.keys() else None
        if groupby is not None:
            groups = np.array([key in groups for key in adata.obs[groupby]])
        else: raise ValueError('groupby attribute not valid.')
    return groups


def most_common_in_list(lst):
    lst = list(lst)
    return max(set(lst), key=lst.count)


def randomized_velocity(adata, vkey='velocity', add_key='velocity_random'):
    V_rnd = adata.layers[vkey].copy()
    for i in range(V_rnd.shape[1]):
        np.random.shuffle(V_rnd[:, i])
        V_rnd[:, i] = V_rnd[:, i] * np.random.choice(np.array([+1, -1]), size=V_rnd.shape[0])
    adata.layers[add_key] = V_rnd

    from .velocity_graph import velocity_graph
    from .velocity_embedding import velocity_embedding

    velocity_graph(adata, vkey=add_key)
    velocity_embedding(adata, vkey=add_key, autoscale=False)


def extract_int_from_str(array):
    def str_to_int(item):
        num = "".join(filter(str.isdigit, item))
        num = int(num) if len(num) > 0 else -1
        return num
    if isinstance(array, str): nums = str_to_int(array)
    elif len(array) > 1 and isinstance(array[0], str):
        nums = []
        for item in array: nums.append(str_to_int(item))
    else: nums = array
    nums = pd.Categorical(nums) if array.dtype == 'category' else np.array(nums)
    return nums


def strings_to_categoricals(adata):
    """Transform string annotations to categoricals.
    """
    from pandas.api.types import is_string_dtype
    from pandas import Categorical
    for df in [adata.obs, adata.var]:
        string_cols = [key for key in df.columns if is_string_dtype(df[key])]
        for key in string_cols:
            c = df[key]
            c = Categorical(c)
            if len(c.categories) < len(c): df[key] = c


def merge_groups(adata, key, map_groups, key_added=None, map_colors=None):
    strings_to_categoricals(adata)
    if len(map_groups) != len(adata.obs[key].cat.categories):
        map_coarse = {}
        for c in adata.obs[key].cat.categories:
            for group in map_groups:
                if any(cluster == c for cluster in map_groups[group]): map_coarse[c] = group
            if c not in map_coarse: map_coarse[c] = c
        map_groups = map_coarse

    if key_added is None:
        key_added = key + '_coarse'

    from pandas.api.types import CategoricalDtype
    adata.obs[key_added] = adata.obs[key].map(map_groups).astype(CategoricalDtype())
    old_categories = adata.obs[key].cat.categories
    new_categories = adata.obs[key_added].cat.categories

    # map_colors is passed
    if map_colors is not None:
        old_colors = None
        if key + '_colors' in adata.uns:
            old_colors = adata.uns[key + '_colors']
        new_colors = []
        for group in adata.obs[key_added].cat.categories:
            if group in map_colors:
                new_colors.append(map_colors[group])
            elif group in old_categories and old_colors is not None:
                new_colors.append(old_colors[old_categories.get_loc(group)])
            else:
                raise ValueError('You didn\'t specify a color for {}.'.format(group))
        adata.uns[key_added + '_colors'] = new_colors

    # map_colors is not passed
    elif key + '_colors' in adata.uns:
        old_colors = adata.uns[key + '_colors']
        inverse_map_groups = {g: [] for g in new_categories}
        for old_group in old_categories:
            inverse_map_groups[map_groups[old_group]].append(old_group)
        new_colors = []
        for group in new_categories:
            # take the largest of the old groups
            old_group = adata.obs[key][adata.obs[key].isin(
                inverse_map_groups[group])].value_counts().index[0]
            new_colors.append(old_colors[old_categories.get_loc(old_group)])
        adata.uns[key_added + '_colors'] = new_colors


def cutoff_small_velocities(adata, vkey='velocity', key_added='velocity_cut', frac_of_max=.5, use_raw=False):
    x = adata.layers['spliced'] if use_raw else adata.layers['Ms']
    y = adata.layers['unspliced'] if use_raw else adata.layers['Mu']

    x_max = x.max(0).A[0] if issparse(x) else x.max(0)
    y_max = y.max(0).A[0] if issparse(y) else y.max(0)

    xy_norm = x / np.clip(x_max, 1e-3, None) + y / np.clip(y_max, 1e-3, None)
    W = xy_norm >= np.percentile(xy_norm, 98, axis=0) * frac_of_max

    adata.layers[key_added] = csr_matrix(W).multiply(adata.layers[vkey]).tocsr()

    from .velocity_graph import velocity_graph
    from .velocity_embedding import velocity_embedding

    velocity_graph(adata, vkey=key_added, approx=True)
    velocity_embedding(adata, vkey=key_added)

