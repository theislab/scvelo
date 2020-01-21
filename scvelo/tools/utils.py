from scipy.sparse import csr_matrix, issparse
import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter("ignore")

from ..preprocessing.neighbors import compute_connectivities_umap


def mean(x, axis=0):
    return x.mean(axis).A1 if issparse(x) else x.mean(axis)


def make_dense(X):
    XA = X.A if issparse(X) and X.ndim == 2 else X.A1 if issparse(X) else X
    if XA.ndim == 2:
        XA = XA[0] if XA.shape[0] == 1 else XA[:, 0] if XA.shape[1] == 1 else XA
    return np.array(XA)


def sum_obs(A):
    """summation over axis 0 (obs) equivalent to np.sum(A, 0)
    """
    return A.sum(0).A1 if issparse(A) else np.einsum('ij -> j', A)


def sum_var(A):
    """summation over axis 1 (var) equivalent to np.sum(A, 1)
    """
    return A.sum(1).A1 if issparse(A) else np.sum(A, axis=1)


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


def get_indices(dist, n_neighbors=None, mode_neighbors='distances'):
    D = dist.copy()
    D.data += 1e-6

    n_counts = sum_var(D > 0)
    n_neighbors = n_counts.min() if n_neighbors is None else min(n_counts.min(), n_neighbors)
    rows = np.where(n_counts > n_neighbors)[0]
    cumsum_neighs = np.insert(n_counts.cumsum(), 0, 0)
    dat = D.data
    for row in rows:
        n0, n1 = cumsum_neighs[row], cumsum_neighs[row + 1]
        rm_idx = n0 + dat[n0:n1].argsort()[n_neighbors:]
        dat[rm_idx] = 0
    D.eliminate_zeros()

    D.data -= 1e-6
    if mode_neighbors == 'distances':
        indices = D.indices.reshape((-1, n_neighbors))
    elif mode_neighbors == 'connectivities':
        knn_indices = D.indices.reshape((-1, n_neighbors))
        knn_distances = D.data.reshape((-1, n_neighbors))
        _ ,conn = compute_connectivities_umap(knn_indices, knn_distances, D.shape[0], n_neighbors)
        indices = get_indices_from_csr(conn)
    return indices, D


def get_indices_from_csr(conn):
    # extracts indices from connectivity matrix, pads with nans
    ixs = np.ones((conn.shape[0], np.max((conn > 0).sum(1)))) * np.nan
    for i in range(ixs.shape[0]):
        cell_indices = conn[i, :].indices
        ixs[i, :len(cell_indices)] = cell_indices
    return ixs


def get_iterative_indices(indices, index, n_recurse_neighbors=2, max_neighs=None):
    def iterate_indices(indices, index, n_recurse_neighbors):
        ix = indices[iterate_indices(indices, index, n_recurse_neighbors - 1)] if n_recurse_neighbors > 1 else indices[index]
        if np.isnan(ix).any(): ix = ix[~np.isnan(ix)]
        return ix.astype(int)

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
    lst = [item for item in lst if item is not np.nan and item is not 'nan']
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
    from pandas.api.types import is_string_dtype, is_integer_dtype, is_bool_dtype
    from pandas import Categorical

    def is_valid_dtype(values):
        return is_string_dtype(values) or is_integer_dtype(values) or is_bool_dtype(values)

    df = adata.obs
    df_keys = [key for key in df.columns if is_valid_dtype(df[key])]
    for key in df_keys:
        c = df[key]
        c = Categorical(c)
        if 1 < len(c.categories) < min(len(c), 100): df[key] = c

    df = adata.var
    df_keys = [key for key in df.columns if is_string_dtype(df[key])]
    for key in df_keys:
        c = df[key]
        c = Categorical(c)
        if 1 < len(c.categories) < min(len(c), 100): df[key] = c


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


def make_unique_list(key, allow_array=False):
    from pandas import unique, Index
    if isinstance(key, Index): key = key.tolist()
    is_list = isinstance(key, (list, tuple, np.record)) if allow_array else isinstance(key, (list, tuple, np.ndarray, np.record))
    is_list_of_str = is_list and all(isinstance(item, str) for item in key)
    return unique(key) if is_list_of_str else key if is_list and len(key) < 20 else [key]


def test_bimodality(x, bins=30, kde=True, plot=False):
    from scipy.stats import gaussian_kde, norm

    grid = np.linspace(np.min(x), np.percentile(x, 99), bins)
    kde_grid = gaussian_kde(x)(grid) if kde else np.histogram(x, bins=grid, density=True)[0]

    idx = int(bins / 2) - 2
    idx += np.argmin(kde_grid[idx: idx + 4])

    peak_0 = kde_grid[:idx].argmax()
    peak_1 = kde_grid[idx:].argmax()
    kde_peak = kde_grid[idx:][peak_1]  # min(kde_grid[:idx][peak_0], kde_grid[idx:][peak_1])
    kde_mid = kde_grid[idx:].mean()  # kde_grid[idx]

    t_stat = (kde_peak - kde_mid) / np.clip(np.std(kde_grid) / np.sqrt(bins), 1, None)
    p_val = norm.sf(t_stat)

    grid_0 = grid[:idx]
    grid_1 = grid[idx:]
    means = [(grid_0[peak_0] + grid_0[min(peak_0 +1, len(grid_0) -1)]) / 2,
             (grid_1[peak_1] + grid_1[min(peak_1 +1, len(grid_1) -1)]) / 2]

    if plot:
        color = 'grey'
        if kde:
            pl.plot(grid, kde_grid, color=color)
            pl.fill_between(grid, 0, kde_grid, alpha=.4, color=color)
        else:
            pl.hist(x, bins=grid, alpha=.4, density=True, color=color);
        pl.axvline(means[0], color=color)
        pl.axvline(means[1], color=color)
        pl.axhline(kde_mid, alpha=.2, linestyle='--', color=color)
        pl.show()

    return t_stat, p_val, means   # ~ t_test (reject unimodality if t_stat > 3)


def random_subsample(adata, fraction=.1, return_subset=False, copy=False):
    adata_sub = adata.copy() if copy else adata
    p, n = fraction, adata.n_obs
    subset = np.random.choice([True, False], size=adata.n_obs, p=[p, 1 - p])
    adata_sub._inplace_subset_obs(subset)

    return adata_sub if copy else subset if return_subset else None


def get_duplicates(array):
    from collections import Counter
    return np.array([item for (item, count) in Counter(array).items() if count > 1])


def corrcoef(x, y, mode='pearsons'):
    from scipy.stats import pearsonr, spearmanr
    corr, _ = spearmanr(x, y) if mode is 'spearmans' else pearsonr(x, y)
    return corr


def vcorrcoef(X, y):
    Xm = np.array(X - (np.nanmean(X, -1)[:, None] if X.ndim > 1 else np.nanmean(X, -1)))
    ym = np.array(y - (np.nanmean(y, -1)[:, None] if y.ndim > 1 else np.nanmean(y, -1)))
    corr = np.nansum(Xm * ym, -1) / np.sqrt(np.nansum(Xm ** 2, -1) * np.nansum(ym ** 2, -1))
    return corr


def isin(x, y):
    return np.array(pd.DataFrame(x).isin(y)).flatten()


def indices_to_bool(indices, n):
    return isin(np.arange(n), indices)


def convolve(adata, x):
    from ..preprocessing.neighbors import get_connectivities
    conn = get_connectivities(adata)
    if isinstance(x, str) and x in adata.layers.keys():
        x = adata.layers[x]
    idx_valid = ~np.isnan(x.sum(0))
    Y = np.ones(x.shape) * np.nan
    Y[:, idx_valid] = conn.dot(x[:, idx_valid])
    return Y
