from .preprocessing.utils import show_proportions, cleanup, set_initial_size, get_initial_size
from .preprocessing.moments import get_connectivities, second_order_moments, second_order_moments_u

from .tools.utils import prod_sum_obs, prod_sum_var, norm, vector_norm, R_squared, cosine_correlation, normalize, \
    scale, get_indices, get_iterative_indices, groups_to_bool, randomized_velocity, extract_int_from_str
from .tools.rank_velocity_genes import get_mean_var
from .tools.run import convert_to_adata, convert_to_loom
from .tools.optimization import leastsq_NxN, leastsq_generalized, maximum_likelihood
from .tools.velocity_confidence import random_subsample
from .tools.velocity_graph import vals_to_csr

from .plotting.utils import is_categorical, clip, interpret_colorkey
from .plotting.velocity_embedding_grid import compute_velocity_on_grid

from .read_load import clean_obs_names, merge


def merge_groups(adata, key, map_groups, key_added=None, map_colors=None):
    adata._sanitize()
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
