from .preprocessing.utils import show_proportions, cleanup, set_initial_size, get_initial_size
from .preprocessing.moments import get_connectivities, second_order_moments

from .tools.utils import prod_sum_obs, prod_sum_var, norm, vector_norm, R_squared, \
    cosine_correlation, normalize, scale, get_indices, get_iterative_indices, groups_to_bool, randomized_velocity
from .tools.rank_velocity_genes import get_mean_var
from .tools.run import convert_to_adata, convert_to_loom
from .tools.optimization import leastsq_NxN, leastsq_generalized, maximum_likelihood
from .tools.velocity_confidence import random_subsample
from .tools.velocity_graph import vals_to_csr

from .plotting.utils import is_categorical, clip, interpret_colorkey
from .plotting.velocity_embedding_grid import compute_velocity_on_grid

from .read_load import clean_obs_names, merge

from scanpy.utils import merge_groups