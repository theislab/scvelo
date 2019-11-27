from .. import settings
from .. import logging as logg
from .utils import strings_to_categoricals

from scipy.sparse import issparse
import numpy as np


def get_mean_var(X, ignore_zeros=False, perc=None):
    data = X.data if issparse(X) else X
    mask_nans = np.isnan(data) | np.isinf(data) | np.isneginf(data)

    n_nonzeros = (X != 0).sum(0)
    n_counts = n_nonzeros if ignore_zeros else X.shape[0]

    if mask_nans.sum() > 0:
        if issparse(X):
            data[np.isnan(data) | np.isinf(data) | np.isneginf(data)] = 0
            n_nans = n_nonzeros - (X != 0).sum(0)
        else:
            X[mask_nans] = 0
            n_nans = mask_nans.sum(0)
        n_counts -= n_nans

    if perc is not None:
        if np.size(perc) < 2: perc = [perc, 100] if perc < 50 else [0, perc]
        lb, ub = np.percentile(data, perc)
        data = np.clip(data, lb, ub)

    mean = (X.sum(0) / n_counts).A1 if issparse(X) else X.sum(0) / n_counts
    mean_sq = (X.multiply(X).sum(0) / n_counts).A1 if issparse(X) else np.multiply(X, X).sum(0) / n_counts
    n_cells = np.clip(X.shape[0], 2, None)  # to avoid division by zero
    var = (mean_sq - mean ** 2) * (n_cells / (n_cells - 1))

    mean = np.nan_to_num(mean)
    var = np.nan_to_num(var)
    return mean, var


def select_groups(adata, groups='all', key='louvain'):
    """Get subset of groups in adata.obs[key].
    """
    strings_to_categoricals(adata)
    if isinstance(groups, list) and isinstance(groups[0], int): groups = [str(n) for n in groups]
    categories = adata.obs[key].cat.categories
    groups_masks = np.array([categories[i] == adata.obs[key].values for i, name in enumerate(categories)])
    if groups == 'all':
        groups = categories.values
    else:
        groups_ids = [np.where(categories.values == name)[0][0] for name in groups]
        groups_masks = groups_masks[groups_ids]
        groups = categories[groups_ids].values
    return groups, groups_masks


def velocity_clusters(data, vkey='velocity', match_with='clusters', sort_by='dpt_pseudotime', resolution=None,
                      min_likelihood=None, copy=False):
    """Computes velocity clusters via louvain on velocities.

    .. code:: python

        scv.tl.velocity_clusters(adata)
        scv.pl.scatter(adata, color='velocity_clusters')

    .. image:: https://user-images.githubusercontent.com/31883718/69625627-484dc480-1047-11ea-847f-6607a3430427.png
       :width: 600px


    Arguments
    ----------
    data : :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Key of velocities computed in `tl.velocity`
    match_with : `int`, optional (default: 100)
        The number of genes that appear in the returned tables.
    match_with: `str` (default: `'clusters'`)
        Match the names of the velocity clusters with the names of this key (.obs).
    sort_by: `str` or `None` (default: `'dpt_pseudotime'`)
        Sort velocity clusters by this key (.obs).
    resolution: `str` or `None` (default: `None`)
        Resolution for louvain modularity.
    min_likelihood: `float` between `0` and `1` or `None` (default: `None`)
        Only rank velocity of genes with a likelihood higher than min_likelihood.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.

    Returns
    -------
    Returns or updates `data` with the attributes
    velocity_clusters : `.obs`
        Velocity clusters obtained from applying louvain modularity on velocity expression.
    """
    adata = data.copy() if copy else data

    logg.info('computing velocity clusters', r=True)

    tmp_filter = ~np.isnan(adata.layers[vkey].sum(0))
    if vkey + '_genes' in adata.var.keys():
        tmp_filter &= np.array(adata.var[vkey + '_genes'].values, dtype=bool)

    if 'unspliced' in adata.layers.keys():
        n_counts = (adata.layers['unspliced'] > 0).sum(0)
        n_counts = n_counts.A1 if issparse(adata.layers['unspliced']) else n_counts
        min_counts = min(50, np.percentile(n_counts, 50))
        tmp_filter &= np.ravel(n_counts > min_counts)

    if 'r2' in adata.var.keys():
        r2 = adata.var.velocity_r2
        min_r2 = np.percentile(r2[r2 > 0], 50)
        tmp_filter &= (r2 > min_r2)

    if 'dispersions_norm' in adata.var.keys():
        dispersions = adata.var.dispersions_norm
        min_dispersion = np.percentile(dispersions, 20)
        tmp_filter &= (dispersions > min_dispersion)

    if 'fit_likelihood' in adata.var.keys():
        l = adata.var['fit_likelihood']
        min_likelihood = .1 if min_likelihood is None else min_likelihood
        tmp_filter &= (l > min_likelihood)

    from .. import AnnData
    vdata = AnnData(adata.layers[vkey][:, tmp_filter])
    vdata.obs = adata.obs.copy()
    vdata.var = adata.var[tmp_filter].copy()

    if 'highly_variable' in vdata.var.keys():
        vdata.var['highly_variable'] = np.array(vdata.var['highly_variable'], dtype=bool)

    import scanpy as sc
    logg.switch_verbosity('off', module='scanpy')
    sc.pp.pca(vdata, n_comps=20, svd_solver='arpack')
    sc.pp.neighbors(vdata, n_pcs=20)
    sc.tl.louvain(vdata, resolution=resolution)
    logg.switch_verbosity('on', module='scanpy')

    if sort_by in vdata.obs.keys():
        vc = vdata.obs['louvain']
        vc_cats = vc.cat.categories
        mean_times = [np.mean(vdata.obs[sort_by][vc == cat]) for cat in vc_cats]
        vdata.obs['louvain'].cat.reorder_categories(vc_cats[np.argsort(mean_times)], inplace=True)
        vdata.obs['louvain'].cat.categories = np.arange(len(vc_cats))

    if isinstance(match_with, str) and match_with in adata.obs.keys():
        from .utils import most_common_in_list
        vc = vdata.obs['louvain']
        for i, cat in enumerate(vc.cat.categories):
            cells_in_cat = np.where(vc == cat)[0]
            new_cat = most_common_in_list(adata.obs[match_with][cells_in_cat])
            vc = vc.cat.rename_categories({cat: str(new_cat) + ' ' + str(cat)})
        vdata.obs['louvain'] = vc

    adata.obs[vkey + '_clusters'] = vdata.obs['louvain'].copy()

    del vdata

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + vkey + '_clusters\', clusters based on modularity on velocity field (adata.obs)')

    return adata if copy else None


def rank_velocity_genes(data, vkey='velocity', n_genes=10, groupby=None, match_with=None, resolution=None,
                        min_counts=None, min_r2=None, min_dispersion=None, min_likelihood=None, copy=False):
    """Rank genes for velocity characterizing groups.

    This applies a differential expression test (Welch t-test with overestimated variance to be conservative) on
    velocity expression, to find genes in a cluster that show dynamics that is transcriptionally regulated differently
    compared to all other clusters (e.g. induction in that cluster and homeostasis in remaining population).
    If no clusters are given, it priorly computes velocity clusters by applying louvain modularity on velocity expression.

    .. code:: python

        scv.tl.rank_velocity_genes(adata, groupby='clusters')
        scv.pl.scatter(adata, basis=adata.uns['rank_velocity_genes']['names']['Beta'][:3])
        pd.DataFrame(adata.uns['rank_velocity_genes']['names']).head()

    .. image:: https://user-images.githubusercontent.com/31883718/69626017-11c47980-1048-11ea-89f4-df3769df5ad5.png
       :width: 600px

    .. image:: https://user-images.githubusercontent.com/31883718/69626572-30774000-1049-11ea-871f-e8a30c42f10e.png
       :width: 600px

    Arguments
    ----------
    data : :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Key of velocities computed in `tl.velocity`
    n_genes : `int`, optional (default: 100)
        The number of genes that appear in the returned tables.
    groupby: `str`, `list` or `np.ndarray` (default: `None`)
        Key of observations grouping to consider.
    match_with: `str` or `None` (default: `None`)
        adata.obs key to separatively rank velocities on.
    resolution: `str` or `None` (default: `None`)
        Resolution for louvain modularity.
    min_counts: `float` (default: None)
        Minimum count of genes for consideration.
    min_r2: `float` (default: None)
        Minimum r2 value of genes for consideration.
    min_dispersion: `float` (default: None)
        Minimum dispersion norm value of genes for consideration.
    min_likelihood: `float` between `0` and `1` or `None` (default: `None`)
        Only rank velocity of genes with a likelihood higher than min_likelihood.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.

    Returns
    -------
    Returns or updates `data` with the attributes
    rank_velocity_genes : `.uns`
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    velocity_score : `.var`
        Storing the score for each gene for each group. Ordered according to scores.
    """
    adata = data.copy() if copy else data

    if groupby is None or groupby is 'velocity_clusters':
        velocity_clusters(adata, vkey=vkey, match_with=match_with, resolution=resolution, min_likelihood=min_likelihood)
        groupby = vkey + '_clusters'

    logg.info('ranking velocity genes', r=True)

    tmp_filter = ~np.isnan(adata.layers[vkey].sum(0))
    if vkey + '_genes' in adata.var.keys():
        tmp_filter &= np.array(adata.var[vkey + '_genes'].values, dtype=bool)

    if 'unspliced' in adata.layers.keys():
        n_counts = (adata.layers['unspliced'] > 0).sum(0)
        n_counts = n_counts.A1 if issparse(adata.layers['unspliced']) else n_counts
        min_counts = min(50, np.percentile(n_counts, 50)) if min_counts is None else min_counts
        tmp_filter &= np.ravel(n_counts > min_counts)

    if vkey + '_r2' in adata.var.keys():
        r2 = adata.var[vkey + '_r2']
        min_r2 = .1 if min_r2 is None else min_r2  # np.percentile(r2[r2 > 0], 50)
        tmp_filter &= (r2 > min_r2)

    if 'dispersions_norm' in adata.var.keys():
        dispersions = adata.var.dispersions_norm
        min_dispersion = 0 if min_dispersion is None else min_dispersion  # np.percentile(dispersions, 20)
        tmp_filter &= (dispersions > min_dispersion)

    if 'fit_likelihood' in adata.var.keys():
        l = adata.var['fit_likelihood']
        min_likelihood = .1 if min_likelihood is None else min_likelihood
        tmp_filter &= (l > min_likelihood)

    X = adata[:, tmp_filter].layers[vkey]
    groups, groups_masks = select_groups(adata, key=groupby)

    n_groups = groups_masks.shape[0]
    sizes = groups_masks.sum(1)

    mean, var = np.zeros((n_groups, X.shape[1])), np.zeros((n_groups, X.shape[1]))
    for i, mask in enumerate(groups_masks): mean[i], var[i] = get_mean_var(X[mask])

    # test each against the union of all other groups
    rankings_gene_names, rankings_gene_scores, indices = [], [], []
    for i in range(n_groups):
        mask_rest = ~groups_masks[i]
        mean_rest, var_rest = get_mean_var(X[mask_rest])
        size_rest = sizes[i]  # else mask_rest.sum() if method == 't-test'

        scores = (mean[i] - mean_rest) / np.sqrt(var[i] / sizes[i] + var_rest / size_rest)
        scores = np.nan_to_num(scores)

        # equivalent to but much faster than np.argsort(scores)[-10:]
        if n_genes > X.shape[1]: n_genes = X.shape[1]
        idx = np.argpartition(scores, -n_genes)[-n_genes:]
        idx = idx[np.argsort(scores[idx])[::-1]]

        rankings_gene_names.append(adata[:, tmp_filter].var_names[idx].values)
        rankings_gene_scores.append(scores[idx])

    rankings_gene_names = np.array([list(n) for n in rankings_gene_names])
    rankings_gene_scores = np.array([list(n) for n in rankings_gene_scores])

    all_names = rankings_gene_names.T.flatten()
    all_scores = rankings_gene_scores.T.flatten()
    vscore = np.zeros(adata.n_vars, dtype=np.int)
    for i, name in enumerate(adata.var_names):
        if name in all_names: vscore[i] = all_scores[np.where(name == all_names)[0][0]]
    adata.var['velocity_score'] = vscore

    key = 'rank_velocity_genes'
    if key not in adata.uns.keys(): adata.uns[key] = {}

    adata.uns[key] = \
        {'names': np.rec.fromarrays([n for n in rankings_gene_names], dtype=[(str(rn), 'U50') for rn in groups]),
         'scores': np.rec.fromarrays([n.round(2) for n in rankings_gene_scores], dtype=[(str(rn), 'float32') for rn in groups]),
         'params': {'groupby': groupby, 'reference': 'rest', 'method': 't-test_overestim_var', 'use_raw': True}}

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added \n'
        '    \'' + key + '\', sorted scores by group ids (adata.uns)')

    return adata if copy else None
