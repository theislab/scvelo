import numpy as np
from scipy.sparse import issparse

from scvelo import logging as logg
from scvelo import settings
from .utils import strings_to_categoricals, vcorrcoef
from .velocity_pseudotime import velocity_pseudotime


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
        if np.size(perc) < 2:
            perc = [perc, 100] if perc < 50 else [0, perc]
        lb, ub = np.percentile(data, perc)
        data = np.clip(data, lb, ub)

    if issparse(X):
        mean = (X.sum(0) / n_counts).A1
        mean_sq = (X.multiply(X).sum(0) / n_counts).A1
    else:
        mean = X.sum(0) / n_counts
        mean_sq = np.multiply(X, X).sum(0) / n_counts
    n_cells = np.clip(X.shape[0], 2, None)  # to avoid division by zero
    var = (mean_sq - mean ** 2) * (n_cells / (n_cells - 1))

    mean = np.nan_to_num(mean)
    var = np.nan_to_num(var)
    return mean, var


def select_groups(adata, groups="all", key="louvain"):
    """Get subset of groups in adata.obs[key]."""
    strings_to_categoricals(adata)
    if isinstance(groups, list) and isinstance(groups[0], int):
        groups = [f"{n}" for n in groups]
    categories = adata.obs[key].cat.categories
    groups_masks = np.array(
        [categories[i] == adata.obs[key].values for i, name in enumerate(categories)]
    )
    if groups == "all":
        groups = categories.values
    else:
        groups_ids = [categories.get_loc(name) for name in groups]
        groups_masks = groups_masks[groups_ids]
        groups = categories[groups_ids].values
    return groups, groups_masks


def velocity_clusters(
    data,
    vkey="velocity",
    match_with="clusters",
    sort_by="velocity_pseudotime",
    resolution=None,
    min_likelihood=None,
    copy=False,
):
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
    resolution: `float` (default: 0.7)
        Resolution for louvain modularity.
    min_likelihood: `float` between `0` and `1` or `None` (default: `None`)
        Only rank velocity of genes with a likelihood higher than min_likelihood.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.

    Returns
    -------
    velocity_clusters : `.obs`
        Clusters obtained from applying louvain modularity on velocity expression.
    """  # noqa E501

    adata = data.copy() if copy else data

    logg.info("computing velocity clusters", r=True)

    tmp_filter = ~np.isnan(adata.layers[vkey].sum(0))
    if f"{vkey}_genes" in adata.var.keys():
        tmp_filter &= np.array(adata.var[f"{vkey}_genes"].values, dtype=bool)

    if "unspliced" in adata.layers.keys():
        n_counts = (adata.layers["unspliced"] > 0).sum(0)
        n_counts = n_counts.A1 if issparse(adata.layers["unspliced"]) else n_counts
        min_counts = min(50, np.percentile(n_counts, 50))
        tmp_filter &= np.ravel(n_counts > min_counts)

    if "r2" in adata.var.keys():
        r2 = adata.var.velocity_r2
        min_r2 = np.percentile(r2[r2 > 0], 50)
        tmp_filter &= r2 > min_r2

    if "dispersions_norm" in adata.var.keys():
        dispersions = adata.var.dispersions_norm
        min_dispersion = np.percentile(dispersions, 20)
        tmp_filter &= dispersions > min_dispersion

    if "fit_likelihood" in adata.var.keys() and min_likelihood is not None:
        tmp_filter &= adata.var["fit_likelihood"] > min_likelihood

    from anndata import AnnData

    vdata = AnnData(adata.layers[vkey][:, tmp_filter])
    vdata.obs = adata.obs.copy()
    vdata.var = adata.var[tmp_filter].copy()

    if "highly_variable" in vdata.var.keys():
        vdata.var["highly_variable"] = np.array(
            vdata.var["highly_variable"], dtype=bool
        )

    import scanpy as sc

    logg.switch_verbosity("off", module="scanpy")
    sc.pp.pca(vdata, n_comps=20, svd_solver="arpack")
    sc.pp.neighbors(vdata, n_pcs=20)
    sc.tl.louvain(vdata, resolution=0.7 if resolution is None else resolution)
    logg.switch_verbosity("on", module="scanpy")

    if sort_by == "velocity_pseudotime" and sort_by not in adata.obs.keys():
        velocity_pseudotime(adata, vkey=vkey)
    if sort_by in vdata.obs.keys():
        vc = vdata.obs["louvain"]
        vc_cats = vc.cat.categories
        mean_times = [np.mean(vdata.obs[sort_by][vc == cat]) for cat in vc_cats]
        vdata.obs["louvain"].cat.reorder_categories(
            vc_cats[np.argsort(mean_times)], inplace=True
        )

    if isinstance(match_with, str) and match_with in adata.obs.keys():
        from .utils import most_common_in_list

        vc = vdata.obs["louvain"]
        cats_nums = {cat: 0 for cat in adata.obs[match_with].cat.categories}
        for cat in vc.cat.categories:
            cells_in_cat = np.where(vc == cat)[0]
            new_cat = most_common_in_list(adata.obs[match_with][cells_in_cat])
            cats_nums[new_cat] += 1
            vc = vc.cat.rename_categories({cat: f"{new_cat} ({cats_nums[new_cat]})"})
        vdata.obs["louvain"] = vc
    else:
        vdata.obs["louvain"].cat.categories = np.arange(
            len(vdata.obs["louvain"].cat.categories)
        )
    adata.obs[f"{vkey}_clusters"] = vdata.obs["louvain"].copy()

    del vdata

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        f"    '{vkey}_clusters', "
        f"clusters based on louvain modularity on velocity vector field (adata.obs)"
    )

    return adata if copy else None


def rank_velocity_genes(
    data,
    vkey="velocity",
    n_genes=100,
    groupby=None,
    match_with=None,
    resolution=None,
    min_counts=None,
    min_r2=None,
    min_corr=None,
    min_dispersion=None,
    min_likelihood=None,
    copy=False,
):
    """Rank genes for velocity characterizing groups.

    This applies a differential expression test (Welch t-test with overestimated
    variance to be conservative) on velocity expression, to find genes in a cluster that
    show dynamics that is transcriptionally regulated differently compared to all other
    clusters (e.g. induction in that cluster and homeostasis in remaining population).
    If no clusters are given, it priorly computes velocity clusters by applying louvain
    modularity on velocity expression.

    .. code:: python

        scv.tl.rank_velocity_genes(adata, groupby='clusters')
        scv.pl.scatter(
            adata, basis=adata.uns['rank_velocity_genes']['names']['Beta'][:3]
        )
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
    min_corr: `float` (default: None)
        Minimum Spearmans correlation coefficient between spliced and unspliced.
    min_dispersion: `float` (default: None)
        Minimum dispersion norm value of genes for consideration.
    min_likelihood: `float` between `0` and `1` or `None` (default: `None`)
        Only rank velocity of genes with a likelihood higher than min_likelihood.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to data.

    Returns
    -------
    rank_velocity_genes : `.uns`
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    velocity_score : `.var`
        Storing the score for each gene for each group. Ordered according to scores.
    """  # noqa E501

    adata = data.copy() if copy else data

    if groupby is None or groupby == "velocity_clusters":
        velocity_clusters(
            adata,
            vkey=vkey,
            match_with=match_with,
            resolution=resolution,
            min_likelihood=min_likelihood,
        )
        groupby = f"{vkey}_clusters"

    logg.info("ranking velocity genes", r=True)

    if "spearmans_score" not in adata.var.keys():
        corr = vcorrcoef(
            np.array(adata.layers["Ms"]).T,
            np.array(adata.layers["Mu"].T),
            mode="spearmans",
        )
        adata.var["spearmans_score"] = np.clip(corr, 0, None)

    tmp_filter = ~np.isnan(adata.layers[vkey].sum(0))
    if f"{vkey}_genes" in adata.var.keys():
        tmp_filter &= np.array(adata.var[f"{vkey}_genes"].values, dtype=bool)

    if "unspliced" in adata.layers.keys():
        n_counts = (adata.layers["unspliced"] > 0).sum(0)
        n_counts = n_counts.A1 if issparse(adata.layers["unspliced"]) else n_counts
        min_counts = (
            min(50, np.percentile(n_counts, 50)) if min_counts is None else min_counts
        )
        tmp_filter &= np.ravel(n_counts > min_counts)

    if f"{vkey}_r2" in adata.var.keys():
        r2 = adata.var[f"{vkey}_r2"]
        min_r2 = 0.1 if min_r2 is None else min_r2  # np.percentile(r2[r2 > 0], 50)
        tmp_filter &= r2 > min_r2

    if "spearmans_score" in adata.var.keys():
        corr = adata.var["spearmans_score"]
        min_corr = (
            0.1 if min_corr is None else min_corr
        )  # np.percentile(r2[r2 > 0], 50)
        tmp_filter &= corr > min_corr

    if "dispersions_norm" in adata.var.keys():
        dispersions = adata.var.dispersions_norm
        min_dispersion = 0 if min_dispersion is None else min_dispersion
        tmp_filter &= dispersions > min_dispersion

    if "fit_likelihood" in adata.var.keys():
        fit_likelihood = adata.var["fit_likelihood"]
        min_likelihood = 0.1 if min_likelihood is None else min_likelihood
        tmp_filter &= fit_likelihood > min_likelihood

    X = adata[:, tmp_filter].layers[vkey]
    groups, groups_masks = select_groups(adata, key=groupby)

    n_groups = groups_masks.shape[0]
    sizes = groups_masks.sum(1)

    mean, var = np.zeros((n_groups, X.shape[1])), np.zeros((n_groups, X.shape[1]))
    for i, mask in enumerate(groups_masks):
        mean[i], var[i] = get_mean_var(X[mask])

    # test each against the union of all other groups
    rankings_gene_names, rankings_gene_scores = [], []
    for i in range(n_groups):
        mask_rest = ~groups_masks[i]
        mean_rest, var_rest = get_mean_var(X[mask_rest])
        size_rest = sizes[i]  # else mask_rest.sum() if method == 't-test'

        scores = (mean[i] - mean_rest) / np.sqrt(
            var[i] / sizes[i] + var_rest / size_rest
        )
        scores = np.nan_to_num(scores)

        # equivalent to but much faster than np.argsort(scores)[-10:]
        if n_genes > X.shape[1]:
            n_genes = X.shape[1]
        idx = np.argpartition(scores, -n_genes)[-n_genes:]
        idx = idx[np.argsort(scores[idx])[::-1]]

        rankings_gene_names.append(adata[:, tmp_filter].var_names[idx].values)
        rankings_gene_scores.append(scores[idx])

    rankings_gene_names = np.array([list(n) for n in rankings_gene_names])
    rankings_gene_scores = np.array([list(n) for n in rankings_gene_scores])

    all_names = rankings_gene_names.T.flatten()
    all_scores = rankings_gene_scores.T.flatten()
    vscore = np.zeros(adata.n_vars, dtype=int)
    for i, name in enumerate(adata.var_names):
        if name in all_names:
            vscore[i] = all_scores[np.where(name == all_names)[0][0]]
    adata.var["velocity_score"] = vscore

    key = "rank_velocity_genes"
    if key not in adata.uns.keys():
        adata.uns[key] = {}

    adata.uns[key] = {
        "names": np.rec.fromarrays(
            [n for n in rankings_gene_names], dtype=[(f"{rn}", "U50") for rn in groups]
        ),
        "scores": np.rec.fromarrays(
            [n.round(2) for n in rankings_gene_scores],
            dtype=[(f"{rn}", "float32") for rn in groups],
        ),
        "params": {
            "groupby": groupby,
            "reference": "rest",
            "method": "t-test_overestim_var",
            "use_raw": True,
        },
    }

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        f"    '{key}', sorted scores by group ids (adata.uns) \n"
        "    'spearmans_score', spearmans correlation scores (adata.var)"
    )

    return adata if copy else None
