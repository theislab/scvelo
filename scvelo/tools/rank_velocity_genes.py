"""Rank genes according to unspliced/spliced correlation and differential expression.
"""

import numpy as np
from scipy.sparse import issparse
from ..logging import logg, settings
from scanpy.utils import select_groups


def get_mean_var(X):
    mean = X.mean(axis=0).A1 if issparse(X) else X.mean(axis=0)
    mean_sq = X.multiply(X).mean(axis=0).A1 if issparse(X) else np.multiply(X, X).mean(axis=0)
    var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    return mean, var


def rank_velocity_genes(adata, groupby=None, groups='all', n_genes=10, method='t-test_overestim_var', use_raw=True, copy=False):
    """Rank genes for characterizing groups.
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groupby : `str`
        The key of the observations grouping to consider.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    groups : `str`, `list`, optional (default: `'all'`)
        Subset of groups, e.g. `['g1', 'g2', 'g3']`, to which comparison shall
        be restricted. If not passed, a ranking will be generated for all
        groups.
    n_genes : `int`, optional (default: 100)
        The number of genes that appear in the returned tables.
    method : {'t-test', 't-test_overestim_var'} (default: 't-test_overestim_var')
        If 't-test' uses t-test, if 't-test_overestim_var' overestimates variance of each group.

    Returns
    -------
    Updates `adata` with the following fields.
    names : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    scores : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the score for each
        gene for each group. Ordered according to scores.

    """
    logg.info('ranking velocity genes', r=True)
    if method not in {'t-test', 't-test_overestim_var'}: raise ValueError('Method not available.')

    adata = adata.copy() if copy else adata

    counts = adata.layers['unspliced'] > 0
    n_counts = counts.sum(0).A1 if issparse(counts) else counts.sum(0)
    filter_counts = n_counts > np.percentile(n_counts, 70)

    filter_r2 = adata.var.velocity_r2 > np.percentile(adata.var.velocity_r2[adata.var.velocity_r2 > 0], 70)
    idx_sub = filter_counts & filter_r2

    adata_sub = adata[:, idx_sub]
    if groupby is None:
        return np.array(adata_sub.var_names[np.argsort(n_counts[idx_sub])[::-1][:n_genes]])
    else:
        X = adata_sub.raw.X if adata_sub.raw is not None and use_raw else adata_sub.X
        if n_genes > X.shape[1]: n_genes = X.shape[1]

        if isinstance(groups, list) and isinstance(groups[0], int): groups = [str(n) for n in groups]
        groups, groups_masks = select_groups(adata_sub, groups, groupby)
        n_groups = groups_masks.shape[0]
        sizes = groups_masks.sum(1)

        mean, var = np.zeros((n_groups, X.shape[1])), np.zeros((n_groups, X.shape[1]))
        for i, mask in enumerate(groups_masks): mean[i], var[i] = get_mean_var(X[mask])

        # test each against the union of all other groups
        rankings_gene_names, rankings_gene_scores, rankings_gene_r2 = [], [], []
        for i in range(n_groups):
            mask_rest = ~groups_masks[i]
            mean_rest, var_rest = get_mean_var(X[mask_rest])
            size_rest = mask_rest.sum() if method == 't-test' else sizes[i]

            scores = (mean[i] - mean_rest) / np.sqrt(var[i] / sizes[i] + var_rest / size_rest)
            scores[np.isnan(scores)] = 0

            # equivalent to but much faster than np.argsort(scores)[-10:]
            idx = np.argpartition(scores, -n_genes)[-n_genes:]
            idx = idx[np.argsort(scores[idx])[::-1]]

            rankings_gene_names.append(adata_sub.var_names[idx].values)
            rankings_gene_scores.append(scores[idx])
            rankings_gene_r2.append(adata_sub.var['velocity_r2'][idx])

        rankings_gene_names = np.array([list(n) for n in rankings_gene_names])
        rankings_gene_scores = np.array([list(n) for n in rankings_gene_scores])
        rankings_gene_gammas = np.array([list(n) for n in rankings_gene_r2])

        #rankings_gene_names = np.rec.fromarrays([n for n in rankings_gene_names], dtype=[(rn, 'U50') for rn in groups])
        #rankings_gene_scores = np.rec.fromarrays([n for n in rankings_gene_scores], dtype=[(rn, 'float32') for rn in groups])
        #rankings_gene_gammas = np.rec.fromarrays([n for n in rankings_gene_gammas], dtype=[(rn, 'float32') for rn in groups])

        key = 'rank_velocity_genes'
        if key not in adata.uns.keys(): adata.uns[key] = {}
        adata.uns[key] = {'groups': groups, 'names': rankings_gene_names,
                          'scores': rankings_gene_scores.round(0), 'r2': rankings_gene_gammas.round(2)}

        logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
        logg.hint(
            'added to `.uns`\n'
            '    \'' + key + '\', sorted np.recarray to be indexed by group ids\n')

        return adata if copy else None
