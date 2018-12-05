from .. import logging as logg

import numpy as np
from scipy.sparse import issparse, csr_matrix
from scanpy.api.pp import log1p


def show_proportions(adata):
    """Fraction of spliced/unspliced/ambiguous abundances

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.

    Returns
    -------
    Prints the fractions of abundances.
    """
    layers_keys = [key for key in ['spliced', 'unspliced', 'ambiguous'] if key in adata.layers.keys()]
    tot_mol_cell_layers = [adata.layers[key].sum(1) for key in layers_keys]

    mean_abundances = np.round(
        [np.mean(tot_mol_cell / np.sum(tot_mol_cell_layers, 0)) for tot_mol_cell in tot_mol_cell_layers], 2)

    print('Abundance of ' + str(layers_keys) + ': ' + str(mean_abundances))


def cleanup(data, clean='layers', keep=None, copy=False):
    """Deletes attributes not needed.

    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    clean: `str` or list of `str` (default: `layers`)
        Which attributes to consider for freeing memory.
    keep: `str` or list of `str` (default: `['spliced', unspliced']`)
        Which attributes to keep.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Returns or updates `adata` with selection of attributes kept.
    """
    adata = data.copy() if copy else data

    keep = list() if keep is None else list(keep)
    keep.extend(['spliced', 'unspliced', 'Ms', 'Mu', 'clusters', 'neighbors'])

    if any(['obs' in clean, 'all' in clean]):
        for key in list(adata.obs.keys()):
            if key not in keep: del adata.obs[key]

    if any(['var' in clean, 'all' in clean]):
        for key in list(adata.var.keys()):
            if key not in keep: del adata.var[key]

    if any(['uns' in clean, 'all' in clean]):
        for key in list(adata.uns.keys()):
            if key not in keep: del adata.uns[key]

    if any(['layers' in clean, 'all' in clean]):
        for key in list(adata.layers.keys()):  # remove layers that are not needed
            if key not in keep: del adata.layers[key]

    return adata if copy else None


def X_is_logarithmized(adata):
    idx = np.random.choice(adata.n_obs, 300)
    X, S = adata.X[idx], adata.layers['spliced'][idx].copy()
    if issparse(S) or issparse(X):
        S, S_log = csr_matrix(S), csr_matrix(X)
        S.data = 1 / np.log1p(S.data)
        W = S_log.multiply(S)
        W.data /= np.mean(W.data)
        val = max(np.max(W.data), np.max(1 / W.data))
    else:
        S = 1 / np.log1p(S)
        S[np.isinf(S)] = 0

        W = X * S
        W /= np.mean(W[W > 0])
        W_inv = 1 / W
        W_inv[np.isinf(W_inv)] = 0

        val = max(np.max(W), np.max(W_inv))

    return val < 3  # approx. 3 * standard deviation


def get_size(adata, layer):
    X = adata.layers[layer]
    return X.sum(1).A1.copy() if issparse(X) else X.sum(1).copy()


def set_initial_size(adata, layers={'spliced', 'unspliced'}):
    if all([layer in adata.layers.keys() for layer in layers]):
        layers = [layer for layer in layers if 'initial_size_' + layer not in adata.obs.keys()]
        total_size = 0
        for layer in layers:
            size = get_size(adata, layer)
            adata.obs['initial_size_' + layer] = size
            total_size += size
        if 'initial_size' not in adata.obs.keys(): adata.obs['initial_size'] = total_size


def get_initial_size(adata, layer, by_total_size=None):
    if layer not in {'spliced', 'unspliced'}: return None
    else: return adata.obs['initial_size'].copy() if by_total_size and 'initial_size' in adata.obs.keys() else \
        adata.obs['initial_size_' + layer].copy() if 'initial_size_' + layer in adata.obs.keys() else get_size(adata, layer)


def filter(X, min_counts=None, min_cells=None, max_counts=None, max_cells=None):
    counts = X if (min_counts is not None or max_counts is not None) else X > 0
    counts = counts.sum(0).A1 if issparse(counts) else counts.sum(0)
    lb = min_counts if min_counts is not None else min_cells if min_cells is not None else -np.inf
    ub = max_counts if max_counts is not None else max_cells if max_cells is not None else np.inf
    return (lb <= counts) & (counts <= ub), counts


def filter_genes(data, min_counts=None, min_cells=None, max_counts=None, max_cells=None,
                 min_counts_u=None,  min_cells_u=None, max_counts_u=None, max_cells_u=None, copy=False):
    """Filter genes based on number of cells or counts.
    Keep genes that have at least `min_counts` counts or are expressed in at
    least `min_cells` cells or have at most `max_counts` counts or are expressed
    in at most `max_cells` cells.
    Only provide one of the optional parameters `min_counts`, `min_cells`,
    `max_counts`, `max_cells` per call.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.spmatrix`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts : `int`, optional (default: `None`)
        Minimum number of counts required for a gene to pass filtering.
    min_cells : `int`, optional (default: `None`)
        Minimum number of cells expressed required for a gene to pass filtering.
    max_counts : `int`, optional (default: `None`)
        Maximum number of counts required for a gene to pass filtering.
    max_cells : `int`, optional (default: `None`)
        Maximum number of cells expressed required for a gene to pass filtering.
    min_counts_u : `int`, optional (default: `None`)
        Minimum number of unspliced counts required for a gene to pass filtering.
    min_cells_u : `int`, optional (default: `None`)
        Minimum number of unspliced cells expressed required for a gene to pass filtering.
    max_counts_u : `int`, optional (default: `None`)
        Maximum number of unspliced counts required for a gene to pass filtering.
    max_cells_u : `int`, optional (default: `None`)
        Maximum number of unspliced cells expressed required for a gene to pass filtering.
    copy : `bool`, optional (default: `False`)
        Determines whether a copy is returned.

    Returns
    -------
    Filters the object and adds `n_counts` to `adata.var`.
    """
    adata = data.copy() if copy else data

    # set initial cell sizes before filtering
    set_initial_size(adata)

    layers = [layer for layer in ['spliced', 'unspliced'] if layer in adata.layers.keys()]
    for layer in layers:
        X = adata.layers[layer]
        if layer == 'spliced':
            _min_counts, _min_cells, _max_counts, _max_cells = min_counts, min_cells, max_counts, max_cells
        else:
            _min_counts, _min_cells, _max_counts, _max_cells = min_counts_u, min_cells_u, max_counts_u, max_cells_u

        if layer in adata.layers.keys():
            gene_subset = np.ones(adata.n_vars, dtype=bool)
            if _min_counts is not None or _max_counts is not None:
                gene_subset, _ = filter(X, min_counts=_min_counts, max_counts=_max_counts)
                adata._inplace_subset_var(gene_subset)

            if _min_cells is not None or _max_cells is not None:
                gene_subset, _ = filter(X, min_cells=_min_cells, max_cells=_max_cells)
                adata._inplace_subset_var(gene_subset)

            s = np.sum(~gene_subset)
            if s > 0:
                logg.info('Filtered out {} genes that are detected'.format(s), end=' ')
                if _min_cells is not None or _min_counts is not None:
                    logg.info('in less than', str(_min_cells) + ' cells (' + str(layer) + ').' if _min_counts is None
                              else str(_min_counts) + ' counts (' + str(layer) + ').', no_indent=True)
                if max_cells is not None or max_counts is not None:
                    logg.info('in more than ', str(_max_cells) + ' cells(' + str(layer) + ').' if _max_counts is None
                              else str(_max_counts) + ' counts (' + str(layer) + ').', no_indent=True)

    return adata if copy else None


def filter_genes_dispersion(data, flavor='seurat', min_disp=None, max_disp=None, min_mean=None, max_mean=None,
                            n_bins=20, n_top_genes=None, log=True, copy=False):
    """Extract highly variable genes.
    The normalized dispersion is obtained by scaling with the mean and standard
    deviation of the dispersions for genes falling into a given bin for mean
    expression of genes. This means that for each bin of mean expression, highly
    variable genes are selected.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    flavor : {'seurat', 'cell_ranger', 'svr'}, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion. If choosing
        'seurat', this expects non-logarithmized data - the logarithm of mean
        and dispersion is taken internally when `log` is at its default value
        `True`. For 'cell_ranger', this is usually called for logarithmized data
        - in this case you should set `log` to `False`. In their default
        workflows, Seurat passes the cutoffs whereas Cell Ranger passes
        `n_top_genes`.
    min_mean=0.0125, max_mean=3, min_disp=0.5, max_disp=`None` : `float`, optional
        If `n_top_genes` unequals `None`, these cutoffs for the means and the
        normalized dispersions are ignored.
    n_bins : `int` (default: 20)
        Number of bins for binning the mean gene expression. Normalization is
        done with respect to each bin. If just a single gene falls into a bin,
        the normalized dispersion is artificially set to 1. You'll be informed
        about this if you set `settings.verbosity = 4`.
    n_top_genes : `int` or `None` (default: `None`)
        Number of highly-variable genes to keep.
    log : `bool`, optional (default: `True`)
        Use the logarithm of the mean to variance ratio.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    If an AnnData `adata` is passed, returns or updates `adata` depending on \
    `copy`. It filters the `adata` and adds the annotations
    """
    adata = data.copy() if copy else data
    set_initial_size(adata)
    if n_top_genes is not None and adata.n_vars < n_top_genes:
        logg.info('Skip filtering by dispersion since number of variables are less than `n_top_genes`')
    else:
        if flavor == 'svr':
            mu = adata.X.mean(0).A1 if issparse(adata.X) else adata.X.mean(0)
            sigma = np.sqrt(adata.X.multiply(adata.X).mean(0).A1 - mu ** 2) if issparse(adata.X) else adata.X.std(0)
            log_mu = np.log2(mu)
            log_cv = np.log2(sigma / mu)

            from sklearn.svm import SVR
            clf = SVR(gamma=150. / len(mu))
            clf.fit(log_mu[:, None], log_cv)
            score = log_cv - clf.predict(log_mu[:, None])
            nth_score = np.sort(score)[::-1][n_top_genes]
            adata._inplace_subset_var(score >= nth_score)
        else:
            from scanpy.api.pp import filter_genes_dispersion
            filter_genes_dispersion(adata, flavor=flavor, min_disp=min_disp, max_disp=max_disp, min_mean=min_mean,
                                    max_mean=max_mean, n_bins=n_bins, n_top_genes=n_top_genes, log=log)
    return adata if copy else None


def counts_per_cell_quantile(X, max_proportion_per_cell=.05, counts_per_cell=None):
    if counts_per_cell is None:
        counts_per_cell = X.sum(1).A1 if issparse(X) else X.sum(1)
    gene_subset = np.all(X <= counts_per_cell[:, None] * max_proportion_per_cell, axis=0)
    if issparse(X): gene_subset = gene_subset.A1
    return X[:, gene_subset].sum(1).A1 if issparse(X) else X[:, gene_subset].sum(1)


def not_yet_normalized(X):
    return np.allclose((X.data[:10] if issparse(X) else X[0]) % 1, 0, atol=1e-3)


def normalize_layers(data, layers=['spliced', 'unspliced'], counts_per_cell_after=None, max_proportion_per_cell=None,
                     by_total_size=None, enforce=False, copy=False):
    """Normalize by total counts to median.
    """
    adata = data.copy() if copy else data
    from scanpy.api.pp import normalize_per_cell

    for layer in layers:
        if not_yet_normalized(adata.layers[layer]) or enforce:
            counts_per_cell = get_initial_size(adata, layer, by_total_size)
            if max_proportion_per_cell is not None and (0 < max_proportion_per_cell < 1):
                counts_per_cell = counts_per_cell_quantile(adata.X, max_proportion_per_cell, counts_per_cell)
            counts_per_cell += counts_per_cell == 0
            adata.layers[layer] = normalize_per_cell(adata.layers[layer], counts_per_cell_after, counts_per_cell, copy=True)
    return adata if copy else None


def normalize_per_cell(data, counts_per_cell_after=None, counts_per_cell=None, key_n_counts=None,
                       max_proportion_per_cell=None, layers=['spliced', 'unspliced'], enforce=False, copy=False):
    """Normalize each cell by total counts over all genes.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    key_n_counts : `str`, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    max_proportion_per_cell : `int` (default: `None`)
        Exclude genes counts that account for more than a specific proportion of cell size, e.g. 0.05.
    layers : `str` or `list` (default: `{'spliced', 'unspliced'}`)
        Keys for layers to be also considered for normalization.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Returns or updates `adata` with normalized version of the original `adata.X`, depending on `copy`.
    """
    adata = data.copy() if copy else data
    from scanpy.api.pp import normalize_per_cell

    if max_proportion_per_cell is not None and (0 < max_proportion_per_cell < 1):
        counts_per_cell = counts_per_cell_quantile(adata.X, max_proportion_per_cell)

    if not_yet_normalized(adata.X) or enforce:
        normalize_per_cell(adata, counts_per_cell_after, counts_per_cell, key_n_counts)
        add_msg_str = 'X and '
    else:
        add_msg_str = ''

    layers = [layers] if isinstance(layers, str) else [layer for layer in layers if layer in adata.layers.keys()]
    if all([not_yet_normalized(adata.layers[layer]) for layer in layers]) or enforce:
        normalize_layers(adata, layers, counts_per_cell_after, max_proportion_per_cell)
        logg.info('Normalized ' + add_msg_str + 'spliced/unspliced count data.')
    else:
        logg.info('Looks like it\'s already normalized.'
                  'If you want to (re-)normalize your data, use `scv.pp.normalize_per_cell(adata, enforce=True)`.')

    return adata if copy else None


def filter_and_normalize(data, min_counts=None, min_counts_u=None, min_cells=None, min_cells_u=None, n_top_genes=None,
                         flavor='seurat', log=True, copy=False):
    """Filtering, normalization and log transform

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    Runs the following steps

    .. code:: python

        scv.pp.filter_genes(adata)
        scv.pp.normalize_per_cell(adata)
        if n_top_genes is not None:
            scv.pp.filter_genes_dispersion(adata)
        if log:
            scv.pp.log1p(adata)


    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    min_counts: `int` (default: `None`)
        Minimum number of counts required for a gene to pass filtering (spliced).
    min_counts_u: `int` (default: `None`)
        Minimum number of counts required for a gene to pass filtering (unspliced).
    min_cells: `int` (default: `None`)
        Minimum number of cells expressed required for a gene to pass filtering (spliced).
    min_cells_u: `int` (default: `None`)
        Minimum number of cells expressed required for a gene to pass filtering (unspliced).
    n_top_genes: `int` (default: `None`)
        Number of genes to keep.
    flavor: {'seurat', 'cell_ranger', 'svr'}, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion. If choosing 'seurat', this expects non-logarithmized data.
    log: `bool` (default: `True`)
        Take logarithm.
    copy: `bool` (default: `False`)
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    adata = data.copy() if copy else data

    if 'spliced' in adata.layers.keys() and 'unspliced' in adata.layers.keys():
        X_not_yet_processed = np.all(adata.X.data[:100] == adata.layers['spliced'].data[:100])
    else:
        raise ValueError('Could not find spliced / unspliced counts.')

    filter_genes(adata, min_counts=min_counts, min_counts_u=min_counts_u, min_cells=min_cells, min_cells_u=min_cells_u)
    normalize_per_cell(adata)
    if n_top_genes is not None: filter_genes_dispersion(adata, n_top_genes=n_top_genes, flavor=flavor)

    if log and X_not_yet_processed:
        log1p(adata)
        logg.info('Logarithmized X.')
    elif log: logg.info('Did not modify X as it looks preprocessed already.')
    elif X_not_yet_processed: logg.info('Consider logarithmizing adata.X with `scv.pp.log1p` for better results.')

    return adata if copy else None


def recipe_velocity(adata, min_counts=3, min_counts_u=3, n_top_genes=None, n_pcs=30, n_neighbors=30, log=True, copy=False):
    """Runs pp.filter_and_normalize() and pp.moments()
    """
    from .moments import moments
    filter_and_normalize(adata, min_counts=min_counts, min_counts_u=min_counts_u, n_top_genes=n_top_genes, log=log)
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata if copy else None
