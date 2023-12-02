import warnings

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.utils import sparsefuncs

from scvelo import logging as logg
from scvelo.core import get_initial_size, get_size, multiply, set_initial_size, sum


def _filter(X, min_counts=None, min_cells=None, max_counts=None, max_cells=None):
    counts = (
        sum(X, axis=0)
        if (min_counts is not None or max_counts is not None)
        else sum(X > 0, axis=0)
    )
    lb = (
        min_counts
        if min_counts is not None
        else min_cells
        if min_cells is not None
        else -np.inf
    )
    ub = (
        max_counts
        if max_counts is not None
        else max_cells
        if max_cells is not None
        else np.inf
    )
    return (lb <= counts) & (counts <= ub), counts


def filter_genes(
    data,
    min_counts=None,
    min_cells=None,
    max_counts=None,
    max_cells=None,
    min_counts_u=None,
    min_cells_u=None,
    max_counts_u=None,
    max_cells_u=None,
    min_shared_counts=None,
    min_shared_cells=None,
    retain_genes=None,
    copy=False,
):
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
        Minimum number of unspliced cells expressed required to pass filtering.
    max_counts_u : `int`, optional (default: `None`)
        Maximum number of unspliced counts required for a gene to pass filtering.
    max_cells_u : `int`, optional (default: `None`)
        Maximum number of unspliced cells expressed required to pass filtering.
    min_shared_counts: `int`, optional (default: `None`)
        Minimum number of counts (both unspliced and spliced) required for a gene.
    min_shared_cells: `int`, optional (default: `None`)
        Minimum number of cells required to be expressed (both unspliced and spliced).
    retain_genes: `list`, optional (default: `None`)
        List of gene names to be retained independent of thresholds.
    copy : `bool`, optional (default: `False`)
        Determines whether a copy is returned.

    Returns
    -------
    Filters the object and adds `n_counts` to `adata.var`.
    """
    adata = data.copy() if copy else data

    # set initial cell sizes before filtering
    set_initial_size(adata)

    layers = [
        layer for layer in ["spliced", "unspliced"] if layer in adata.layers.keys()
    ]
    if min_shared_counts is not None or min_shared_cells is not None:
        layers.extend(["shared"])

    for layer in layers:
        if layer == "spliced":
            _min_counts, _min_cells, _max_counts, _max_cells = (
                min_counts,
                min_cells,
                max_counts,
                max_cells,
            )
        elif layer == "unspliced":
            _min_counts, _min_cells, _max_counts, _max_cells = (
                min_counts_u,
                min_cells_u,
                max_counts_u,
                max_cells_u,
            )
        else:  # shared counts/cells
            _min_counts, _min_cells, _max_counts, _max_cells = (
                min_shared_counts,
                min_shared_cells,
                None,
                None,
            )

        if layer in adata.layers.keys():
            X = adata.layers[layer]
        else:  # shared counts/cells
            Xs, Xu = adata.layers["spliced"], adata.layers["unspliced"]

            nonzeros = multiply(Xs > 0, Xu > 0)
            X = multiply(nonzeros, Xs) + multiply(nonzeros, Xu)

        gene_subset = np.ones(adata.n_vars, dtype=bool)

        if _min_counts is not None or _max_counts is not None:
            gene_subset &= _filter(X, min_counts=_min_counts, max_counts=_max_counts)[0]

        if _min_cells is not None or _max_cells is not None:
            gene_subset &= _filter(X, min_cells=_min_cells, max_cells=_max_cells)[0]

        if retain_genes is not None:
            if isinstance(retain_genes, str):
                retain_genes = [retain_genes]
            gene_subset |= adata.var_names.isin(retain_genes)

        adata._inplace_subset_var(gene_subset)

        s = np.sum(~gene_subset)
        if s > 0:
            logg.info(f"Filtered out {s} genes that are detected", end=" ")
            if _min_cells is not None or _min_counts is not None:
                logg.info(
                    f"in less than {_min_cells} cells ({layer})."
                    if _min_counts is None
                    else f"{_min_counts} counts ({layer}).",
                    no_indent=True,
                )
            if max_cells is not None or max_counts is not None:
                logg.info(
                    f"in more than {_max_cells} cells ({layer})."
                    if _max_counts is None
                    else f"{_max_counts} counts ({layer}).",
                    no_indent=True,
                )

    return adata if copy else None


# TODO: Add docstrings
def get_mean_var(X, ignore_zeros=False, perc=None):
    """TODO."""
    data = X.data if issparse(X) else X
    mask_nans = np.isnan(data) | np.isinf(data) | np.isneginf(data)

    if issparse(X):
        n_nonzeros = X.getnnz(axis=0)
    else:
        n_nonzeros = (X != 0).sum(axis=0)

    if ignore_zeros:
        n_counts = n_nonzeros
    else:
        n_counts = X.shape[0]

    if mask_nans.sum() > 0:
        if issparse(X):
            data[mask_nans] = 0
            n_nans = (n_nonzeros - (X != 0).sum(0)).A1
        else:
            X[mask_nans] = 0
            n_nans = mask_nans.sum(0)
        n_counts -= n_nans

    if perc is not None:
        if np.size(perc) < 2:
            perc = [perc, 100] if perc < 50 else [0, perc]
        lb, ub = np.percentile(data, perc)
        if issparse(X):
            X.data = np.clip(data, lb, ub)
        else:
            X = np.clip(data, lb, ub)

    if issparse(X):
        mean = (X.sum(0) / n_counts).A1
        mean_sq = (X.multiply(X).sum(0) / n_counts).A1
    else:
        mean = X.sum(0) / n_counts
        mean_sq = np.multiply(X, X).sum(0) / n_counts

    n_counts = np.clip(n_counts, 2, None)  # to avoid division by zero
    var = (mean_sq - mean**2) * (n_counts / (n_counts - 1))

    mean = np.nan_to_num(mean)
    var = np.nan_to_num(var)

    return mean, var


# TODO: Finish docstrings
def materialize_as_ndarray(key):
    """Convert distributed arrays to ndarrays."""
    if isinstance(key, (list, tuple)):
        return tuple(np.asarray(arr) for arr in key)
    return np.asarray(key)


def filter_genes_dispersion(
    data,
    flavor="seurat",
    min_disp=None,
    max_disp=None,
    min_mean=None,
    max_mean=None,
    n_bins=20,
    n_top_genes=None,
    retain_genes=None,
    log=True,
    subset=True,
    copy=False,
):
    """Extract highly variable genes.

    Expects non-logarithmized data.
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
    retain_genes: `list`, optional (default: `None`)
        List of gene names to be retained independent of thresholds.
    log : `bool`, optional (default: `True`)
        Use the logarithm of the mean to variance ratio.
    subset : `bool`, optional (default: `True`)
        Keep highly-variable genes only (if True) else write a bool
        array for highly-variable genes while keeping all genes.
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

    mean, var = materialize_as_ndarray(get_mean_var(adata.X))

    if n_top_genes is not None and adata.n_vars < n_top_genes:
        logg.info(
            "Skip filtering by dispersion since number "
            "of variables are less than `n_top_genes`."
        )
    else:
        if flavor == "svr":
            from sklearn.svm import SVR

            log_mu = np.log2(mean)
            log_cv = np.log2(np.sqrt(var) / mean)
            clf = SVR(gamma=150.0 / len(mean))
            clf.fit(log_mu[:, None], log_cv)
            score = log_cv - clf.predict(log_mu[:, None])
            nth_score = np.sort(score)[::-1][n_top_genes - 1]
            adata.var["highly_variable"] = score >= nth_score

        else:
            cut_disp = [min_disp, max_disp, min_mean, max_mean]
            if n_top_genes is not None and not all(x is None for x in cut_disp):
                logg.info("If you pass `n_top_genes`, all cutoffs are ignored.")
            if min_disp is None:
                min_disp = 0.5
            if max_disp is None:
                max_disp = np.inf
            if min_mean is None:
                min_mean = 0.0125
            if max_mean is None:
                max_mean = 3

            mean[mean == 0] = 1e-12  # set entries equal to zero to small value
            dispersion = var / mean
            if log:  # logarithmized mean as in Seurat
                dispersion[dispersion == 0] = np.nan
                dispersion = np.log(dispersion)
                mean = np.log1p(mean)

            # all of the following quantities are "per-gene" here
            df = pd.DataFrame()
            df["mean"], df["dispersion"] = mean, dispersion

            if flavor == "seurat":
                df["mean_bin"] = pd.cut(df["mean"], bins=n_bins)
                disp_grouped = df.groupby("mean_bin")["dispersion"]
                disp_mean_bin = disp_grouped.mean()
                disp_std_bin = disp_grouped.std(ddof=1)

                # retrieve genes that have nan std (i.e. single gene fell in one bin)
                # and implicitly set them to have a normalized disperion of 1
                one_gene_per_bin = disp_std_bin.isnull()

                disp_std_bin[one_gene_per_bin] = disp_mean_bin[one_gene_per_bin].values
                disp_mean_bin[one_gene_per_bin] = 0

                # normalized dispersion
                mu = disp_mean_bin[df["mean_bin"].values].values
                std = disp_std_bin[df["mean_bin"].values].values
                df["dispersion_norm"] = ((df["dispersion"] - mu) / std).fillna(0)
            elif flavor == "cell_ranger":
                from statsmodels import robust

                cut = np.percentile(df["mean"], np.arange(10, 105, 5))
                df["mean_bin"] = pd.cut(df["mean"], np.r_[-np.inf, cut, np.inf])
                disp_grouped = df.groupby("mean_bin")["dispersion"]
                disp_median_bin = disp_grouped.median()
                with warnings.catch_warnings():  # ignore warning: "Mean of empty slice"
                    warnings.simplefilter("ignore")
                    disp_mad_bin = disp_grouped.apply(robust.mad)
                mu = disp_median_bin[df["mean_bin"].values].values
                std = disp_mad_bin[df["mean_bin"].values].values
                df["dispersion_norm"] = (np.abs(df["dispersion"] - mu) / std).fillna(0)
            else:
                raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')
            dispersion_norm = df["dispersion_norm"].values
            if n_top_genes is not None:
                cut_off = df["dispersion_norm"].nlargest(n_top_genes).values[-1]
                gene_subset = df["dispersion_norm"].values >= cut_off
            else:
                gene_subset = np.logical_and.reduce(
                    (
                        mean > min_mean,
                        mean < max_mean,
                        dispersion_norm > min_disp,
                        dispersion_norm < max_disp,
                    )
                )

            adata.var["means"] = df["mean"].values
            adata.var["dispersions"] = df["dispersion"].values
            adata.var["dispersions_norm"] = df["dispersion_norm"].values
            adata.var["highly_variable"] = gene_subset

        if subset:
            gene_subset = adata.var["highly_variable"]
            if retain_genes is not None:
                if isinstance(retain_genes, str):
                    retain_genes = [retain_genes]
                gene_subset = gene_subset | adata.var_names.isin(retain_genes)
            adata._inplace_subset_var(gene_subset)

        logg.info(f"Extracted {np.sum(gene_subset)} highly variable genes.")
    return adata if copy else None


# TODO: Add docstrings
def csr_vcorrcoef(X, y):
    """TODO."""
    mu_x = np.ravel(np.mean(X, axis=-1))
    mu_y = np.ravel(np.mean(y, axis=-1))
    nom = X.dot(y) - X.dot(np.repeat(mu_y, len(y))) - mu_x * np.sum(y - mu_y)

    if X.ndim == 1:
        n_features = len(X)
    else:
        n_features = X.shape[1]

    denom_x = (
        np.ravel(np.sum(X.multiply(X), axis=-1))
        if issparse(X)
        else np.sum(X * X, axis=-1)
    )
    denom_x = denom_x - 2 * np.ravel(np.sum(X, axis=-1)) * mu_x + n_features * mu_x**2
    denom_y = (
        np.ravel(np.sum(y * y, axis=-1))
        - 2 * (np.ravel(np.sum(y, axis=-1)) * mu_y)
        + n_features * mu_y**2
    )
    return nom / np.sqrt(denom_x * denom_y)


# TODO: Add docstrings
def counts_per_cell_quantile(X, max_proportion_per_cell=0.05, counts_per_cell=None):
    """TODO."""
    if counts_per_cell is None:
        counts_per_cell = sum(X, axis=1)
    gene_subset = np.all(
        X <= counts_per_cell[:, None] * max_proportion_per_cell, axis=0
    )
    if issparse(X):
        gene_subset = gene_subset.A1
    return sum(X[:, gene_subset], axis=1)


# TODO: Add docstrings
def not_yet_normalized(X):
    """TODO."""
    return np.allclose(np.ravel(X[:5].data if issparse(X) else X[:5]) % 1, 0, atol=1e-3)


# TODO: Add docstrings
def check_if_valid_dtype(adata, layer="X"):
    """TODO."""
    X = adata.X if layer == "X" else adata.layers[layer]
    if "int" in X.dtype.name:
        if layer == "X":
            adata.X = adata.X.astype(np.float32)
        elif layer in adata.layers.keys():
            adata.layers[layer] = adata.layers[layer].astype(np.float32)


def normalize_per_cell(
    data,
    counts_per_cell_after=None,
    counts_per_cell=None,
    key_n_counts=None,
    max_proportion_per_cell=None,
    use_initial_size=True,
    layers=None,
    enforce=None,
    copy=False,
):
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
        Exclude genes counts that account for more than
        a specific proportion of cell size, e.g. 0.05.
    use_initial_size : `bool` (default: `True`)
        Whether to use initial cell sizes oder actual cell sizes.
    layers : `str` or `list` (default: `['spliced', 'unspliced']`)
        Keys for layers to be also considered for normalization.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Returns or updates `adata` with normalized counts.
    """
    adata = data.copy() if copy else data
    if layers is None:
        layers = ["spliced", "unspliced"]
    elif layers == "all":
        layers = adata.layers.keys()
    elif isinstance(layers, str):
        layers = [layers]
    layers = ["X"] + [layer for layer in layers if layer in adata.layers.keys()]
    modified_layers = []

    if isinstance(counts_per_cell, str):
        if counts_per_cell not in adata.obs.keys():
            set_initial_size(adata, layers)
        counts_per_cell = (
            adata.obs[counts_per_cell].values
            if counts_per_cell in adata.obs.keys()
            else None
        )

    for layer in layers:
        check_if_valid_dtype(adata, layer)
        X = adata.X if layer == "X" else adata.layers[layer]

        if not_yet_normalized(X) or enforce:
            counts = (
                counts_per_cell
                if counts_per_cell is not None
                else get_initial_size(adata, layer)
                if use_initial_size
                else get_size(adata, layer)
            )
            if max_proportion_per_cell is not None and (
                0 < max_proportion_per_cell < 1
            ):
                counts = counts_per_cell_quantile(X, max_proportion_per_cell, counts)
            # equivalent to sc.pp.normalize_per_cell(X, counts_per_cell_after, counts)
            counts_after = (
                np.median(counts)
                if counts_per_cell_after is None
                else counts_per_cell_after
            )

            counts_after += counts_after == 0
            counts = counts / counts_after
            counts += counts == 0  # to avoid division by zero

            if issparse(X):
                sparsefuncs.inplace_row_scale(X, 1 / counts)
            else:
                X /= np.array(counts[:, None])
            modified_layers.append(layer)
            if (
                layer == "X"
                and "gene_count_corr" not in adata.var.keys()
                and X.shape[-1] > 3e3
            ):
                # TODO: Proper handling of exception. Check why try-except is needed in
                # the first place.
                try:
                    adata.var["gene_count_corr"] = np.round(
                        csr_vcorrcoef(X.T, np.ravel((X > 0).sum(1))), 4
                    )
                except ValueError:
                    pass
        else:
            logg.warn(
                f"Did not normalize {layer} as it looks processed already. "
                "To enforce normalization, set `enforce=True`."
            )

    adata.obs["n_counts" if key_n_counts is None else key_n_counts] = get_size(adata)
    if len(modified_layers) > 0:
        logg.info("Normalized count data:", f"{', '.join(modified_layers)}.")

    return adata if copy else None


def log1p(data, copy=False):
    r"""Logarithmize the data matrix.

    Computes :math:`X = \log(X + 1)`, where :math:`log` denotes the natural logarithm.

    Parameters
    ----------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    copy: `bool` (default: `False`)
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    warnings.warn(
        "`log1p` is deprecated since scVelo v0.3.0 and will be removed in a "
        "future version. Please use `log1p` from `scanpy.pp` instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    from scanpy.pp import log1p as scanpy_log1p

    res = scanpy_log1p(data, copy=copy)
    return res if copy else None


def filter_and_normalize(
    data,
    min_counts=None,
    min_counts_u=None,
    min_cells=None,
    min_cells_u=None,
    min_shared_counts=None,
    min_shared_cells=None,
    n_top_genes=None,
    retain_genes=None,
    subset_highly_variable=True,
    flavor="seurat",
    log=True,
    layers_normalize=None,
    copy=False,
    **kwargs,
):
    """Filtering, normalization and log transform.

    Expects non-logarithmized data. If using logarithmized data, pass `log=False`.

    Runs the following steps

    .. code:: python

        scv.pp.filter_genes(adata)
        scv.pp.normalize_per_cell(adata)
        if n_top_genes is not None:
            scv.pp.filter_genes_dispersion(adata)
        if log:
            scv.pp.log1p(adata)


    Arguments:
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    min_counts: `int` (default: `None`)
        Minimum number of counts required for a gene to pass filtering (spliced).
    min_counts_u: `int` (default: `None`)
        Minimum number of counts required for a gene to pass filtering (unspliced).
    min_cells: `int` (default: `None`)
        Minimum number of cells expressed required to pass filtering (spliced).
    min_cells_u: `int` (default: `None`)
        Minimum number of cells expressed required to pass filtering (unspliced).
    min_shared_counts: `int`, optional (default: `None`)
        Minimum number of counts (both unspliced and spliced) required for a gene.
    min_shared_cells: `int`, optional (default: `None`)
        Minimum number of cells required to be expressed (both unspliced and spliced).
    n_top_genes: `int` (default: `None`)
        Number of genes to keep.
    retain_genes: `list`, optional (default: `None`)
        List of gene names to be retained independent of thresholds.
    subset_highly_variable: `bool` (default: True)
        Whether to subset highly variable genes or to store in .var['highly_variable'].
    flavor: {'seurat', 'cell_ranger', 'svr'}, optional (default: 'seurat')
        Choose the flavor for computing normalized dispersion.
        If choosing 'seurat', this expects non-logarithmized data.
    log: `bool` (default: `True`)
        Take logarithm.
    layers_normalize: list of `str` (default: None)
        List of layers to be normalized.
        If set to None, the layers {'X', 'spliced', 'unspliced'} are considered for
        normalization upon testing whether they have already been normalized
        (by checking type of entries: int -> unprocessed, float -> processed).
    copy: `bool` (default: `False`)
        Return a copy of `adata` instead of updating it.
    **kwargs:
        Keyword arguments passed to pp.normalize_per_cell (e.g. counts_per_cell).

    Returns
    -------
    Returns or updates `adata` depending on `copy`.
    """
    adata = data.copy() if copy else data

    if "spliced" not in adata.layers.keys() or "unspliced" not in adata.layers.keys():
        logg.warn("Could not find spliced / unspliced counts.")

    filter_genes(
        adata,
        min_counts=min_counts,
        min_counts_u=min_counts_u,
        min_cells=min_cells,
        min_cells_u=min_cells_u,
        min_shared_counts=min_shared_counts,
        min_shared_cells=min_shared_cells,
        retain_genes=retain_genes,
    )

    if layers_normalize is not None and "enforce" not in kwargs:
        kwargs["enforce"] = True
    normalize_per_cell(adata, layers=layers_normalize, **kwargs)

    if n_top_genes is not None:
        filter_genes_dispersion(
            adata,
            n_top_genes=n_top_genes,
            retain_genes=retain_genes,
            flavor=flavor,
            subset=subset_highly_variable,
        )

    if log:
        log1p(adata)
        logg.info("Logarithmized X.")

    return adata if copy else None


# TODO: Finish docstrings
def recipe_velocity(
    adata,
    min_counts=3,
    min_counts_u=3,
    n_top_genes=None,
    n_pcs=30,
    n_neighbors=30,
    log=True,
    copy=False,
):
    """Runs pp.filter_and_normalize() and pp.moments()."""
    from .moments import moments

    filter_and_normalize(
        adata,
        min_counts=min_counts,
        min_counts_u=min_counts_u,
        n_top_genes=n_top_genes,
        log=log,
    )
    moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata if copy else None
