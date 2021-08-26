import numpy as np
import pandas as pd
from scipy.sparse import issparse

from scvelo import logging as logg
from .utils import (
    interpret_colorkey,
    is_categorical,
    savefig_or_show,
    set_colors_for_categorical_obs,
    strings_to_categoricals,
    to_list,
)


def heatmap(
    adata,
    var_names,
    sortby="latent_time",
    layer="Ms",
    color_map="viridis",
    col_color=None,
    palette="viridis",
    n_convolve=30,
    standard_scale=0,
    sort=True,
    colorbar=None,
    col_cluster=False,
    row_cluster=False,
    context=None,
    font_scale=None,
    figsize=(8, 4),
    show=None,
    save=None,
    **kwargs,
):
    """\
    Plot time series for genes as heatmap.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names: `str`,  list of `str`
        Names of variables to use for the plot.
    sortby: `str` (default: `'latent_time'`)
        Observation key to extract time data from.
    layer: `str` (default: `'Ms'`)
        Layer key to extract count data from.
    color_map: `str` (default: `'viridis'`)
        String denoting matplotlib color map.
    col_color: `str` or list of `str` (default: `None`)
        String denoting matplotlib color map to use along the columns.
    palette: list of `str` (default: `'viridis'`)
        Colors to use for plotting groups (categorical annotation).
    n_convolve: `int` or `None` (default: `30`)
        If `int` is given, data is smoothed by convolution
        along the x-axis with kernel size n_convolve.
    standard_scale : `int` or `None` (default: `0`)
        Either 0 (rows) or 1 (columns). Whether or not to standardize that dimension
        (each row or column), subtract minimum and divide each by its maximum.
    sort: `bool` (default: `True`)
        Wether to sort the expression values given by xkey.
    colorbar: `bool` or `None` (default: `None`)
        Whether to show colorbar.
    {row,col}_cluster : `bool` or `None`
        If True, cluster the {rows, columns}.
    context : `None`, or one of {paper, notebook, talk, poster}
        A dictionary of parameters or the name of a preconfigured set.
    font_scale : float, optional
        Scaling factor to scale the size of the font elements.
    figsize: tuple (default: `(8,4)`)
        Figure size.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the default
        filename. Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    kwargs:
        Arguments passed to seaborns clustermap,
        e.g., set `yticklabels=True` to display all gene names in all rows.

    Returns
    -------
    If `show==False` a `matplotlib.Axis`
    """

    import seaborn as sns

    var_names = [name for name in var_names if name in adata.var_names]

    tkey, xkey = kwargs.pop("tkey", sortby), kwargs.pop("xkey", layer)
    time = adata.obs[tkey].values
    time = time[np.isfinite(time)]

    X = (
        adata[:, var_names].layers[xkey]
        if xkey in adata.layers.keys()
        else adata[:, var_names].X
    )
    if issparse(X):
        X = X.A
    df = pd.DataFrame(X[np.argsort(time)], columns=var_names)

    if n_convolve is not None:
        weights = np.ones(n_convolve) / n_convolve
        for gene in var_names:
            try:
                df[gene] = np.convolve(df[gene].values, weights, mode="same")
            except Exception:
                pass  # e.g. all-zero counts or nans cannot be convolved

    if sort:
        max_sort = np.argsort(np.argmax(df.values, axis=0))
        df = pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort])
    strings_to_categoricals(adata)

    if col_color is not None:
        col_colors = to_list(col_color)
        col_color = []
        for _, col in enumerate(col_colors):
            if not is_categorical(adata, col):
                obs_col = adata.obs[col]
                cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
                adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
                col += "_categorical"
                set_colors_for_categorical_obs(adata, col, palette)
            col_color.append(interpret_colorkey(adata, col)[np.argsort(time)])

    if "dendrogram_ratio" not in kwargs:
        kwargs["dendrogram_ratio"] = (
            0.1 if row_cluster else 0,
            0.2 if col_cluster else 0,
        )
    if "cbar_pos" not in kwargs or not colorbar:
        kwargs["cbar_pos"] = None

    kwargs.update(
        dict(
            col_colors=col_color,
            col_cluster=col_cluster,
            row_cluster=row_cluster,
            cmap=color_map,
            xticklabels=False,
            standard_scale=standard_scale,
            figsize=figsize,
        )
    )

    args = {}
    if font_scale is not None:
        args = {"font_scale": font_scale}
        context = context or "notebook"

    with sns.plotting_context(context=context, **args):
        try:
            cm = sns.clustermap(df.T, **kwargs)
        except Exception:
            logg.warn("Please upgrade seaborn with `pip install -U seaborn`.")
            kwargs.pop("dendrogram_ratio")
            kwargs.pop("cbar_pos")
            cm = sns.clustermap(df.T, **kwargs)

    savefig_or_show("heatmap", save=save, show=show)
    if show is False:
        return cm
