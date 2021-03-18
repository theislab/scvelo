import numpy as np
import pandas as pd
from pandas import isnull, unique
from scipy.sparse import issparse

import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib.colors import ColorConverter

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
    figsize=(8, 4),
    font_scale=None,
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
    {row,col}_cluster : bool, optional
        If True, cluster the {rows, columns}.
    figsize: tuple (default: `(7,5)`)
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

    if font_scale is not None:
        sns.set(font_scale=font_scale)
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


def heatmap_deprecated(
    adata,
    var_names,
    groups=None,
    groupby=None,
    annotations=None,
    use_raw=False,
    layers=None,
    color_map=None,
    color_map_anno=None,
    colorbar=True,
    row_width=None,
    xlabel=None,
    title=None,
    figsize=None,
    dpi=None,
    show=True,
    save=None,
    ax=None,
    **kwargs,
):

    """\
    Plot pseudotimeseries for genes as heatmap.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names: `str`,  list of `str`
        Names of variables to use for the plot.
    groups: `str`,  list of `str` or `None` (default: `None`)
        Groups selected to plot. Must be an element of adata.obs[groupby].
    groupby: `str` or `None` (default: `None`)
        Key in adata.obs. Indicates how to group the plot.
    annotations: `str`,  list of `str` or `None` (default: `None`)
        Key in adata.obs. Annotations are plotted in the last row.
    use_raw: `bool` (default: `False`)
        If true, moments are used instead of raw data.
    layers: `str`,  list of `str` or `None` (default: `['X']`)
        Selected layers.
    color_map: `str`,  list of `str` or `None` (default: `None`)
        String denoting matplotlib color map for the heat map.
        There must be one list entry for each layer.
    color_map_anno: `str`,  list of `str` or `None` (default: `None`)
        String denoting matplotlib color map for the annotations.
        There must be one list entry for each annotation.
    colorbar: `bool` (default: `True`)
        If True, a colormap for each layer is added on the right bottom corner.
    row_width: `float` (default: `None`)
        Constant width of all rows.
    xlabel:
        Label for the x-axis.
    title: `str` or `None` (default: `None`)
        Main plot title.
    figsize: tuple (default: `(7,5)`)
        Figure size.
    dpi: `int` (default: 80)
        Figure dpi.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure.
        Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    ax: `matplotlib.Axes`, optional (default: `None`)
        A matplotlib axes object. Only works if plotting a single component.

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """

    # catch
    if "velocity_pseudotime" not in adata.obs.keys():
        raise ValueError(
            "A function requires computation of the pseudotime"
            "for ordering at single-cell resolution"
        )
    if layers is None:
        layers = ["X"]
    if annotations is None:
        annotations = []
    if isinstance(var_names, str):
        var_names = [var_names]
    if len(var_names) == 0:
        var_names = np.arange(adata.X.shape[1])
    if var_names.ndim == 2:
        var_names = var_names[:, 0]
    var_names = [name for name in var_names if name in adata.var_names]
    if len(var_names) == 0:
        raise ValueError(
            "The specified var_names are all not" "contained in the adata.var_names."
        )

    if layers is None:
        layers = ["X"]
    if isinstance(layers, str):
        layers = [layers]
    layers = [layer for layer in layers if layer in adata.layers.keys() or layer == "X"]
    if len(layers) == 0:
        raise ValueError(
            "The selected layers are not contained" "in adata.layers.keys()."
        )
    if not use_raw:
        layers = np.array(layers)
        if "X" in layers:
            layers[np.array([layer == "X" for layer in layers])] = "Ms"
        if "spliced" in layers:
            layers[np.array([layer == "spliced" for layer in layers])] = "Ms"
        if "unspliced" in layers:
            layers[np.array([layer == "unspliced" for layer in layers])] = "Ms"
        layers = list(layers)
    if "Ms" in layers and "Ms" not in adata.layers.keys():
        raise ValueError(
            "Moments have to be computed before" "using this plot function."
        )
    if "Mu" in layers and "Mu" not in adata.layers.keys():
        raise ValueError(
            "Moments have to be computed before" "using this plot function."
        )
    layers = unique(layers)

    # Number of rows to plot
    tot_len = len(var_names) * len(layers) + len(annotations)

    # init main figure
    figsize = rcParams["figure.figsize"] if figsize is None else figsize
    if row_width is not None:
        figsize[1] = row_width * tot_len
    ax = pl.figure(figsize=figsize, dpi=dpi).gca() if ax is None else ax
    ax.set_yticks([])
    ax.set_xticks([])

    # groups bar
    ax_bounds = ax.get_position().bounds
    if groupby is not None:
        # catch
        if groupby not in adata.obs_keys():
            raise ValueError(
                "The selected groupby is not contained" "in adata.obs_keys()."
            )
        if groups is None:  # Then use everything of that obs
            groups = unique(adata.obs.clusters.values)

        imlist = []

        for igroup, group in enumerate(groups):
            for ivar, var in enumerate(var_names):
                for ilayer, layer in enumerate(layers):
                    groups_axis = pl.axes(
                        [
                            ax_bounds[0] + igroup * ax_bounds[2] / len(groups),
                            ax_bounds[1]
                            + ax_bounds[3]
                            * (tot_len - ivar * len(layers) - ilayer - 1)
                            / tot_len,
                            ax_bounds[2] / len(groups),
                            (ax_bounds[3] - ax_bounds[3] / tot_len * len(annotations))
                            / (len(var_names) * len(layers)),
                        ]
                    )

                    # Get data to fill and reshape
                    dat = adata[:, var]

                    idx_group = [adata.obs[groupby] == group]
                    idx_group = np.array(idx_group[0].tolist())

                    idx_pt = np.array(adata.obs.velocity_pseudotime).argsort()
                    idx_pt = idx_pt[
                        np.array(
                            isnull(np.array(dat.obs.velocity_pseudotime)[idx_pt])
                            is False
                        )
                    ]

                    if layer == "X":
                        laydat = dat.X
                    else:
                        laydat = dat.layers[layer]

                    t1, t3 = idx_group, idx_pt
                    t1 = t1[t3]
                    # idx_var = [vn in var_names for vn in adata.var_names]
                    # t2 = idx_var
                    # laydat = laydat[:, t2]  # select vars
                    laydat = laydat[t3]
                    laydat = laydat[t1]  # select ordered groups

                    if issparse(laydat):
                        laydat = laydat.A

                    # transpose X for ordering in direction var_names: up->downwards
                    laydat = laydat.T[::-1]
                    laydat = laydat.reshape((1, len(laydat)))  # ensure 1dimty

                    # plot
                    im = groups_axis.imshow(
                        laydat,
                        aspect="auto",
                        interpolation="nearest",
                        cmap=color_map[ilayer],
                    )

                    # Frames
                    if ilayer == 0:
                        groups_axis.spines["bottom"].set_visible(False)
                    elif ilayer == len(layer) - 1:
                        groups_axis.spines["top"].set_visible(False)
                    else:
                        groups_axis.spines["top"].set_visible(False)
                        groups_axis.spines["bottom"].set_visible(False)

                    # Further visuals
                    if igroup == 0:
                        if colorbar:
                            if len(layers) % 2 == 0:
                                if ilayer == len(layers) / 2 - 1:
                                    pl.yticks([0.5], [var])
                                else:
                                    groups_axis.set_yticks([])
                            else:
                                if ilayer == (len(layers) - 1) / 2:
                                    pl.yticks([0], [var])
                                else:
                                    groups_axis.set_yticks([])
                        else:
                            pl.yticks([0], [f"{layer} {var}"])
                    else:
                        groups_axis.set_yticks([])

                    groups_axis.set_xticks([])
                    if ilayer == 0 and ivar == 0:
                        groups_axis.set_title(f"{group}")
                    groups_axis.grid(False)

                    # handle needed as mappable for colorbar
                    if igroup == len(groups) - 1:
                        imlist.append(im)

            # further annotations for each group
            if annotations is not None:
                for ianno, anno in enumerate(annotations):
                    anno_axis = pl.axes(
                        [
                            ax_bounds[0] + igroup * ax_bounds[2] / len(groups),
                            ax_bounds[1]
                            + ax_bounds[3] / tot_len * (len(annotations) - ianno - 1),
                            ax_bounds[2] / len(groups),
                            ax_bounds[3] / tot_len,
                        ]
                    )
                    if is_categorical(adata, anno):
                        colo = interpret_colorkey(adata, anno)[t3][t1]
                        colo.reshape(1, len(colo))
                        mapper = np.vectorize(ColorConverter.to_rgb)
                        a = mapper(colo)
                        a = np.array(a).T
                        Y = a.reshape(1, len(colo), 3)
                    else:
                        Y = np.array(interpret_colorkey(adata, anno))[t3][t1]
                        Y = Y.reshape(1, len(Y))
                    _ = anno_axis.imshow(
                        Y, aspect="auto", interpolation="nearest", cmap=color_map_anno
                    )
                    if igroup == 0:
                        anno_axis.set_yticklabels(
                            ["", anno, ""]
                        )  # , fontsize=ytick_fontsize)
                        anno_axis.tick_params(axis="both", which="both", length=0)
                    else:
                        anno_axis.set_yticklabels([])
                        anno_axis.set_yticks([])
                    anno_axis.set_xticks([])
                    anno_axis.set_xticklabels([])
                    anno_axis.grid(False)
                    pl.ylim([0.5, -0.5])  # center ticks

    else:  # groupby is False
        imlist = []
        for ivar, var in enumerate(var_names):
            for ilayer, layer in enumerate(layers):
                ax_bounds = ax.get_position().bounds
                groups_axis = pl.axes(
                    [
                        ax_bounds[0],
                        ax_bounds[1]
                        + ax_bounds[3]
                        * (tot_len - ivar * len(layers) - ilayer - 1)
                        / tot_len,
                        ax_bounds[2],
                        (ax_bounds[3] - ax_bounds[3] / tot_len * len(annotations))
                        / (len(var_names) * len(layers)),
                    ]
                )
                # Get data to fill
                dat = adata[:, var]
                idx = np.array(dat.obs.velocity_pseudotime).argsort()
                idx = idx[
                    np.array(
                        isnull(np.array(dat.obs.velocity_pseudotime)[idx]) is False
                    )
                ]

                if layer == "X":
                    laydat = dat.X
                else:
                    laydat = dat.layers[layer]
                laydat = laydat[idx]
                if issparse(laydat):
                    laydat = laydat.A

                # transpose X for ordering in direction var_names: up->downwards
                laydat = laydat.T[::-1]
                laydat = laydat.reshape((1, len(laydat)))

                # plot
                im = groups_axis.imshow(
                    laydat,
                    aspect="auto",
                    interpolation="nearest",
                    cmap=color_map[ilayer],
                )
                imlist.append(im)

                # Frames
                if ilayer == 0:
                    groups_axis.spines["bottom"].set_visible(False)
                elif ilayer == len(layer) - 1:
                    groups_axis.spines["top"].set_visible(False)
                else:
                    groups_axis.spines["top"].set_visible(False)
                    groups_axis.spines["bottom"].set_visible(False)

                # Further visuals
                groups_axis.set_xticks([])
                groups_axis.grid(False)
                pl.ylim([0.5, -0.5])  # center
                if colorbar:
                    if len(layers) % 2 == 0:
                        if ilayer == len(layers) / 2 - 1:
                            pl.yticks([0.5], [var])
                        else:
                            groups_axis.set_yticks([])
                    else:
                        if ilayer == (len(layers) - 1) / 2:
                            pl.yticks([0], [var])
                        else:
                            groups_axis.set_yticks([])
                else:
                    pl.yticks([0], [f"{layer} {var}"])

        # further annotations bars
        if annotations is not None:
            for ianno, anno in enumerate(annotations):
                anno_axis = pl.axes(
                    [
                        ax_bounds[0],
                        ax_bounds[1]
                        + ax_bounds[3] / tot_len * (len(annotations) - ianno - 1),
                        ax_bounds[2],
                        ax_bounds[3] / tot_len,
                    ]
                )
                dat = adata[:, var_names]
                if is_categorical(dat, anno):
                    colo = interpret_colorkey(dat, anno)[idx]
                    colo.reshape(1, len(colo))
                    mapper = np.vectorize(ColorConverter.to_rgb)
                    a = mapper(colo)
                    a = np.array(a).T
                    Y = a.reshape(1, len(idx), 3)
                else:
                    Y = np.array(interpret_colorkey(dat, anno)[idx]).reshape(
                        1, len(idx)
                    )
                _ = anno_axis.imshow(
                    Y, aspect="auto", interpolation="nearest", cmap=color_map_anno
                )

                anno_axis.set_yticklabels(["", anno, ""])  # , fontsize=ytick_fontsize)
                anno_axis.tick_params(axis="both", which="both", length=0)
                anno_axis.grid(False)
                anno_axis.set_xticks([])
                anno_axis.set_xticklabels([])
                pl.ylim([-0.5, +0.5])

    # Colorbar
    if colorbar:
        if len(layers) > 1:
            # I must admit, this part is chaotic
            for ilayer, layer in enumerate(layers):
                w = 0.015 * 10 / figsize[0]  # 0.02 * ax_bounds[2]
                x = ax_bounds[0] + ax_bounds[2] * 0.99 + 1.5 * w + w * 1.2 * ilayer
                y = ax_bounds[1]
                h = ax_bounds[3] * 0.3
                cbaxes = pl.axes([x, y, w, h])
                cb = pl.colorbar(mappable=imlist[ilayer], cax=cbaxes)
                pl.text(
                    x - 40 * w,
                    y + h * 4,
                    layer,
                    rotation=45,
                    horizontalalignment="left",
                    verticalalignment="bottom",
                )
                if ilayer == len(layers) - 1:
                    ext = abs(cb.vmin - cb.vmax)
                    cb.set_ticks([cb.vmin + 0.07 * ext, cb.vmax - 0.07 * ext])
                    cb.ax.set_yticklabels(["Low", "High"])  # vertical colorbar
                else:
                    cb.set_ticks([])
        else:
            cbaxes = pl.axes(
                [
                    ax_bounds[0] + ax_bounds[2] + 0.01,
                    ax_bounds[1],
                    0.02,
                    ax_bounds[3] * 0.3,
                ]
            )
            cb = pl.colorbar(mappable=im, cax=cbaxes)
            cb.set_ticks([cb.vmin, cb.vmax])
            cb.ax.set_yticklabels(["Low", "High"])

    if xlabel is None:
        xlabel = "velocity" + " " + "pseudotime"
    if title is not None:
        ax.set_title(title, pad=30)
    if len(annotations) == 0:
        ax.set_xlabel(xlabel)
        ax.xaxis.labelpad = 20

    # set_label(xlabel, None, fontsize, basis)
    # set_title(title, None, None, fontsize)
    # update_axes(ax, fontsize)

    savefig_or_show("heatmap", dpi=dpi, save=save, show=show)
    if not show:
        return ax
