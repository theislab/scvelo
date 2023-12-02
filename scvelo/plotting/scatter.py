from inspect import signature

import numpy as np
import pandas as pd
from pandas import unique

import matplotlib.pyplot as pl
from matplotlib.colors import is_color_like

from anndata import AnnData

from scvelo import logging as logg
from scvelo import settings
from scvelo.preprocessing.neighbors import get_connectivities
from .docs import doc_params, doc_scatter
from .utils import (
    default_basis,
    default_color,
    default_color_map,
    default_legend_loc,
    default_size,
    default_xkey,
    default_ykey,
    get_ax,
    get_components,
    get_figure_params,
    get_kwargs,
    get_obs_vector,
    get_value_counts,
    gets_vals_from_color_gradients,
    groups_to_bool,
    interpret_colorkey,
    is_categorical,
    is_int,
    is_list,
    is_list_of_int,
    is_list_of_list,
    is_list_of_str,
    make_dense,
    plot_density,
    plot_linfit,
    plot_outline,
    plot_polyfit,
    plot_rug,
    plot_velocity_fits,
    rgb_custom_colormap,
    savefig_or_show,
    set_colorbar,
    set_colors_for_categorical_obs,
    set_label,
    set_legend,
    set_margin,
    set_title,
    to_list,
    to_val,
    to_valid_bases_list,
    update_axes,
)


@doc_params(scatter=doc_scatter)
def scatter(
    adata=None,
    basis=None,
    x=None,
    y=None,
    vkey=None,
    color=None,
    use_raw=None,
    layer=None,
    color_map=None,
    colorbar=None,
    palette=None,
    size=None,
    alpha=None,
    linewidth=None,
    linecolor=None,
    perc=None,
    groups=None,
    sort_order=True,
    components=None,
    projection=None,
    legend_loc=None,
    legend_loc_lines=None,
    legend_fontsize=None,
    legend_fontweight=None,
    legend_fontoutline=None,
    legend_align_text=None,
    xlabel=None,
    ylabel=None,
    title=None,
    fontsize=None,
    figsize=None,
    xlim=None,
    ylim=None,
    add_density=None,
    add_assignments=None,
    add_linfit=None,
    add_polyfit=None,
    add_rug=None,
    add_text=None,
    add_text_pos=None,
    add_margin=None,
    add_outline=None,
    outline_width=None,
    outline_color=None,
    n_convolve=None,
    smooth=None,
    normalize_data=None,
    rescale_color=None,
    color_gradients=None,
    dpi=None,
    frameon=None,
    zorder=None,
    ncols=None,
    nrows=None,
    wspace=None,
    hspace=None,
    show=None,
    save=None,
    ax=None,
    **kwargs,
):
    """Scatter plot along observations or variables axes.

    Arguments:
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate
    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate
    {scatter}

    Returns
    -------
    If `show==False` a `matplotlib.Axis`
    """
    if adata is None and (x is not None and y is not None):
        adata = AnnData(np.stack([x, y]).T)

    # restore old conventions
    add_assignments = kwargs.pop("show_assignments", add_assignments)
    add_linfit = kwargs.pop("show_linear_fit", add_linfit)
    add_polyfit = kwargs.pop("show_polyfit", add_polyfit)
    add_density = kwargs.pop("show_density", add_density)
    add_rug = kwargs.pop("rug", add_rug)
    basis = kwargs.pop("var_names", basis)

    # keys for figures (fkeys) and multiple plots (mkeys)
    fkeys = ["adata", "show", "save", "groups", "ncols", "nrows", "wspace", "hspace"]
    fkeys += ["add_margin", "ax", "kwargs"]
    mkeys = ["color", "layer", "basis", "components", "x", "y", "xlabel", "ylabel"]
    mkeys += ["title", "color_map", "add_text"]
    scatter_kwargs = {"show": False, "save": False}
    for key in signature(scatter).parameters:
        if key not in mkeys + fkeys:
            scatter_kwargs[key] = eval(key)
    mkwargs = {}
    for key in mkeys:  # mkwargs[key] = key for key in mkeys
        mkwargs[key] = eval("{0}[0] if is_list({0}) else {0}".format(key))

    # use c & color and cmap & color_map interchangeably,
    # and plot each group separately if groups is 'all'
    if "c" in kwargs:
        color = kwargs.pop("c")
    if "cmap" in kwargs:
        color_map = kwargs.pop("cmap")
    if "rasterized" not in kwargs:
        kwargs["rasterized"] = settings._vector_friendly
    if isinstance(color_map, (list, tuple)) and all(
        [is_color_like(c) or c == "transparent" for c in color_map]
    ):
        color_map = rgb_custom_colormap(colors=color_map)
    if isinstance(groups, str) and groups == "all":
        if color is None:
            color = default_color(adata)
        if is_categorical(adata, color):
            vc = adata.obs[color].value_counts()
            groups = [[c] for c in vc[vc > 0].index]
    if isinstance(add_text, (list, tuple, np.ndarray, np.record)):
        add_text = list(np.array(add_text, dtype=str))

    # create list of each mkey and check if all bases are valid.
    color = to_list(color, max_len=None)
    layer, components = to_list(layer), to_list(components)
    x, y, basis = to_list(x), to_list(y), to_valid_bases_list(adata, basis)

    # get multikey (with more than one element)
    multikeys = eval(f"[{','.join(mkeys)}]")
    if is_list_of_list(groups):
        multikeys.append(groups)
    key_lengths = np.array([len(key) if is_list(key) else 1 for key in multikeys])
    multikey = (
        multikeys[np.where(key_lengths > 1)[0][0]] if np.max(key_lengths) > 1 else None
    )

    # gridspec frame for plotting multiple keys (mkeys: list or tuple)
    if multikey is not None:
        if np.sum(key_lengths > 1) == 1 and is_list_of_str(multikey):
            multikey = unique(multikey)  # take unique set if no more than one multikey
        if len(multikey) > 20:
            raise ValueError("Please restrict the passed list to max 20 elements.")
        if ax is not None:
            logg.warn("Cannot specify `ax` when plotting multiple panels.")
        if is_list(title):
            title *= int(np.ceil(len(multikey) / len(title)))
        if nrows is None:
            ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
            nrows = int(np.ceil(len(multikey) / ncols))
        else:
            ncols = int(np.ceil(len(multikey) / nrows))
        if not frameon or frameon == "artist":
            lloc, llines = "legend_loc", "legend_loc_lines"
            if lloc in scatter_kwargs and scatter_kwargs[lloc] is None:
                scatter_kwargs[lloc] = "none"
            if llines in scatter_kwargs and scatter_kwargs[llines] is None:
                scatter_kwargs[llines] = "none"

        grid_figsize, dpi = get_figure_params(figsize, dpi, ncols)
        grid_figsize = (grid_figsize[0] * ncols, grid_figsize[1] * nrows)
        fig = pl.figure(None, grid_figsize, dpi=dpi)
        hspace = 0.3 if hspace is None else hspace
        gspec = pl.GridSpec(nrows, ncols, fig, hspace=hspace, wspace=wspace)

        ax = []
        for i, gs in enumerate(gspec):
            if i < len(multikey):
                g = groups[i * (len(groups) > i)] if is_list_of_list(groups) else groups
                multi_kwargs = {"groups": g}
                for key in mkeys:  # multi_kwargs[key] = key[i] if is multikey else key
                    multi_kwargs[key] = eval(
                        "{0}[i * (len({0}) > i)] if is_list({0}) else {0}".format(key)
                    )
                ax.append(
                    scatter(
                        adata,
                        ax=pl.subplot(gs),
                        **multi_kwargs,
                        **scatter_kwargs,
                        **kwargs,
                    )
                )

        if not frameon and isinstance(ylabel, str):
            set_label(xlabel, ylabel, fontsize, ax=ax[0], fontweight="bold")
        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False:
            return ax

    else:
        # make sure that there are no more lists, e.g. ['clusters'] becomes 'clusters'
        color_map = to_val(color_map)
        color, layer, basis = to_val(color), to_val(layer), to_val(basis)
        x, y, components = to_val(x), to_val(y), to_val(components)
        xlabel, ylabel, title = to_val(xlabel), to_val(ylabel), to_val(title)

        # multiple plots within one ax for comma-separated y or layers (string).

        if any([isinstance(key, str) and "," in key for key in [y, layer]]):
            # comma split
            y, layer, color = (
                [k.strip() for k in key.split(",")]
                if isinstance(key, str) and "," in key
                else to_list(key)
                for key in [y, layer, color]
            )
            multikey = y if len(y) > 1 else layer if len(layer) > 1 else None

            if multikey is not None:
                for i in multikey:
                    ax = scatter(
                        adata,
                        x=x,
                        y=y[i * (len(y) > i)],
                        color=color[i * (len(color) > i)],
                        layer=layer[i * (len(layer) > i)],
                        basis=basis,
                        components=components,
                        groups=groups,
                        xlabel=xlabel,
                        ylabel="expression" if ylabel is None else ylabel,
                        color_map=color_map,
                        title=y[i * (len(y) > i)] if title is None else title,
                        ax=ax,
                        **scatter_kwargs,
                    )
                if legend_loc is None:
                    legend_loc = "best"
                if legend_loc and legend_loc != "none":
                    multikey = [key.replace("Ms", "spliced") for key in multikey]
                    multikey = [key.replace("Mu", "unspliced") for key in multikey]
                    ax.legend(multikey, fontsize=legend_fontsize, loc=legend_loc)

                savefig_or_show(dpi=dpi, save=save, show=show)
                if show is False:
                    return ax

        elif color_gradients is not None and color_gradients is not False:
            vals, names, color, scatter_kwargs = gets_vals_from_color_gradients(
                adata, color, **scatter_kwargs
            )
            cols = zip(adata.obs[color].cat.categories, adata.uns[f"{color}_colors"])
            c_colors = {cat: col for (cat, col) in cols}
            mkwargs.pop("color")
            ax = scatter(
                adata,
                color="grey",
                ax=ax,
                **mkwargs,
                **get_kwargs(scatter_kwargs, {"alpha": 0.05}),
            )  # background
            ax = scatter(
                adata,
                color=color,
                ax=ax,
                **mkwargs,
                **get_kwargs(scatter_kwargs, {"s": 0}),
            )  # set legend
            sorted_idx = np.argsort(vals, 1)[:, ::-1][:, :2]
            for id0 in range(len(names)):
                for id1 in range(id0 + 1, len(names)):
                    cmap = rgb_custom_colormap(
                        [c_colors[names[id0]], "white", c_colors[names[id1]]],
                        alpha=[1, 0, 1],
                    )
                    mkwargs.update({"color_map": cmap})
                    c_vals = np.array(vals[:, id1] - vals[:, id0]).flatten()
                    c_bool = np.array([id0 in c and id1 in c for c in sorted_idx])
                    if np.sum(c_bool) > 1:
                        _adata = adata[c_bool] if np.sum(~c_bool) > 0 else adata
                        mkwargs["color"] = c_vals[c_bool]
                        ax = scatter(
                            _adata, ax=ax, **mkwargs, **scatter_kwargs, **kwargs
                        )
            savefig_or_show(dpi=dpi, save=save, show=show)
            if show is False:
                return ax

        # actual scatter plot
        else:
            # set color, color_map, edgecolor, basis, linewidth, frameon, use_raw
            if color is None:
                color = default_color(adata, add_outline)
            if "cmap" not in kwargs:
                kwargs["cmap"] = (
                    default_color_map(adata, color) if color_map is None else color_map
                )
            if "s" not in kwargs:
                kwargs["s"] = default_size(adata) if size is None else size
            if "edgecolor" not in kwargs:
                kwargs["edgecolor"] = "none"
            is_embedding = ((x is None) | (y is None)) and basis not in adata.var_names
            if basis is None and is_embedding:
                basis = default_basis(adata)
            if linewidth is None:
                linewidth = 1
            if frameon is None:
                frameon = True if not is_embedding else settings._frameon
            if isinstance(groups, str):
                groups = [groups]
            if use_raw is None and basis not in adata.var_names:
                use_raw = layer is None and adata.raw is not None

            ax, show = get_ax(ax, show, figsize, dpi, projection)

            # phase portrait: get x and y from .layers (e.g. spliced vs. unspliced)
            if basis in adata.var_names:
                if title is None:
                    title = basis
                if x is None and y is None:
                    x = default_xkey(adata, use_raw=use_raw)
                    y = default_ykey(adata, use_raw=use_raw)
                elif x is None or y is None:
                    raise ValueError("Both x and y have to specified.")
                if isinstance(x, str) and isinstance(y, str):
                    layers_keys = list(adata.layers.keys()) + ["X"]
                    if any([key not in layers_keys for key in [x, y]]):
                        raise ValueError("Could not find x or y in layers.")

                    if xlabel is None:
                        xlabel = x
                    if ylabel is None:
                        ylabel = y

                    x = get_obs_vector(adata, basis, layer=x, use_raw=use_raw)
                    y = get_obs_vector(adata, basis, layer=y, use_raw=use_raw)

                if legend_loc is None:
                    legend_loc = "none"

                if use_raw and perc is not None:
                    ub = np.percentile(x, 99.9 if not isinstance(perc, int) else perc)
                    ax.set_xlim(right=ub * 1.05)
                    ub = np.percentile(y, 99.9 if not isinstance(perc, int) else perc)
                    ax.set_ylim(top=ub * 1.05)

                # velocity model fits (full dynamics and steady-state ratios)
                if any(["gamma" in key or "alpha" in key for key in adata.var.keys()]):
                    plot_velocity_fits(
                        adata,
                        basis,
                        vkey,
                        use_raw,
                        linewidth,
                        linecolor,
                        legend_loc_lines,
                        legend_fontsize,
                        add_assignments,
                        ax=ax,
                    )

            # embedding: set x and y to embedding coordinates
            elif is_embedding:
                X_emb = adata.obsm[f"X_{basis}"][:, get_components(components, basis)]
                x, y = X_emb[:, 0], X_emb[:, 1]
                # todo: 3d plotting
                # z = X_emb[:, 2] if projection == "3d" and X_emb.shape[1] > 2 else None

            elif isinstance(x, str) and isinstance(y, str):
                var_names = (
                    adata.raw.var_names
                    if use_raw and adata.raw is not None
                    else adata.var_names
                )
                if layer is None:
                    layer = default_xkey(adata, use_raw=use_raw)
                x_keys = list(adata.obs.keys()) + list(adata.layers.keys())
                is_timeseries = y in var_names and x in x_keys
                if xlabel is None:
                    xlabel = x
                if ylabel is None:
                    ylabel = layer if is_timeseries else y
                if title is None:
                    title = y if is_timeseries else color
                if legend_loc is None:
                    legend_loc = "none"

                # gene trend: x and y as gene along obs/layers (e.g. pseudotime)
                if is_timeseries:
                    x = (
                        adata.obs[x]
                        if x in adata.obs.keys()
                        else adata.obs_vector(y, layer=x)
                    )
                    y = get_obs_vector(adata, basis=y, layer=layer, use_raw=use_raw)
                # get x and y from var_names, var or obs
                else:
                    if x in var_names and y in var_names:
                        if layer in adata.layers.keys():
                            x = adata.obs_vector(x, layer=layer)
                            y = adata.obs_vector(y, layer=layer)
                        else:
                            data = adata.raw if use_raw else adata
                            x, y = data.obs_vector(x), data.obs_vector(y)
                    elif x in adata.var.keys() and y in adata.var.keys():
                        x, y = adata.var[x], adata.var[y]
                    elif x in adata.obs.keys() and y in adata.obs.keys():
                        x, y = adata.obs[x], adata.obs[y]
                    elif np.any(
                        [var_key in x or var_key in y for var_key in adata.var.keys()]
                    ):
                        var_keys = [
                            k
                            for k in adata.var.keys()
                            if not isinstance(adata.var[k][0], str)
                        ]
                        var = adata.var[var_keys]
                        x = var.astype(np.float32).eval(x)
                        y = var.astype(np.float32).eval(y)
                    elif np.any(
                        [obs_key in x or obs_key in y for obs_key in adata.obs.keys()]
                    ):
                        obs_keys = [
                            k
                            for k in adata.obs.keys()
                            if not isinstance(adata.obs[k][0], str)
                        ]
                        obs = adata.obs[obs_keys]
                        x = obs.astype(np.float32).eval(x)
                        y = obs.astype(np.float32).eval(y)
                    else:
                        raise ValueError(
                            "x or y is invalid! pass valid observation or a gene name"
                        )

            x, y = make_dense(x).flatten(), make_dense(y).flatten()

            # convolve along x axes (e.g. pseudotime)
            if n_convolve is not None:
                vec_conv = np.ones(n_convolve) / n_convolve
                y[np.argsort(x)] = np.convolve(y[np.argsort(x)], vec_conv, mode="same")

            # if color is set to a cell index, plot that cell on top
            if is_int(color) or is_list_of_int(color) and len(color) != len(x):
                color = np.array(np.isin(np.arange(len(x)), color), dtype=bool)
                size = kwargs["s"] * 2 if np.sum(color) == 1 else kwargs["s"]
                if zorder is None:
                    zorder = 10
                ax.scatter(
                    np.ravel(x[color]),
                    np.ravel(y[color]),
                    s=size,
                    zorder=zorder,
                    color=palette[-1] if palette is not None else "darkblue",
                )
                color = (
                    palette[0] if palette is not None and len(palette) > 1 else "gold"
                )
                zorder -= 1

            # if color is in {'ascending', 'descending'}
            elif isinstance(color, str):
                if color == "ascending":
                    color = np.linspace(0, 1, len(x))
                elif color == "descending":
                    color = np.linspace(1, 0, len(x))

            # set palette if categorical color vals
            if is_categorical(adata, color):
                set_colors_for_categorical_obs(adata, color, palette)

            # set color
            if (
                basis in adata.var_names
                and isinstance(color, str)
                and color in adata.layers.keys()
            ):
                # phase portrait: color=basis, layer=color
                c = interpret_colorkey(adata, basis, color, perc, use_raw)
            else:
                # embedding, gene trend etc.
                c = interpret_colorkey(adata, color, layer, perc, use_raw)

            if c is not None and not isinstance(c, str) and not isinstance(c[0], str):
                # smooth color values across neighbors and rescale
                if smooth and len(c) == adata.n_obs:
                    n_neighbors = None if isinstance(smooth, bool) else smooth
                    c = get_connectivities(adata, n_neighbors=n_neighbors).dot(c)
                # rescale color values to min and max acc. to rescale_color tuple
                if rescale_color is not None:
                    # TODO: Use type checking instead
                    try:
                        c += rescale_color[0] - np.nanmin(c)
                        c *= rescale_color[1] / np.nanmax(c)
                    except TypeError:
                        logg.warn("Could not rescale colors. Pass a tuple, e.g. [0,1].")

            # set vmid to 0 if color values obtained from velocity expression
            if not np.any([v in kwargs for v in ["vmin", "vmid", "vmax"]]) and np.any(
                [
                    isinstance(v, str)
                    and "time" not in v
                    and (v.endswith("velocity") or v.endswith("transition"))
                    for v in [color, layer]
                ]
            ):
                kwargs["vmid"] = 0

            # introduce vmid by setting vmin and vmax accordingly
            if "vmid" in kwargs:
                vmid = kwargs.pop("vmid")
                if vmid is not None:
                    if not (isinstance(c, str) or isinstance(c[0], str)):
                        lb, ub = np.min(c), np.max(c)
                        crange = max(np.abs(vmid - lb), np.abs(ub - vmid))
                        kwargs.update({"vmin": vmid - crange, "vmax": vmid + crange})

            x, y = np.ravel(x), np.ravel(y)
            if len(x) != len(y):
                raise ValueError("x or y do not share the same dimension.")

            if normalize_data:
                x = (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))
                y = (y - np.nanmin(x)) / (np.nanmax(y) - np.nanmin(y))

            if not isinstance(c, str):
                c = np.ravel(c) if len(np.ravel(c)) == len(x) else c

            # store original order of color values
            color_array, scatter_array = c, np.stack([x, y]).T

            # set color to grey for NAN values and for cells that are not in groups
            if (
                groups is not None
                or is_categorical(adata, color)
                and np.any(pd.isnull(adata.obs[color]))
            ):
                if isinstance(groups, (list, tuple, np.record)):
                    groups = unique(groups)
                zorder = 0 if zorder is None else zorder
                pop_keys = ["groups", "add_linfit", "add_polyfit", "add_density"]
                _ = [scatter_kwargs.pop(key, None) for key in pop_keys]
                ax = scatter(
                    adata,
                    x=x,
                    y=y,
                    basis=basis,
                    layer=layer,
                    color="lightgrey",
                    ax=ax,
                    **scatter_kwargs,
                )
                if groups is not None and len(groups) == 1:
                    if (
                        isinstance(groups[0], str)
                        and groups[0] in adata.var.keys()
                        and basis in adata.var_names
                    ):
                        groups = f"{adata[:, basis].var[groups[0]][0]}"
                idx = groups_to_bool(adata, groups, color)
                if idx is not None:
                    if np.sum(idx) > 0:  # if any group to be highlighted
                        x, y = x[idx], y[idx]
                        if not isinstance(c, str) and len(c) == adata.n_obs:
                            c = c[idx]
                        if isinstance(kwargs["s"], np.ndarray):
                            kwargs["s"] = np.array(kwargs["s"])[idx]
                        if (
                            title is None
                            and groups is not None
                            and len(groups) == 1
                            and isinstance(groups[0], str)
                        ):
                            title = groups[0]
                    else:  # if nothing to be highlighted
                        add_linfit, add_polyfit, add_density = None, None, None
            else:
                idx = None

            if not isinstance(c, str) and len(c) != len(x):
                c = "grey"
                if not isinstance(color, str) or color != default_color(adata):
                    logg.warn("Invalid color key. Using grey instead.")

            # check if higher value points should be plotted on top
            if not isinstance(c, str) and len(c) == len(x):
                order = None
                if sort_order and not is_categorical(adata, color):
                    order = np.argsort(c)
                elif not sort_order and is_categorical(adata, color):
                    counts = get_value_counts(
                        adata[idx] if idx is not None else adata, color
                    )
                    np.random.seed(0)
                    nums, p = np.arange(0, len(x)), counts / np.sum(counts)
                    order = np.random.choice(nums, len(x), replace=False, p=p)
                if order is not None:
                    x, y, c = x[order], y[order], c[order]
                    if isinstance(kwargs["s"], np.ndarray):  # sort sizes if array-type
                        kwargs["s"] = np.array(kwargs["s"])[order]

            marker = kwargs.pop("marker", ".")
            smp = ax.scatter(
                x, y, c=c, alpha=alpha, marker=marker, zorder=zorder, **kwargs
            )

            outline_dtypes = (list, tuple, np.ndarray, int, np.int_, str)
            if isinstance(add_outline, outline_dtypes) or add_outline:
                if isinstance(add_outline, (list, tuple, np.record)):
                    add_outline = unique(add_outline)
                if (
                    add_outline is not True
                    and isinstance(add_outline, (int, np.int_))
                    or is_list_of_int(add_outline)
                    and len(add_outline) != len(x)
                ):
                    add_outline = np.isin(np.arange(len(x)), add_outline)
                    add_outline = np.array(add_outline, dtype=bool)
                    if outline_width is None:
                        outline_width = (0.6, 0.3)
                if isinstance(add_outline, str):
                    if add_outline in adata.var.keys() and basis in adata.var_names:
                        add_outline = f"{adata[:, basis].var[add_outline][0]}"
                idx = groups_to_bool(adata, add_outline, color)
                if idx is not None and np.sum(idx) > 0:  # if anything to be outlined
                    zorder = 2 if zorder is None else zorder + 2
                    if kwargs["s"] is not None:
                        kwargs["s"] *= 1.2
                    # restore order of values
                    x, y = scatter_array[:, 0][idx], scatter_array[:, 1][idx]
                    c = color_array
                    if not isinstance(c, str) and len(c) == adata.n_obs:
                        c = c[idx]
                    if isinstance(kwargs["s"], np.ndarray):
                        kwargs["s"] = np.array(kwargs["s"])[idx]
                    if isinstance(c, np.ndarray) and not isinstance(c[0], str):
                        if "vmid" not in kwargs and "vmin" not in kwargs:
                            kwargs["vmin"] = np.min(color_array)
                        if "vmid" not in kwargs and "vmax" not in kwargs:
                            kwargs["vmax"] = np.max(color_array)
                    ax.scatter(
                        x, y, c=c, alpha=alpha, marker=".", zorder=zorder, **kwargs
                    )
                if idx is None or np.sum(idx) > 0:  # if all or anything to be outlined
                    plot_outline(
                        x, y, kwargs, outline_width, outline_color, zorder, ax=ax
                    )
                if idx is not None and np.sum(idx) == 0:  # if nothing to be outlined
                    add_linfit, add_polyfit, add_density = None, None, None

            # set legend if categorical categorical color vals
            if is_categorical(adata, color) and len(scatter_array) == adata.n_obs:
                legend_loc = default_legend_loc(adata, color, legend_loc)
                g_bool = groups_to_bool(adata, add_outline, color)
                if not (add_outline is None or g_bool is None):
                    groups = add_outline
                set_legend(
                    adata,
                    ax,
                    color,
                    legend_loc,
                    scatter_array,
                    legend_fontweight,
                    legend_fontsize,
                    legend_fontoutline,
                    legend_align_text,
                    groups,
                )
            if add_density:
                plot_density(x, y, add_density, ax=ax)

            if add_linfit:
                if add_linfit is True and basis in adata.var_names:
                    add_linfit = "no_intercept"  # without intercept
                plot_linfit(
                    x,
                    y,
                    add_linfit,
                    legend_loc != "none",
                    linecolor,
                    linewidth,
                    fontsize,
                    ax=ax,
                )

            if add_polyfit:
                if add_polyfit is True and basis in adata.var_names:
                    add_polyfit = "no_intercept"  # without intercept
                plot_polyfit(
                    x,
                    y,
                    add_polyfit,
                    legend_loc != "none",
                    linecolor,
                    linewidth,
                    fontsize,
                    ax=ax,
                )

            if add_rug:
                rug_color = add_rug if isinstance(add_rug, str) else color
                rug_color = np.ravel(interpret_colorkey(adata, rug_color))
                plot_rug(np.ravel(x), color=rug_color, ax=ax)

            if add_text:
                if add_text_pos is None:
                    add_text_pos = [0.05, 0.95]
                ax.text(
                    add_text_pos[0],
                    add_text_pos[1],
                    f"{add_text}",
                    ha="left",
                    va="top",
                    fontsize=fontsize,
                    transform=ax.transAxes,
                    bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.2},
                )

            set_label(xlabel, ylabel, fontsize, basis, ax=ax)
            set_title(title, layer, color, fontsize, ax=ax)
            update_axes(ax, xlim, ylim, fontsize, is_embedding, frameon, figsize)
            if add_margin:
                set_margin(ax, x, y, add_margin)
            if colorbar is not False:
                if not isinstance(c, str) and not is_categorical(adata, color):
                    labelsize = fontsize * 0.75 if fontsize is not None else None
                    set_colorbar(smp, ax=ax, labelsize=labelsize)

            savefig_or_show(dpi=dpi, save=save, show=show)
            if show is False:
                return ax


def _wraps_plot_scatter(wrapper):
    annots_orig = {
        k: v for k, v in wrapper.__annotations__.items() if k not in {"adata", "kwargs"}
    }
    annots_scatter = {k: v for k, v in scatter.__annotations__.items() if k != "basis"}
    wrapper.__annotations__ = {**annots_scatter, **annots_orig}
    wrapper.__wrapped__ = scatter
    return wrapper


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def trimap(adata, **kwargs):
    """Scatter plot in trimap basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="trimap", **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def umap(adata, **kwargs):
    """Scatter plot in UMAP basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="umap", **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def tsne(adata, **kwargs):
    """Scatter plot in tsne basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="tsne", **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def diffmap(adata, **kwargs):
    """Scatter plot in diffmap basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="diffmap", **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def phate(adata, **kwargs):
    """Scatter plot in phate basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="phate", **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def draw_graph(adata, layout=None, **kwargs):
    """Scatter plot in draw_graph basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    if layout is None:
        layout = f"{adata.uns['draw_graph']['params']['layout']}"
    basis = f"draw_graph_{layout}"
    if f"X_{basis}" not in adata.obsm_keys():
        raise ValueError(f"Could not find draw_graph_{layout} in adata.obs.")
    return scatter(adata, basis=basis, **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def pca(adata, **kwargs):
    """Scatter plot in pca basis.

    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis="pca", **kwargs)
