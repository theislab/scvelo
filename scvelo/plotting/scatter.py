from .. import settings
from .. import logging as logg
from .. import AnnData
from .docs import doc_scatter, doc_params

from .utils import is_categorical, is_list, is_list_of_str, is_list_of_list, to_list, to_valid_bases_list, to_val
from .utils import default_basis, default_color, default_size, default_color_map, default_legend_loc, default_xkey, default_ykey
from .utils import unique, make_dense, get_components, get_connectivities, groups_to_bool, interpret_colorkey, get_obs_vector
from .utils import update_axes, set_label, set_title, set_colorbar, _set_colors_for_categorical_obs, _add_legend
from .utils import plot_linfit, plot_polyfit, plot_density, plot_outline, plot_rug, plot_velocity_fits, savefig_or_show


from inspect import signature
from matplotlib import rcParams, patheffects
from matplotlib.gridspec import SubplotSpec
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd


@doc_params(scatter=doc_scatter)
def scatter(adata=None, x=None, y=None, basis=None, vkey=None, color=None, use_raw=None, layer=None, color_map=None,
            colorbar=None, palette=None, size=None, alpha=None, linewidth=None, linecolor=None, perc=None, groups=None,
            sort_order=True, components=None, projection=None, legend_loc=None, legend_fontsize=None, legend_fontweight=None,
            xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None, xlim=None, ylim=None, show_density=None,
            show_assignments=None, show_linear_fit=None, show_polyfit=None, rug=None, add_outline=None,
            outline_width=None, outline_color=None, n_convolve=None, smooth=None, rescale_color=None, dpi=None,
            frameon=None, zorder=None, ncols=None, wspace=None, hspace=None, show=None, save=None, ax=None, **kwargs):
    """\
    Scatter plot along observations or variables axes.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate
    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    {scatter}

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """
    adata = AnnData(np.stack([x, y]).T) if adata is None and (x is not None and y is not None) else adata

    # keys for figures (fkeys) and multiple plots (mkeys)
    fkeys = ['adata', 'show', 'save', 'groups', 'figsize', 'dpi', 'ncols', 'wspace', 'hspace', 'ax', 'kwargs']
    mkeys = ['color', 'layer', 'basis', 'components', 'x', 'y', 'xlabel', 'ylabel', 'title', 'color_map']
    scatter_kwargs = {'show': False, 'save': False}
    for key in signature(scatter).parameters:
        if key not in mkeys + fkeys: scatter_kwargs[key] = eval(key)

    # use c & color and cmap & color_map interchangeably, and plot each group separately if groups is 'all'
    if 'c' in kwargs: color = kwargs.pop('c')
    if 'cmap' in kwargs: color_map = kwargs.pop('cmap')
    if groups is 'all':
        if color is None:  color = default_color(adata)
        if is_categorical(adata, color): groups = [[c] for c in adata.obs[color].cat.categories]

    # create list of each mkey (won't be needed in the future) and check if all bases are valid.
    color, layer, x, y, components = to_list(color), to_list(layer), to_list(x), to_list(y), to_list(components)
    basis = to_valid_bases_list(adata, basis)

    # get multikey (with more than one element)
    multikeys = eval('[' + ','.join(mkeys) + ']')
    if is_list_of_list(groups): multikeys.append(groups)
    key_lengths = np.array([len(key) if is_list(key) else 1 for key in multikeys])
    multikey = multikeys[np.where(key_lengths > 1)[0][0]] if np.max(key_lengths) > 1 else None

    # gridspec frame for plotting multiple colors, layers and bases, xs, ys, components and groups (lists or tuples)
    if multikey is not None:
        if np.sum(key_lengths > 1) == 1 and is_list_of_str(multikey):
            multikey = unique(multikey)  # take unique set if no more than one multikey
        if len(multikey) > 20:
            raise ValueError('Please restrict the passed list to no more than 20 elements.')
        if ax is not None:
            logg.warn("Cannot specify `ax` when plotting multiple panels.")
        if is_list(title):
            title *= int(np.ceil(len(multikey) / len(title)))

        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        fig = pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi)
        gspec = pl.GridSpec(nrows, ncols, fig, hspace=0.25 if hspace is None else hspace, wspace=wspace)

        ax = []
        for i, gs in enumerate(gspec):
            if i < len(multikey):
                multi_kwargs = {'groups': groups[i * (len(groups) > i)] if is_list_of_list(groups) else groups}
                for key in mkeys:  # multi_kwargs[key] = key[i] if is multikey (list with more than 1 element) else key
                    multi_kwargs[key] = eval('{0}[i * (len({0}) > i)] if is_list({0}) else {0}'.format(key))
                ax.append(scatter(adata, ax=pl.subplot(gs), **multi_kwargs, **scatter_kwargs, **kwargs))

        savefig_or_show(dpi=dpi, save=save, show=show)
        if show is False: return ax

    else:
        # make sure that there are no more lists, e.g. input ['clusters'] becomes 'clusters'
        color_map = to_val(color_map)
        color, layer, basis, components = to_val(color), to_val(layer), to_val(basis), to_val(components)
        x, y, xlabel, ylabel, title = to_val(x), to_val(y), to_val(xlabel), to_val(ylabel), to_val(title)

        # multiple plots within one ax for comma-separated y or layers (string).

        if any([isinstance(key, str) and ',' in key for key in [y, layer]]):
            y, layer, color = [[k.strip() for k in key.split(',')] if isinstance(key, str) and ',' in key
                               else to_list(key) for key in [y, layer, color]]  # comma split
            multikey = y if len(y) > 1 else layer if len(layer) > 1 else None

            if multikey is not None:
                for i, mi in enumerate(multikey):
                    ax = scatter(adata, x=x, y=y[i * (len(y) > i)], color=color[i * (len(color) > i)],
                                 layer=layer[i * (len(layer) > i)], basis=basis, components=components, groups=groups,
                                 xlabel=xlabel, ylabel='expression' if ylabel is None else ylabel, color_map=color_map,
                                 title=y[i * (len(y) > i)] if title is None else title, ax=ax, **scatter_kwargs)
                if legend_loc is not False and legend_loc is not 'none':
                    multikey = [key.replace('Mu', 'unspliced').replace('Ms', 'spliced') for key in multikey]
                    ax.legend(multikey, fontsize=legend_fontsize, loc='best' if legend_loc is None else legend_loc)

                savefig_or_show(dpi=dpi, save=save, show=show)
                if show is False: return ax

        # actual scatter plot
        else:
            # set color, color_map, edgecolor, basis, linewidth, frameon, use_raw, projection
            if color is None:
                color = default_color(adata)
            if 'cmap' not in kwargs:
                kwargs['cmap'] = default_color_map(adata, color) if color_map is None else color_map
            if 's' not in kwargs:
                kwargs['s'] = default_size(adata) if size is None else size
            if 'edgecolor' not in kwargs:
                kwargs['edgecolor'] = 'none'
            is_embedding = ((x is None) | (y is None)) and basis not in adata.var_names
            if basis is None and is_embedding: basis = default_basis(adata)
            if linewidth is None: linewidth = 1
            if linecolor is None: linecolor = 'k'
            if frameon is None: frameon = True if not is_embedding else settings._frameon
            if isinstance(groups, str): groups = [groups]
            if use_raw is None and basis not in adata.var_names:
                use_raw = layer is None and adata.raw is not None
            if projection is '3d':
                from mpl_toolkits.mplot3d import Axes3D
            else:
                projection = None

            if ax is None:
                ax = pl.figure(None, figsize, dpi=dpi).gca(projection=projection)
                if show is None: show = True
            elif isinstance(ax, SubplotSpec):
                geo = ax.get_geometry()
                if show is None: show = geo[-1] + 1 == geo[0] * geo[1]
                ax = pl.subplot(ax)

            # phase portrait: get x and y from .layers (e.g. spliced vs. unspliced) when basis is in var_names
            if basis in adata.var_names:
                if title is None: title = basis
                if x is None and y is None:
                    x, y = default_xkey(adata, use_raw=use_raw), default_ykey(adata, use_raw=use_raw)
                elif x is None or y is None:
                    raise ValueError('Both x and y have to specified.')
                if any([key not in list(adata.layers.keys()) + ['X'] for key in [x, y]]):
                    raise ValueError('Could not find x or y in layers.')

                if xlabel is None: xlabel = x
                if ylabel is None: ylabel = y

                x = get_obs_vector(adata, basis, layer=x, use_raw=use_raw)
                y = get_obs_vector(adata, basis, layer=y, use_raw=use_raw)

                if use_raw and perc is not None:
                    ax.set_xlim(right=np.percentile(x, 99.9 if not isinstance(perc, int) else perc) * 1.05)
                    ax.set_ylim(top=np.percentile(y, 99.9 if not isinstance(perc, int) else perc) * 1.05)

                # velocity model fits (full dynamics and steady-state ratios)
                if any(['gamma' in key or 'alpha' in key for key in adata.var.keys()]):
                    plot_velocity_fits(adata, basis, vkey, use_raw, linewidth, linecolor, legend_loc, legend_fontsize,
                                       show_assignments, ax=ax)

            # embedding: set x and y to embedding coordinates
            elif is_embedding:
                X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)]
                x, y = X_emb[:, 0], X_emb[:, 1]
                z = X_emb[:, 2] if projection is '3d' and X_emb.shape[1] > 2 else None

                # set legend if categorical color vals in embedding
                if is_categorical(adata, color):
                    _set_colors_for_categorical_obs(adata, color, palette)
                    legend_loc = default_legend_loc(adata, color, legend_loc)
                    _add_legend(adata, ax, color, legend_loc, X_emb, legend_fontweight, legend_fontsize,
                                [patheffects.withStroke(linewidth=True, foreground='w')], groups)

            elif isinstance(x, str) and isinstance(y, str):
                if layer is None:
                    layer = default_xkey(adata, use_raw=use_raw)
                is_timeseries = y in adata.var_names and x in list(adata.obs.keys()) + list(adata.layers.keys())
                if xlabel is None: xlabel = x
                if ylabel is None: ylabel = layer if is_timeseries else y
                if title is None: title = y if is_timeseries else color

                # gene trend: get x and y as gene (var_names) along obs/layers (e.g. pseudotime)
                if is_timeseries:
                    x = adata.obs[x] if x in adata.obs.keys() else adata.obs_vector(y, layer=x)
                    y = get_obs_vector(adata, basis=y, layer=layer, use_raw=use_raw)
                # get x and y from var_names, var or obs
                else:
                    if x in adata.var_names and y in adata.var_names:
                        if layer in adata.layers.keys():
                            x = adata.obs_vector(x, layer=layer)
                            y = adata.obs_vector(y, layer=layer)
                        else:
                            x = adata.raw.obs_vector(x) if use_raw else adata.obs_vector(x)
                            y = adata.raw.obs_vector(y) if use_raw else adata.obs_vector(y)
                    elif x in adata.var.keys() and y in adata.var.keys():
                        x, y = adata.var[x], adata.var[y]
                    elif x in adata.obs.keys() and y in adata.obs.keys():
                        x, y = adata.obs[x], adata.obs[y]
                    elif np.any([var_key in x or var_key in y for var_key in adata.var.keys()]):
                        var = adata.var[[key for key in adata.var.keys() if not isinstance(adata.var[key][0], str)]]
                        x = var.astype(np.float32).eval(x)
                        y = var.astype(np.float32).eval(y)
                    elif np.any([obs_key in x or obs_key in y for obs_key in adata.obs.keys()]):
                        obs = adata.obs[[key for key in adata.obs.keys() if not isinstance(adata.obs[key][0], str)]]
                        x = obs.astype(np.float32).eval(x)
                        y = obs.astype(np.float32).eval(y)
                    else:
                        raise ValueError('x or y is invalid! pass valid observation annotation or a gene name')

            x, y = make_dense(x).flatten(), make_dense(y).flatten()

            # convolve along x axes (e.g. pseudotime)
            if n_convolve is not None:
                y[np.argsort(x)] = np.convolve(y[np.argsort(x)], np.ones(n_convolve) / n_convolve, mode='same')

            # if color is set to a cell index, plot that cell on top
            if isinstance(color, int):
                color = np.array(np.arange(len(x)) == color, dtype=bool)
                if zorder is None: zorder = 10
                ax.scatter(np.ravel(x[color]), np.ravel(y[color]), color='darkblue', s=kwargs['s'] * 2, zorder=zorder)
                zorder -= 1

            # set color
            if basis in adata.var_names and isinstance(color, str) and color in adata.layers.keys():
                c = interpret_colorkey(adata, basis, color, perc, use_raw)  # phase portrait: color=basis, layer=color
            else:
                c = interpret_colorkey(adata, color, layer, perc, use_raw)  # embedding, gene trend etc.

            if c is not None and not isinstance(c, str) and not isinstance(c[0], str):
                # smooth color values across neighbors and rescale
                if smooth and len(c) == adata.n_obs:
                    c = get_connectivities(adata, n_neighbors=(None if isinstance(smooth, bool) else smooth)).dot(c)
                # rescale color values to min and max acc. to rescale_color tuple
                if rescale_color is not None:
                    if len(rescale_color) != 2:
                        raise ValueError('rescale_color has to be a tuple with two values, e.g. [0,1].')
                    c += rescale_color[0] - np.min(c)
                    c *= rescale_color[1] / np.max(c)

            # check if higher value points should be plotted on top
            if sort_order and not isinstance(c, str) and not is_categorical(adata, color):
                order = np.argsort(c)
                x, y, c = x[order], y[order], c[order]
                # sort order of size if given as vector
                if isinstance(kwargs['s'], np.ndarray):
                    kwargs['s'] = np.array(kwargs['s'])[order]

            if "vmid" not in kwargs and "vmin" not in kwargs and layer is not None and 'velocity' in layer:
                kwargs['vmid'] = 0  # set vmid to 0 if color values obtained from velocity expression

            # introduce vmid by setting vmin and vmax accordingly
            if "vmid" in kwargs:
                vmid = kwargs.pop("vmid")
                if not isinstance(c, str) and not isinstance(c[0], str) and vmid is not None:
                    lb, ub = np.min(c), np.max(c)
                    crange = max(np.abs(vmid - lb), np.abs(ub - vmid))
                    kwargs.update({"vmin": vmid - crange, "vmax": vmid + crange})

            # set color to grey for NAN values and for cells that are not in groups
            if groups is not None or np.any(pd.isnull(c)):
                zorder = 0 if zorder is None else zorder
                ax = scatter(adata, x=x, y=y, basis=basis, layer=layer, color='lightgrey', ax=ax, groups=None, **scatter_kwargs)
                idx = groups_to_bool(adata, groups, color)
                x, y = x[idx], y[idx]
                if not isinstance(c, str) and len(c) == adata.n_obs:
                    c = c[idx]
                if title is None and groups is not None and len(groups) == 1 and isinstance(groups[0], str):
                    title = groups[0]
                zorder += 1

            x, y = np.ravel(x), np.ravel(y)
            if len(x) != len(y):
                raise ValueError('x or y do not share the same dimension.')

            if not isinstance(c, str):
                c = np.ravel(c)
                if len(c) != len(x):
                    c = 'grey'
                    if color is not default_color(adata):
                        logg.warn('Invalid color key. Using grey instead.')

            smp = ax.scatter(x, y, c=c, alpha=alpha, marker='.', zorder=zorder, **kwargs)

            if add_outline:
                plot_outline(x, y, kwargs, outline_width, outline_color, zorder, ax=ax)

            if show_density:
                plot_density(x, y, show_density, ax=ax)

            if show_linear_fit or show_linear_fit is 0:
                plot_linfit(x, y, show_linear_fit, legend_loc is not 'none', linecolor, linewidth, fontsize, ax=ax)

            if show_polyfit:
                plot_polyfit(x, y, show_polyfit, legend_loc is not 'none', linecolor, linewidth, fontsize, ax=ax)

            if rug:
                plot_rug(np.ravel(x), color=np.ravel(interpret_colorkey(adata, rug)), ax=ax)

            set_label(xlabel, ylabel, fontsize, basis, ax=ax)
            set_title(title, layer, color, fontsize, ax=ax)
            update_axes(ax, xlim, ylim, fontsize, is_embedding, frameon)

            if colorbar is not False and not isinstance(c, str) and not is_categorical(adata, color):
                set_colorbar(smp, ax=ax, labelsize=fontsize * .75 if fontsize is not None else None)

            savefig_or_show(dpi=dpi, save=save, show=show)
            if not show: return ax


def gridspec(nrows=1, ncols=2, figsize=None, dpi=None):
    if figsize is None: figsize = rcParams['figure.figsize']
    gs = pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))
    return gs


def _wraps_plot_scatter(wrapper):
    annots_orig = {
        k: v for k, v in wrapper.__annotations__.items()
        if k not in {'adata', 'kwargs'}
    }
    annots_scatter = {
        k: v for k, v in scatter.__annotations__.items()
        if k != 'basis'
    }
    wrapper.__annotations__ = {**annots_scatter, **annots_orig}
    wrapper.__wrapped__ = scatter
    return wrapper


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def trimap(adata, **kwargs):
    """\
    Scatter plot in trimap basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='trimap', **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def umap(adata, **kwargs):
    """\
    Scatter plot in UMAP basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='umap', **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def tsne(adata, **kwargs):
    """\
    Scatter plot in tsne basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='tsne', **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def diffmap(adata, **kwargs):
    """\
    Scatter plot in diffmap basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='diffmap', **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def phate(adata, **kwargs):
    """\
    Scatter plot in phate basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='phate', **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def draw_graph(adata, layout=None, **kwargs):
    """\
    Scatter plot in draw_graph basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    if layout is None:
        layout = str(adata.uns['draw_graph']['params']['layout'])
    basis = 'draw_graph_' + layout
    if 'X_' + basis not in adata.obsm_keys():
        raise ValueError('Did not find {} in adata.obs. Did you compute layout {}?'
                         .format('draw_graph_' + layout, layout))
    return scatter(adata, basis=basis, **kwargs)


@_wraps_plot_scatter
@doc_params(scatter=doc_scatter)
def pca(adata, **kwargs):
    """\
    Scatter plot in pca basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='pca', **kwargs)
