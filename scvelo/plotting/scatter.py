from .. import settings
from .. import logging as logg
from .. import AnnData
from .docs import doc_scatter, doc_params
from .utils import make_dense, is_categorical, update_axes, set_label, set_title, interpret_colorkey, set_colorbar, \
    default_basis, default_color, default_size, default_color_map, get_components, savefig_or_show, make_unique_list, \
    plot_linear_fit, plot_density, default_legend_loc, make_unique_valid_list, rugplot, groups_to_bool, \
    _set_colors_for_categorical_obs, _add_legend, get_connectivities, plot_outline

from matplotlib import rcParams, patheffects
import matplotlib.pyplot as pl
from scipy.stats import pearsonr
import numpy as np
import pandas as pd


@doc_params(scatter=doc_scatter)
def scatter(adata=None, x=None, y=None, basis=None, vkey=None, color=None, use_raw=None, layer=None, color_map=None,
            colorbar=None, palette=None, size=None, alpha=None, linewidth=None, perc=None, sort_order=True, groups=None,
            components=None, projection='2d', legend_loc=None, legend_fontsize=None, legend_fontweight=None,
            xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None, xlim=None, ylim=None, show_density=None,
            show_assignments=None, show_linear_fit=None, show_polyfit=None, rug=None, add_outline=False,
            outline_width=None, outline_color=None, n_convolve=None, smooth=None, rescale_color=None, dpi=None,
            frameon=None, zorder=None, ncols=None, wspace=None, hspace=None, show=True, save=None, ax=None, **kwargs):
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
    scatter_kwargs = {'vkey': vkey, 'use_raw': use_raw, 'color_map': color_map, 'colorbar': colorbar, 'palette': palette,
                      'size': size, 'alpha': alpha, 'linewidth': linewidth, 'perc': perc, 'sort_order': sort_order,
                      'projection': projection, 'legend_loc': legend_loc, 'legend_fontsize': legend_fontsize,
                      'legend_fontweight': legend_fontweight, 'xlabel': xlabel, 'fontsize': fontsize, 'xlim': xlim,
                      'ylim': ylim, 'show_density': show_density, 'show_assignments': show_assignments,
                      'show_linear_fit': show_linear_fit, 'show_polyfit': show_polyfit, 'rug': rug,
                      'add_outline': add_outline, 'outline_color': outline_color, 'n_convolve': n_convolve, 'smooth': smooth,
                      'rescale_color': rescale_color, 'frameon': frameon, 'zorder': zorder, 'show': False, 'save': False}

    adata = AnnData(np.stack([x, y]).T) if adata is None and (x is not None and y is not None) else adata
    # multiple colors, layers and bases (string)
    if 'c' in kwargs: color = kwargs.pop('c')
    colors = make_unique_list(color, allow_array=True)
    xs, ys = make_unique_list(x, allow_array=True), make_unique_list(y, allow_array=True)
    layers, bases = make_unique_list(layer), make_unique_valid_list(adata, basis)
    components = make_unique_list(components)
    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 \
        else bases if len(bases) > 1 else xs if len(xs) > 1 else ys if len(ys) > 1 \
        else components if len(components) > 1 else None
    if multikey is not None:
        if ax is not None: logg.warn("Cannot specify `ax` when plotting multiple panels.")
        if isinstance(title, (list, tuple)): title *= int(np.ceil(len(multikey) / len(title)))
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        ax = []
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi),
                            hspace=hspace, wspace=wspace)):
            if i < len(multikey):
                ax.append(scatter(adata, ax=pl.subplot(gs), ylabel=ylabel, groups=groups,
                                  x=xs[i] if len(xs) > 1 else x, y=ys[i] if len(ys) > 1 else y,
                                  color=colors[i] if len(colors) > 1 else color,
                                  layer=layers[i] if len(layers) > 1 else layer,
                                  basis=bases[i] if len(bases) > 1 else basis,
                                  components=components[i] if len(components) > 1 else components,
                                  title=title[i] if isinstance(title, (list, tuple)) else title,
                                  **scatter_kwargs, **kwargs))
        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax

    else:
        color, layer, basis, components = colors[0], layers[0], bases[0], components[0]

        # comma-separated y or layers (string)
        ys = [yi.strip() for yi in y.split(',')] if isinstance(y, str) and ',' in y else [y]
        layers = [li.strip() for li in layer.split(',')] if isinstance(layer, str) and ',' in layer else [layer]
        multikey = ys if len(ys) > 1 else layers if len(layers) > 1 else None
        if multikey is not None:
            colors = [ci.strip() for ci in color.split(',')] if isinstance(color, str) and ',' in color else [color]
            for i, mi in enumerate(multikey):
                ax = scatter(adata, x=x, ax=ax, basis=basis, title=y if title is None else title, groups=groups,
                             y=ys[i] if len(ys) > 1 else y, ylabel='expression' if ylabel is None else ylabel,
                             layer=layers[i] if len(layers) > 1 else layer,
                             color=colors[i] if len(colors) > 1 else color, **scatter_kwargs)
            if legend_loc is not False and legend_loc is not 'none':
                multikey = [key.replace('Mu', 'unspliced').replace('Ms', 'spliced') for key in multikey]
                ax.legend(multikey, fontsize=legend_fontsize, loc='best' if legend_loc is None else legend_loc)
            savefig_or_show(dpi=dpi, save=save, show=show)
            if not show: return ax

        # perform regular plot
        else:
            # set color, color_map, edgecolor, basis, linewidth, frameon, use_raw, dim
            if color is None:  color = default_color(adata)
            if 'cmap' not in kwargs:
                kwargs['cmap'] = default_color_map(adata, color) if color_map is None else color_map
            if 's' not in kwargs:
                kwargs['s'] = default_size(adata) if size is None else size
            if 'edgecolor' not in kwargs:
                kwargs['edgecolor'] = 'none'
            is_embedding = ((x is None) | (y is None)) and basis not in adata.var_names
            if basis is None and is_embedding: basis = default_basis(adata)
            if linewidth is None: linewidth = 1
            if frameon is None: frameon = True if not is_embedding else settings._frameon
            if isinstance(groups, str): groups = [groups]
            if use_raw is None and basis not in adata.var_names:
                use_raw = layer is None and adata.raw is not None
            dim = 3 if '3' in projection else 2
            if dim == 3: from mpl_toolkits.mplot3d import Axes3D

            if ax is None:
                ax = pl.figure(None, figsize, dpi=dpi).gca(projection='3d')\
                     if dim == 3 else pl.figure(None, figsize, dpi=dpi).gca()

            # set legend
            if is_categorical(adata, color) and is_embedding:
                _set_colors_for_categorical_obs(adata, color, palette)
                legend_loc = default_legend_loc(adata, color, legend_loc)
                legend_fontweight = 'bold' if legend_fontweight is None else legend_fontweight
                _add_legend(adata, ax, color, legend_loc, adata.obsm['X_' + basis][:, :dim], legend_fontweight,
                            legend_fontsize, [patheffects.withStroke(linewidth=True, foreground='w')], groups, False)

            # phase portrait: get x and y from .layers (e.g. spliced vs. unspliced) when basis is in var_names
            if basis in adata.var_names:
                if title is None: title = basis
                if x is None:
                    x = 'spliced' if 'spliced' in adata.layers.keys() and (use_raw or 'Ms' not in adata.layers.keys())\
                        else 'Ms' if 'Ms' in adata.layers.keys() else 'X'
                if y is None:
                    y = 'unspliced' if 'unspliced' in adata.layers.keys() and (use_raw or 'Mu' not in adata.layers.keys())\
                        else 'Mu' if 'Mu' in adata.layers.keys() else 'X'
                if xlabel is None: xlabel = x
                if ylabel is None: ylabel = y

                if not (x in adata.layers.keys() or x is 'X') and not (y in adata.layers.keys() or y is 'X'):
                    raise ValueError('Could not find x or y in layers.')
                x = adata[:, basis].layers[x] if x in adata.layers.keys() \
                    else adata.raw.obs_vector(basis) if use_raw else adata.obs_vector(basis)
                y = adata[:, basis].layers[y] if y in adata.layers.keys() \
                    else adata.raw.obs_vector(basis) if use_raw else adata.obs_vector(basis)

            # embedding: set x and y to embedding coordinates
            elif is_embedding:
                X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)]
                x, y = X_emb[:, 0], X_emb[:, 1]
                z = X_emb[:, 2] if dim == 3 and X_emb.shape[1] > 2 else None

            elif isinstance(x, str) and isinstance(y, str):
                if layer is None:
                    layer = 'spliced' if 'spliced' in adata.layers.keys() and (use_raw or 'Ms' not in adata.layers.keys())\
                        else 'Ms' if 'Ms' in adata.layers.keys() else 'X'

                # gene trend: get x and y as gene (var_names) along obs/layers (e.g. pseudotime)
                if y in adata.var_names and (x in adata.obs.keys() or x in adata.layers.keys()):
                    if xlabel is None: xlabel = x
                    if ylabel is None: ylabel = layer
                    if title is None: title = y
                    x = adata.obs[x] if x in adata.obs.keys() else adata[:, y].layers[x]
                    y = adata[:, y].layers[layer] if layer in adata.layers.keys() \
                        else adata.raw.obs_vector(y) if use_raw else adata.obs_vector(y)
                # get x and y from var_names, var or obs
                else:
                    if xlabel is None: xlabel = x
                    if ylabel is None: ylabel = y
                    if title is None: title = color

                    if x in adata.var_names and y in adata.var_names:
                        x = adata[:, x].layers[layer] if layer in adata.layers.keys() else adata.raw.obs_vector(x) if use_raw else adata.obs_vector(x)
                        y = adata[:, y].layers[layer] if layer in adata.layers.keys() else adata.raw.obs_vector(y) if use_raw else adata.obs_vector(y)
                    elif x in adata.var.keys() and y in adata.var.keys():
                        x, y = adata.var[x], adata.var[y]
                        if colors[0] is None: color = 'grey'  # since default_color (clusters) does not match with .var
                    elif x in adata.obs.keys() and y in adata.obs.keys():
                        x, y = adata.obs[x], adata.obs[y]
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

            # smooth color values across neighbors and rescale
            if smooth and len(c) == adata.n_obs:
                c = get_connectivities(adata, n_neighbors=(None if isinstance(smooth, bool) else smooth)).dot(c)
            if rescale_color is not None:
                c += rescale_color[0] - np.min(c)
                c *= rescale_color[1] / np.max(c)

            # check if higher value points should be plotted on top
            if sort_order and not isinstance(c, str) and not is_categorical(adata, color):
                order = np.argsort(c)
                c = c[order]
                x = x[order]
                y = y[order]

                # check if 'size' is given as a vector and reorder it.
                if isinstance(kwargs['s'], np.ndarray):
                    kwargs['s'] = np.array(kwargs['s'])[order]

            # adjust coloring to ignore extreme outliers since these layers are not logarithmized
            if layer is not None and any(l in layer for l in ['spliced', 'unspliced', 'Ms', 'Mu', 'velocity']) \
                    and isinstance(color, str) and color in adata.var_names:
                ub = np.percentile(np.abs(c), 98)
                if "vmax" not in kwargs:
                    kwargs.update({"vmax": ub})
                if "vmin" not in kwargs and 'velocity' in layer:
                    kwargs.update({"vmin": -ub})

            # introduce vmid by setting vmin and vmax accordingly
            if "vmid" in kwargs:
                if not isinstance(c, str) and not isinstance(c[0], str):
                    vmid, lb, ub = kwargs["vmid"], np.min(c), np.max(c)
                    crange = min(np.abs(vmid - lb), np.abs(ub - vmid))
                    kwargs.update({"vmin": vmid - crange, "vmax": vmid + crange})
                kwargs.pop("vmid")

            # set color to grey for NAN values and for cells that are not in groups
            if groups is not None or np.any(pd.isnull(c)):
                zorder = 0 if zorder is None else zorder
                ax = scatter(adata, x=x, y=y, basis=basis, layer=layer, color='lightgrey', ax=ax, groups=None, **scatter_kwargs)
                idx = groups_to_bool(adata, groups, color)
                x, y = x[idx], y[idx]
                if not isinstance(c, str) and len(c) == adata.n_obs:
                    c = c[idx]
                zorder += 1

            # set color to grey for NAN values and for cells that are not in groups
            if basis in adata.var_names:
                if use_raw is None: use_raw = 'Ms' not in adata.layers.keys()
                lines, fits = plot_linear_fit(adata, basis, vkey, 'spliced' if use_raw else 'Ms', linewidth, ax=ax)
                from .simulation import show_full_dynamics
                if 'true_alpha' in adata.var.keys() and (vkey is not None and 'true_dynamics' in vkey):
                    line, fit = show_full_dynamics(adata, basis, 'true', use_raw, linewidth, ax=ax)
                    fits.append(fit)
                    lines.append(line)
                if 'fit_alpha' in adata.var.keys() and (vkey is None or 'dynamics' in vkey):
                    line, fit = show_full_dynamics(adata, basis, 'fit', use_raw, linewidth, show_assignments=show_assignments, ax=ax)
                    fits.append(fit)
                    lines.append(line)
                if len(fits) > 0 and legend_loc is not False and legend_loc is not 'none':
                    ax.legend(handles=lines, labels=fits, fontsize=legend_fontsize,
                              loc='lower right' if legend_loc is None else legend_loc)
                if use_raw and perc is not None:
                    ax.set_xlim(right=np.percentile(x, 99.9 if not isinstance(perc, int) else perc) * 1.05)
                    ax.set_ylim(top=np.percentile(y, 99.9 if not isinstance(perc, int) else perc) * 1.05)

            x, y, c = np.ravel(x), np.ravel(y), np.ravel(c)
            if not isinstance(c, str) and len(c) != len(x): c = 'grey'

            if len(x) != len(y):
                raise ValueError('x or y do not share the same dimension.')
            elif not isinstance(c, str) and len(x) != len(c):
                raise ValueError('color and x do not share the same dimensions.')

            smp = ax.scatter(x, y, c=c, alpha=alpha, marker='.', zorder=zorder, **kwargs)

            if add_outline:
                plot_outline(x, y, kwargs, outline_width, outline_color, zorder, ax=ax)

            if show_density:
                plot_density(x, y, color=show_density if isinstance(show_density, str) else 'grey', ax=ax)

            if show_linear_fit:
                idx_valid = ~np.isnan(x + y)
                x, y = x[idx_valid], y[idx_valid]
                xnew = np.linspace(np.min(x), np.max(x) * 1.02)
                ax.plot(xnew, xnew * (x * y).sum() / (x ** 2).sum(), linewidth=linewidth,
                        color=show_linear_fit if isinstance(show_linear_fit, str) else 'grey')
                corr, _ = pearsonr(x, y)
                if legend_loc is not 'none':
                    ax.text(.05, .95, r'$\rho = $' + str(np.round(corr, 2)), ha='left', va='top', fontsize=fontsize,
                            transform=ax.transAxes, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))
            if show_polyfit:
                idx_valid = ~np.isnan(x + y)
                x, y = x[idx_valid], y[idx_valid]
                fit = np.polyfit(x, y, deg=2 if isinstance(show_polyfit, (str, bool)) else show_polyfit)
                f = np.poly1d(fit)
                xnew = np.linspace(np.min(x), np.max(x), num=100)
                color = show_polyfit if isinstance(show_polyfit, str) else c if isinstance(c, str) else 'grey'
                ax.plot(xnew, f(xnew), color=color, linewidth=linewidth)

            if rug:
                ax = rugplot(np.ravel(x), color=np.ravel(interpret_colorkey(adata, rug)), ax=ax)

            set_label(xlabel, ylabel, fontsize, basis, ax=ax)
            set_title(title, layer, color, fontsize, ax=ax)
            ax = update_axes(ax, xlim, ylim, fontsize, is_embedding, frameon)

            if colorbar is not False and not isinstance(c, str) and not is_categorical(adata, color):
                set_colorbar(smp, ax=ax, labelsize=fontsize * .75 if fontsize is not None else None)

            savefig_or_show(dpi=dpi, save=save, show=show)
            if not show: return ax


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
    Scatter plot in UMAP basis.
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
    Scatter plot in UMAP basis.
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
def draw_graph(adata, layout=None, **kwargs):
    """\
    Scatter plot in UMAP basis.
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
    Scatter plot in UMAP basis.
    Parameters
    ----------
    {scatter}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
    return scatter(adata, basis='pca', **kwargs)
