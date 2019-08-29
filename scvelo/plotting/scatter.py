from .. import settings
from .. import AnnData
from .utils import make_dense, is_categorical, update_axes, set_label, set_title, interpret_colorkey, set_colorbar, \
    default_basis, default_color, default_size, default_color_map, get_components, savefig_or_show, make_unique_list, \
    plot_linear_fit, plot_density, default_legend_loc, make_unique_valid_list, rugplot
from .docs import doc_scatter, doc_params

from matplotlib import rcParams
import matplotlib.pyplot as pl
from scipy.stats import pearsonr
import numpy as np
import pandas as pd


@doc_params(scatter=doc_scatter)
def scatter(adata=None, x=None, y=None, basis=None, vkey=None, color=None, use_raw=None, layer=None, color_map=None,
            colorbar=None, palette=None, size=None, alpha=None, linewidth=None, perc=None, sort_order=True, groups=None,
            components=None, projection='2d', legend_loc=None, legend_fontsize=None, legend_fontweight=None,
            right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None,
            xlim=None, ylim=None, show_density=None, show_assignments=None, show_linear_fit=None, show_polyfit=None,
            show_rug=None, n_convolve=None, dpi=None, frameon=None, show=True, save=None, ax=None, zorder=None, ncols=None, **kwargs):
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
    {scatter}

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """
    scatter_kwargs = {"use_raw": use_raw, "sort_order": sort_order, "alpha": alpha, "components": components,
                      "projection": projection, "groups": groups, "palette": palette, "legend_fontsize": legend_fontsize,
                      "legend_fontweight": legend_fontweight, "right_margin": right_margin, "left_margin": left_margin,
                      "show": False, "save": False}

    ext_kwargs = {'size': size, 'linewidth': linewidth, 'xlabel': xlabel, 'vkey': vkey, 'color_map': color_map,
                  'colorbar': colorbar, 'perc': perc, 'frameon': frameon, 'zorder': zorder, 'legend_loc': legend_loc,
                  'fontsize': fontsize, 'xlim': xlim, 'ylim': ylim, 'n_convolve': n_convolve,
                  'show_density': show_density, 'show_assignments': show_assignments, 'show_rug': show_rug,
                  'show_linear_fit': show_linear_fit, 'show_polyfit': show_polyfit}

    adata = AnnData(np.stack([x, y]).T) if adata is None and (x is not None and y is not None) else adata

    # multiple colors, layers and bases (string)
    colors, ys = make_unique_list(color, allow_array=True), make_unique_list(y, allow_array=True)
    layers, bases = make_unique_list(layer), make_unique_valid_list(adata, basis)
    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 \
        else bases if len(bases) > 1 else ys if len(ys) > 1 else None
    if multikey is not None:
        if isinstance(title, (list, tuple)): title *= int(np.ceil(len(multikey) / len(title)))
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        ax = []
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))):
            if i < len(multikey):
                ax.append(scatter(adata, x=x, ax=pl.subplot(gs), ylabel=ylabel, y=ys[i] if len(ys) > 1 else y,
                                  color=colors[i] if len(colors) > 1 else color,
                                  layer=layers[i] if len(layers) > 1 else layer,
                                  basis=bases[i] if len(bases) > 1 else basis,
                                  title=title[i] if isinstance(title, (list, tuple)) else title,
                                  **scatter_kwargs, **ext_kwargs, **kwargs))
        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax

    else:
        color, layer, basis = colors[0], layers[0], bases[0]

        # comma-separated y or layers (string)
        ys = [yi.strip() for yi in y.split(',')] if isinstance(y, str) and ',' in y else [y]
        layers = [li.strip() for li in layer.split(',')] if isinstance(layer, str) and ',' in layer else [layer]
        multikey = ys if len(ys) > 1 else layers if len(layers) > 1 else None
        if multikey is not None:
            colors = [ci.strip() for ci in color.split(',')] if isinstance(color, str) and ',' in color else [color]
            for i, mi in enumerate(multikey):
                ax = scatter(adata, x=x, ax=ax, basis=basis, title=y if title is None else title,
                             y=ys[i] if len(ys) > 1 else y, ylabel='expression' if ylabel is None else ylabel,
                             layer=layers[i] if len(layers) > 1 else layer,
                             color=colors[i] if len(colors) > 1 else color, **scatter_kwargs, **ext_kwargs)
            if legend_loc is not False and legend_loc is not 'none':
                multikey = [key.replace('Mu', 'unspliced').replace('Ms', 'spliced') for key in multikey]
                ax.legend(multikey, fontsize=legend_fontsize, loc='best' if legend_loc is None else legend_loc)
            savefig_or_show(dpi=dpi, save=save, show=show)
            if not show: return ax

        # perform regular plot
        else:
            if color is None:  color = default_color(adata)
            if color_map is None: color_map = default_color_map(adata, color)

            is_embedding = ((x is None) | (y is None)) and basis not in adata.var_names
            if basis is None and is_embedding: basis = default_basis(adata)
            if size is None: size = default_size(adata)
            if linewidth is None: linewidth = 1
            if frameon is None: frameon = True if not is_embedding else settings._frameon
            if use_raw is None and 'Ms' not in adata.layers.keys(): use_raw = True

            if projection == '3d':
                from mpl_toolkits.mplot3d import Axes3D
                if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca(projection=projection)
            else:
                if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca()

            if is_categorical(adata, color) and is_embedding:
                from scanpy.api.pl import scatter as scatter_
                legend_loc = default_legend_loc(adata, color, legend_loc)
                ax = scatter_(adata, basis=basis, color=color, color_map=color_map, size=size, frameon=frameon, ax=ax,
                              title=title, legend_loc=legend_loc, **scatter_kwargs, **kwargs)

            else:
                if basis in adata.var_names:
                    xkey, ykey = ('spliced', 'unspliced') if use_raw or 'Ms' not in adata.layers.keys() else ('Ms', 'Mu')
                    x, y = adata[:, basis].layers[xkey], adata[:, basis].layers[ykey]
                    if xlabel is None: xlabel = 'spliced'
                    if ylabel is None: ylabel = 'unspliced'
                    if title is None: title = basis

                elif is_embedding:
                    X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)]
                    x, y = X_emb[:, 0], X_emb[:, 1]

                elif isinstance(x, str) and isinstance(y, str):
                    if xlabel is None: xlabel = x
                    if ylabel is None: ylabel = 'expression' if y in adata.var_names and x not in adata.var_names else y
                    if title is None: title = y if y in adata.var_names and x not in adata.var_names else color
                    if layer is None: layer = 'spliced' if use_raw else 'Ms'

                    if x in adata.var_names and y in adata.var_names:
                        x = adata[:, x].layers[layer] if layer in adata.layers.keys() else adata[:, x].X
                        y = adata[:, y].layers[layer] if layer in adata.layers.keys() else adata[:, y].X
                    elif x in adata.var.keys() and y in adata.var.keys():
                        x, y = adata.var[x], adata.var[y]
                        if colors[0] is None: color = 'grey'
                    elif y in adata.var_names:
                        if x in adata.obs.keys(): x = adata.obs[x]
                        elif x in adata.layers.keys(): x = adata[:, y].layers[x]
                        y = adata[:, y].layers[layer] if layer in adata.layers.keys() else adata[:, y].X
                    elif x in adata.obs.keys() and y in adata.obs.keys():
                        x, y = adata.obs[x], adata.obs[y]

                x, y = make_dense(x).flatten(), make_dense(y).flatten()

                if n_convolve is not None:
                    y[np.argsort(x)] = np.convolve(y[np.argsort(x)], np.ones(n_convolve) / n_convolve, mode='same')
                    #y = get_temporal_connectivities(adata, x, n_convolve=n_convolve).dot(y)

                if isinstance(color, int): color = pd.Categorical(np.array(np.arange(len(x)) == color, dtype=bool))
                if basis in adata.var_names and isinstance(color, str) and color in adata.layers.keys():
                    c = interpret_colorkey(adata, basis, color, perc)
                else:
                    c = interpret_colorkey(adata, color, layer, perc)

                if layer is not None and any(l in layer for l in ['spliced', 'Ms', 'Mu', 'velocity']) \
                        and isinstance(color, str) and color in adata.var_names:
                    ub = np.percentile(np.abs(c), 98)
                    if "vmax" not in kwargs:
                        kwargs.update({"vmax": ub})
                    if "vmin" not in kwargs and 'velocity' in layer:
                        kwargs.update({"vmin": -ub})

                if "vmid" in kwargs:
                    if not isinstance(c, str) and not isinstance(c[0], str):
                        vmid, lb, ub = kwargs["vmid"], np.min(c), np.max(c)
                        crange = min(np.abs(vmid - lb), np.abs(ub - vmid))
                        kwargs.update({"vmin": vmid - crange, "vmax": vmid + crange})
                    kwargs.pop("vmid")

                if groups is not None or np.any(pd.isnull(c)):
                    zorder = 0 if zorder is None else zorder
                    ax = scatter(adata, basis=basis, color='lightgrey', ax=ax, zorder=zorder, **scatter_kwargs)
                    zorder += 1

                if basis in adata.var_names:
                    fits = plot_linear_fit(adata, basis, vkey, xkey, linewidth)
                    from .simulation import show_full_dynamics
                    if 'true_alpha' in adata.var.keys() and (vkey is None or 'true_dynamics' in vkey):
                        fit = show_full_dynamics(adata, basis, 'true', use_raw, linewidth)
                        fits.append(fit)
                    if 'fit_alpha' in adata.var.keys() and (vkey is None or 'dynamic' in vkey):
                        fit = show_full_dynamics(adata, basis, 'fit', use_raw, linewidth, show_assignments=show_assignments)
                        fits.append(fit)
                    if len(fits) > 0 and legend_loc is not False and legend_loc is not 'none':
                        ax.legend(fits, fontsize=legend_fontsize, loc='lower right' if legend_loc is None else legend_loc)
                    if use_raw and perc is not None:
                        ax.set_xlim(right=np.percentile(x, 99.9 if not isinstance(perc, int) else perc) * 1.05)
                        ax.set_ylim(top=np.percentile(y, 99.9 if not isinstance(perc, int) else perc) * 1.05)

                smp = ax.scatter(np.ravel(x), np.ravel(y), c=np.ravel(c), cmap=color_map, s=size, alpha=alpha,
                                 edgecolors='none', marker='.', zorder=zorder, **kwargs)

                if show_density:
                    plot_density(x, y, color=show_density if isinstance(show_density, str) else 'grey', ax=ax)

                if show_linear_fit:
                    xnew = np.linspace(np.min(x), np.max(x) * 1.02)
                    ax.plot(xnew, xnew * (x * y).sum() / (x ** 2).sum(),
                            color=show_linear_fit if isinstance(show_linear_fit, str) else 'grey')
                    corr, _ = pearsonr(x, y)
                    ax.text(.05, .95, r'$\rho = $' + str(np.round(corr, 2)),
                            ha='left', va='top', transform=ax.transAxes, fontsize=fontsize,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))
                if show_polyfit:
                    fit = np.polyfit(x, y, deg=2 if isinstance(show_polyfit, (str, bool)) else show_polyfit)
                    f = np.poly1d(fit)
                    xnew = np.linspace(np.min(x), np.max(x), num=100)
                    color = show_polyfit if isinstance(show_polyfit, str) else c if isinstance(c, str) else 'grey'
                    ax.plot(xnew, f(xnew), color=color, linewidth=linewidth)

                if show_rug:
                    ax = rugplot(np.ravel(x), color=np.ravel(interpret_colorkey(adata, show_rug)), ax=ax)

                set_label(xlabel, ylabel, fontsize, basis, ax=ax)
                set_title(title, layer, color, fontsize, ax=ax)
                ax = update_axes(ax, xlim, ylim, fontsize, is_embedding, frameon)

                if colorbar is None and not isinstance(c, str): colorbar = True
                if colorbar and not is_categorical(adata, color):
                    set_colorbar(smp, ax=ax, labelsize=fontsize * .75 if fontsize is not None else None)

            savefig_or_show(dpi=dpi, save=save, show=show)
            if not show: return ax
