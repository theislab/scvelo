from .. import settings
from .. import logging as logg
from . import palettes

import os
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import is_color_like, ListedColormap, to_rgb, cnames
from matplotlib import rcParams
from scipy.sparse import issparse
from cycler import Cycler, cycler


def make_dense(X):
    return X.A if issparse(X) and X.ndim == 2 else X.A1 if issparse(X) else X


def strings_to_categoricals(adata):
    """Transform string annotations to categoricals.
    """
    if not adata._isview:
        from pandas.api.types import is_string_dtype
        from pandas import Categorical
        for df in [adata.obs, adata.var]:
            string_cols = [key for key in df.columns if is_string_dtype(df[key])]
            for key in string_cols:
                c = df[key]
                c = Categorical(c)
                if len(c.categories) < len(c): df[key] = c


def is_categorical(adata, c):
    from pandas.api.types import is_categorical as cat
    strings_to_categoricals(adata)
    str_not_var = isinstance(c, str) and c not in adata.var_names
    return str_not_var and (c in adata.obs.keys() and cat(adata.obs[c]) or is_color_like(c))


def n_categories(adata, c):
    from pandas.api.types import is_categorical as cat
    return len(adata.obs[c].cat.categories) if (c in adata.obs.keys() and cat(adata.obs[c])) else 0


def default_basis(adata):
    keys = [key for key in ['pca', 'tsne', 'umap'] if 'X_' + key in adata.obsm.keys()]
    if not keys:
        raise ValueError('No basis specified.')
    return keys[-1] if len(keys) > 0 else None


def check_basis(adata, basis):
    if basis in adata.obsm.keys() and 'X_' + basis not in adata.obsm.keys():
        adata.obsm['X_' + basis] = adata.obsm[basis]
        logg.info('Renamed', '\'' + basis + '\'', 'to convention', '\'X_' + basis + '\' (adata.obsm).')


def get_basis(adata, basis):
    if isinstance(basis, str) and basis.startswith('X_'):
        basis = basis[2:]
    check_basis(adata, basis)
    return basis


def make_unique_list(key, allow_array=False):
    from pandas import unique, Index
    if isinstance(key, Index): key = key.tolist()
    is_list = isinstance(key, (list, tuple, np.record)) if allow_array else isinstance(key, (list, tuple, np.ndarray, np.record))
    is_list_of_str = is_list and all(isinstance(item, str) for item in key)
    return unique(key) if is_list_of_str else key if is_list and len(key) < 20 else [key]


def make_unique_valid_list(adata, keys):
    keys = make_unique_list(keys)
    if all(isinstance(item, str) for item in keys):
        for i, key in enumerate(keys):
            if key.startswith('X_'):
                keys[i] = key = key[2:]
            check_basis(adata, key)
        valid_keys = np.hstack([adata.obs.keys(), adata.var.keys(), adata.varm.keys(), adata.obsm.keys(),
                                [key[2:] for key in adata.obsm.keys()], list(adata.layers.keys())])
        keys_ = keys
        keys = [key for key in keys if key in valid_keys or key in adata.var_names]
        keys_ = [key for key in keys_ if key not in keys]
        if len(keys_) > 0:
            logg.warn(', '.join(keys_), 'not found.')
    return keys


def update_axes(ax, xlim=None, ylim=None, fontsize=None, is_embedding=False, frameon=None):
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    frameon = settings._frameon if frameon is None else frameon
    if frameon:
        if is_embedding:
            ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            labelsize = int(fontsize * .75) if fontsize is not None else None
            ax.tick_params(axis='both', which='major', labelsize=labelsize)
    else:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        ax.set_frame_on(False)

    return ax


def set_label(xlabel, ylabel, fontsize=None, basis=None):
    if isinstance(xlabel, str) and isinstance(ylabel, str):
        pl.xlabel(xlabel, fontsize=fontsize)
        pl.ylabel(ylabel, fontsize=fontsize)
    elif basis is not None:
        component_name = ('DC' if 'diffmap' in basis else 'tSNE' if basis == 'tsne' else 'UMAP' if basis == 'umap'
        else 'PC' if basis == 'pca' else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis else basis)
        pl.xlabel(component_name + '1')
        pl.ylabel(component_name + '2')


def set_title(title, layer=None, color=None, fontsize=None):
    if isinstance(title, str):
        title = title.replace('_', ' ')
        pl.title(title, fontsize=fontsize)
    elif isinstance(layer, str) and isinstance(color, str):
        title = (color + ' ' + layer).replace('_', ' ')
        pl.title(title, fontsize=fontsize)
    elif isinstance(color, str):
        title = color.replace('_', ' ')
        pl.title(title, fontsize=fontsize)


def set_frame(ax, frameon):
    frameon = settings._frameon if frameon is None else frameon
    if not frameon:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_frame_on(False)
    return ax


def default_size(adata):
    adjusted, classic = 1.2e5 / adata.n_obs, 20
    return np.mean([adjusted, classic]) if settings._rcParams_style == 'scvelo' else adjusted


def default_color(adata):
    return 'clusters' if 'clusters' in adata.obs.keys() else 'louvain' if 'louvain' in adata.obs.keys() else 'grey'


def default_color_map(adata, c):
    if isinstance(c, str) and c in adata.obs.keys() and not is_categorical(adata, c): c = adata.obs[c]
    obs_unitint = len(np.array(c).flatten()) == adata.n_obs and \
                  (np.min(c) == 0 or np.min(c) is False) and (np.max(c) == 1 or np.max(c) is True)
    return 'viridis_r' if obs_unitint else None


def default_legend_loc(adata, color, legend_loc):
    if legend_loc is False:
        legend_loc = 'none'
    elif legend_loc is None:
        legend_loc = 'upper right' if n_categories(adata, color) <= 4 else 'on data'
    return legend_loc


def clip(c, perc):
    if np.size(perc) < 2: perc = [perc, 100] if perc < 50 else [0, perc]
    lb, ub = np.percentile(c, perc)
    return np.clip(c, lb, ub)


def get_colors(adata, c):
    if is_color_like(c):
        return c
    else:
        if c+'_colors' not in adata.uns.keys():
            palette = default_palette(None)
            palette = adjust_palette(palette, length=len(adata.obs[c].cat.categories))
            adata.uns[c + '_colors'] = palette[:len(adata.obs[c].cat.categories)].by_key()['color']
        cluster_ix = adata.obs[c].cat.codes.values
        return np.array([adata.uns[c + '_colors'][cluster_ix[i]] for i in range(adata.n_obs)])


def interpret_colorkey(adata, c=None, layer=None, perc=None):
    if c is None: c = default_color(adata)
    if is_categorical(adata, c): c = get_colors(adata, c)
    elif isinstance(c, str):
        if c in adata.obs.keys():  # color by observation key
            c = adata.obs[c]
        elif c in adata.var_names:  # color by var in specific layer
            c = adata[:, c].layers[layer] if layer in adata.layers.keys() else adata[:, c].X
            c = c.A.flatten() if issparse(c) else c
        else:
            raise ValueError('color key is invalid! pass valid observation annotation or a gene name')
        if perc is not None: c = clip(c, perc=perc)
    elif len(np.array(c).flatten()) == adata.n_obs:  # continuous coloring
        c = np.array(c).flatten()
        if perc is not None: c = clip(c, perc=perc)
    else:
        raise ValueError('color key is invalid! pass valid observation annotation or a gene name')
    return c


def get_components(components=None, basis=None, projection=None):
    if components is None: components = '1,2,3' if projection == '3d' else '1,2'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1
    if 'diffmap' in basis or 'vmap' in basis: components += 1
    return components


def set_colorbar(ax, orientation='vertical'):
    cb = pl.colorbar(orientation=orientation, cax=inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0))
    cb.set_alpha(1)
    cb.draw_all()
    cb.locator = MaxNLocator(nbins=3, integer=True)
    cb.update_ticks()


def default_arrow(size):
    if isinstance(size, (list, tuple)) and len(size) == 3:
        head_l, head_w, ax_l = size
    elif type(size) == int or type(size) == float:
        head_l, head_w, ax_l = 12 * size, 10 * size, 8 * size
    else:
        head_l, head_w, ax_l = 12, 10, 8
    return head_l, head_w, ax_l


def savefig_or_show(writekey=None, show=None, dpi=None, ext=None, save=None):
    if isinstance(save, str):
        # check whether `save` contains a figure extension
        if ext is None:
            for try_ext in ['.svg', '.pdf', '.png']:
                if save.endswith(try_ext):
                    ext = try_ext[1:]
                    save = save.replace(try_ext, '')
                    break
        # append it
        writekey = (writekey + '_' if writekey is not None and len(writekey) > 0 else '') + save
        save = True
    save = settings.autosave if save is None else save
    show = settings.autoshow if show is None else show

    if save:
        if dpi is None:
            # we need this as in notebooks, the internal figures are also influenced by 'savefig.dpi' this...
            if not isinstance(rcParams['savefig.dpi'], str) and rcParams['savefig.dpi'] < 150:
                if settings._low_resolution_warning:
                    logg.warn(
                        'You are using a low resolution (dpi<150) for saving figures.\n'
                        'Consider running `set_figure_params(dpi_save=...)`, which will '
                        'adjust `matplotlib.rcParams[\'savefig.dpi\']`')
                    settings._low_resolution_warning = False
            else:
                dpi = rcParams['savefig.dpi']
        if not os.path.exists(settings.figdir): os.makedirs(settings.figdir)
        if ext is None: ext = settings.file_format_figs
        filename = settings.figdir + f'{settings.plot_prefix}{writekey}{settings.plot_suffix}.{ext}'
        logg.msg('saving figure to file', filename, v=1)
        pl.savefig(filename, dpi=dpi, bbox_inches='tight')

    if show: pl.show()
    if save: pl.close()  # clear figure


def default_palette(palette=None):
    if palette is None: return rcParams['axes.prop_cycle']
    elif not isinstance(palette, Cycler): return cycler(color=palette)
    else: return palette


def adjust_palette(palette, length):
    islist = False
    if isinstance(palette, list):
        islist = True
    if ((islist and len(palette) < length)
       or (not isinstance(palette, list) and len(palette.by_key()['color']) < length)):
        if length <= 28:
            palette = palettes.default_26
        elif length <= len(palettes.default_64):  # 103 colors
            palette = palettes.default_64
        else:
            palette = ['grey' for i in range(length)]
            logg.info('more than 103 colors would be required, initializing as \'grey\'')
        return palette if islist else cycler(color=palette)
    elif islist:
        return palette
    elif not isinstance(palette, Cycler):
        return cycler(color=palette)
    else:
        return palette


def rgb_custom_colormap(colors=['royalblue', 'white', 'forestgreen'], alpha=None, N=256):
    """Creates a custom colormap with the given colors. Colors can be given as names or as rgb values.

    Arguments
    ---------
    colors: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    alpha: `list`, `np.ndarray` or `None` (default: `None`)
        Alpha of the colors. Must be same length as colors.
    N: `int` (default: `256`)
        y coordinate

    Returns
    -------
        A ListedColormap
    """
    c = []
    for color in colors:
        if type(color) is str:
            c.append(to_rgb(cnames[color]))

    vals = np.ones((N, 4))
    ints = len(c) - 1
    n = int(N / ints)

    alpha = np.ones(len(c)) if alpha is None else alpha

    for j in range(ints):
        for i in range(3):
            vals[n * j:n * (j + 1), i] = np.linspace(c[j][i], c[j + 1][i], n)
        vals[n * j:n * (j + 1), -1] = np.linspace(alpha[j], alpha[j + 1], n)
    return ListedColormap(vals)


def plot_linear_fit(adata, basis, vkey, xkey, linewidth=1):
    xnew = np.linspace(0, np.percentile(make_dense(adata[:, basis].layers[xkey]), 98))
    vkeys = adata.layers.keys() if vkey is None else make_unique_list(vkey)
    fits = [fit for fit in vkeys if all(['velocity' in fit, fit + '_gamma' in adata.var.keys()])]
    for i, fit in enumerate(fits):
        linestyle = '--' if 'variance_' + fit in adata.layers.keys() else '-'
        gamma = adata[:, basis].var[fit + '_gamma'].values if fit + '_gamma' in adata.var.keys() else 1
        beta = adata[:, basis].var[fit + '_beta'].values if fit + '_beta' in adata.var.keys() else 1
        offset = adata[:, basis].var[fit + '_offset'].values if fit + '_offset' in adata.var.keys() else 0
        pl.plot(xnew, gamma / beta * xnew + offset / beta, linestyle=linestyle, linewidth=linewidth,
                c='k' if i == 0 else None)
    return fits


def plot_density(x, y=None, eval_pts=50, scale=10, alpha=.3):
    from scipy.stats import gaussian_kde as kde

    # density plot along x-coordinate
    xnew = np.linspace(min(x), max(x), eval_pts)
    vals = kde(x)(xnew)

    if y is not None:
        offset = max(y) / scale
        scale_s = offset / np.max(vals)
        offset *= 1.3
        vals = vals * scale_s - offset
    else:
        offset = 0

    pl.plot(xnew, vals, color='grey')
    pl.fill_between(xnew, -offset, vals, alpha=alpha, color='grey')
    pl.ylim(-offset)

    # density plot along y-coordinate
    if y is not None:
        ynew = np.linspace(min(y), max(y), eval_pts)
        vals = kde(y)(ynew)

        offset = max(x) / scale
        scale_u = offset / np.max(vals)
        offset *= 1.3
        vals = vals * scale_u - offset

        pl.plot(vals, ynew, color='grey')
        pl.fill_betweenx(ynew, -offset, vals, alpha=alpha, color='grey')
        pl.xlim(-offset)


def hist(arrays, alpha=.5, bins=50, colors=None, labels=None, kde=False, xlabel=None, ylabel=None, xlim=None, ylim=None,
         cutoff=None, xscale=None, yscale=None, fontsize=None, legend_fontsize=None, figsize=None, ax=None, dpi=None, show=True):
    fig = pl.figure(None, figsize, dpi=dpi) if ax is None else ax
    ax = pl.subplot()
    update_axes(ax, xlim, ylim, fontsize, frameon=True)

    arrays = arrays if isinstance(arrays, (list, tuple)) or arrays.ndim > 1 else [arrays]

    palette = default_palette(None).by_key()['color'][::-1]
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    masked_arrays = np.ma.masked_invalid(arrays)
    bmin, bmax = masked_arrays.min(), masked_arrays.max()
    bins = np.arange(bmin, bmax + (bmax - bmin) / bins, (bmax - bmin) / bins)
    if xscale is 'log':
        bins = bins[bins > 0]
        bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

    if cutoff is not None:
        bins = bins[(bins > cutoff[0]) & (bins < cutoff[1])] if isinstance(cutoff, list) else bins[bins < cutoff]

    for i, x in enumerate(arrays):
        if kde:
            from scipy.stats import gaussian_kde
            kde_bins = gaussian_kde(x)(bins)
            pl.plot(bins, kde_bins, color='grey')
            pl.fill_between(bins, 0, kde_bins, alpha=.4, color='grey')
        else:
            pl.hist(x[np.isfinite(x)], bins=bins, alpha=alpha, color=colors[i],
                    label=labels[i] if labels is not None else None)

    set_label(xlabel if xlabel is not None else '', ylabel if xlabel is not None else '', fontsize=fontsize)

    if labels is not None:
        pl.legend(fontsize=legend_fontsize)

    if xscale is not None: pl.xscale(xscale)
    if yscale is not None: pl.yscale(yscale)

    if not show:
        return ax
    else:
        pl.show()


def plot(arrays, normalize=False, colors=None, labels=None, xlabel=None, ylabel=None, xscale=None, yscale=None, ax=None,
         figsize=None, dpi=None, show=True):
    ax = pl.figure(None, figsize, dpi=dpi) if ax is None else ax
    arrays = np.array(arrays)
    arrays = arrays if isinstance(arrays, (list, tuple)) or arrays.ndim > 1 else [arrays]

    palette = default_palette(None).by_key()['color'][::-1]
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    for i, array in enumerate(arrays):
        X = array[np.isfinite(array)]
        X = X / np.max(X) if normalize else X
        pl.plot(X, color=colors[i], label=labels[i] if labels is not None else None)

    pl.xlabel(xlabel if xlabel is not None else '')
    pl.ylabel(ylabel if xlabel is not None else '')
    if labels is not None: pl.legend()
    if xscale is not None: pl.xscale(xscale)
    if yscale is not None: pl.yscale(yscale)

    if not show:
        return ax
    else:
        pl.show()


def fraction_timeseries(adata, xkey='clusters', tkey='dpt_pseudotime', bins=30, legend_loc='best', title=None,
                        fontsize=None, ax=None, figsize=None, dpi=None, xlabel=None, ylabel=None, show=True):
    t = np.linspace(0, 1 + 1 / bins, bins)
    types = np.unique(adata.obs[xkey].values)

    y = []
    for i in range(bins - 1):
        mask = np.all([adata.obs[tkey].values <= t[i + 1], adata.obs[tkey].values > t[i]], axis=0)
        x = list(adata[mask].obs[xkey].values)
        y.append([])
        for name in types:
            occur = x.count(name)
            y[-1].append(occur)
        y[-1] /= np.sum(y[-1])
    y = np.array(y).T

    ax = pl.figure(figsize=figsize, dpi=dpi) if ax is None else ax

    pl.stackplot(t[:-1], y, baseline='zero', labels=types,
                 colors=adata.uns['clusters_colors'] if 'clusters_colors' in adata.uns.keys() else None,
                 edgecolor='white')

    pl.legend(types, loc=legend_loc)
    if title is not None:
        pl.title(title, fontsize=fontsize)
    pl.xlabel(tkey if xlabel is None else xlabel, fontsize=fontsize)
    pl.ylabel(xkey + ' fractions' if ylabel is None else ylabel, fontsize=fontsize)
    pl.xlim(adata.obs[tkey].values.min(), adata.obs[tkey].values.max())
    pl.ylim(0, 1)

    if not show:
        return ax
    else:
        pl.show()


# def phase(adata, var=None, x=None, y=None, color='louvain', fits='all', xlabel='spliced', ylabel='unspliced',
#           fontsize=None, show=True, ax=None, **kwargs):
#     if isinstance(var, str) and (var in adata.var_names):
#         if (x is None) or (y is None):
#             ix = np.where(adata.var_names == var)[0][0]
#             x, y = adata.layers['Ms'][:, ix], adata.layers['Mu'][:, ix]
#     else:
#         ValueError('var not found in adata.var_names.')
#
#     ax = scatter(adata, x=x, y=y, color=color, frameon=True, title=var, xlabel=xlabel, ylabel=ylabel, ax=ax, **kwargs)
#
#     xnew = np.linspace(0, x.max() * 1.02)
#     fits = adata.layers.keys() if fits == 'all' else fits
#     fits = [fit for fit in fits if 'velocity' in fit]
#     for fit in fits:
#         linestyle = '--' if 'stochastic' in fit else '-'
#         pl.plot(xnew, adata.var[fit+'_gamma'][ix] / adata.var[fit+'_beta'][ix] * xnew
#                 + adata.var[fit+'_offset'][ix] / adata.var[fit+'_beta'][ix], c='k', linestyle=linestyle)
#
#     if show: pl.show()
#     else: return ax