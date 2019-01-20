from .. import settings
from .. import logging as logg
from . import palettes

import os
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import is_color_like
from matplotlib import rcParams
from scipy.sparse import issparse
from cycler import Cycler, cycler


def strings_to_categoricals(adata):
    """Transform string annotations to categoricals.
    """
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
    return isinstance(c, str) and (c in adata.obs.keys() and cat(adata.obs[c]) or is_color_like(c))


def default_basis(adata):
    keys = [key for key in ['pca', 'tsne', 'umap'] if 'X_' + key in adata.obsm.keys()]
    return keys[-1] if len(keys) > 0 else None


def make_unique_list(key, allow_array=False):
    from pandas import unique
    is_list = isinstance(key, (list, tuple, np.record)) if allow_array else isinstance(key, (list, tuple, np.ndarray, np.record))
    is_list_of_str = is_list and all(isinstance(item, str) for item in key)
    return unique(key) if is_list_of_str else key if is_list and len(key) < 20 else [key]


def update_axes(ax, fontsize, is_embedding, frameon):
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
    return 'viridis_r' if isinstance(c, str) and c in adata.obs.keys() and not is_categorical(adata, c)\
                          and adata.obs[c].min() == 0 and adata.obs[c].max() == 1 else None


def clip(c, perc=None):
    if isinstance(perc, int): perc = [perc, 100] if perc < 50 else [0, perc]
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
        cluster_ix = adata.obs[c].cat.codes
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


def savefig_or_show(writekey, show=None, dpi=None, ext=None, save=None):
    if isinstance(save, str):
        # check whether `save` contains a figure extension
        if ext is None:
            for try_ext in ['.svg', '.pdf', '.png']:
                if save.endswith(try_ext):
                    ext = try_ext[1:]
                    save = save.replace(try_ext, '')
                    break
        # append it
        writekey = 'velocity_' + (writekey + '_' if len(writekey) > 0 else '') + save
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
        if settings.figdir[-1] != '/': settings.figdir += '/'
        if ext is None: ext = settings.file_format_figs
        filename = settings.figdir + writekey + settings.plot_suffix + '.' + ext
        # output the following msg at warning level; it's really important for the user
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


def hist(arrays, alpha=.5, bins=None, colors=None, labels=None, xlabel=None, ylabel=None, ax=None, figsize=None, dpi=None):
    ax = pl.figure(None, figsize, dpi=dpi) if ax is None else ax
    arrays = arrays if isinstance(arrays, (list, tuple)) else [arrays]

    palette = default_palette(None)[::3][:len(arrays)].by_key()['color']
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    for i, array in enumerate(arrays):
        pl.hist(array[np.isfinite(array)], bins=bins, alpha=alpha, color=colors[i], label=labels[i] if labels is not None else None)
    pl.legend()
    pl.xlabel(xlabel if xlabel is not None else '')
    pl.ylabel(ylabel if xlabel is not None else '')
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