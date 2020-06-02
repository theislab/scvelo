from .. import settings
from .. import logging as logg
from .. import AnnData
from ..preprocessing.moments import get_connectivities
from ..tools.utils import strings_to_categoricals
from . import palettes

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import is_color_like, ListedColormap, to_rgb, cnames
from matplotlib.collections import LineCollection
from matplotlib.gridspec import SubplotSpec
from matplotlib import patheffects
import matplotlib.transforms as tx
from matplotlib import rcParams
from pandas import unique, Index
from scipy.sparse import issparse
from scipy.stats import pearsonr
from cycler import Cycler, cycler
from collections import abc


"""helper functions"""


def make_dense(X):
    XA = X.A if issparse(X) and X.ndim == 2 else X.A1 if issparse(X) or isinstance(X, np.matrix) else X
    return np.array(XA)


def is_view(adata):
    return adata.is_view if hasattr(adata, 'is_view') \
        else adata.isview if hasattr(adata, 'isview') \
        else adata._isview if hasattr(adata, '_isview') else True


def is_categorical(adata, c=None):
    from pandas.api.types import is_categorical as cat
    if c is None: return cat(adata)
    if not is_view(adata): strings_to_categoricals(adata)
    return isinstance(c, str) and c in adata.obs.keys() and cat(adata.obs[c])


def is_list(key):
    return isinstance(key, (list, tuple, np.record))


def is_list_of_str(key, max_len=None):
    if max_len is not None:
        return isinstance(key, (list, tuple, np.record, np.ndarray)) \
               and len(key) < max_len \
               and all(isinstance(item, str) for item in key)
    else:
        return isinstance(key, (list, tuple, np.record)) \
               and all(isinstance(item, str) for item in key)


def is_list_of_list(lst):
    return lst is not None and any(isinstance(l, list) for l in lst)


def is_list_of_int(lst):
    return isinstance(lst, (list, tuple, np.ndarray, np.record)) \
           and all(isinstance(item, (int, np.integer)) for item in lst)


def to_list(key, max_len=20):
    if isinstance(key, Index) or is_list_of_str(key, max_len):
        key = list(key)
    return key if is_list(key) and len(key) < max_len else [key]


def to_val(key):
    return key[0] if isinstance(key, (list, tuple)) and len(key) == 1 else key


def get_figure_params(figsize, dpi=None, ncols=1):
    figsize = rcParams['figure.figsize'] if figsize is None else figsize
    dpi = rcParams['figure.dpi'] if dpi is None else dpi
    if settings.presenter_view and figsize[0] * ncols * (dpi / 80) > 12:
        figscale = 12 / (figsize[0] * ncols)
        figsize = (figsize[0] * figscale, figsize[1] * figscale)
        dpi = min(80, dpi)
    return figsize, dpi


def get_ax(ax, show=None, figsize=None, dpi=None, projection=None):
    figsize, _ = get_figure_params(figsize)
    if ax is None:
        ax = pl.figure(None, figsize, dpi=dpi).gca(projection='3d' if projection == '3d' else None)
        if show is None: show = True
    elif isinstance(ax, SubplotSpec):
        geo = ax.get_geometry()
        if show is None: show = geo[-1] + 1 == geo[0] * geo[1]
        ax = pl.subplot(ax)
    return ax, show


def get_kwargs(kwargs, dict_new_kwargs):
    kwargs = kwargs.copy()
    kwargs.update(dict_new_kwargs)
    return kwargs


def check_basis(adata, basis):
    if basis in adata.obsm.keys() and 'X_' + basis not in adata.obsm.keys():
        adata.obsm['X_' + basis] = adata.obsm[basis]
        logg.info('Renamed', '\'' + basis + '\'', 'to convention', '\'X_' + basis + '\' (adata.obsm).')


def get_basis(adata, basis):
    if isinstance(basis, str) and basis.startswith('X_'):
        basis = basis[2:]
    check_basis(adata, basis)
    return basis


def to_valid_bases_list(adata, keys):
    if isinstance(keys, pd.DataFrame):
        keys = keys.index
    if not isinstance(keys, str):
        keys = list(np.ravel(keys))
    keys = to_list(keys, max_len=np.inf)
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
            msg_embedding = ''
            if len(keys_) == 1 and keys_[0] in {'diffmap', 'umap', 'tsne'}:
                msg_embedding = f"You need to run `scv.tl.{keys_[0]}` first."
            logg.warn(', '.join(keys_), 'not found.', msg_embedding)
    return keys


def get_components(components=None, basis=None, projection=None):
    if components is None: components = '1,2,3' if projection == '3d' else '1,2'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1
    if 'diffmap' in basis or 'vmap' in basis: components += 1
    return components


def get_obs_vector(adata, basis, layer=None, use_raw=None):
    return adata.obs_vector(basis, layer=layer) if layer in adata.layers.keys() \
        else adata.raw.obs_vector(basis) if use_raw else adata.obs_vector(basis)


def get_value_counts(adata, color):
    value_counts = adata.obs[color].value_counts()
    probs = np.array(adata.obs[color])
    for cat in value_counts.index:
        probs[probs == cat] = value_counts[cat]
    return np.array(probs, dtype=np.float32)


def get_groups(adata, groups, groupby=None):
    if not isinstance(groupby, str) or groupby not in adata.obs.keys():
        groupby = 'clusters' if 'clusters' in adata.obs.keys() \
            else 'louvain' if 'louvain' in adata.obs.keys() else None
    if groups is True:
        return None, groupby
    if groups is not None and not isinstance(groups, str) and len(groups) == 1:
        groups = groups[0]
    if isinstance(groups, str):
        cats = ['']
        if is_categorical(adata, groupby):
            cats = adata.obs[groupby].cat.categories
        if ':' in groups and not np.any([':' in cat for cat in cats]):
            groupby, groups = groups.split(':')
            groups = groups.strip()
        if ',' in groups and not np.any([',' in cat for cat in cats]):
            groups = [a.strip() for a in groups.split(',')]
    if isinstance(groups, str):
        groups = [groups]
    return groups, groupby


def groups_to_bool(adata, groups, groupby=None):
    groups, groupby = get_groups(adata, groups, groupby)
    if isinstance(groups, (list, tuple, np.ndarray, np.record)):
        if groupby is not None and isinstance(groups[0], str):
            groups = np.array([key in groups for key in adata.obs[groupby]])

    if groupby is not None and groupby in adata.obs.keys():
        c = adata.obs[groupby]
        if np.any(pd.isnull(c)):
            valid = np.array(~pd.isnull(c))
            groups = valid if groups is None or len(groups) != len(c) else groups & valid

    groups = np.ravel(groups) if isinstance(groups, (list, tuple, np.ndarray, np.record)) else None
    return groups


def gets_vals_from_color_gradients(adata, color=None, **scatter_kwargs):
    color_gradients = scatter_kwargs.pop('color_gradients')
    scatter_kwargs.update({'color_gradients': None})
    if 'colorbar' not in scatter_kwargs or scatter_kwargs['colorbar'] is None:
        scatter_kwargs.update({'colorbar': False})
    if 's' not in scatter_kwargs:
        scatter_kwargs['s'] = default_size(adata) if scatter_kwargs['size'] is None else scatter_kwargs['size']
    if not any([v in scatter_kwargs for v in ['vmin', 'vmax', 'vmid']]):
        scatter_kwargs['vmid'] = 0
    if isinstance(color_gradients, str) and color_gradients in adata.obsm.keys():
        if color is None: color = color_gradients
        color_gradients = adata.obsm[color_gradients]
    elif isinstance(color_gradients, (list, tuple)) and color_gradients[0] in adata.obs.keys():
        color_gradients = pd.DataFrame(np.stack([adata.obs[c] for c in color_gradients]).T, columns=color_gradients)
    if color is None:
        color = 'clusters_gradients'
    palette = scatter_kwargs.pop('palette')
    if palette is None and hasattr(color_gradients, 'colors'):
        palette = list(color_gradients.colors)

    pd_colgrad = pd.DataFrame(color_gradients)
    vals = np.clip(pd_colgrad.values, 0, None)
    names = color_gradients.names if hasattr(color_gradients, 'names') else pd_colgrad.columns

    adata.obs[color] = pd.Categorical([str(names[i]) for i in np.argmax(vals, 1)], categories=names)
    set_colors_for_categorical_obs(adata, color, palette)

    return vals, names, color, scatter_kwargs


"""get default parameters"""


def default_basis(adata, **kwargs):
    if 'x' in kwargs and 'y' in kwargs:
        keys, x, y = ['embedding'], kwargs.pop('x'), kwargs.pop('y')
        adata.obsm['X_embedding'] = np.stack([x, y]).T
        if 'velocity_embedding' in adata.obsm.keys(): del adata.obsm['velocity_embedding']
    else:
        keys = [key for key in ['pca', 'tsne', 'umap'] if 'X_' + key in adata.obsm.keys()]
    if not keys:
        raise ValueError('No basis specified.')
    return keys[-1] if len(keys) > 0 else None


def default_size(adata):
    adjusted, classic = 1.2e5 / adata.n_obs, 20
    return np.mean([adjusted, classic]) if settings._rcParams_style == 'scvelo' else adjusted


def default_color(adata, add_outline=None):
    if isinstance(add_outline, str) and add_outline in adata.var.keys() \
            and 'recover_dynamics' in adata.uns.keys() and add_outline in adata.uns['recover_dynamics']:
        return adata.uns['recover_dynamics'][add_outline]
    return 'clusters' if 'clusters' in adata.obs.keys() else 'louvain' if 'louvain' in adata.obs.keys() else 'grey'


def default_color_map(adata, c):
    cmap = None
    if isinstance(c, str) and c in adata.obs.keys() and not is_categorical(adata, c): c = adata.obs[c]
    elif isinstance(c, int): cmap = 'viridis_r'
    if len(np.array(c).flatten()) == adata.n_obs:
        try:
            if np.min(c) in [-1, 0, False] and np.max(c) in [1, True]: cmap = 'viridis_r'
        except:
            cmap = None
    return cmap


def default_legend_loc(adata, color, legend_loc):
    n_categories = len(adata.obs[color].cat.categories) if is_categorical(adata, color) else 0
    if legend_loc is False:
        legend_loc = 'none'
    elif legend_loc is None:
        legend_loc = 'upper right' if n_categories <= 4 else 'on data'
    return legend_loc


def default_xkey(adata, use_raw):
    use_raw = 'spliced' in adata.layers.keys() and (use_raw or 'Ms' not in adata.layers.keys())
    return 'spliced' if use_raw else 'Ms' if 'Ms' in adata.layers.keys() else 'X'


def default_ykey(adata, use_raw):
    use_raw = 'unspliced' in adata.layers.keys() and (use_raw or 'Mu' not in adata.layers.keys())
    return 'unspliced' if use_raw else 'Mu' if 'Mu' in adata.layers.keys() else None


def default_arrow(size):
    if isinstance(size, (list, tuple)) and len(size) == 3:
        head_l, head_w, ax_l = size
    elif type(size) == int or type(size) == float:
        head_l, head_w, ax_l = 12 * size, 10 * size, 8 * size
    else:
        head_l, head_w, ax_l = 12, 10, 8
    return head_l, head_w, ax_l


"""set axes parameters (ticks, frame, labels, title, """


def update_axes(ax, xlim=None, ylim=None, fontsize=None, is_embedding=False, frameon=None, figsize=None):
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    frameon = settings._frameon if frameon is None else frameon
    if isinstance(frameon, str) and frameon == 'artist':
        set_artist_frame(ax, figsize=figsize)
    elif frameon:
        if is_embedding:
            ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            labelsize = int(fontsize * .75) if fontsize is not None else None
            ax.tick_params(axis='both', which='major', labelsize=labelsize)
        if isinstance(frameon, str) and frameon != 'full':
            frameon = 'bl' if frameon == 'half' else frameon
            bf, lf, tf, rf = [f in frameon for f in ['bottom', 'left', 'top', 'right']]
            if not np.any([bf, lf, tf, rf]):
                bf, lf, tf, rf = [f in frameon for f in ['b', 'l', 't', 'r']]
            ax.spines['top'].set_visible(tf)
            ax.spines['right'].set_visible(rf)
            if not bf:
                ax.set_xlabel('')
                ax.spines['bottom'].set_visible(False)
            if not lf:
                ax.set_ylabel('')
                ax.spines['left'].set_visible(False)
            ax.tick_params(which='both', bottom=bf, left=lf, top=tf, right=rf, labelbottom=bf, labelleft=lf)

    else:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)
        ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        ax.set_frame_on(False)

    if rcParams['savefig.transparent']:
        ax.patch.set_alpha(0)


def set_artist_frame(ax, length=.2, figsize=None):
    ax.tick_params(axis='both', which='both', labelbottom=False, labelleft=False,
                   bottom=False, left=False, top=False, right=False)
    for side in ['bottom', 'right', 'top', 'left']:
        ax.spines[side].set_visible(False)
    figsize = rcParams['figure.figsize'] if figsize is None else figsize
    aspect_ratio = figsize[0] / figsize[1]
    ax.xaxis.set_label_coords(length * .45, -0.035)
    ax.yaxis.set_label_coords(-0.0175, length * aspect_ratio * .45)
    ax.xaxis.label.set_size(ax.xaxis.label.get_size() / 1.2)
    ax.yaxis.label.set_size(ax.yaxis.label.get_size() / 1.2)

    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDirectionArrows as arrows
    arr = arrows(ax.transAxes, 'none', 'none', loc=3, pad=-1, back_length=0, aspect_ratio=aspect_ratio,
                 arrow_props={'ec': 'none', 'fc': 'k'}, fontsize=0, length=length,
                 text_props={'ec': 'k', 'fc': 'k', 'lw': .1})
    ax.add_artist(arr)


def set_label(xlabel, ylabel, fontsize=None, basis=None, ax=None, **kwargs):
    labels, labels_new = np.array(['Ms', 'Mu', 'X']), np.array(['spliced', 'unspliced', 'expression'])
    if xlabel in labels: xlabel = labels_new[xlabel == labels][0]
    if ylabel in labels: ylabel = labels_new[ylabel == labels][0]
    if ax is None: ax = pl.gca()
    if basis is not None:
        component_name = ('DC' if 'diffmap' in basis else 'tSNE' if basis == 'tsne' else 'UMAP' if basis == 'umap'
        else 'PC' if basis == 'pca' else basis.replace('draw_graph_', '').upper() if 'draw_graph' in basis else basis)
        ax.set_xlabel(component_name + '1', fontsize=fontsize, **kwargs)
        ax.set_ylabel(component_name + '2', fontsize=fontsize, **kwargs)
    if isinstance(xlabel, str):
        ax.set_xlabel(xlabel.replace('_', ' '), fontsize=fontsize, **kwargs)
    if isinstance(ylabel, str):
        ax.set_ylabel(ylabel.replace('_', ' '), fontsize=fontsize,
                      rotation=0 if ylabel.startswith('$') or len(ylabel) == 1 else 90, **kwargs)


def set_title(title, layer=None, color=None, fontsize=None, ax=None):
    if ax is None: ax = pl.gca()
    color = color if isinstance(color, str) and not is_color_like(color) else None
    title = title if isinstance(title, str) and not is_color_like(title) else None
    if isinstance(title, str):
        title = title.replace('_', ' ')
    elif isinstance(layer, str) and isinstance(color, str):
        title = (color + ' ' + layer).replace('_', ' ')
    elif isinstance(color, str):
        title = color.replace('_', ' ')
    else:
        title = ''
    ax.set_title(title, fontsize=fontsize)


def set_frame(ax, frameon):
    frameon = settings._frameon if frameon is None else frameon
    if not frameon:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_frame_on(False)


def set_legend(adata, ax, value_to_plot, legend_loc, scatter_array, legend_fontweight, legend_fontsize,
                legend_fontoutline, groups):
    """
    Adds a legend to the given ax with categorial data.
    """
    # add legend
    if legend_fontoutline is None:
        legend_fontoutline = 1
    obs_vals = adata.obs[value_to_plot]
    obs_vals.cat.categories = obs_vals.cat.categories.astype(str)
    color_keys = adata.uns[value_to_plot + '_colors']
    if isinstance(color_keys, dict):
        color_keys = np.array([color_keys[c] for c in obs_vals.cat.categories])
    valid_cats = np.where(obs_vals.value_counts()[obs_vals.cat.categories] > 0)[0]
    categories = np.array(obs_vals.cat.categories)[valid_cats]
    colors = np.array(color_keys)[valid_cats]

    if groups is not None:
        groups, groupby = get_groups(adata, groups, value_to_plot)
        # only label groups with the respective color
        groups = [g for g in groups if g in categories]
        colors = [colors[list(categories).index(x)] for x in groups]
        categories = groups

    if legend_loc == 'on data':
        legend_fontweight = 'bold' if legend_fontweight is None else legend_fontweight
        # identify centroids to put labels
        texts = []
        for ilabel, label in enumerate(categories):
            x_pos, y_pos = np.nanmedian(scatter_array[obs_vals == label, :], axis=0)
            if isinstance(label, str): label = label.replace('_', ' ')
            text = ax.text(x_pos, y_pos, label, weight=legend_fontweight, verticalalignment='center',
                           horizontalalignment='center', fontsize=legend_fontsize,
                           path_effects=[patheffects.withStroke(linewidth=legend_fontoutline, foreground='w')])
            texts.append(text)

        # todo: adjust text positions to minimize overlaps, e.g. using https://github.com/Phlya/adjustText
        # from adjustText import adjust_text
        # adjust_text(texts, ax=ax)

    else:
        for idx, label in enumerate(categories):
            if isinstance(label, str): label = label.replace('_', ' ')
            ax.scatter([], [], c=[colors[idx]], label=label)
        ncol = (1 if len(categories) <= 14 else 2 if len(categories) <= 30 else 3)
        if legend_loc == 'upper right':
            ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1), fontsize=legend_fontsize, ncol=ncol)
        elif legend_loc == 'lower right':
            ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(1, 0), fontsize=legend_fontsize, ncol=ncol)
        elif 'right' in legend_loc:  # 'right', 'center right', 'right margin'
            ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, .5), fontsize=legend_fontsize, ncol=ncol)
        elif legend_loc != 'none':
            ax.legend(frameon=False, loc=legend_loc, fontsize=legend_fontsize, ncol=ncol)


"""get color values"""


def clip(c, perc):
    if np.size(perc) < 2: perc = [perc, 100] if perc < 50 else [0, perc]
    lb, ub = np.percentile(c, perc)
    return np.clip(c, lb, ub)


def get_colors(adata, c):
    if is_color_like(c):
        return c
    else:
        if c + '_colors' not in adata.uns.keys():
            palette = default_palette(None)
            palette = adjust_palette(palette, length=len(adata.obs[c].cat.categories))
            adata.uns[c + '_colors'] = palette[:len(adata.obs[c].cat.categories)].by_key()['color']
        cluster_ix = adata.obs[c].values if isinstance(adata.uns[c + '_colors'], dict) else adata.obs[c].cat.codes.values
        return np.array([adata.uns[c + '_colors'][cluster_ix[i]]
                         if cluster_ix[i] != -1 else 'lightgrey' for i in range(adata.n_obs)])


def interpret_colorkey(adata, c=None, layer=None, perc=None, use_raw=None):
    if c is None: c = default_color(adata)
    if issparse(c): c = make_dense(c).flatten()
    if is_categorical(adata, c): c = get_colors(adata, c)
    elif isinstance(c, str):
        if c in adata.obs.keys():  # color by observation key
            c = adata.obs[c]
        elif c in adata.var_names or (use_raw and adata.raw is not None and c in adata.raw.var_names):  # by gene
            if layer in adata.layers.keys():
                if perc is None and any(l in layer for l in ['spliced', 'unspliced', 'Ms', 'Mu', 'velocity']):
                    perc = [1, 99]  # clip values to ignore extreme outliers since these layers are not logarithmized
                c = adata.obs_vector(c, layer=layer)
            elif layer is not None and np.any([l in layer or 'X' in layer for l in adata.layers.keys()]):
                l_array = np.hstack([adata.obs_vector(c, layer=l)[:, None] for l in adata.layers.keys()])
                l_array = pd.DataFrame(l_array, columns=adata.layers.keys())
                l_array.insert(0, 'X', adata.obs_vector(c))
                c = np.array(l_array.astype(np.float32).eval(layer))
            else:
                if layer is not None and layer != 'X':
                    logg.warn(layer, 'not found. Using .X instead.')
                if adata.raw is None and use_raw:
                    raise ValueError("`use_raw` is set to True but AnnData object does not have raw. Please check.")
                c = adata.raw.obs_vector(c) if use_raw else adata.obs_vector(c)
            c = c.A.flatten() if issparse(c) else c
        elif c in adata.var.keys():  # color by observation key
            c = adata.var[c]
        elif np.any([var_key in c for var_key in adata.var.keys()]):
            var = adata.var[[key for key in adata.var.keys() if not isinstance(adata.var[key][0], str)]]
            c = var.astype(np.float32).eval(c)
        elif np.any([obs_key in c for obs_key in adata.obs.keys()]):
            obs = adata.obs[[key for key in adata.obs.keys() if not isinstance(adata.obs[key][0], str)]]
            c = obs.astype(np.float32).eval(c)
        elif not is_color_like(c):
            raise ValueError('color key is invalid! pass valid observation annotation or a gene name')
        if not isinstance(c, str) and perc is not None: c = clip(c, perc=perc)
    else:
        c = np.array(c).flatten()
        if perc is not None: c = clip(c, perc=perc)
    return c


# colors in addition to matplotlib's colors
additional_colors = {
    'gold2': '#eec900', 'firebrick3': '#cd2626', 'khaki2': '#eee685', 'slategray3': '#9fb6cd', 'palegreen3': '#7ccd7c',
    'tomato2': '#ee5c42', 'grey80': '#cccccc', 'grey90': '#e5e5e5', 'wheat4': '#8b7e66', 'grey65': '#a6a6a6',
    'grey10': '#1a1a1a', 'grey20': '#333333', 'grey50': '#7f7f7f', 'grey30': '#4d4d4d', 'grey40': '#666666',
    'antiquewhite2': '#eedfcc', 'grey77': '#c4c4c4', 'snow4': '#8b8989', 'chartreuse3': '#66cd00', 'yellow4': '#8b8b00',
    'darkolivegreen2': '#bcee68', 'olivedrab3': '#9acd32', 'azure3': '#c1cdcd', 'violetred': '#d02090',
    'mediumpurple3': '#8968cd', 'purple4': '#551a8b', 'seagreen4': '#2e8b57', 'lightblue3': '#9ac0cd',
    'orchid3': '#b452cd', 'indianred 3': '#cd5555', 'grey60': '#999999', 'mediumorchid1': '#e066ff', 'plum3': '#cd96cd',
    'palevioletred3': '#cd6889'}


# adapted from scanpy
def set_colors_for_categorical_obs(adata, value_to_plot, palette=None):
    """
    Sets the adata.uns[value_to_plot + '_colors'] according to the given palette or default colors.
    Parameters
    ----------
    adata
        annData object
    value_to_plot
        name of a valid categorical observation
    palette
        Palette should be either a valid :func:`~matplotlib.pyplot.colormaps` string,
        a sequence of colors (in a format that can be understood by matplotlib,
        eg. RGB, RGBS, hex, or a cycler object with key='color'
    Returns
    -------
    None
    """
    from matplotlib.colors import to_hex
    color_key = f"{value_to_plot}_colors"
    valid = True
    categories = adata.obs[value_to_plot].cat.categories
    length = len(categories)

    if isinstance(palette, str) and 'default' in palette:
        palette = palettes.default_26 if length <= 28 else palettes.default_64
    if isinstance(palette, str) and palette in adata.uns:
        palette = adata.uns[palette].values() if isinstance(adata.uns[palette], dict) else adata.uns[palette]

    if palette is None and color_key in adata.uns:
        # Check if colors already exist in adata.uns and if they are a valid palette
        _palette = []
        color_keys = adata.uns[color_key].values() if isinstance(adata.uns[color_key], dict) else adata.uns[color_key]
        for color in color_keys:
            if not is_color_like(color):
                # check if the color is a valid R color and translate it to a valid hex color value
                if color in additional_colors:
                    color = additional_colors[color]
                else:
                    logg.warn(
                        f"The following color value found in adata.uns['{value_to_plot}_colors'] "
                        f"is not valid: '{color}'. Default colors will be used instead."
                    )
                    valid = False
                    break
        if len(adata.uns[color_key]) < len(adata.obs[value_to_plot].cat.categories):
            valid = False
    elif palette is not None:
        # check is palette given is a valid matplotlib colormap
        if isinstance(palette, str) and palette in pl.colormaps():
            # this creates a palette from a colormap. E.g. 'Accent, Dark2, tab20'
            cmap = pl.get_cmap(palette)
            colors_list = [to_hex(x) for x in cmap(np.linspace(0, 1, len(categories)))]

        else:
            # check if palette is as dict and convert it to an ordered list
            if isinstance(palette, dict):
                palette = [palette[c] for c in categories]
            # check if palette is a list and convert it to a cycler
            if isinstance(palette, abc.Sequence):
                if len(palette) < len(categories):
                    logg.warn(
                        "Length of palette colors is smaller than the number of "
                        f"categories (palette length: {len(palette)}, "
                        f"categories length: {len(categories)}. "
                        "Some categories will have the same color.")
                # check that colors are valid
                _color_list = []
                for color in palette:
                    if not is_color_like(color):
                        # check if the color is a valid R color and translate it to a valid hex color value
                        if color in additional_colors:
                            color = additional_colors[color]
                        else:
                            logg.warn(
                                f"The following color value is not valid: '{color}'. "
                                f"Default colors will be used instead.")
                            valid = False
                            break
                    _color_list.append(color)
                palette = cycler(color=_color_list)

            if not isinstance(palette, Cycler) or 'color' not in palette.keys:
                logg.warn(
                    "Please check that the value of 'palette' is a valid "
                    "matplotlib colormap string (eg. Set2), a  list of color names "
                    "or a cycler with a 'color' key. Default colors will be used instead"
                )
                valid = False

            if valid:
                cc = palette()
                colors_list = [to_hex(next(cc)['color']) for x in range(len(categories))]
        if valid:
            adata.uns[value_to_plot + '_colors'] = colors_list
    else:
        # No valid palette exists or was given
        valid = False

    # Set to defaults:
    if not valid:
        # check if default matplotlib palette has enough colors
        if len(rcParams['axes.prop_cycle'].by_key()['color']) >= length:
            cc = rcParams['axes.prop_cycle']()
            palette = [next(cc)['color'] for _ in range(length)]
        # Else fall back to default palettes
        else:
            if length <= 28:
                palette = palettes.default_26
            elif length <= len(palettes.default_64):  # 103 colors
                palette = palettes.default_64
            else:
                palette = ['grey' for _ in range(length)]
                logg.info(
                    f'the obs value {value_to_plot!r} has more than 103 categories. Uniform '
                    "'grey' color will be used for all categories."
                )

        adata.uns[value_to_plot + '_colors'] = palette[:length]


def set_colorbar(smp, ax, orientation='vertical', labelsize=None):
    cb = pl.colorbar(smp, orientation=orientation, cax=inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0))
    cb.set_alpha(1)
    cb.ax.tick_params(labelsize=labelsize)
    cb.draw_all()
    cb.locator = MaxNLocator(nbins=3, integer=True)
    cb.update_ticks()


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
    if 'transparent' in colors:
        alpha = [1 if i != 'transparent' else 0 for i in colors] if alpha is None else alpha
        colors = [i if i != 'transparent' else 'white' for i in colors]

    for color in colors:
        if isinstance(color, str):
            color = to_rgb(color if color.startswith('#') else cnames[color])
            c.append(color)
    alpha = np.ones(len(c)) if alpha is None else alpha

    vals = np.ones((N, 4))
    ints = len(c) - 1
    n = int(N / ints)

    for j in range(ints):
        for i in range(3):
            vals[n * j:n * (j + 1), i] = np.linspace(c[j][i], c[j + 1][i], n)
        vals[n * j:n * (j + 1), -1] = np.linspace(alpha[j], alpha[j + 1], n)
    return ListedColormap(vals)


"""save figure"""


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
        if len(settings.figdir) > 0:
            if settings.figdir[-1] != '/': settings.figdir += '/'
            if not os.path.exists(settings.figdir): os.makedirs(settings.figdir)
        if ext is None: ext = settings.file_format_figs
        try:
            filename = settings.figdir + f'{settings.plot_prefix}{writekey}{settings.plot_suffix}.{ext}'
            pl.savefig(filename, dpi=dpi, bbox_inches='tight')
        except:  # save as .png if it cannot be saved as .pdf (e.g. streamplot under specific settings)
            logg.msg(f'figure cannot be saved as {ext}, using png instead.', v=1)
            filename = settings.figdir + f'{settings.plot_prefix}{writekey}{settings.plot_suffix}.png'
            pl.savefig(filename, dpi=dpi, bbox_inches='tight')
        logg.msg('saving figure to file', filename, v=1)

    if show: pl.show()
    if save: pl.close()  # clear figure


"""additional plots (linear fit, density, outline, rug)"""


def plot_linfit(x, y, add_linfit=True, add_legend=True, color=None, linewidth=None, fontsize=None, ax=None):
    if ax is None: ax = pl.gca()
    idx_valid = ~np.isnan(x + y)
    x, y = x[idx_valid], y[idx_valid]

    mu_x, mu_y = (0, 0) if isinstance(add_linfit, str) and 'no_intercept' in add_linfit else (np.mean(x), np.mean(y))
    slope = (np.mean(x * y) - mu_x * mu_y) / (np.mean(x ** 2) - mu_x ** 2)
    offset = mu_y - slope * mu_x

    if isinstance(add_linfit, str) and 'intercept' in add_linfit:
        str_list = [a.strip() for a in add_linfit.split(',')]
        add_linfit = [a for a in str_list if 'intercept' not in a][0] if ',' in add_linfit else None

    color = add_linfit if isinstance(add_linfit, str) else color if isinstance(color, str) else 'grey'
    xnew = np.linspace(np.min(x), np.max(x) * 1.02)
    ax.plot(xnew, offset + xnew * slope, linewidth=linewidth, color=color)
    if add_legend:
        ax.text(.05, .95, r'$\rho = $' + str(np.round(pearsonr(x, y)[0], 2)), ha='left', va='top', fontsize=fontsize,
                transform=ax.transAxes, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))


def plot_polyfit(x, y, add_polyfit=True, add_legend=True, color=None, linewidth=None, fontsize=None, ax=None):
    if ax is None: ax = pl.gca()
    idx_valid = ~np.isnan(x + y)
    x, y = x[idx_valid], y[idx_valid]
    lb, ub = np.percentile(x, [.5, 99.5])
    idx = (x > lb) & (x < ub)
    x, y = x[idx], y[idx]
    if isinstance(add_polyfit, str) and 'no_intercept' in add_polyfit:
        x, y = np.hstack([np.zeros(len(x)), x]), np.hstack([np.zeros(len(x)), y])
    fit = np.polyfit(x, y, deg=2 if isinstance(add_polyfit, (str, bool)) else add_polyfit)
    f = np.poly1d(fit)

    if isinstance(add_polyfit, str) and 'intercept' in add_polyfit:
        str_list = [a.strip() for a in add_polyfit.split(',')]
        add_polyfit = [a for a in str_list if 'intercept' not in a][0] if ',' in add_polyfit else None
    color = add_polyfit if isinstance(add_polyfit, str) else color if isinstance(color, str) else 'grey'

    xnew = np.linspace(np.min(x), np.max(x), num=100)
    ax.plot(xnew, f(xnew), color=color, linewidth=linewidth)

    if add_legend:
        R2 = np.sum((f(x) - np.mean(y)) ** 2) / np.sum((y - np.mean(y)) ** 2)
        ax.text(.05, .95, r'$R^2 = $' + str(np.round(R2, 2)), ha='left', va='top', fontsize=fontsize,
                transform=ax.transAxes, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))


def plot_vlines(adata, basis, vkey, xkey, linewidth=1, linecolor=None, ax=None):
    if ax is None: ax = pl.gca()
    xnew = np.linspace(0, np.percentile(make_dense(adata[:, basis].layers[xkey]), 98))
    fits = [fit for fit in make_unique_list(vkey) if all(['velocity' in fit, fit + '_gamma' in adata.var.keys()])]
    linecolor, lines = to_list(linecolor), []
    for i, fit in enumerate(fits):
        linestyle = '--' if 'variance_' + fit in adata.layers.keys() else '-'
        gamma = adata[:, basis].var[fit + '_gamma'].values if fit + '_gamma' in adata.var.keys() else 1
        beta = adata[:, basis].var[fit + '_beta'].values if fit + '_beta' in adata.var.keys() else 1
        offset = adata[:, basis].var[fit + '_offset'].values if fit + '_offset' in adata.var.keys() else 0
        line, = ax.plot(xnew, gamma / beta * xnew + offset / beta, linestyle=linestyle, linewidth=linewidth,
                        c=linecolor[i] if len(linecolor) > i and linecolor[i] is not None else 'k' if i == 0 else None)
        lines.append(line)
        fits[i] = 'steady-state ratio ({})'.format(fit) if len(fits) > 1 else 'steady-state ratio'
    return lines, fits


def plot_velocity_fits(adata, basis, vkey=None, use_raw=None, linewidth=None, linecolor=None, legend_loc=None,
                       legend_fontsize=None, show_assignments=None, ax=None):
    if ax is None: ax = pl.gca()
    if use_raw is None: use_raw = 'Ms' not in adata.layers.keys()

    # linear fits
    if vkey is None: vkey = 'dynamics' if 'fit_alpha' in adata.var.keys() else 'velocity'
    lines, fits = plot_vlines(adata, basis, vkey, 'spliced' if use_raw else 'Ms', linewidth, linecolor, ax=ax)

    # full dynamic fits
    from .simulation import show_full_dynamics
    if 'true_alpha' in adata.var.keys() and (vkey is not None and 'true_dynamics' in vkey):
        line, fit = show_full_dynamics(adata, basis, 'true', use_raw, linewidth, ax=ax)
        fits.append(fit); lines.append(line)
    if 'fit_alpha' in adata.var.keys() and (vkey is None or 'dynamics' in vkey):
        line, fit = show_full_dynamics(adata, basis, 'fit', use_raw, linewidth, show_assignments=show_assignments, ax=ax)
        fits.append(fit); lines.append(line)

    if legend_loc in {True, None, 'bottom right'}:
        legend_loc = 'lower right'
    if len(fits) > 0 and legend_loc and legend_loc != 'none':
        ax.legend(handles=lines, labels=fits, fontsize=legend_fontsize, loc=legend_loc)


def plot_density(x, y=None, add_density=True, eval_pts=50, scale=10, alpha=.3, color='grey', ax=None):
    from scipy.stats import gaussian_kde as kde
    if ax is None: ax = pl.gca()
    color = add_density if isinstance(add_density, str) else color if isinstance(color, str) else 'grey'

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

    ax.plot(xnew, vals, color=color)
    ax.fill_between(xnew, -offset, vals, alpha=alpha, color=color)
    ax.set_ylim(-offset)

    # density plot along y-coordinate
    if y is not None:
        ynew = np.linspace(min(y), max(y), eval_pts)
        vals = kde(y)(ynew)

        offset = max(x) / scale
        scale_u = offset / np.max(vals)
        offset *= 1.3
        vals = vals * scale_u - offset

        ax.plot(vals, ynew, color=color)
        ax.fill_betweenx(ynew, -offset, vals, alpha=alpha, color=color)
        ax.set_xlim(-offset)


def plot_outline(x, y, kwargs, outline_width=None, outline_color=None, zorder=None, ax=None):
    # Adapted from scanpy. The default outline is a black edge followed by a thin white edged added around connected
    # clusters. Three overlapping scatter plots are drawn: First black dots with slightly larger size, then, white dots
    # a bit smaller, but still larger than the final dots. Then the final dots are drawn with some transparency.
    if ax is None: ax = pl.gca()

    bg_width, gp_width = (0.3, 0.05) if outline_width is None else outline_width
    bg_color, gp_color = ('black', 'white') if outline_color is None else outline_color

    s = kwargs.pop('s')
    point = np.sqrt(s)

    gp_size = (2 * (point * gp_width) + point) ** 2
    bg_size = (2 * (point * bg_width) + np.sqrt(gp_size)) ** 2

    kwargs['edgecolor'] = 'none'
    zord = 0 if zorder is None else zorder
    if 'rasterized' not in kwargs: kwargs['rasterized'] = settings._vector_friendly
    ax.scatter(x, y, s=bg_size, marker=".", c=bg_color, zorder=zord - 2, **kwargs)
    ax.scatter(x, y, s=gp_size, marker=".", c=gp_color, zorder=zord - 1, **kwargs)
    # restore size
    kwargs['s'] = s


def plot_rug(x, height=.03, color=None, ax=None, **kwargs):
    if ax is None: ax = pl.gca()
    x = np.asarray(x)

    transform = tx.blended_transform_factory(ax.transData, ax.transAxes)
    line_segs = np.column_stack([np.repeat(x, 2), np.tile([0, height], len(x))]).reshape([len(x), 2, 2])
    ax.add_collection(LineCollection(line_segs, transform=transform, color=color, **kwargs))
    ax.autoscale_view(scalex=True, scaley=False)


"""for velocity_embedding"""


def velocity_embedding_changed(adata, basis, vkey):
    if 'X_' + basis not in adata.obsm.keys(): changed = False
    else:
        changed = vkey + '_' + basis not in adata.obsm_keys()
        if vkey + '_params' in adata.uns.keys():
            sett = adata.uns[vkey + '_params']
            changed |= 'embeddings' not in sett or basis not in sett['embeddings']
    return changed


"""additional plots (linear fit, density, outline, rug)"""


def hist(arrays, alpha=.5, bins=50, color=None, colors=None, labels=None, hist=None, kde=None, bw_method=None,
         xlabel=None, ylabel=None, xlim=None, ylim=None, cutoff=None, xscale=None, yscale=None, xticks=None, yticks=None,
         fontsize=None, legend_fontsize=None, figsize=None, normed=None, perc=None, exclude_zeros=None,
         axvline=None, axhline=None, pdf=None, ax=None, dpi=None, show=True, **kwargs):
    """\
    Plot a histogram.

    Arguments
    ---------
    arrays: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    alpha: `list`, `np.ndarray` or `None` (default: `None`)
        Alpha of the colors. Must be same length as colors.
    bins : `int` or `sequence` (default: `50`)
        If an integer is given, ``bins + 1`` bin edges are calculated and
        returned, consistent with `numpy.histogram`.
        If `bins` is a sequence, gives bin edges, including left edge of
        first bin and right edge of last bin.  In this case, `bins` is
        returned unmodified.
    colors: : `list` or `array` (default `['royalblue', 'white', 'forestgreen']`)
        List of colors, either as names or rgb values.
    labels : str or None (default: `None`)
            String, or sequence of strings to label clusters.
    hist: `bool` or `None` (default: `None`)
        Whether to show histogram.
    kde: `bool` or `None` (default: `None`)
        Whether to use kernel density estimation on data.
    bw_method : `str`, `scalar` or `callable`, (default: `None`)
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a
        scalar, this will be used directly as `kde.factor`.  If a callable,
        it should take a `gaussian_kde` instance as only parameter and
        return a scalar.  If None (default), nothing happens; the current
        `kde.covariance_factor` method is kept.
    xlabel: `str` (default: `None`)
        Label of x-axis.
    ylabel: `str` (default: `None`)
        Label of y-axis.
    xlim: tuple, e.g. [0,1] or `None` (default: `None`)
        Restrict x-limits of the axis.
    ylim: tuple, e.g. [0,1] or `None` (default: `None`)
        Restrict y-limits of the axis.
    cutoff: tuple, e.g. [0,1] or `float` or `None` (default: `None`)
        Bins will be cut off below and above the cutoff values.
    xscale: `log` or `None` (default: `None`)
        Scale of the x-axis.
    yscale: `log` or `None` (default: `None`)
        Scale of the y-axis.
    fontsize: `float` (default: `None`)
        Label font size.
    legend_fontsize: `int` (default: `None`)
        Legend font size.
    figsize: tuple (default: `(7,5)`)
        Figure size.
    normed: `bool` or `None` (default: `None`)
        Whether to normalize data.
    perc: tuple, e.g. [2,98] (default: `None`)
        Specify percentile for continuous coloring.
    exclude_zeros: `bool` or `None` (default: `None`)
        Whether to exclude zeros in data for the kde and hist plot.
    axvline: `float` or `None` (default: `None`)
        Plot a vertical line at the specified x-value.
    axhline `float` or `None` (default: `None`)
        Plot a horizontal line at the specified y-value.
    pdf: `str` or `None` (default: `None`)
        probability density function to be fitted. E.g. 'norm', 't', 'chi', 'beta', 'gamma', 'laplace' etc.
    ax: `matplotlib.Axes`, optional (default: `None`)
        A matplotlib axes object. Only works if plotting a single component.
    dpi: `int` (default: 80)
        Figure dpi.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """

    if ax is None:
        fig, ax = pl.subplots(figsize=figsize, dpi=dpi)

    arrays = [np.ravel(array) for array in arrays] if isinstance(arrays, (list, tuple)) or arrays.ndim > 1 else [arrays]
    if normed is None: normed = kde or pdf
    if hist is None: hist = not kde

    palette = default_palette(None).by_key()['color'][::-1]
    colors = to_list(colors) if color is None else to_list(color)
    colors = palette if colors is None or len(colors) < len(arrays) else colors

    masked_arrays = np.ma.masked_invalid(np.hstack(arrays)).compressed()
    bmin, bmax = masked_arrays.min(), masked_arrays.max()
    if xlim is not None:
        bmin, bmax = max(bmin, xlim[0]), min(bmax, xlim[1])
    elif perc is not None:
        if np.size(perc) < 2: perc = [perc, 100] if perc < 50 else [0, perc]
        bmin, bmax = np.nanpercentile(masked_arrays, perc)
    bins = np.arange(bmin, bmax + (bmax - bmin) / bins, (bmax - bmin) / bins)

    if xscale == 'log':
        bins = bins[bins > 0]
        bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

    if cutoff is not None:
        bins = bins[(bins > cutoff[0]) & (bins < cutoff[1])] if isinstance(cutoff, list) else bins[bins < cutoff]

    if isinstance(labels, str):
        labels = [labels]

    for i, x in enumerate(arrays):
        x_vals = np.array(x[np.isfinite(x)])
        if exclude_zeros: x_vals = np.array(x_vals[x_vals != 0])
        if kde:
            from scipy.stats import gaussian_kde
            kde_bins = gaussian_kde(x_vals, bw_method=bw_method)(bins)
            if not normed:
                kde_bins *= (bins[1] - bins[0]) * len(x_vals)
            ax.plot(bins, kde_bins, color=colors[i])
            ax.fill_between(bins, 0, kde_bins, alpha=.4, color=colors[i],
                            label=labels[i] if labels is not None else None)
            ylim = np.min(kde_bins) if ylim is None else ylim
        if hist:
            kwargs.update({'color': colors[i], 'label': labels[i] if labels is not None else None})
            ax.hist(x_vals, bins=bins, alpha=alpha, density=normed, **kwargs)
            #try:
            #    ax.hist(x_vals, bins=bins, alpha=alpha, normed=normed, **kwargs)
            #except:
            #    ax.hist(x_vals, bins=bins, alpha=alpha, **kwargs)
    set_label(xlabel if xlabel is not None else '', ylabel if xlabel is not None else '', fontsize=fontsize, ax=ax)

    if labels is not None:
        ax.legend(fontsize=legend_fontsize)

    if axvline: ax.axvline(axvline)
    if axhline: ax.axhline(axhline)

    if xscale is not None: ax.set_xscale(xscale)
    if yscale is not None: ax.set_yscale(yscale)

    update_axes(ax, xlim, ylim, fontsize, frameon=True)

    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    @pl.FuncFormatter
    def log_fmt(x, pos):
        return r'${%.2g}$' % (x)

    if xscale == 'log':
        if xticks is None:
            ticks = [a for a in [10**a for a in np.linspace(-10, 10, 21)] if bmin < a < bmax]
            ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(log_fmt)
        ax.minorticks_off()

    pdf = [pdf] if isinstance(pdf, str) else pdf
    if pdf is not None:
        fits = []
        for i, pd in enumerate(pdf):
            from scipy import stats
            xt = ax.get_xticks()
            xmin, xmax = min(xt), max(xt)
            lnspc = np.linspace(xmin, xmax, len(bins))

            if '(' in pd:  # used passed parameters
                args, pd = eval(pd[pd.rfind('('):]), pd[:pd.rfind('(')]
            else:  # fit parameters
                args = eval('stats.' + pd + '.fit(x_vals)')
            pd_vals = eval('stats.' + pd + '.pdf(lnspc, *args)')
            logg.info('Fitting', pd, np.round(args, 4), '.')
            fit = ax.plot(lnspc, pd_vals, label=pd, color=colors[i])
            fits.extend(fit)
        ax.legend(handles=fits, labels=pdf, fontsize=legend_fontsize)

    if rcParams['savefig.transparent']:
        ax.patch.set_alpha(0)

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


"""deprecated"""


def make_unique_list(key, allow_array=False):
    if isinstance(key, (Index, abc.KeysView)): key = list(key)
    is_list = isinstance(key, (list, tuple, np.record)) if allow_array else isinstance(key, (list, tuple, np.ndarray, np.record))
    is_list_of_str = is_list and all(isinstance(item, str) for item in key)
    return key if is_list_of_str else key if is_list and len(key) < 20 else [key]


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


def get_temporal_connectivities(adata, tkey, n_convolve=30):
    from ..tools.velocity_graph import vals_to_csr
    from ..tools.utils import normalize, get_indices

    #c_idx = get_indices(get_connectivities(adata, recurse_neighbors=True))[0]
    #c_idx = np.hstack([c_idx, np.linspace(0, len(c_idx) - 1, len(c_idx), dtype=int)[:, None]])

    rows, cols, vals, n_obs, n_convolve = [], [], [], len(tkey), int(n_convolve / 2)
    for i in range(n_obs):
        i_max = None if i + n_convolve >= n_obs else i + n_convolve
        i_min = np.max([0, i - n_convolve])
        t_idx = np.argsort(tkey)[i_min: i_max]  # temporal neighbourhood
        #t_idx = np.intersect1d(t_idx, c_idx[i])
        rows.extend(np.ones(len(t_idx), dtype=int) * np.argsort(tkey)[i])
        cols.extend(t_idx)
        vals.extend(np.ones(len(t_idx), dtype=int))

    c_conn = get_connectivities(adata, recurse_neighbors=True)
    t_conn = vals_to_csr(vals, rows, cols, shape=(n_obs, n_obs))
    return normalize(t_conn)  # normalize(t_conn.multiply(c_conn))



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
