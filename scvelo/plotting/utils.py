import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.colors import rgb2hex
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def get_colors(adata, basis):
    if basis+'_colors' in adata.uns.keys():
        clusters_colors = adata.uns[basis + '_colors']
    else:
        colormap = np.vstack((pl.cm.tab20b(np.linspace(0., 1, 20))[::2], pl.cm.tab20c(np.linspace(0, 1, 20))[1::2],
                              pl.cm.tab20b(np.linspace(0., 1, 20))[1::2], pl.cm.tab20c(np.linspace(0, 1, 20))[0::2]))
        clusters_colors = [rgb2hex(c) for c in colormap[:len(adata.obs[basis])]]
    cluster_ix = adata.obs[basis].cat.codes
    return np.array([clusters_colors[cluster_ix[i]] for i in range(adata.X.shape[0])])


def bound(c, perc=None):
    lb, ub = np.percentile(c, perc)
    return np.clip(c, lb, ub)


def interpret_colorkey(adata, c=None, layer=None, perc=None):
    if (c == 'louvain' or c == 'clusters') and adata.obs[c].dtype.name != 'category':
        adata.obs[c] = pd.Categorical(adata.obs[c])
    if isinstance(c, str):
        if c in adata.obs.keys():  # color by observation key
            c = get_colors(adata, basis=c) if adata.obs[c].dtype.name == 'category' else adata.obs[c]
        elif c in adata.var_names:  # color by var in specific layer
            c = adata[:, c].layers[layer] if layer in adata.layers.keys() else adata[:, c].X
        c = c if perc is None else bound(c, perc=perc)
    elif c is None:  # color by cluster or louvain or grey if no color is specified
        c = get_colors(adata, 'clusters') if 'clusters' in adata.obs.keys() \
            else get_colors(adata, 'louvain') if 'louvain' in adata.obs.keys() else 'grey'
    elif isinstance(c, np.ndarray) and len(c.flatten()) == len(adata.obs):  # continuous coloring
        c = c.flatten() if perc is None else bound(c.flatten(), perc=perc)
    else: raise ValueError('color key is invalid! pass valid observation annotation or a gene name')
    return c


def get_components(components=None):
    if components is None: components = '1,2'
    if isinstance(components, str): components = components.split(',')
    return np.array(components).astype(int) - 1


def quiver_autoscale(X, Y, U, V):
    Q = pl.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=None)
    Q._init()
    pl.clf()
    return Q.scale


def plot_colorbar(ax, orientation='vertical'):
    cb = pl.colorbar(orientation=orientation, cax=inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0))
    cb.locator = (MaxNLocator(nbins=3))
    cb.update_ticks()


def plot_legend(ax, legend_loc=None, legend_fontsize=None, right_margin=None):
    if legend_loc.startswith('on data') and legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']
    elif legend_fontsize is None:
        legend_fontsize = rcParams['legend.fontsize']
    if right_margin is None and len(categoricals) > 0:
        if legend_loc == 'right margin': right_margin = 0.5

    legend = None
    if legend_loc.startswith('on data'):
        if legend_fontweight is None:
            legend_fontweight = 'bold'
        for name, pos in centroids.items():
            ax.text(pos[0], pos[1], name, weight=legend_fontweight, verticalalignment='center',
                    horizontalalignment='center', fontsize=legend_fontsize)
        all_pos = np.zeros((len(adata.obs[key].cat.categories), 2))
        for iname, name in enumerate(adata.obs[key].cat.categories):
            if name in centroids:
                all_pos[iname] = centroids[name]
            else:
                all_pos[iname] = [np.nan, np.nan]
        utils._tmp_cluster_pos = all_pos
        if legend_loc == 'on data export':
            filename = settings.writedir + 'pos.csv'
            logg.msg('exporting label positions to {}'.format(filename), v=1)
            if settings.writedir != '' and not os.path.exists(settings.writedir):
                os.makedirs(settings.writedir)
            np.savetxt(filename, all_pos, delimiter=',')
    elif legend_loc == 'right margin':
        legend = ax.legend(
            frameon=False, loc='center left',
            bbox_to_anchor=(1, 0.5),
            ncol=(1 if len(adata.obs[key].cat.categories) <= 14
                  else 2 if len(adata.obs[key].cat.categories) <= 30 else 3),
            fontsize=legend_fontsize)



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