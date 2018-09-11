import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


zeileis_26 = [
    "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
    "#f3e1eb", "#f6c4e1", "#f79cd4", '#7f7f7f', "#c7c7c7", "#1CE6FF", "#336600"]


def get_colors(adata, basis):
    clusters_colors = adata.uns[basis + '_colors'] if basis+'_colors' in adata.uns.keys() else zeileis_26
    cluster_ix = adata.obs[basis].cat.codes
    return np.array([clusters_colors[cluster_ix[i]] for i in range(adata.n_obs)])


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