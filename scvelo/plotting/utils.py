import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import rgb2hex
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def get_colors(adata, basis):
    cluster_ix = adata.obs[basis].cat.codes
    if basis+'_colors' in adata.uns.keys():
        clusters_colors = adata.uns[basis + '_colors']
    else:
        colormap = np.vstack((pl.cm.tab20b(np.linspace(0., 1, 20))[::2], pl.cm.tab20c(np.linspace(0, 1, 20))[1::2],
                              pl.cm.tab20b(np.linspace(0., 1, 20))[1::2], pl.cm.tab20c(np.linspace(0, 1, 20))[0::2]))
        clusters_colors = [rgb2hex(c) for c in colormap[:len(adata.obs[basis])]]
    return np.array([clusters_colors[cluster_ix[i]] for i in range(adata.X.shape[0])])


def interpret_colorkey(adata, c=None, layer=None, perc=None):
    if perc is None: perc = [2, 98]
    if isinstance(c, str):
        if c in adata.obs.keys():  # color by observation key
            if adata.obs[c].dtype.name == 'category':
                c = get_colors(adata, basis=c)
            else:
                c = adata.obs[c]
                lb, ub = np.percentile(c, perc)
                c = np.clip(c, lb, ub)
        elif c in adata.var_names:  # color by var in specific layer
            c = adata[:, c].layers[layer] if layer in adata.layers.keys() else adata[:, c].X
            lb, ub = np.percentile(c, perc)
            c = np.clip(c, lb, ub)
    elif c is None:  # color by cluster or louvain or grey if no color is specified
        c = get_colors(adata, 'clusters') if 'clusters' in adata.obs_keys() \
            else get_colors(adata, 'louvain') if 'louvain' in adata.obs_keys() else 'grey'
    else:  # continuous coloring
        lb, ub = np.percentile(c, perc)
        c = np.clip(c, lb, ub)
    return c


def scatter(adata, x=None, y=None, basis='umap', layer=None, color=None, xlabel=None, ylabel=None, color_map=None,
            size=None, alpha=1, fontsize=None, frameon=False, title=None, show=True, colorbar=False, save=None, ax=None, **kwargs):
    """Scatter plot along observations or variables axes.
    Color the plot using annotations of observations (`.obs`), variables (`.var`) or expression of genes (`.var_names`).

    Arguments
    ---------
    adata: `AnnData`
        Annotated data matrix.

    basis: `str` (default='tsne')
        plots embedding obsm['X_' + basis]

    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate

    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate

    color : `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes

    Returns
    -------
        If `show==False` a `matplotlib.Axis`
    """
    if (x is None) | (y is None):
            X_emb = adata.obsm['X_' + basis][:, :2]
            x, y = X_emb[:, 0], X_emb[:, 1]

    pl.scatter(x, y, c=interpret_colorkey(adata, color, layer), cmap=color_map, s=size, alpha=alpha, **kwargs)

    if isinstance(title, str):
        pl.title(title, fontsize=fontsize)

    if isinstance(xlabel, str) or isinstance(ylabel, str):
        pl.xlabel(xlabel, fontsize=fontsize)
        pl.ylabel(ylabel, fontsize=fontsize)

    if not frameon: pl.axis('off')

    if isinstance(save, str):
        pl.savefig(save)

    if ax is not None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        labelsize = int(fontsize * .75) if fontsize is not None else None
        ax.tick_params(axis='both', which='major', labelsize=labelsize)

    if colorbar and (ax is not None):
        cb = pl.colorbar(orientation='vertical', cax= inset_axes(ax, width="2%", height="30%", loc=4, borderpad=0))
        cb.locator = (MaxNLocator(nbins=3))
        cb.update_ticks()

    if show: pl.show()
    else: return ax


def phase(adata, var=None, x=None, y=None, color='louvain', fits='all', xlabel='spliced', ylabel='unspliced',
          fontsize=None, show=True, ax=None, **kwargs):
    if isinstance(var, str) and (var in adata.var_names):
        if (x is None) or (y is None):
            ix = np.where(adata.var_names == var)[0][0]
            x, y = adata.layers['Ms'][:, ix], adata.layers['Mu'][:, ix]
    else:
        ValueError('var not found in adata.var_names.')

    ax = scatter(adata, x=x, y=y, color=color, frameon=True, title=var, xlabel=xlabel, ylabel=ylabel, ax=ax, **kwargs)

    xnew = np.linspace(0, x.max() * 1.02)
    fits = adata.layers.keys() if fits == 'all' else fits
    fits = [fit for fit in fits if 'velocity' in fit]
    for fit in fits:
        linestyle = '--' if 'stochastic' in fit else '-'
        pl.plot(xnew, adata.var[fit+'_gamma'][ix] / adata.var[fit+'_beta'][ix] * xnew
                + adata.var[fit+'_offset'][ix] / adata.var[fit+'_beta'][ix], c='k', linestyle=linestyle)

    if show: pl.show()
    else: return ax