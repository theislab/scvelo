from ..preprocessing.moments import second_order_moments
from ..tools.rank_velocity_genes import rank_velocity_genes
from .scatter import scatter
from .utils import savefig, default_basis, default_size

import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as pl
from scipy.sparse import issparse


def velocity(adata, var_names=None, basis=None, groupby=None, groups=None, mode=None, fits='all', layers='all',
             color=None, color_map='RdBu_r', perc=[2,98], use_raw=False, size=None, alpha=.5, fontsize=None,
             figsize=None, dpi=None, show=True, save=None, ax=None, ncols=None, **kwargs):
    """Phase and velocity plot for set of genes.

    The phase plot shows spliced against unspliced expressions with steady-state fit.
    Further the embedding is shown colored by velocity and expression.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names: `str` or list of `str` (default: `None`)
        Which variables to show.
    basis: `str` (default: `'umap'`)
        Key for embedding coordinates.
    mode: `'stochastic'` or `None` (default: `None`)
        Whether to show show covariability phase portrait.
    fits: `str` or list of `str` (default: `'all'`)
        Which steady-state estimates to show.
    layers: `str` or list of `str` (default: `'all'`)
        Which layers to show.
    color: `str`,  list of `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes
    color_map: `str` (default: `matplotlib.rcParams['image.cmap']`)
        String denoting matplotlib color map.
    perc: tuple, e.g. [2,98] (default: `None`)
        Specify percentile for continuous coloring.
    size: `float` (default: 5)
        Point size.
    alpha: `float` (default: 1)
        Set blending - 0 transparent to 1 opaque.
    fontsize: `float` (default: `None`)
        Label font size.
    figsize: tuple (default: `(7,5)`)
        Figure size.
    dpi: `int` (default: 80)
        Figure dpi.
    show: `bool`, optional (default: `None`)
        Show the plot, do not return axis.
    save: `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure. A string is appended to the default filename.
        Infer the filetype if ending on {'.pdf', '.png', '.svg'}.
    ax: `matplotlib.Axes`, optional (default: `None`)
        A matplotlib axes object. Only works if plotting a single component.

    """
    basis = default_basis(adata) if basis is None else basis

    if isinstance(groupby, str) and groupby in adata.obs.keys():
        if 'rank_velocity_genes' not in adata.uns.keys() or adata.uns['rank_velocity_genes']['params']['groupby'] != groupby:
            rank_velocity_genes(adata, vkey='velocity', n_genes=10, groupby=groupby)
        names = np.array(adata.uns['rank_velocity_genes']['names'].tolist())
        if groups is None:
            var_names = names[:, 0]
        else:
            groups = [groups] if isinstance(groups, str) else groups
            idx = np.array([any([g in group for g in groups]) for group in adata.obs[groupby].cat.categories])
            var_names = np.hstack(names[idx, :int(10 / idx.sum())])
    elif var_names is not None:
        var_names = [var_names] if isinstance(var_names, str) else [var for var in var_names if var in adata.var_names]
    else:
        raise ValueError('No var_names or groups specified.')
    var_names = pd.unique(var_names)

    (skey, ukey) = ('spliced', 'unspliced') if use_raw else ('Ms', 'Mu')
    layers = ['velocity', skey, 'variance_velocity'] if layers == 'all' else layers
    layers = [layer for layer in layers if layer in adata.layers.keys()]

    fits = adata.layers.keys() if fits == 'all' else fits
    fits = [fit for fit in fits if all(['velocity' in fit, fit + '_gamma' in adata.var.keys()])]
    stochastic_fits = [fit for fit in fits if 'variance_' + fit in adata.layers.keys()]

    nplts = (1 + len(layers) + (mode == 'stochastic') * 2)
    ncols = 1 if ncols is None else ncols
    nrows = int(np.ceil(len(var_names) / ncols))
    ncols = int(ncols * nplts)
    figsize = rcParams['figure.figsize'] if figsize is None else figsize
    #n_row, n_col = len(var_names), (1 + len(layers) + (mode == 'stochastic')*2)
    ax = pl.figure(figsize=(figsize[0] * ncols / 2, figsize[1] * nrows / 2), dpi=dpi) if ax is None else ax
    gs = pl.GridSpec(nrows, ncols, wspace=0.3, hspace=0.5)

    size = default_size(adata) / 2 if size is None else size  # since fontsize is halved in width and height
    fontsize = rcParams['font.size'] if fontsize is None else fontsize
    for v, var in enumerate(var_names):
        _adata = adata[:, var]
        s, u = _adata.layers[skey], _adata.layers[ukey]
        if issparse(s): s, u = s.A, u.A

        # spliced/unspliced phase portrait with steady-state estimate
        ax = pl.subplot(gs[v * nplts])
        scatter(adata, basis=var, color=color, frameon=True, title=var, size=size, alpha=alpha, fontsize=fontsize,
                xlabel='spliced', ylabel='unspliced', show=False, ax=ax, save=False, **kwargs)

        xnew = np.linspace(0, s.max() * 1.02)
        for fit in fits:
            linestyle = '--' if fit in stochastic_fits else '-'
            gamma = _adata.var[fit + '_gamma'].values if fit + '_gamma' in adata.var.keys() else 1
            beta = _adata.var[fit + '_beta'].values if fit + '_beta' in adata.var.keys() else 1
            offset = _adata.var[fit + '_offset'].values if fit + '_offset' in adata.var.keys() else 0
            pl.plot(xnew, gamma / beta * xnew + offset / beta, c='k', linestyle=linestyle)
        if v == len(var_names)-1: pl.legend(fits, loc='lower right', prop={'size': .5*fontsize})

        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis='both', which='major', labelsize=.7*fontsize)

        # velocity and expression plots
        for l, layer in enumerate(layers):
            ax = pl.subplot(gs[v * nplts + l + 1])
            title = 'expression' if layer == skey else layer
            scatter(adata, basis=basis, color=var, layer=layer, color_map=color_map, title=title, perc=perc,
                    fontsize=fontsize, size=size, alpha=alpha, frameon=False, show=False, ax=ax, save=False, **kwargs)

        if mode == 'stochastic':
            ss, us = second_order_moments(_adata)
            ss, us = ss.flatten(), us.flatten()
            fit = stochastic_fits[0]

            ax = pl.subplot(gs[v * nplts + len(layers) + 1])
            offset = _adata.var[fit + '_offset'] if fit + '_offset' in adata.var.keys() else 0
            beta = _adata.var[fit + '_beta'] if fit + '_beta' in adata.var.keys() else 1
            x = 2 * (ss - s**2) - s
            y = 2 * (us - u * s) + u + 2 * s * offset / beta

            scatter(adata, x=x, y=y, color=color, title=var, fontsize=40/ncols, size=size, perc=perc, show=False,
                    xlabel=r'2 $\Sigma_s - \langle s \rangle$', ylabel=r'2 $\Sigma_{us} + \langle u \rangle$',
                    frameon=True, ax=ax, save=False, **kwargs)

            xnew = np.linspace(x.min(), x.max() * 1.02)
            for fit in stochastic_fits:
                gamma = _adata.var[fit + '_gamma'].values if fit + '_gamma' in adata.var.keys() else 1
                beta = _adata.var[fit + '_beta'].values if fit + '_beta' in adata.var.keys() else 1
                offset2 = _adata.var[fit + '_offset2'].values if fit + '_offset2' in adata.var.keys() else 0

                pl.plot(xnew, gamma / beta * xnew + offset2 / beta, c='k', linestyle='--')
            if v == len(var_names) - 1: pl.legend(fits, loc='lower right', prop={'size': 34/ncols})

    if isinstance(save, str): savefig('', dpi=dpi, save=save, show=show)

    if show: pl.show()
    else: return ax
