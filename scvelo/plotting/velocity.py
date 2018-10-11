from ..preprocessing.moments import second_order_moments
from ..tools import rank_velocity_genes
from .scatter import scatter
from .utils import savefig
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as pl
from scipy.sparse import issparse


def velocity(adata, var_names=None, basis='umap', mode=None, fits='all', layers='all', use_raw=False, color=None, color_map='RdBu_r',
             perc=[2,98], size=.2, alpha=.5, fontsize=8, figsize=(3,2), dpi=120, show=True, save=None, ax=None, **kwargs):
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
    var_names = [var_names] if isinstance(var_names, str) else var_names if var_names is not None \
        else rank_velocity_genes(adata, n_genes=4)
    var_names = pd.unique([var for var in var_names if var in adata.var_names[adata.var.velocity_genes]])

    (skey, ukey) = ('spliced', 'unspliced') if use_raw else ('Ms', 'Mu')
    layers = ['velocity', skey, 'variance_velocity'] if layers == 'all' else layers
    layers = [layer for layer in layers if layer in adata.layers.keys()]

    fits = adata.layers.keys() if fits == 'all' else fits
    fits = [fit for fit in fits if all(['velocity' in fit, fit + '_gamma' in adata.var.keys()])]

    n_row, n_col = len(var_names), (1 + len(layers) + (mode == 'stochastic')*2)

    ax = pl.figure(figsize=(figsize[0]*n_col, figsize[1]*n_row), dpi=dpi) if ax is None else ax
    gs = pl.GridSpec(n_row, n_col, wspace=0.3, hspace=0.5)

    for v, var in enumerate(var_names):
        ix = np.where(adata.var_names == var)[0][0]
        s, u = adata.layers[skey][:, ix], adata.layers[ukey][:, ix]
        if issparse(s): s, u = s.A, u.A

        # spliced/unspliced phase portrait with steady-state estimate
        ax = pl.subplot(gs[v * n_col])
        scatter(adata, x=s, y=u, color=color, frameon=True, title=var, xlabel='spliced', ylabel='unspliced',
                show=False, save=False, ax=ax, fontsize=fontsize, size=size, alpha=alpha, **kwargs)

        xnew = np.linspace(0, s.max() * 1.02)
        for fit in fits:
            linestyle = '--' if 'stochastic' in fit else '-'
            pl.plot(xnew, adata.var[fit + '_gamma'][ix] / adata.var[fit + '_beta'][ix] * xnew
                    + adata.var[fit + '_offset'][ix] / adata.var[fit + '_beta'][ix], c='k', linestyle=linestyle)
        if v == len(var_names)-1: pl.legend(fits, loc='lower right', prop={'size': .5*fontsize})

        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis='both', which='major', labelsize=.7*fontsize)

        # velocity and expression plots
        for l, layer in enumerate(layers):
            ax = pl.subplot(gs[v*n_col + l + 1])
            title = 'expression' if layer == skey else layer
            scatter(adata, basis=basis, color=var, layer=layer, color_map=color_map, title=title,
                    perc=perc, fontsize=fontsize, size=size, alpha=alpha, show=False, ax=ax, save=False, **kwargs)

        if mode == 'stochastic':
            ss, us = second_order_moments(adata)

            ax = pl.subplot(gs[v*n_col + len(layers) + 1])
            x = 2 * (ss - s**2) - s
            y = 2 * (us - u * s) + u + 2 * s * \
                adata.var['stochastic_velocity_offset'][ix] / adata.var['stochastic_velocity_beta'][ix]

            scatter(adata, x=x, y=y, color=color, title=var, fontsize=40/n_col, show=False, ax=ax, save=False,
                    perc=perc, xlabel=r'2 $\Sigma_s - \langle s \rangle$', ylabel=r'2 $\Sigma_{us} + \langle u \rangle$', **kwargs)

            xnew = np.linspace(x.min(), x.max() * 1.02)
            fits = adata.layers.keys() if fits == 'all' else fits
            fits = [fit for fit in fits if 'velocity' in fit]
            for fit in fits:
                linestyle = '--' if 'stochastic' in fit else '-'
                pl.plot(xnew, adata.var[fit + '_gamma'][ix] / adata.var[fit + '_beta'][ix] * xnew +
                        adata.var[fit + '_offset2'][ix] / adata.var[fit + '_beta'][ix], c='k', linestyle=linestyle)
            if v == len(var_names) - 1: pl.legend(fits, loc='lower right', prop={'size': 34/n_col})

    if isinstance(save, str): savefig('', dpi=dpi, save=save, show=show)

    if show: pl.show()
    else: return ax
