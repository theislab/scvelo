from .utils import *
from ..preprocessing.moments import second_order_moments


def velocity(adata, var_names=None, basis='umap', mode='deterministic', fits='all', layers='all', color=None,
             fontsize=8, color_map='RdBu_r', size=.2, alpha=.5, ax=None, **kwargs):
    """Phase and velocity plot for set of genes.

    The phase plot shows pliced against unspliced expressions with steady-state fit.
    Further the embedding is shown colored by velocity and expression.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    var_names: `str` or list of `str` (default: `None`)
        Which variables to show.
    basis: `str` (default: `'umap'`)
        Key for embedding coordinates.
    """
    var_names = [var_names] if isinstance(var_names, str) else var_names
    var_names = [var for var in var_names if var in adata.var_names] if isinstance(var_names, list) \
        else adata.var_names[(adata.layers['spliced'] > 0).sum(0).argsort()[::-1][:4]]

    layers = ['velocity', 'Ms', 'variance_velocity'] if layers == 'all' else layers
    layers = [layer for layer in layers if layer in adata.layers.keys()]

    fits = adata.layers.keys() if fits == 'all' else fits
    fits = [fit for fit in fits if all(['velocity' in fit, fit + '_gamma' in adata.var.keys()])]

    n_row, n_col = len(var_names), (1 + len(layers) + (mode == 'stochastic')*2)

    ax = pl.figure(figsize=(3*n_col, 2*n_row), dpi=120) if ax is None else ax
    gs = pl.GridSpec(n_row, n_col, wspace=0.3, hspace=0.5)

    for v, var in enumerate(var_names):
        ix = np.where(adata.var_names == var)[0][0]
        s, u = adata.layers['Ms'][:, ix], adata.layers['Mu'][:, ix]

        # spliced/unspliced phase portrait with steady-state estimate
        ax = pl.subplot(gs[v * n_col])
        scatter(adata, x=s, y=u, color=color, frameon=True, title=var, xlabel='spliced', ylabel='unspliced',
                show=False, ax=ax, fontsize=fontsize, size=size, alpha=alpha, **kwargs)

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
            title = 'expression' if layer == 'Ms' else layer
            scatter(adata, basis=basis, color=var, layer=layer, color_map=color_map, title=title,
                    fontsize=fontsize, size=size, alpha=alpha, show=False, ax=ax, **kwargs)

        if mode == 'stochastic' is not None:
            ss, us = second_order_moments(adata)

            ax = pl.subplot(gs[v*n_col + len(layers) + 1])
            x = 2 * (ss - s**2) - s
            y = 2 * (us - u * s) + u + 2 * s * \
                adata.var['stochastic_velocity_offset'][ix] / adata.var['stochastic_velocity_beta'][ix]

            scatter(adata, x=x, y=y, color=color, title=var, fontsize=40/n_col, show=False, ax=ax,
                    xlabel=r'2 $\Sigma_s - \langle s \rangle$', ylabel=r'2 $\Sigma_{us} + \langle u \rangle$', **kwargs)

            xnew = np.linspace(x.min(), x.max() * 1.02)
            fits = adata.layers.keys() if fits == 'all' else fits
            fits = [fit for fit in fits if 'velocity' in fit]
            for fit in fits:
                linestyle = '--' if 'stochastic' in fit else '-'
                pl.plot(xnew, adata.var[fit + '_gamma'][ix] / adata.var[fit + '_beta'][ix] * xnew +
                        adata.var[fit + '_offset2'][ix] / adata.var[fit + '_beta'][ix], c='k', linestyle=linestyle)
            if v == len(var_names) - 1: pl.legend(fits, loc='lower right', prop={'size': 34/n_col})