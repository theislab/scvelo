from ..tools.velocity_embedding import velocity_embedding as compute_velocity_embedding
from .utils import interpret_colorkey, default_basis, default_size, get_components, savefig
from .scatter import scatter
from .docs import doc_scatter, doc_params

from matplotlib import rcParams
from matplotlib.colors import is_color_like
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd


@doc_params(scatter=doc_scatter)
def velocity_embedding(adata, basis=None, vkey='velocity', density=1, scale=1, X=None, V=None, color=None, use_raw=None,
                       layer=None, color_map=None, colorbar=False, palette=None, size=None, alpha=.2, perc=None,
                       sort_order=True, groups=None, components=None, projection='2d', legend_loc='none', legend_fontsize=None,
                       legend_fontweight=None, right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None,
                       fontsize=None, figsize=None, dpi=None, frameon=None, show=True, save=None, ax=None, **kwargs):
    """\
    Scatter plot of velocities on the embedding

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    x: `str`, `np.ndarray` or `None` (default: `None`)
        x coordinate
    y: `str`, `np.ndarray` or `None` (default: `None`)
        y coordinate
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    density: `float` (default: 1)
        Amount of velocities to show - 0 none to 1 all
    scale: `float` (default: 1)
        Length of velocities in the embedding.
    {scatter}

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    basis = default_basis(adata) if basis is None else basis
    colors = pd.unique(color) if isinstance(color, (list, tuple, np.record)) else [color]
    layers = pd.unique(layer) if isinstance(layer, (list, tuple, np.ndarray, np.record)) else [layer]
    vkeys = pd.unique(vkey) if isinstance(vkey, (list, tuple, np.ndarray, np.record)) else [vkey]
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys() and V is None:
            compute_velocity_embedding(adata, basis=basis, vkey=key)

    if len(colors) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, basis=basis, vkey=vkey, layer=layer, density=density, scale=scale,
                               perc=perc, color=colors[i], use_raw=use_raw, sort_order=sort_order, X=X, V=V,
                               alpha=alpha, groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                               legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                               color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                               left_margin=left_margin, size=size, title=title, show=False, figsize=figsize, dpi=dpi,
                               save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                               fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(layers) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, basis=basis, vkey=vkey, layer=layers[i], density=density, scale=scale,
                               perc=None, color=color, use_raw=use_raw, sort_order=sort_order, X=X, V=V, alpha=alpha,
                               groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                               legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                               color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                               left_margin=left_margin, size=size, title=title, show=False, figsize=figsize, dpi=dpi,
                               save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                               fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(vkeys) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(vkeys), pl.figure(None, (figsize[0] * len(vkeys), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, basis=basis, vkey=vkeys[i], layer=layer, density=density, scale=scale,
                               perc=None, color=color, use_raw=use_raw, sort_order=sort_order, X=X, V=V, alpha=alpha,
                               groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                               legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                               color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                               left_margin=left_margin, size=size, title=title, show=False, figsize=figsize, dpi=dpi,
                               save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel, colorbar=colorbar,
                               fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    else:
        ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        color, layer, vkey = colors[0], layers[0], vkeys[0]
        ix_choice = np.random.choice(adata.n_obs, size=int(density * adata.n_obs), replace=False)

        x, y = None if X is None else X[:, 0], None if X is None else X[:, 1]
        X = adata.obsm['X_' + basis][:, get_components(components, basis)][ix_choice] if X is None else X[:, :2][ix_choice]
        V = adata.obsm[vkey + '_' + basis][:, get_components(components, basis)][ix_choice] if V is None else V[:, :2][ix_choice]

        quiver_kwargs = {"scale": scale, "cmap": color_map, "angles": 'xy', "scale_units": 'xy', "width": .0005,
                         "edgecolors": 'k', "headlength": 12, "headwidth": 10, "headaxislength": 8, "linewidth": .1}
        quiver_kwargs.update(kwargs)

        c = interpret_colorkey(adata, color, layer, perc)
        c = c[ix_choice] if len(c) == adata.n_obs else c

        if is_color_like(c[0]): pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=c, zorder=3, **quiver_kwargs)
        else: pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], c, zorder=3, **quiver_kwargs)

        size = default_size(adata)/2 if size is None else size
        ax = scatter(adata, x=x, y=y, basis=basis, layer=layer, color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                     perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                     colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=ax,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax
