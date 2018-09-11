from ..tools.velocity_embedding import velocity_embedding as tl_velocity_embedding
from .utils import interpret_colorkey, quiver_autoscale, get_components, savefig
from .scatter import scatter
from matplotlib.colors import is_color_like
import matplotlib.pyplot as pl
import numpy as np


def velocity_embedding(adata, basis='umap', vbasis='velocity', layer=None, density=1, scale=1, autoscale=True,
                       perc=None, color=None, use_raw=True, sort_order=True, alpha=.2, groups=None, components=None,
                       projection='2d', legend_loc='none', legend_fontsize=None, legend_fontweight=None,
                       color_map=None, palette=None, frameon=False, right_margin=None, left_margin=None,
                       size=1, title=None, show=True, figsize=(14,10), dpi=150, save=None, ax=None,
                       xlabel=None, ylabel=None, colorbar=False, fontsize=None, **kwargs):
    """Scatter plot with velocities along `.obs` or `.var` axes.
    Color the plot using annotations of observations (`.obs`), variables (`.var`) or expression of genes (`.var_names`).

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    basis: `str` (default: `'umap'`)
        Key for embedding coordinates.
    vbasis: `str` (default: `'velocity'`)
        Key for velocity embedding coordinates.
    color : `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    vkey = vbasis + '_' + basis
    if vkey not in adata.obsm_keys(): tl_velocity_embedding(adata, basis=basis, vkey=vkey)

    ix_choice = np.random.choice(adata.n_obs, size=int(density * adata.n_obs), replace=False)
    X = adata.obsm['X_' + basis][:, get_components(components)][ix_choice]
    V = adata.obsm[vkey][ix_choice]
    if autoscale: scale *= 2 * quiver_autoscale(X[:, 0], X[:, 1], V[:, 0], V[:, 1])

    colors = color if isinstance(color, (list, tuple)) else [color]
    layers = layer if isinstance(layer, (list, tuple)) else [layer]

    if len(colors) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, basis=basis, vbasis=vbasis, layer=layer, density=density, scale=scale,
                               autoscale=False, perc=perc, color=colors[i], use_raw=use_raw, sort_order=sort_order,
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
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, basis=basis, vbasis=vbasis, layer=layers[i], density=density, scale=scale, autoscale=False,
                               perc=None, color=color, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
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
        ix_choice = np.random.choice(adata.n_obs, size=int(density * adata.n_obs), replace=False)
        X = adata.obsm['X_' + basis][:, get_components(components)][ix_choice]
        V = adata.obsm[vkey][ix_choice]

        if color_map is None: color_map = 'viridis_r' if (color == 'root' or color == 'end') else 'RdBu_r'
        _kwargs = {"scale": scale, "cmap": color_map, "angles": 'xy', "scale_units": 'xy', "width": .0005,
                   "edgecolors": 'k', "headwidth": 9, "headlength": 10, "headaxislength": 6, "linewidth": .25}
        _kwargs.update(kwargs)

        if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca()

        C = interpret_colorkey(adata, color, layer, perc)[ix_choice]
        if is_color_like(C[0]): pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=C, zorder=1, **_kwargs)
        else: pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], C, zorder=1, **_kwargs)

        ax = scatter(adata, basis=basis, layer=layer, color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                     perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                     colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=ax, zorder=0,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin, ** kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)

        if show: pl.show()
        else: return ax
