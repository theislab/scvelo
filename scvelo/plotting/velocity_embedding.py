from ..tools.velocity_embedding import velocity_embedding as compute_velocity_embedding
from ..tools.utils import groups_to_bool
from .utils import interpret_colorkey, default_basis, default_size, get_components, savefig, default_color
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

    scatter_kwargs = {"basis": basis, "perc": perc, "use_raw": use_raw, "sort_order": sort_order, "alpha": alpha,
                      "components": components, "projection": projection, "legend_loc": legend_loc, "groups": groups,
                      "legend_fontsize": legend_fontsize, "legend_fontweight": legend_fontweight, "palette": palette,
                      "color_map": color_map, "frameon": frameon, "title": title, "xlabel": xlabel, "ylabel": ylabel,
                      "right_margin": right_margin, "left_margin": left_margin, "colorbar": colorbar, "dpi": dpi,
                      "fontsize": fontsize, "show": False, "save": None}

    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 else vkeys if len(vkeys) > 1 else None
    if multikey is not None:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(multikey), pl.figure(None, (figsize[0]*len(multikey), figsize[1]), dpi=dpi))):
            velocity_embedding(adata, X=X, V=V, density=density, scale=scale, size=size, ax=pl.subplot(gs),
                               color=colors[i] if len(colors) > 1 else color,
                               layer=layers[i] if len(layers) > 1 else layer,
                               vkey=vkeys[i] if len(vkeys) > 1 else vkey, **scatter_kwargs, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show:
            pl.show()
        else:
            return ax

    else:
        if projection == '3d':
            from mpl_toolkits.mplot3d import Axes3D
            ax = pl.figure(None, figsize, dpi=dpi).gca(projection=projection) if ax is None else ax
        else:
            ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        color, layer, vkey = colors[0], layers[0], vkeys[0]
        color = default_color(adata) if color is None else color
        _adata = adata[groups_to_bool(adata, groups, groupby=color)] if groups is not None and color in adata.obs.keys() else adata

        ix_choice = np.random.choice(_adata.n_obs, size=int(density * _adata.n_obs), replace=False)

        x, y = None if X is None else X[:, 0], None if X is None else X[:, 1]
        X = _adata.obsm['X_' + basis][:, get_components(components, basis)][ix_choice] if X is None else X[:, :2][ix_choice]
        V = _adata.obsm[vkey + '_' + basis][:, get_components(components, basis)][ix_choice] if V is None else V[:, :2][ix_choice]

        quiver_kwargs = {"scale": scale, "cmap": color_map, "angles": 'xy', "scale_units": 'xy', "width": .0005,
                         "edgecolors": 'k', "headlength": 12, "headwidth": 10, "headaxislength": 8, "linewidth": .1}
        quiver_kwargs.update(kwargs)

        c = interpret_colorkey(_adata, color, layer, perc)
        c = c[ix_choice] if len(c) == _adata.n_obs else c

        if is_color_like(c[0]): pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=c, zorder=3, **quiver_kwargs)
        else: pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], c, zorder=3, **quiver_kwargs)

        size = default_size(adata) / 2 if size is None else size,
        ax = scatter(adata, x=x, y=y, layer=layer, color=color, size=size, ax=ax, **scatter_kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax
