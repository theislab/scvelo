from ..tools.velocity_embedding import velocity_embedding as compute_velocity_embedding
from ..tools.utils import groups_to_bool
from .utils import interpret_colorkey, default_basis, default_size, get_components, savefig_or_show, default_color, \
    default_arrow, make_unique_list, make_unique_valid_list, get_basis, default_color_map, velocity_embedding_changed
from .scatter import scatter
from .docs import doc_scatter, doc_params

from matplotlib import rcParams
from matplotlib.colors import is_color_like
import matplotlib.pyplot as pl
import numpy as np


@doc_params(scatter=doc_scatter)
def velocity_embedding(adata, basis=None, vkey='velocity', density=None, arrow_size=None, arrow_length=None, scale=None,
                       X=None, V=None, recompute=None, color=None, use_raw=None, layer=None, color_map=None, colorbar=True,
                       palette=None, size=None, alpha=.2, perc=None, sort_order=True, groups=None, components=None,
                       projection='2d', legend_loc='none', legend_fontsize=None, legend_fontweight=None,
                       xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None, dpi=None, frameon=None,
                       show=True, save=None, ax=None, ncols=None, **kwargs):
    """\
    Scatter plot of velocities on the embedding.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` or `None` (default: `None`)
        Key for annotations of observations/cells or variables/genes.
    density: `float` (default: 1)
        Amount of velocities to show - 0 none to 1 all
    arrow_size: `float` or 3-tuple for headlength, headwidth and headaxislength (default: 1)
        Size of arrows.
    arrow_length: `float` (default: 1)
        Length of arrows.
    scale: `float` (default: 1)
        Length of velocities in the embedding.
    {scatter}

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    #fkeys = ['adata', 'show', 'save', 'groups', 'figsize', 'dpi', 'ncols', 'wspace', 'hspace', 'ax', 'kwargs']

    vkey = [key for key in adata.layers.keys() if 'velocity' in key and '_u' not in key] if vkey is 'all' else vkey
    layers, vkeys, colors = make_unique_list(layer), make_unique_list(vkey), make_unique_list(color, allow_array=True)
    bases = [default_basis(adata) if basis is None else basis for basis in make_unique_valid_list(adata, basis)]

    if V is None:
        for key in vkeys:
            for basis in bases:
                if recompute or velocity_embedding_changed(adata, basis=basis, vkey=key):
                    compute_velocity_embedding(adata, basis=basis, vkey=key)

    scatter_kwargs = {"perc": perc, "use_raw": use_raw, "sort_order": sort_order, "alpha": alpha,
                      "components": components, "projection": projection, "legend_loc": legend_loc, "groups": groups,
                      "legend_fontsize": legend_fontsize, "legend_fontweight": legend_fontweight, "palette": palette,
                      "color_map": color_map, "frameon": frameon, "xlabel": xlabel, "ylabel": ylabel,
                      "colorbar": colorbar, "dpi": dpi, "fontsize": fontsize, "show": False, "save": False}

    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 \
        else vkeys if len(vkeys) > 1 else bases if len(bases) > 1 else None
    if multikey is not None:
        if title is None: title = list(multikey)
        elif isinstance(title, (list, tuple)): title *= int(np.ceil(len(multikey) / len(title)))
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        ax = []
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))):
            if i < len(multikey):
                ax.append(velocity_embedding(adata, density=density, scale=scale, size=size, ax=pl.subplot(gs),
                                             arrow_size=arrow_size, arrow_length=arrow_length,
                                             basis=bases[i] if len(bases) > 1 else basis,
                                             color=colors[i] if len(colors) > 1 else color,
                                             layer=layers[i] if len(layers) > 1 else layer,
                                             vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                                             title=title[i] if isinstance(title, (list, tuple)) else title,
                                             **scatter_kwargs, **kwargs))
        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax

    else:
        if projection == '3d':
            from mpl_toolkits.mplot3d import Axes3D
            ax = pl.figure(None, figsize, dpi=dpi).gca(projection=projection) if ax is None else ax
        else:
            ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        color, layer, vkey, basis = colors[0], layers[0], vkeys[0], bases[0]
        color = default_color(adata) if color is None else color
        color_map = default_color_map(adata, color) if color_map is None else color_map
        size = default_size(adata) / 2 if size is None else size
        if use_raw is None and 'Ms' not in adata.layers.keys(): use_raw = True
        _adata = adata[groups_to_bool(adata, groups, groupby=color)] if groups is not None and color in adata.obs.keys() else adata

        quiver_kwargs = {"scale": scale, "cmap": color_map, "angles": 'xy',
                         "scale_units": 'xy', "edgecolors": 'k', "linewidth": .1, "width": None}
        if basis in adata.var_names:
            x = adata[:, basis].layers['spliced'] if use_raw else adata[:, basis].layers['Ms']
            y = adata[:, basis].layers['unspliced'] if use_raw else adata[:, basis].layers['Mu']
            dx = adata[:, basis].layers[vkey]
            dy = adata[:, basis].layers[vkey + '_u'] if vkey + '_u' in adata.layers.keys() else np.zeros(adata.n_obs)
            X = np.stack([np.ravel(x), np.ravel(y)]).T
            V = np.stack([np.ravel(dx), np.ravel(dy)]).T
        else:
            x = None if X is None else X[:, 0]
            y = None if X is None else X[:, 1]
            X = _adata.obsm['X_' + basis][:, get_components(components, basis, projection)] if X is None else X[:, :2]
            V = _adata.obsm[vkey + '_' + basis][:, get_components(components, basis, projection)] if V is None else V[:, :2]

            hl, hw, hal = default_arrow(arrow_size)
            scale = 1 / arrow_length if arrow_length is not None else scale if scale is not None else 1
            quiver_kwargs.update({"scale": scale, "width": .0005, "headlength": hl, "headwidth": hw, "headaxislength": hal})

        for arg in list(kwargs):
            if arg in quiver_kwargs: quiver_kwargs.update({arg: kwargs[arg]})
            else: scatter_kwargs.update({arg: kwargs[arg]})

        if basis in adata.var_names and isinstance(color, str) and color in adata.layers.keys():
            c = interpret_colorkey(_adata, basis, color, perc)
        else:
            c = interpret_colorkey(_adata, color, layer, perc)

        if density is not None and 0 < density < 1:
            ix_choice = np.random.choice(_adata.n_obs, size=int(density * _adata.n_obs), replace=False)
            c = c[ix_choice] if len(c) == _adata.n_obs else c
            X = X[ix_choice]
            V = V[ix_choice]

        if projection == '3d' and X.shape[1] > 2 and V.shape[1] > 2:
            V, size = V / scale / 5, size / 10
            x0, x1, x2, v0, v1, v2 = X[:, 0], X[:, 1], X[:, 2], V[:, 0], V[:, 1], V[:, 2]
            quiver3d_kwargs = {"zorder": 3, "linewidth": .5, "arrow_length_ratio": .3}
            c = list(c) + [element for element in list(c) for _ in range(2)]
            if is_color_like(c[0]): ax.quiver(x0, x1, x2, v0, v1, v2, color=c, **quiver3d_kwargs)
            else: ax.quiver(x0, x1, x2, v0, v1, v2, c, **quiver3d_kwargs)
        else:
            if is_color_like(c[0]): ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], color=c, zorder=3, **quiver_kwargs)
            else: ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], c, zorder=3, **quiver_kwargs)

        ax = scatter(adata, basis=basis, x=x, y=y, vkey=vkey, layer=layer, color=color, size=size, title=title, ax=ax,
                     zorder=0, **scatter_kwargs)

        savefig_or_show(dpi=dpi, save=save, show=show)
        if not show: return ax
