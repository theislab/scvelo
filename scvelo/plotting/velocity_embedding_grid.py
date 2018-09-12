from ..tools.velocity_embedding import quiver_autoscale, velocity_embedding
from .utils import get_components, savefig
from .scatter import scatter
from .docs import doc_scatter, doc_params

from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal
import matplotlib.pyplot as pl
import numpy as np


def compute_velocity_on_grid(X_emb, V_emb, density=1, smooth=0.5, n_neighbors=None, min_mass=.5):
    # prepare grid
    n_obs, n_dim = X_emb.shape

    grs = []
    for dim_i in range(n_dim):
        m, M = np.min(X_emb[:, dim_i]), np.max(X_emb[:, dim_i])
        m = m - .01 * np.abs(M - m)
        M = M + .01 * np.abs(M - m)
        gr = np.linspace(m, M, 70 * density)
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    X_grid = np.vstack([i.flat for i in meshes_tuple]).T

    # estimate grid velocities
    if n_neighbors is None: n_neighbors = int(n_obs/50)
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=-1)
    nn.fit(X_emb)
    dists, neighs = nn.kneighbors(X_grid)

    std = np.mean([(g[1] - g[0]) for g in grs])
    weight = normal.pdf(loc=0, scale=smooth * std, x=dists)
    p_mass = weight.sum(1)

    V_grid = (V_emb[neighs] * weight[:, :, None]).sum(1) / np.maximum(1, p_mass)[:, None]
    X_grid, V_grid = X_grid[p_mass > min_mass], V_grid[p_mass > min_mass]

    V_grid /= 3.5 * quiver_autoscale(X_grid[:, 0], X_grid[:, 1], V_grid[:, 0], V_grid[:, 1])

    return X_grid, V_grid


@doc_params(scatter=doc_scatter)
def velocity_embedding_grid(adata, basis='umap', vkey='velocity', density=1, scale=1, min_mass=.5, smooth=.5,
                            n_neighbors=None, X=None, V=None, principal_curve=False, color=None, use_raw=None, layer=None,
                            color_map=None, colorbar=False, palette=None, size=5, alpha=.2, perc=None, sort_order=True,
                            groups=None, components=None, projection='2d', legend_loc='none', legend_fontsize=None,
                            legend_fontweight=None, right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None,
                            fontsize=None, figsize=(7,5), dpi=100, frameon=False, show=True, save=None, ax=None, **kwargs):
    """\
    Scatter plot of velocities for the grid points on the embedding

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
    min_mass: `float` (default: 0.5)
        Minimum threshold for mass to be shown.
    smooth: `float` (default: 0.5)
        Multiplication factor for scale in Gaussian kernel around grid point.
    n_neighbors: `int` (default: None)
        Number of neighbors to consider around grid point.
    X: `np.ndarray` (default: None)
        embedding grid point coordinates
    V: `np.ndarray` (default: None)
        embedding grid velocity coordinates
    {scatter}

    Returns
    -------
        `matplotlib.Axis` if `show==False`
    """
    colors = color if isinstance(color, (list, tuple)) else [color]
    layers = layer if isinstance(layer, (list, tuple)) else [layer]
    vkeys = vkey if isinstance(vkey, (list, tuple)) else [vkey]
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys(): velocity_embedding(adata, basis=basis, vkey=key)

    if X is None and V is None:
        X, V = compute_velocity_on_grid(X_emb=adata.obsm['X_' + basis][:, get_components(components)],
                                        V_emb=adata.obsm[vkeys[0] + '_' + basis], density=density, smooth=smooth,
                                        n_neighbors=n_neighbors, min_mass=min_mass)

    if len(colors) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkey, layer=layer, density=density, scale=scale, X=X, V=V,
                                    perc=perc, color=colors[i], min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(layers) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkey, layer=layers[i], density=density, scale=scale, X=X, V=V,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(vkeys) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(vkeys), pl.figure(None, (figsize[0] * len(vkeys), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkeys[i], layer=layer, density=density, scale=scale, X=X, V=V,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    else:
        _kwargs = {"scale": scale, "angles": 'xy', "scale_units": 'xy', "width": .001, "color": 'black',
                   "edgecolors": 'k', "headwidth": 4.5, "headlength": 5, "headaxislength": 3, "linewidth": .2}
        _kwargs.update(kwargs)

        if color_map is None: color_map = 'viridis_r' if (color == 'root' or color == 'end') else 'RdBu_r'
        if ax is None: ax = pl.figure(None, figsize, dpi=dpi).gca()

        pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], **_kwargs, zorder=1)

        if principal_curve:
            curve = adata.uns['principal_curve']['projections']
            pl.plot(curve[:, 0], curve[:, 1], c="w", lw=6, zorder=2)
            pl.plot(curve[:, 0], curve[:, 1], c="k", lw=3, zorder=3)

        ax = scatter(adata, basis=basis, layer=layer, color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                     perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                     colorbar=colorbar, components=components, figsize=(7, 5), dpi=80, save=None, ax=ax,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)

        if show: pl.show()
        else: return ax
