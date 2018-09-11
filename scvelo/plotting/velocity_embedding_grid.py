from ..tools.velocity_embedding import velocity_embedding as tl_velocity_embedding
from .utils import quiver_autoscale, get_components, savefig
from .scatter import scatter
from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal
import matplotlib.pyplot as pl
import numpy as np


def compute_velocity_on_grid(X_emb, V_emb, density=1, smooth=0.5, n_neighbors=None, min_mass=.5):
    """Computes the velocities for the grid points on the embedding

    Arguments
    ---------
    X_emb: np.ndarray
        embedding coordinates

    V_emb: np.ndarray
        embedded single cell velocity (obtained with 'velocity_embedding')

    density: float, default=1

    smooth: float, default=.5

    n_neighbors: int

    min_mass: float, default=.5

    Returns
    -------
    grid_coord: np.array
        grid point coordinates
    grid_velocity: np.array
        velocity for each grid point in the embedding
    """
    # prepare grid
    n_obs, n_dim = X_emb.shape
    steps = (int(np.sqrt(n_obs) * density), int(np.sqrt(n_obs) * density))

    grs = []
    for dim_i in range(n_dim):
        m, M = np.min(X_emb[:, dim_i]), np.max(X_emb[:, dim_i])
        m = m - .01 * np.abs(M - m)
        M = M + .01 * np.abs(M - m)
        gr = np.linspace(m, M, steps[dim_i])
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    grid_coord = np.vstack([i.flat for i in meshes_tuple]).T

    # estimate grid velocities
    if n_neighbors is None: n_neighbors = int(n_obs/50)
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=-1)
    nn.fit(X_emb)
    dists, neighs = nn.kneighbors(grid_coord)

    std = np.mean([(g[1] - g[0]) for g in grs])
    weight = normal.pdf(loc=0, scale=smooth * std, x=dists)
    p_mass = weight.sum(1)

    grid_velocity = (V_emb[neighs] * weight[:, :, None]).sum(1) / np.maximum(1, p_mass)[:, None]
    grid_coord, grid_velocity = grid_coord[p_mass > min_mass], grid_velocity[p_mass > min_mass]

    return grid_coord, grid_velocity


def velocity_embedding_grid(adata, basis='umap', vbasis='velocity', layer=None, density=1, scale=1, autoscale=True,
                            color=None, perc=None, min_mass=.5, smooth=.5, n_neighbors=None, principal_curve=False,
                            use_raw=True, sort_order=True, alpha=.2, groups=None, components=None, projection='2d',
                            legend_loc='none', legend_fontsize=None, legend_fontweight=None,
                            color_map=None, palette=None, frameon=False, right_margin=None, left_margin=None,
                            size=None, title=None, show=True, figsize=(14,10), dpi=150,
                            xlabel=None, ylabel=None, colorbar=False, fontsize=None, save=None, ax=None, **kwargs):
    """Scatter plot with grid velocities along `.obs` or `.var` axes.
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
    if vkey not in adata.obsm_keys():
        tl_velocity_embedding(adata, basis=basis, vkey=vkey)

    X, V = compute_velocity_on_grid(X_emb=adata.obsm['X_' + basis][:, get_components(components)],
                                    V_emb=adata.obsm[vkey], density=density, smooth=smooth,
                                    n_neighbors=n_neighbors, min_mass=min_mass)
    if autoscale: scale *= 3.5 * quiver_autoscale(X[:, 0], X[:, 1], V[:, 0], V[:, 1])

    colors = color if isinstance(color, (list, tuple)) else [color]
    layers = layer if isinstance(layer, (list, tuple)) else [layer]

    if len(colors) > 1:
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vbasis=vbasis, layer=layer, density=density, scale=scale, autoscale=False,
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
            velocity_embedding_grid(adata, basis=basis, vbasis=vbasis, layer=layers[i], density=density, scale=scale, autoscale=False,
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
        X, V = compute_velocity_on_grid(X_emb=adata.obsm['X_' + basis][:, get_components(components)],
                                        V_emb=adata.obsm[vkey], density=density, smooth=smooth,
                                        n_neighbors=n_neighbors, min_mass=min_mass)

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
                     colorbar=colorbar, components=components, figsize=(7, 5), dpi=80, save=None, ax=ax, zorder=0,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)

        if show: pl.show()
        else: return ax
