from ..tools.velocity_embedding import quiver_autoscale, velocity_embedding
from .utils import default_basis, default_size, get_components, savefig
from .scatter import scatter
from .docs import doc_scatter, doc_params

from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal
from scipy.interpolate import griddata
from matplotlib import rcParams
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd


def compute_velocity_on_grid(X_emb, V_emb, density=1, smooth=0.5, n_neighbors=None, min_mass=None, autoscale=True):
    # prepare grid
    n_obs, n_dim = X_emb.shape

    grs = []
    for dim_i in range(n_dim):
        m, M = np.min(X_emb[:, dim_i]), np.max(X_emb[:, dim_i])
        m = m - .01 * np.abs(M - m)
        M = M + .01 * np.abs(M - m)
        gr = np.linspace(m, M, 50 * density)
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    X_grid = np.vstack([i.flat for i in meshes_tuple]).T

    # estimate grid velocities
    if n_neighbors is None: n_neighbors = int(n_obs/50)
    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=-1)
    nn.fit(X_emb)
    dists, neighs = nn.kneighbors(X_grid)

    scale = np.mean([(g[1] - g[0]) for g in grs]) * smooth
    weight = normal.pdf(x=dists, scale=scale)
    p_mass = weight.sum(1)

    V_grid = (V_emb[neighs] * weight[:, :, None]).sum(1) / np.maximum(1, p_mass)[:, None]
    if min_mass is None: min_mass = np.clip(np.percentile(p_mass, 99) / 100, 1e-2, 1)
    X_grid, V_grid = X_grid[p_mass > min_mass], V_grid[p_mass > min_mass]

    if autoscale: V_grid /= 3 * quiver_autoscale(X_grid, V_grid)

    return X_grid, V_grid


@doc_params(scatter=doc_scatter)
def velocity_embedding_grid(adata, basis=None, vkey='velocity', density=1, scale=1, smooth=.5, min_mass=None, autoscale=True,
                            n_neighbors=None, X=None, V=None, X_grid=None, V_grid=None, principal_curve=False, color=None, use_raw=None, layer=None,
                            color_map=None, colorbar=False, palette=None, size=None, alpha=.2, perc=None, sort_order=True,
                            groups=None, components=None, projection='2d', legend_loc='none', legend_fontsize=None,
                            legend_fontweight=None, right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None,
                            fontsize=None, figsize=None, dpi=None, frameon=None, show=True, save=None, ax=None, **kwargs):
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
    basis = default_basis(adata) if basis is None else basis
    colors = pd.unique(color) if isinstance(color, (list, tuple, np.record)) else [color]
    layers = pd.unique(layer) if isinstance(layer, (list, tuple, np.ndarray, np.record)) else [layer]
    vkeys = pd.unique(vkey) if isinstance(vkey, (list, tuple, np.ndarray, np.record)) else [vkey]
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys() and V is None:
            velocity_embedding(adata, basis=basis, vkey=key)
    color, layer, vkey = colors[0], layers[0], vkeys[0]

    if X_grid is None or V_grid is None:
        X_emb  = adata.obsm['X_' + basis][:, get_components(components, basis)] if X is None else X[:, :2]
        V_emb = adata.obsm[vkey + '_' + basis][:, get_components(components, basis)] if V is None else V[:, :2]
        X_grid, V_grid = compute_velocity_on_grid(X_emb=X_emb, V_emb=V_emb, density=density, autoscale=autoscale,
                                                  smooth=smooth, n_neighbors=n_neighbors, min_mass=min_mass)

    if len(colors) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkey, layer=layer, density=density, scale=scale, X_grid=X_grid, V_grid=V_grid,
                                    perc=perc, color=colors[i], min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, autoscale=autoscale,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(layers) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkey, layer=layers[i], density=density, scale=scale, X_grid=X_grid, V_grid=V_grid,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, autoscale=autoscale,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    elif len(vkeys) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(vkeys), pl.figure(None, (figsize[0] * len(vkeys), figsize[1]), dpi=dpi))):
            velocity_embedding_grid(adata, basis=basis, vkey=vkeys[i], layer=layer, density=density, scale=scale,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    principal_curve=principal_curve, use_raw=use_raw, sort_order=sort_order, alpha=alpha,
                                    groups=groups, components=components, projection=projection, legend_loc=legend_loc,
                                    legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, autoscale=autoscale,
                                    color_map=color_map, palette=palette, frameon=frameon, right_margin=right_margin,
                                    left_margin=left_margin, size=size, title=title, show=False, figsize=figsize,
                                    dpi=dpi, save=None, ax=pl.subplot(gs), xlabel=xlabel, ylabel=ylabel,
                                    colorbar=colorbar, fontsize=fontsize, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax

    else:
        ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        quiver_kwargs = {"scale": scale, "angles": 'xy', "scale_units": 'xy', "width": .001, "color": 'black',
                   "edgecolors": 'k', "headwidth": 4.5, "headlength": 5, "headaxislength": 3, "linewidth": .2}
        quiver_kwargs.update(kwargs)
        pl.quiver(X_grid[:, 0], X_grid[:, 1], V_grid[:, 0], V_grid[:, 1], **quiver_kwargs, zorder=3)

        if principal_curve:
            curve = adata.uns['principal_curve']['projections']
            pl.plot(curve[:, 0], curve[:, 1], c="w", lw=6, zorder=4)
            pl.plot(curve[:, 0], curve[:, 1], c="k", lw=3, zorder=5)

        size = default_size(adata) if size is None else size
        ax = scatter(adata, basis=basis, layer=layer, color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                     perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                     colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=ax,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax


def compute_grid_for_stream(X_grid, V_grid):
    x, y, u, v = X_grid[:, 0], X_grid[:, 1], V_grid[:, 0], V_grid[:, 1]

    pts = np.vstack((y, x)).T
    vals = np.vstack((v, u)).T

    n_gridpoints = 50
    x = np.linspace(x.min(), x.max(), n_gridpoints)
    y = np.linspace(y.min(), y.max(), n_gridpoints)

    # an (nx * ny, 2) array of x,y coordinates to interpolate at
    ipts = np.vstack(a.ravel() for a in np.meshgrid(x, y)[::-1]).T

    # an (nx * ny, 2) array of interpolated u, v values
    v, u = griddata(pts, vals, ipts, method='cubic').T
    u.shape = v.shape = (n_gridpoints, n_gridpoints)

    return x, y, u, v


@doc_params(scatter=doc_scatter)
def velocity_embedding_stream(adata, basis=None, vkey='velocity', density=1, scale=1, smooth=.5, min_mass=None,
                              n_neighbors=None, X=None, V=None, X_grid=None, V_grid=None, linewidth=1, color=None, use_raw=None, layer=None,
                              color_map=None, colorbar=False, palette=None, size=None, alpha=.2, perc=None, sort_order=True,
                              groups=None, components=None, projection='2d', legend_loc='none', legend_fontsize=None,
                              legend_fontweight=None, right_margin=None, left_margin=None, xlabel=None, ylabel=None, title=None,
                              fontsize=None, figsize=None, dpi=None, frameon=None, show=True, save=None, ax=None, **kwargs):
    """\
    Stream plot of velocities computed from the grid points on the embedding

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
    basis = default_basis(adata) if basis is None else basis
    colors = pd.unique(color) if isinstance(color, (list, tuple, np.record)) else [color]
    layers = pd.unique(layer) if isinstance(layer, (list, tuple, np.ndarray, np.record)) else [layer]
    vkeys = pd.unique(vkey) if isinstance(vkey, (list, tuple, np.ndarray, np.record)) else [vkey]
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys() and V is None:
            velocity_embedding(adata, basis=basis, vkey=key)
    color, layer, vkey = colors[0], layers[0], vkeys[0]

    if X_grid is None or V_grid is None:
        X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)] if X is None else X[:, :2]
        V_emb = adata.obsm[vkey + '_' + basis][:, get_components(components, basis)] if V is None else V[:, :2]
        X_grid, V_grid = compute_velocity_on_grid(X_emb=X_emb, V_emb=V_emb, density=1, smooth=smooth,
                                                  n_neighbors=n_neighbors, min_mass=-np.inf, autoscale=False)
    if len(colors) > 1:
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(colors), pl.figure(None, (figsize[0] * len(colors), figsize[1]), dpi=dpi))):
            velocity_embedding_stream(adata, basis=basis, vkey=vkey, layer=layer, density=density, scale=scale, X_grid=X_grid, V_grid=V_grid,
                                    perc=perc, color=colors[i], min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    use_raw=use_raw, sort_order=sort_order, alpha=alpha, linewidth=linewidth,
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
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(layers), pl.figure(None, (figsize[0] * len(layers), figsize[1]), dpi=dpi))):
            velocity_embedding_stream(adata, basis=basis, vkey=vkey, layer=layers[i], density=density, scale=scale, X_grid=X_grid, V_grid=V_grid,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    use_raw=use_raw, sort_order=sort_order, alpha=alpha, linewidth=linewidth,
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
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(pl.GridSpec(1, len(vkeys), pl.figure(None, (figsize[0] * len(vkeys), figsize[1]), dpi=dpi))):
            velocity_embedding_stream(adata, basis=basis, vkey=vkeys[i], layer=layer, density=density, scale=scale, X_grid=X_grid, V_grid=V_grid,
                                    perc=perc, color=color, min_mass=min_mass, smooth=smooth, n_neighbors=n_neighbors,
                                    use_raw=use_raw, sort_order=sort_order, alpha=alpha, linewidth=linewidth,
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
        ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        x, y, u, v = compute_grid_for_stream(X_grid, V_grid)
        lw = np.sqrt(u ** 2 + v ** 2)
        u[lw.reshape(u.shape) < 1e-5] = np.nan

        stream_kwargs = {"linewidth": 8 * lw * linewidth, "density": 2 * density}
        stream_kwargs.update(kwargs)
        pl.streamplot(x, y, u, v, color='grey', zorder=3, **stream_kwargs)

        size = default_size(adata) if size is None else size
        ax = scatter(adata, basis=basis, layer=layer, color=color, xlabel=xlabel, ylabel=ylabel, color_map=color_map,
                     perc=perc, size=size, alpha=alpha, fontsize=fontsize, frameon=frameon, title=title, show=False,
                     colorbar=colorbar, components=components, figsize=figsize, dpi=dpi, save=None, ax=ax,
                     use_raw=use_raw, sort_order=sort_order, groups=groups, projection=projection,
                     legend_loc=legend_loc, legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
                     palette=palette, right_margin=right_margin, left_margin=left_margin, **kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax
