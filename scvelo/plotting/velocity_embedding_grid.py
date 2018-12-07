from ..tools.velocity_embedding import quiver_autoscale, velocity_embedding
from .utils import default_basis, default_size, get_components, savefig, default_arrow, make_unique_list
from .scatter import scatter
from .docs import doc_scatter, doc_params

from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal
from matplotlib import rcParams
import matplotlib.pyplot as pl
import numpy as np


def compute_velocity_on_grid(X_emb, V_emb, density=None, smooth=None, n_neighbors=None, min_mass=None, autoscale=True, adjust_for_stream=False):
    # prepare grid
    n_obs, n_dim = X_emb.shape
    density = 1 if density is None else density
    smooth = .5 if smooth is None else smooth

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

    if adjust_for_stream:
        X_grid = np.stack([np.unique(X_grid[:, 0]), np.unique(X_grid[:, 1])])
        ns = int(np.sqrt(len(V_grid[:, 0])))
        V_grid = V_grid.T.reshape(2, ns, ns)

        mass = np.sqrt((V_grid ** 2).sum(0))
        V_grid[0][mass.reshape(V_grid[0].shape) < 1e-5] = np.nan
    else:
        if min_mass is None: min_mass = np.clip(np.percentile(p_mass, 99) / 100, 1e-2, 1)
        X_grid, V_grid = X_grid[p_mass > min_mass], V_grid[p_mass > min_mass]

        if autoscale: V_grid /= 3 * quiver_autoscale(X_grid, V_grid)

    return X_grid, V_grid


@doc_params(scatter=doc_scatter)
def velocity_embedding_grid(adata, basis=None, vkey='velocity', density=None, smooth=None, min_mass=None, arrow_size=None,
                            arrow_length=None, arrow_color=None, scale=None, autoscale=True, n_neighbors=None,
                            X=None, V=None, X_grid=None, V_grid=None, principal_curve=False, color=None, use_raw=None,
                            layer=None, color_map=None, colorbar=False, palette=None, size=None, alpha=.2, perc=None,
                            sort_order=True, groups=None, components=None, projection='2d', legend_loc='none',
                            legend_fontsize=None, legend_fontweight=None, right_margin=None, left_margin=None,
                            xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None, dpi=None, frameon=None,
                            show=True, save=None, ax=None, ncols=None, **kwargs):
    """\
    Scatter plot of velocities on a grid.

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
    arrow_size: `float` or 3-tuple for headlength, headwidth and headaxislength (default: 1)
        Size of arrows.
    arrow_length: `float` (default: 1)
        Length of arrows.
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
    colors, layers, vkeys = make_unique_list(color, allow_array=True), make_unique_list(layer), make_unique_list(vkey)
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys() and V is None:
            velocity_embedding(adata, basis=basis, vkey=key)
    color, layer, vkey = colors[0], layers[0], vkeys[0]

    if X_grid is None or V_grid is None:
        X_emb  = adata.obsm['X_' + basis][:, get_components(components, basis)] if X is None else X[:, :2]
        V_emb = adata.obsm[vkey + '_' + basis][:, get_components(components, basis)] if V is None else V[:, :2]
        X_grid, V_grid = compute_velocity_on_grid(X_emb=X_emb, V_emb=V_emb, density=density, autoscale=autoscale,
                                                  smooth=smooth, n_neighbors=n_neighbors, min_mass=min_mass)

    scatter_kwargs = {"basis": basis, "perc": perc, "use_raw": use_raw, "sort_order": sort_order, "alpha": alpha,
                      "components": components, "projection": projection, "legend_loc": legend_loc, "groups": groups,
                      "legend_fontsize": legend_fontsize, "legend_fontweight": legend_fontweight, "palette": palette,
                      "color_map": color_map, "frameon": frameon, "title": title, "xlabel": xlabel, "ylabel": ylabel,
                      "right_margin": right_margin, "left_margin": left_margin, "colorbar": colorbar, "dpi": dpi,
                      "fontsize": fontsize, "show": False, "save": None}

    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 else vkeys if len(vkeys) > 1 else None
    if multikey is not None:
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))):
            if i < len(multikey):
                velocity_embedding_grid(adata, density=density, scale=scale, size=size, min_mass=min_mass, smooth=smooth,
                                        n_neighbors=n_neighbors, principal_curve=principal_curve, ax=pl.subplot(gs),
                                        arrow_size=arrow_size, arrow_length=arrow_length,
                                        color=colors[i] if len(colors) > 1 else color,
                                        layer=layers[i] if len(layers) > 1 else layer,
                                        vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                                        X_grid=None if len(vkeys) > 1 else X_grid,
                                        V_grid=None if len(vkeys) > 1 else V_grid,
                                        autoscale=False if len(vkeys) > 1 else autoscale, **scatter_kwargs, **kwargs)
        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show:
            pl.show()
        else:
            return ax

    else:
        ax = pl.figure(None, figsize, dpi=dpi).gca() if ax is None else ax

        hl, hw, hal = default_arrow(arrow_size)
        scale = 1 / arrow_length if arrow_length is not None else scale if scale is not None else 1
        quiver_kwargs = {"scale": scale, "angles": 'xy', "scale_units": 'xy', "width": .001,
                         "color": 'grey' if arrow_color is None else arrow_color, "edgecolors": 'k',
                         "headlength": hl/2, "headwidth": hw/2, "headaxislength": hal/2, "linewidth": .2}
        quiver_kwargs.update(kwargs)
        pl.quiver(X_grid[:, 0], X_grid[:, 1], V_grid[:, 0], V_grid[:, 1], **quiver_kwargs, zorder=3)

        if principal_curve:
            curve = adata.uns['principal_curve']['projections']
            pl.plot(curve[:, 0], curve[:, 1], c="w", lw=6, zorder=4)
            pl.plot(curve[:, 0], curve[:, 1], c="k", lw=3, zorder=5)

        size = 4 * default_size(adata) if size is None else size
        ax = scatter(adata, layer=layer, color=color, size=size, ax=ax, zorder=0, **scatter_kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show:
            pl.show()
        else:
            return ax


@doc_params(scatter=doc_scatter)
def velocity_embedding_stream(adata, basis=None, vkey='velocity', density=None, smooth=None, linewidth=None,
                              n_neighbors=None, X=None, V=None, X_grid=None, V_grid=None, color=None, use_raw=None,
                              layer=None, color_map=None, colorbar=False, palette=None, size=None, alpha=.1, perc=None,
                              sort_order=True, groups=None, components=None, projection='2d', legend_loc='none',
                              legend_fontsize=None, legend_fontweight=None, right_margin=None, left_margin=None,
                              xlabel=None, ylabel=None, title=None, fontsize=None, figsize=None, dpi=None, frameon=None,
                              show=True, save=None, ax=None, ncols=None, **kwargs):
    """\
    Stream plot of velocities on the embedding.

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
    smooth: `float` (default: 0.5)
        Multiplication factor for scale in Gaussian kernel around grid point.
    linewidth: `float` (default: 1)
        Line width for streamplot.
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
    colors, layers, vkeys = make_unique_list(color, allow_array=True), make_unique_list(layer), make_unique_list(vkey)
    for key in vkeys:
        if key + '_' + basis not in adata.obsm_keys() and V is None:
            velocity_embedding(adata, basis=basis, vkey=key)
    color, layer, vkey = colors[0], layers[0], vkeys[0]

    if X_grid is None or V_grid is None:
        X_emb = adata.obsm['X_' + basis][:, get_components(components, basis)] if X is None else X[:, :2]
        V_emb = adata.obsm[vkey + '_' + basis][:, get_components(components, basis)] if V is None else V[:, :2]
        X_grid, V_grid = compute_velocity_on_grid(X_emb=X_emb, V_emb=V_emb, density=1, smooth=smooth,
                                                  n_neighbors=n_neighbors, autoscale=False, adjust_for_stream=True)
        lengths = np.sqrt((V_grid ** 2).sum(0))
        linewidth = 1 if linewidth is None else linewidth
        linewidth *= 2 * lengths / lengths[~np.isnan(lengths)].max()

    scatter_kwargs = {"basis": basis, "perc": perc, "use_raw": use_raw, "sort_order": sort_order, "alpha": alpha,
                      "components": components, "projection": projection, "legend_loc": legend_loc, "groups": groups,
                      "legend_fontsize": legend_fontsize, "legend_fontweight": legend_fontweight, "palette": palette,
                      "color_map": color_map, "frameon": frameon, "title": title, "xlabel": xlabel, "ylabel": ylabel,
                      "right_margin": right_margin, "left_margin": left_margin, "colorbar": colorbar, "dpi": dpi,
                      "fontsize": fontsize, "show": False, "save": None}

    multikey = colors if len(colors) > 1 else layers if len(layers) > 1 else vkeys if len(vkeys) > 1 else None
    if multikey is not None:
        ncols = len(multikey) if ncols is None else min(len(multikey), ncols)
        nrows = int(np.ceil(len(multikey) / ncols))
        figsize = rcParams['figure.figsize'] if figsize is None else figsize
        for i, gs in enumerate(
                pl.GridSpec(nrows, ncols, pl.figure(None, (figsize[0] * ncols, figsize[1] * nrows), dpi=dpi))):
            if i < len(multikey):
                velocity_embedding_stream(adata, density=density, size=size, smooth=smooth, n_neighbors=n_neighbors,
                                          linewidth=linewidth, ax=pl.subplot(gs),
                                          color=colors[i] if len(colors) > 1 else color,
                                          layer=layers[i] if len(layers) > 1 else layer,
                                          vkey=vkeys[i] if len(vkeys) > 1 else vkey,
                                          X_grid=None if len(vkeys) > 1 else X_grid,
                                          V_grid=None if len(vkeys) > 1 else V_grid, **scatter_kwargs, **kwargs)
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

        density = 1 if density is None else density
        stream_kwargs = {"linewidth": linewidth, "density": 2 * density}
        stream_kwargs.update(kwargs)
        pl.streamplot(X_grid[0], X_grid[1], V_grid[0], V_grid[1], color='grey', zorder=3, **stream_kwargs)

        size = 4 * default_size(adata) if size is None else size
        ax = scatter(adata, layer=layer, color=color, size=size, ax=ax, zorder=0, **scatter_kwargs)

        if isinstance(save, str): savefig('' if basis is None else basis, dpi=dpi, save=save, show=show)
        if show: pl.show()
        else: return ax
