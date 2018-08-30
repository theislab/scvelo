from scanpy.api.pl import scatter
from sklearn.neighbors import NearestNeighbors
from scipy.stats import norm as normal
import numpy as np
import matplotlib.pyplot as pl


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


def velocity_embedding_grid(adata, basis='umap', vbasis='velocity', density=1, color=None,
                            min_mass=.5, smooth=.5, n_neighbors=None, principal_curve=False,
                            use_raw=True, sort_order=True, alpha=.2, groups=None, components=None, projection='2d',
                            legend_loc='right margin', legend_fontsize=None, legend_fontweight=None,
                            color_map=None, palette=None, frameon=False, right_margin=None, left_margin=None,
                            size=None, title=None, show=True, save=None, ax=None, **kwargs):
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
    if ax is None: ax = pl.figure(None, (14, 10), dpi=120).gca()

    scatter(adata, color=color, use_raw=use_raw, sort_order=sort_order, alpha=alpha, basis=basis,
            groups=groups, components=components, projection=projection, legend_loc=legend_loc,
            legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight, color_map=color_map,
            palette=palette, frameon=frameon, right_margin=right_margin, left_margin=left_margin,
            size=size, title=title, show=False, save=save, ax=ax)

    vbasis += '_' + basis
    if vbasis not in adata.obsm_keys():
        raise ValueError(
            'You need to run `tl.velocity_embedding` first to compute embedded velocity vectors.')

    X_emb = adata.obsm['X_' + basis][:, :2]
    V_emb = adata.obsm[vbasis]

    _kwargs = {"scale": .5, "width": .001, "color": 'black', "edgecolors": 'k', "headwidth": 4.5, "headlength": 5,
                   "headaxislength": 3, "linewidth": .2}
    _kwargs.update(kwargs)

    X, V = compute_velocity_on_grid(X_emb, V_emb, density, smooth, n_neighbors, min_mass)
    pl.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], angles='xy', scale_units='xy', **_kwargs)

    if principal_curve:
        curve = adata.uns['principal_curve']['projections']
        pl.plot(curve[:, 0], curve[:, 1], c="w", lw=6, zorder=1000000)
        pl.plot(curve[:, 0], curve[:, 1], c="k", lw=3, zorder=2000000)

    if show: pl.show()
    else: return ax